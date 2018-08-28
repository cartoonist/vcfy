# coding=utf-8

"""
    vcfy.ksnper
    ~~~~~~~~~~~

    Report the number of SNPs in all k-mers.

    :copyright: (c) 2018 by Ali Ghaffaari.
    :license: MIT, see LICENSE for more details.
"""

import sys
import collections
import csv

import click
import vcf
from Bio import SeqIO
from BitVector import BitVector


def reflen(ref_file, region=None):
    """Get the length of reference genome sequence length.

    Args:
        ref_file : file-like object or str
            The reference genome file or file path.
        region : str
            The region in the reference genome. If it is not provided, the first
            region will be used.
    """
    parser = SeqIO.parse(ref_file, 'fasta')
    if region is None:
        return len(next(parser).seq)
    for r in parser:
        if r.name == region:
            return len(r.seq)
    raise RuntimeError("region not found")


def compute_snpbv(vcf_reader, length):
    """Compute SNP bitvector. The SNP bitvector is a bitvector of the length of
    the reference genome where SNP loci are set.

    Args:
        vcf_reader : vcf.Reader
            VCF reader.
        length : int
            The length of reference genome.

    Return:
        A BitVector.BitVector instance of size `length` where all SNP loci are
        set.
    """
    cnt = BitVector(size=length)
    for record in vcf_reader:
        cnt[record.POS-1] = 1
    return cnt


def ksnpcounts(snpbv, k):
    """Generate number of SNPs in all k-mers in the reference genome.

    Args:
        snpbv : BitVector.BitVector
            The SNP bitvector.
        k : int
            The length of the k-mer.

    Return:
        The number of SNPs in a k-mer. It yields this for all k-mers in the
        reference genome.
    """
    kcount = 0
    for i in range(k):
        kcount += snpbv[i]
    yield kcount
    for i in range(k, len(snpbv)):
        kcount += snpbv[i]
        kcount -= snpbv[i-k]
        yield kcount


def write_csv(output, vcf_file, ref_file, k, region=None, dialect='unix',
              compressed=None, frequency=False):
    """Write CSV file.

    Args:
        output : file-like object
            The output CSV file.
        vcf_file : file-like object or str
            The input VCF file or file path.
        ref_file : file-like object or str
            The input reference genome. If it is set to `None`, the reference
            genome path will be inferred from VCF metadata (header).
        k : int
            The length of the k-mer.
        region : str
            The region in the reference genome.
        dialect : str
            This string specifies the dialect of the output CSV file.
        compressed : bool
            Whether input VCF is compressed or not. It is determined by file
            extension if it is not specified.
        frequency : bool
            If set, instead of reporting SNPs counts of each k-mer, it reports
            the frequency of each counts.
    """
    fieldnames = ['k', 'count', 'freq'] if frequency else ['k', 'count']
    csv_writer = csv.DictWriter(output,
                                fieldnames=fieldnames,
                                dialect=dialect,
                                quoting=csv.QUOTE_NONE)
    csv_writer.writeheader()

    if isinstance(vcf_file, str):
        vcf_reader = vcf.Reader(filename=vcf_file, compressed=compressed)
    else:
        vcf_reader = vcf.Reader(vcf_file, compressed=compressed)
    if ref_file is None:
        ref_file = open(vcf_reader.metadata['reference'], 'r')
    bv = compute_snpbv(vcf_reader, reflen(ref_file, region))
    freqs = collections.defaultdict(int)
    for count in ksnpcounts(bv, k):
        if frequency:
            freqs[count] += 1
        else:
            csv_writer.writerow(dict(k=k, count=count))
    if not frequency:
        return
    for key, value in sorted(freqs.items()):
        csv_writer.writerow(dict(k=k, count=key, freq=value))


@click.command()
@click.argument('vcf', type=str, default="-")
@click.option('-o', '--output', type=click.File('w'), default="-",
              help="Write to this file instead of standard output.")
@click.option('-r', '--reference', type=click.File('r'), default=None,
              help=("Reference genome FASTA file. It will be inferred from VCF "
                    "header, if not specified."))
@click.option('-R', '--region', type=str, default=None,
              help="Use this region from reference genome.")
@click.option('-k', type=int, required=True, help="The value of k.")
@click.option('-c', is_flag=True, default=None,
              help="Set if the input VCF is compressed.")
@click.option('-d', '--dialect', type=click.Choice(csv.list_dialects()),
              default='unix', show_default=True,
              help="Use this CSV dialect.")
@click.option('-f', '--frequency', is_flag=True, default=False,
              help="Report frequency instead of reporting count in each k-mer.")
def cli(**kwargs):
    """Report the number of SNPs in all k-mers. Specify the k and the VCF file,
    it reports number of SNPS occurred in each k-mer.
    """
    stdin_fsock = sys.stdin.buffer if kwargs['c'] else sys.stdin
    write_csv(kwargs.pop('output'),
              kwargs.pop('vcf') if kwargs['vcf'] != '-' else stdin_fsock,
              kwargs.pop('reference'),
              kwargs.pop('k'),
              kwargs.pop('region'),
              kwargs.pop('dialect'),
              kwargs.pop('c'),
              kwargs.pop('frequency'))
