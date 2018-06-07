# coding=utf-8

"""
    vcfy.ksnper
    ~~~~~~~~~~~

    Report the number of SNPs in all k-mers.

    :copyright: (c) 2018 by Ali Ghaffaari.
    :license: MIT, see LICENSE for more details.
"""

import csv

import click
import vcf
from Bio import SeqIO
from BitVector import BitVector


def reflen(ref_file):
    """Get the length of reference genome sequence length.

    Args:
        ref_file : file-like object or str
            The reference genome file or file path.
    """
    return len(next(SeqIO.parse(ref_file, 'fasta')).seq)


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


def write_csv(output, vcf_file, ref_file, k, dialect='unix'):
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
        dialect : str
            This string specifies the dialect of the output CSV file.
    """
    csv_writer = csv.DictWriter(output,
                                fieldnames=['k', 'count'],
                                dialect=dialect,
                                quoting=csv.QUOTE_NONE)
    csv_writer.writeheader()

    vcf_reader = vcf.Reader(vcf_file)
    if ref_file is None:
        ref_file = open(vcf_reader.metadata['reference'], 'r')
    bv = compute_snpbv(vcf_reader, reflen(ref_file))
    for count in ksnpcounts(bv, k):
        csv_writer.writerow(dict(k=k, count=count))


@click.command()
@click.argument('vcf', type=click.File('r'), default="-")
@click.option('-o', '--output', type=click.File('w'), default="-",
              help="Write to this file instead of standard output.")
@click.option('-r', '--reference', type=click.File('r'), default=None,
              help=("Reference genome FASTA file. It will be inferred from VCF "
                    "header, if not specified."))
@click.option('-k', type=int, required=True, help="The value of k.")
@click.option('-d', '--dialect', type=click.Choice(csv.list_dialects()),
              default='unix', show_default=True,
              help="Use this CSV dialect.")
def cli(**kwargs):
    write_csv(kwargs.pop('output'),
              kwargs.pop('vcf'),
              kwargs.pop('reference'),
              kwargs.pop('k'),
              kwargs.pop('dialect'))
