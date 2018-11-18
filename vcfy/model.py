# coding=utf-8

"""
    vcfy.model
    ~~~~~~~~~~

    This module implements the model for generating variants based on the given
    parameters.

    :copyright: (c) 2018 by Ali Ghaffaari.
    :license: MIT, see LICENSE for more details.
"""

from numpy import random, arange
import vcf

from . import util


def rnd_snp(locus, seq):
    """Get random SNP for the given `locus` in the `seq`.

    Parameters:
        locus : int
            The position in the sequence for which a random SNP is generated.
        seq : str
            The reference DNA sequence.

    Return:
        A tuple indicating the REF and ALT alleles.
    """
    ref = seq[locus-1].upper()
    if ref not in util.BASES:
        raise RuntimeError("invalid base character")

    return ref, random.choice([b for b in util.BASES if b != ref])


def rnd_indel(locus, seq, insertion=False, length=1):
    """Get insertion or deletion for the given `locus` in the `seq`.

    Parameters:
        locus : int
            The position in the sequence for which an indel is generated.
        seq : str
            The reference DNA sequence.
        insertion : bool
            Generate and insertion if set; otherwise deletion.
        length : int
            Length of the indel.

    Return:
        A tuple indicating the REF and ALT alleles.
    """
    if insertion:
        alt = ''.join(random.choice(util.BASES) for _ in range(length))
        return util.VCF_MISSING_VALUE, [alt]
    else:
        alt = seq[locus-1:locus-1+length].upper()
        return alt, util.VCF_MISSING_VALUE


def simulate(region, mrate, indrate, extrate, low=None, high=None):
    """Simulate variants for the given region in the one-based range
    [low, high) using the probability model defined by mass function `pmf`.

    Parameters:
        region : Bio.SeqRecord.SeqRecord
            The BioPython's `SeqRecord`-like object of the region containing the
            region ID and its sequence.
        mrate : float
            Base mutation rate.
        indrate : float
            Indel fraction rate. If a locus is decided to be a mutation site, it
            would be indel with this probability.
        extrate : float
            Indel extension rate. If a locus is decided to be a indel site, it
            would be extended with this probability.
        low : int, optional
            The lower bound of the range in which the variants are simulated. It
            is assumed to be 1, if not provided.
        high : int, optional
            If provided, it is one above the upper bound of the range; otherwise
            it is set to the length of the region sequence.
    Return:
        A dictionary containing values for these VCF fields: POS, ID, REF, ALT,
        QUAL, and FILTER.
    """
    low = 1 if low is None else max(1, low)
    high = len(region.seq) + 1 if high is None else min(high, len(region.seq)+1)

    indel_len = 0
    for locus in arange(low, high):
        if indel_len != 0:
            indel_len -= 1
            continue
        if util.toss(mrate):
            if util.toss(indrate):  # indel
                indel_len = 1
                while (util.toss(extrate) and
                       indel_len < util.max_indel_len(high - low)):  # extension
                    indel_len += 1
                insertion = util.toss()
                ref, alt = rnd_indel(locus, region.seq, insertion, indel_len)
                yield dict(POS=locus,
                           ID=util.VCF_MISSING_VALUE,
                           REF=ref,
                           ALT=alt,
                           QUAL=util.VCF_MISSING_VALUE,
                           FILTER=util.VCF_MISSING_VALUE)
                indel_len = 0 if insertion else indel_len - 1
            else:  # substitution
                try:
                    ref, alt = rnd_snp(locus, region.seq)
                except RuntimeError as err:
                    util.warn(err)
                    continue
                yield dict(POS=locus,
                           ID=util.VCF_MISSING_VALUE,
                           REF=ref,
                           ALT=alt,
                           QUAL=util.VCF_MISSING_VALUE,
                           FILTER=util.VCF_MISSING_VALUE)


def generate_vcf(ref, vcf_out, region_id=None, **sim_params):
    """Generate simulated VCF file.

    Parameters:
        ref : str or opened file
            File path or object of the reference genome FASTA file.
        vcf_out : str or opened file
            File path or writable file object corresponding to output VCF file.
        region_id : str
            The ID of the region for which the variants are simulated.
        sim_params: keyword paramters
            The parameters required for simulation passed directly to the
            `simulate` function.
    """
    if region_id is None:
        region = next(util.filter_regions(ref, n=1))
    else:
        region = next(util.filter_regions(ref, include=[region_id]))

    template = vcf.Reader(util.make_template(ref, region, **sim_params))

    if isinstance(vcf_out, str):
        writer = vcf.Writer(open(vcf_out, 'w'), template)
    else:
        writer = vcf.Writer(vcf_out, template)

    tmpl_record = next(template)
    for svar in simulate(region, **sim_params):
        record = util.update_record(tmpl_record, CHROM=region.id, **svar)
        writer.write_record(record)
