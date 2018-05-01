# coding=utf-8

"""
    vcfy.model
    ~~~~~~~~~~

    This module implements the model for generating variants based on the given
    parameters.

    :copyright: (c) 2018 by Ali Ghaffaari.
    :license: MIT, see LICENSE for more details.
"""

import os

from numpy import random, arange
import vcf

from . import util


def rnd_sv(locus, seq):
    """Get random structural variant for the given `locus` in the `seq`.

    Parameters:
        locus : int
            The position in the sequence for which a random SV is generated.
        seq : str
            The reference DNA sequence.

    Return:
        A tuple indicating the REF and ALT alleles.
    """
    ref = seq[locus].upper()
    if ref not in util.BASES:
        raise RuntimeError("invalid base character")

    alts = [b for b in util.BASES if b != ref]
    idx = random.choice(len(alts), 1)[0]
    return ref, alts[idx]


def synthesise(region, num, low=None, high=None, pmf=None):
    """Synthesise variants for the given region in the one-based range
    [low, high) using the probability model defined by mass function `pmf`.

    Parameters:
        region : Bio.SeqRecord.SeqRecord
            The BioPython's `SeqRecord`-like object of the region containing the
            region ID and its sequence.
        num : int
            The number of variants to be synthesised.
        low : int, optional
            The lower bound of the range in which the variants are simulated. It
            is assumed to be 1, if not provided.
        high : int, optional
            If provided, it is one above the upper bound of the range; otherwise
            it is set to the length of the region sequence.
        pmf : 1-D array-like, optional
            The probabilities associated with each entry in [low, high). If not
            given the sample assumes a uniform distribution over all entries in
            [low, high).
    Return:
        A dictionary containing values for these VCF fields: POS, ID, REF, ALT,
        QUAL, and FILTER.
    """
    low = 1 if low is None else max(1, low)
    high = len(region.seq) if high is None else min(high, len(region.seq))

    loci = random.choice(arange(low, high), num, replace=False, p=pmf)
    for locus in sorted(loci):
        ref, alt = rnd_sv(locus, region.seq)
        yield dict(POS=locus,
                   ID=util.VCF_MISSING_VALUE,
                   REF=ref,
                   ALT=alt,
                   QUAL=util.VCF_MISSING_VALUE,
                   FILTER=util.VCF_MISSING_VALUE)


def generate_vcf(ref, vcf_out, region_id=None, **synth_params):
    """Generate synthetic VCF file.

    Parameters:
        ref : str or opened file
            File path or object of the reference genome FASTA file.
        vcf_out : str or opened file
            File path or writable file object corresponding to output VCF file.
        region_id : str
            The ID of the region for which the variants are simulated.
        synth_params: keyword paramters
            The parameters required for simulation passed directly to the
            `synthesise` function.
    """
    if region_id is None:
        region = next(util.filter_regions(ref, n=1))
    else:
        region = next(util.filter_regions(ref, include=[region_id]))

    template_fpath = util.make_template(ref, region, **synth_params)
    template = vcf.Reader(open(template_fpath, 'r'))

    if isinstance(vcf_out, str):
        writer = vcf.Writer(open(vcf_out, 'w'), template)
    else:
        writer = vcf.Writer(vcf_out, template)

    tmpl_record = next(template)
    for svar in synthesise(region, **synth_params):
        record = util.update_record(tmpl_record, CHROM=region.id, **svar)
        writer.write_record(record)

    os.remove(template_fpath)
