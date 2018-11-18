# coding=utf-8

"""
    vcfy.cli
    ~~~~~~~~

    This module implements command-line interface.

    :copyright: (c) 2018 by Ali Ghaffaari.
    :license: MIT, see LICENSE for more details.
"""

import click

from . import model


@click.command()
@click.argument('reference', type=click.File('r'))
@click.option('-o', '--output', default="-", type=click.File('w'),
              help="Write to this file instead of standard output.")
@click.option('-m', '--mutation-rate', type=float, required=True,
              help="Base mutation rate.")
@click.option('-i', '--indel-rate', type=float, default=0.0,
              help="Fraction of indels.")
@click.option('-e', '--indel-extension', type=float, default=0.0,
              help="Probability an indel is extended.")
@click.option('-r', '--region', type=str, default=None,
              help="Region ID (default=first region in the reference)")
@click.option('-l', '--low', type=int, default=None,
              help="Range lower bound (default=first locus in the region)")
@click.option('-h', '--high', type=int, default=None,
              help="Range upper bound (default=last locus in the region)")
def cli(**kwargs):
    """Generate VCF file with simulated variants in specified range [low, high)
    for the given region of the reference genome. In case that the region is not
    specified the first region is used. If no range is provided, it is assumed
    that the variants are scattered throughout the region.

    For more information, consult with the README file.
    """
    model.generate_vcf(
        ref=kwargs.pop('reference'),
        vcf_out=kwargs.pop('output'),
        region_id=kwargs.pop('region'),
        mrate=kwargs.pop('mutation_rate'),
        indrate=kwargs.pop('indel_rate'),
        extrate=kwargs.pop('indel_extension'),
        **kwargs)
