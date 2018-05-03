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
@click.option('-n', '--num', type=int, required=True,
              help="Number of variants to be simulated.")
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
        kwargs.pop('reference'),
        kwargs.pop('output'),
        kwargs.pop('region'),
        **kwargs)
