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
@click.argument('output', type=click.File('w'))
@click.option('-n', '--num', type=int, required=True, help="Number of variants")
@click.option('-r', '--region', type=str, default=None, help="Region ID")
@click.option('-l', '--low', type=int, default=None, help="Range lower bound")
@click.option('-h', '--high', type=int, default=None, help="Range upper bound")
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
