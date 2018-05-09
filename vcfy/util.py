# coding=utf-8

"""
    vcfy.util
    ~~~~~~~~~

    The utility module for helper classes and functions.

    :copyright: (c) by 2018 by Ali Ghaffaari.
    :license: MIT, see LICENSE for more details.
"""

import os
import copy
import datetime
import pkg_resources as pkg_res

from Bio import SeqIO

from . import release


BASES = ['A', 'C', 'G', 'T']
VCF_MISSING_VALUE = '.'
VCF_TEMPLATE_PATH = 'resources/templates/template.vcf'


def make_template(ref, region, **params):
    """Fill wildcards in the templated VCF file and create a ready-to-use
    template.

    Parameters:
        ref : str or opened file
            File path or object of the reference genome FASTA file.
        region : SeqIORecord
            The region FASTA record.

    Return:
        The file path of newly created template VCF file.
    """
    wildcards = dict()
    wildcards['date'] = datetime.datetime.today().strftime('%Y%m%d')
    wildcards['version'] = release.__version__
    wildcards['reference'] = os.path.realpath(ref.name if hasattr(ref, 'name')
                                              else ref)
    wildcards['contigs'] = ('<ID=' + region.id +
                            ',length=' + str(len(region.seq)) + '>')
    opts = ['vcfy']
    opts.append('-m ' + str(params['mrate']))
    opts.append("-r '" + region.id + "'")
    opts.append('-l ' + str(params['low']) if params['low'] else '')
    opts.append('-h ' + str(params['high']) if params['high'] else '')
    wildcards['cmd'] = ' '.join(o for o in opts if o)
    tmpl = open(pkg_res.resource_filename(__name__, VCF_TEMPLATE_PATH), 'r')
    new_tmpl = open('.template.vcf', 'w')
    for line in tmpl:
        new_tmpl.write(line.format(**wildcards))
    return new_tmpl.name


def filter_regions(fasta_file, **kwargs):
    """Yield the SeqIORecord of the regions in the FASTA file filtered based on
    the given parameters. In case of no filter parameters are provided, all
    regions is reported.

    Parameters:
        fasta_file : str or opened file
            File path or object of the input FASTA file.
        include : list of str, optional
            List of region IDs to be *included* to the return dictionary (will
            be overriden by `exclude` if there is any conflict between these two
            lists).
        exclude : list of str, optional
            List of region IDs to be *excluded* from the return dictionary.
        n : int, optional
            If provided, the first `n` number of regions are returned; otherwise
            all regions are be returned if neither `include` nor `exclude` is
            provided. If this parameter is set to negative values, all regions
            are returned.

    Return:
        A dictionary containing DNA sequence of each region with its ID as key.
    """
    include = kwargs.pop('include', [])
    exclude = kwargs.pop('exclude', [])
    num = kwargs.pop('n', -1) if not include and not exclude else -1

    fasta = SeqIO.parse(fasta_file, 'fasta')
    for region in fasta:
        if num == 0:
            break
        if region.id in exclude or (include and region.id not in include):
            continue
        yield region
        num -= 1


def update_record(old, **kwargs):
    """Get a copy of the given VCF record with updated fields with new values.

    Parameters:
        old : vcf.model._Record
            The record to be update.
        kwargs : dict
            The new values for each field provided by keyword parameters.

    Return:
        The updated copy of the record.
    """
    record = copy.deepcopy(old)
    for key, value in kwargs.items():
        setattr(record, key, value)
    return record
