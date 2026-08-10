# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Assess the types of features having expression annotations in chado.

Author(s):
    Gil dos Santos dossantos@morgan.harvard.edu

Usage:
    report_expression_feature_categories.py [-h] [-v VERBOSE] [-c CONFIG]

Example:
    python report_expression_feature_categories.py -v -c /path/to/config.cfg

Notes:
    Every current feature having at least one feature_expression annotation is
    called a "gene product" here, and is assigned zero to many category labels
    that describe how (if at all) that gene product can be traced to a current
    Dmel gene. Categories are NOT mutually exclusive: a gene product having no
    category gets the "UNKNOWN" label. The script writes two files:
    1. A log file documenting the run and reporting bulk counts: counts of gene
       products and expression annotations per category, all pairwise category
       overlaps, gene products tracing to MANY genes, and gene products
       representing aberrant gene expression (the "ABERRANT" category).
    2. A TSV report file listing each gene product and its category labels.
    The categorization SQL follows the category queries in
    "docs/plans/2026-08-07-mapping_transgenic_xprn.md", with these differences:
    1) each query returns gene product feature_ids and related gene uniquenames
    instead of aggregate counts, so that category overlaps can be assessed;
    2) the "gene product has expression annotation" filter is expressed as an
    EXISTS clause rather than a feature_expression join, which is equivalent for
    this purpose but avoids row multiplication; 3) the ABERRANT query also requires
    that the gene be a Dmel gene, which a curation invariant check confirms to be a
    no-op; 4) in the five categories whose gene product is "associated_with" an FBal
    allele, the gene product is identified by its FBtr/FBpp uniquename alone, with no
    "](R|P)A$" name restriction. The name restriction is retained for CURATED and
    ABERRANT, where curation guarantees it holds anyway. Expression
    annotation counts are taken per gene product from a single feature_expression
    query, so annotation counts for any set of gene products (a category, or an
    intersection of categories) are just sums over the gene products in that set.

"""

import argparse
import re
from itertools import combinations
from harvdev_utils.general_functions import (
    generic_FB_tsv_dict, tsv_report_dump
)
from harvdev_utils.psycopg_functions import (
    set_up_db_reading, connect
)

# Global variables for the output file. Header order will match list order below.
REPORT_LABEL = 'expression_feature_categories'
REPORT_TITLE = 'FlyBase Expression Annotation Feature Category Report'
HEADER_LIST = [
    'feature_id',
    'uniquename',
    'name',
    'categories',
    'n_categories',
    'n_expression_annotations',
    'n_dmel_genes',
    'dmel_genes',
]

# The label given to gene products belonging to none of the categories below.
UNKNOWN_LABEL = 'UNKNOWN'

# Categories of gene product having expression annotations.
# Each category has a label, a description, and an SQL query returning two columns:
# 1. the gene product's feature_id;
# 2. the uniquename of the current Dmel gene reached by that category's path (NULL if the path reaches no gene).
CATEGORIES = [
    {
        'label': 'SPLIT_COMBO',
        'description': 'The annotation subject is the result of a split system combination; by design, no Dmel gene is traced for these.',
        'sql': """
            SELECT DISTINCT gp.feature_id, NULL::VARCHAR AS gene_uniquename
            FROM feature gp
            WHERE gp.is_obsolete IS FALSE
              AND gp.uniquename ~ '^FBco[0-9]{7}$'
              AND EXISTS (SELECT 1 FROM feature_expression fe WHERE fe.feature_id = gp.feature_id);
        """,
    },
    {
        'label': 'ENDO',
        'description': 'A gene product directly "associated_with" a Dmel gene.',
        'sql': """
            SELECT DISTINCT gp.feature_id, g.uniquename
            FROM feature gp
            JOIN feature_relationship gpg ON gpg.subject_id = gp.feature_id
              AND gpg.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'associated_with')
            JOIN feature g ON g.feature_id = gpg.object_id
            WHERE gp.is_obsolete IS FALSE
              AND gp.name ~ '-X(R|P)$'
              AND g.is_obsolete IS FALSE
              AND g.uniquename ~ '^FBgn[0-9]{7}$'
              AND g.organism_id = 1
              AND EXISTS (SELECT 1 FROM feature_expression fe WHERE fe.feature_id = gp.feature_id);
        """,
    },
    {
        'label': 'ENDO_ISO',
        'description': 'A gene product "isa" gene product that is directly "associated_with" a Dmel gene.',
        'sql': """
            SELECT DISTINCT gp.feature_id, g.uniquename
            FROM feature gp
            JOIN feature_relationship gpg ON gpg.subject_id = gp.feature_id
              AND gpg.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'isa')
            JOIN feature gp2 ON gp2.feature_id = gpg.object_id
            JOIN feature_relationship gp2g ON gp2g.subject_id = gp2.feature_id
              AND gp2g.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'associated_with')
            JOIN feature g ON g.feature_id = gp2g.object_id
            WHERE gp.is_obsolete IS FALSE
              AND gp2.is_obsolete IS FALSE
              AND gp2.name ~ '-X(R|P)$'
              AND g.is_obsolete IS FALSE
              AND g.uniquename ~ '^FBgn[0-9]{7}$'
              AND g.organism_id = 1
              AND EXISTS (SELECT 1 FROM feature_expression fe WHERE fe.feature_id = gp.feature_id);
        """,
    },
    {
        'label': 'TAGGED',
        'description': 'A gene product is directly "associated_with" a Dmel allele, which is in turn directly "alleleof" a Dmel gene.',
        'sql': """
            SELECT DISTINCT gp.feature_id, g.uniquename
            FROM feature gp
            JOIN feature_relationship gpa ON gpa.subject_id = gp.feature_id
              AND gpa.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'associated_with')
            JOIN feature a ON a.feature_id = gpa.object_id
            JOIN feature_relationship ag ON ag.subject_id = a.feature_id
              AND ag.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'alleleof')
            JOIN feature g ON g.feature_id = ag.object_id
            WHERE gp.is_obsolete IS FALSE
              AND gp.uniquename ~ '^FB(tr|pp)[0-9]{7}$'
              AND a.is_obsolete IS FALSE
              AND a.uniquename ~ '^FBal[0-9]{7}$'
              AND a.organism_id = 1
              AND g.is_obsolete IS FALSE
              AND g.uniquename ~ '^FBgn[0-9]{7}$'
              AND EXISTS (SELECT 1 FROM feature_expression fe WHERE fe.feature_id = gp.feature_id);
        """,
    },
    {
        'label': 'PROMOTER_CHAR_GENE',
        'description': "A transgenic non-Dmel allele has a Dmel gene's regulatory region.",
        'sql': """
            SELECT DISTINCT gp.feature_id, g.uniquename
            FROM feature gp
            JOIN feature_relationship gpa ON gpa.subject_id = gp.feature_id
              AND gpa.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'associated_with')
            JOIN feature a ON a.feature_id = gpa.object_id
            JOIN feature_relationship ag ON ag.subject_id = a.feature_id
              AND ag.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'has_reg_region')
            JOIN feature g ON g.feature_id = ag.object_id
            WHERE gp.is_obsolete IS FALSE
              AND gp.uniquename ~ '^FB(tr|pp)[0-9]{7}$'
              AND a.is_obsolete IS FALSE
              AND a.uniquename ~ '^FBal[0-9]{7}$'
              AND a.organism_id != 1
              AND g.is_obsolete IS FALSE
              AND g.uniquename ~ '^FBgn[0-9]{7}$'
              AND g.organism_id = 1
              AND EXISTS (SELECT 1 FROM feature_expression fe WHERE fe.feature_id = gp.feature_id);
        """,
    },
    {
        'label': 'PROMOTER_CHAR_FBSF',
        'description': 'A transgenic non-Dmel allele has a Dmel regulatory region that is associated with a gene.',
        'sql': """
            SELECT DISTINCT gp.feature_id, g.uniquename
            FROM feature gp
            JOIN feature_relationship gpa ON gpa.subject_id = gp.feature_id
              AND gpa.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'associated_with')
            JOIN feature a ON a.feature_id = gpa.object_id
            JOIN feature_relationship asf ON asf.subject_id = a.feature_id
              AND asf.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'has_reg_region')
            JOIN feature sf ON sf.feature_id = asf.object_id
            JOIN feature_relationship sfg ON sfg.subject_id = sf.feature_id
              AND sfg.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'associated_with')
            JOIN feature g ON g.feature_id = sfg.object_id
            WHERE gp.is_obsolete IS FALSE
              AND gp.uniquename ~ '^FB(tr|pp)[0-9]{7}$'
              AND a.is_obsolete IS FALSE
              AND a.uniquename ~ '^FBal[0-9]{7}$'
              AND a.organism_id != 1
              AND sf.is_obsolete IS FALSE
              AND sf.uniquename ~ '^FBsf[0-9]{10}$'
              AND g.is_obsolete IS FALSE
              AND g.uniquename ~ '^FBgn[0-9]{7}$'
              AND g.organism_id = 1
              AND EXISTS (SELECT 1 FROM feature_expression fe WHERE fe.feature_id = gp.feature_id);
        """,
    },
    {
        'label': 'INS_TRAP_KNOWN',
        'description': 'An insertion traps nearby regulatory elements; there is an indirectly associated Dmel allele and gene.',
        'sql': """
            SELECT DISTINCT gp.feature_id, g.uniquename
            FROM feature gp
            JOIN feature_relationship gpa ON gpa.subject_id = gp.feature_id
              AND gpa.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'associated_with')
            JOIN feature a ON a.feature_id = gpa.object_id
            JOIN feature_relationship ai ON ai.subject_id = a.feature_id
              AND ai.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'associated_with')
            JOIN feature i ON i.feature_id = ai.object_id
            JOIN feature_relationship ia ON ia.object_id = i.feature_id
              AND ia.subject_id != ai.subject_id
              AND ia.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'associated_with')
            JOIN feature a2 ON a2.feature_id = ia.subject_id AND a2.feature_id != a.feature_id
            JOIN feature_relationship a2g ON a2g.subject_id = a2.feature_id
              AND a2g.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'alleleof')
            JOIN feature g ON g.feature_id = a2g.object_id
            WHERE gp.is_obsolete IS FALSE
              AND gp.uniquename ~ '^FB(tr|pp)[0-9]{7}$'
              AND a.is_obsolete IS FALSE
              AND a.uniquename ~ '^FBal[0-9]{7}$'
              AND a.organism_id != 1
              AND i.is_obsolete IS FALSE
              AND i.uniquename ~ '^FBti[0-9]{7}$'
              AND i.organism_id = 1
              AND a2.is_obsolete IS FALSE
              AND a2.uniquename ~ '^FBal[0-9]{7}$'
              AND a2.organism_id = 1
              AND g.is_obsolete IS FALSE
              AND g.uniquename ~ '^FBgn[0-9]{7}$'
              AND g.organism_id = 1
              AND EXISTS (SELECT 1 FROM feature_expression fe WHERE fe.feature_id = gp.feature_id);
        """,
    },
    {
        'label': 'INS_TRAP_UNK',
        'description': 'An insertion traps nearby regulatory elements, but there is no known Dmel allele/gene for that insertion.',
        'sql': """
            SELECT DISTINCT gp.feature_id, NULL::VARCHAR AS gene_uniquename
            FROM feature gp
            JOIN feature_relationship gpa ON gpa.subject_id = gp.feature_id
              AND gpa.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'associated_with')
            JOIN feature a ON a.feature_id = gpa.object_id
            JOIN feature_relationship ai ON ai.subject_id = a.feature_id
              AND ai.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'associated_with')
            JOIN feature i ON i.feature_id = ai.object_id
            WHERE gp.is_obsolete IS FALSE
              AND gp.uniquename ~ '^FB(tr|pp)[0-9]{7}$'
              AND a.is_obsolete IS FALSE
              AND a.uniquename ~ '^FBal[0-9]{7}$'
              AND a.organism_id != 1
              AND i.is_obsolete IS FALSE
              AND i.uniquename ~ '^FBti[0-9]{7}$'
              AND i.organism_id = 1
              AND EXISTS (SELECT 1 FROM feature_expression fe WHERE fe.feature_id = gp.feature_id)
              AND i.feature_id NOT IN
              (
                  SELECT DISTINCT i2.feature_id
                  FROM feature a2
                  JOIN feature_relationship a2i2 ON a2i2.subject_id = a2.feature_id
                    AND a2i2.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'associated_with')
                  JOIN feature i2 ON i2.feature_id = a2i2.object_id
                  WHERE a2.is_obsolete IS FALSE
                    AND a2.uniquename ~ '^FBal[0-9]{7}$'
                    AND a2.organism_id = 1
                    AND i2.is_obsolete IS FALSE
                    AND i2.uniquename ~ '^FBti[0-9]{7}$'
                    AND i2.organism_id = 1
              );
        """,
    },
    {
        'label': 'CURATED',
        'description': 'A curator has attributed the gene product to the expression of a Dmel gene, without flagging that attribution as aberrant.',
        'sql': """
            SELECT DISTINCT gp.feature_id, g.uniquename
            FROM feature gp
            JOIN feature_relationship gpg ON gpg.subject_id = gp.feature_id
              AND gpg.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'attributed_as_expression_of')
            JOIN feature g ON g.feature_id = gpg.object_id
            WHERE gp.is_obsolete IS FALSE
              AND gp.name ~ '](R|P)A$'
              AND g.is_obsolete IS FALSE
              AND g.uniquename ~ '^FBgn[0-9]{7}$'
              AND g.organism_id = 1
              AND NOT EXISTS
              (
                  SELECT 1
                  FROM feature_relationshipprop frp
                  WHERE frp.feature_relationship_id = gpg.feature_relationship_id
                    AND frp.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'is_relative_wildtype')
              )
              AND EXISTS (SELECT 1 FROM feature_expression fe WHERE fe.feature_id = gp.feature_id);
        """,
    },
    {
        'label': 'ABERRANT',
        'description': 'A curator has flagged the gene product as an aberrant representation of the expression of a Dmel gene.',
        'sql': """
            SELECT DISTINCT gp.feature_id, g.uniquename
            FROM feature gp
            JOIN feature_relationship gpg ON gpg.subject_id = gp.feature_id
              AND gpg.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'attributed_as_expression_of')
            JOIN feature g ON g.feature_id = gpg.object_id
            JOIN feature_relationshipprop frp ON frp.feature_relationship_id = gpg.feature_relationship_id
              AND frp.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'is_relative_wildtype')
            WHERE gp.is_obsolete IS FALSE
              AND gp.name ~ '](R|P)A$'
              AND g.is_obsolete IS FALSE
              AND g.uniquename ~ '^FBgn[0-9]{7}$'
              AND g.organism_id = 1
              AND EXISTS (SELECT 1 FROM feature_expression fe WHERE fe.feature_id = gp.feature_id);
        """,
    },
]

# Supplementary counts reported in the log only.
COUNT_QUERIES = [
    {
        'question': 'Obsolete features having expression annotations (excluded from this report): (features, annotations)',
        'sql': """
            SELECT COUNT(DISTINCT f.feature_id), COUNT(DISTINCT fe.feature_expression_id)
            FROM feature f
            JOIN feature_expression fe ON fe.feature_id = f.feature_id
            WHERE f.is_obsolete IS TRUE;
        """,
    },
]

# Checks on the curation invariants that the CURATED and ABERRANT queries above depend upon.
# Each query is expected to return no rows: any row returned means the invariant no longer holds.
INVARIANT_QUERIES = [
    {
        'invariant': 'Every "is_relative_wildtype" feature_relationshipprop has value = "y"',
        'sql': """
            SELECT DISTINCT frp.value
            FROM feature_relationshipprop frp
            WHERE frp.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'is_relative_wildtype')
              AND (frp.value IS DISTINCT FROM 'y');
        """,
    },
    {
        'invariant': 'Every gene product "attributed_as_expression_of" a gene has a "](R|P)A$" name ending',
        'sql': """
            SELECT DISTINCT gp.uniquename, gp.name
            FROM feature gp
            JOIN feature_relationship gpg ON gpg.subject_id = gp.feature_id
              AND gpg.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'attributed_as_expression_of')
            WHERE gp.is_obsolete IS FALSE
              AND gp.name !~ '](R|P)A$'
              AND EXISTS (SELECT 1 FROM feature_expression fe WHERE fe.feature_id = gp.feature_id);
        """,
    },
    {
        'invariant': 'Every current gene "attributed_as_expression_of" by a gene product is a Dmel gene',
        'sql': """
            SELECT DISTINCT g.uniquename, g.name
            FROM feature gp
            JOIN feature_relationship gpg ON gpg.subject_id = gp.feature_id
              AND gpg.type_id IN (SELECT DISTINCT cvterm_id FROM cvterm WHERE name = 'attributed_as_expression_of')
            JOIN feature g ON g.feature_id = gpg.object_id
            WHERE gp.is_obsolete IS FALSE
              AND g.is_obsolete IS FALSE
              AND g.uniquename ~ '^FBgn[0-9]{7}$'
              AND g.organism_id != 1
              AND EXISTS (SELECT 1 FROM feature_expression fe WHERE fe.feature_id = gp.feature_id);
        """,
    },
]

# The universe of features to be categorized: current features having expression annotations.
GENE_PRODUCT_QUERY = """
    SELECT f.feature_id,
           f.uniquename,
           f.name,
           cvt.name,
           o.abbreviation,
           COUNT(DISTINCT fe.feature_expression_id)
    FROM feature f
    JOIN cvterm cvt ON cvt.cvterm_id = f.type_id
    JOIN organism o ON o.organism_id = f.organism_id
    JOIN feature_expression fe ON fe.feature_id = f.feature_id
    WHERE f.is_obsolete IS FALSE
    GROUP BY f.feature_id, f.uniquename, f.name, cvt.name, o.abbreviation;
"""

# Proceed with generic setup.
set_up_dict = set_up_db_reading(REPORT_LABEL)
DATABASE = set_up_dict['database']
DATABASE_RELEASE = set_up_dict['database_release']
OUTPUT_FILENAME = set_up_dict['output_filename']
log = set_up_dict['log']
CONN = set_up_dict['conn']

# Process more input parameters (-c and -v handled by set_up_db_reading() function above).
parser = argparse.ArgumentParser(description='inputs')
# Use parse_known_args(), not parse_args(), to handle args specific to this script (outside of set_up_db_reading()).
args, extra_args = parser.parse_known_args()
log.info('Parsing args specific to this script; ignoring these: {}'.format(extra_args))


# Basic process of the script.
def main():
    """Retrieve, repackage and print out database information."""
    log.info('Started main function.')
    log.info('Assessing feature_expression annotations in database {} (release {}).'.format(DATABASE, DATABASE_RELEASE))
    gene_products = get_gene_products()
    category_members = get_category_members(gene_products)
    report_feature_types(gene_products)
    report_category_counts(gene_products, category_members)
    report_category_overlaps(gene_products, category_members)
    report_gene_mapping(gene_products, category_members)
    report_supplementary_counts()
    report_invariants()
    data_to_export_as_tsv = generic_FB_tsv_dict(REPORT_TITLE, DATABASE)
    data_to_export_as_tsv['metaData']['note'] = get_report_notes()
    data_to_export_as_tsv['data'] = process_gene_products(gene_products)
    tsv_report_dump(data_to_export_as_tsv, OUTPUT_FILENAME, headers=HEADER_LIST)
    CONN.close()
    log.info('Ended main function.')


# BELOW: Functions for retrieval and processing of specific data types.
def get_gene_products():
    """Retrieve current features having expression annotations.

    Returns:
        A dict of gene product dicts, keyed by feature_id.

    """
    log.info('Querying database for current features having expression annotations.')
    ret_gene_products = connect(GENE_PRODUCT_QUERY, 'no_query', CONN)
    gene_products = {}
    for row in ret_gene_products:
        gene_products[row[0]] = {
            'feature_id': row[0],
            'uniquename': row[1],
            'name': row[2],
            'type': row[3],
            'organism': row[4],
            'n_annotations': row[5],
            'categories': [],
            'genes': set(),
        }
    n_annotations = count_annotations(gene_products, gene_products.keys())
    log.info('Found {} current features having {} expression annotations.'.format(len(gene_products), n_annotations))

    return gene_products


def get_category_members(gene_products):
    """Assign category labels and related Dmel genes to gene products.

    Args:
        gene_products (dict): A dict of gene product dicts, keyed by feature_id.

    Returns:
        A dict of sets of gene product feature_ids, keyed by category label.

    """
    category_members = {}
    for category in CATEGORIES:
        label = category['label']
        log.info('Querying database for gene products of category {}.'.format(label))
        ret_category_info = connect(category['sql'], 'no_query', CONN)
        category_members[label] = set()
        unexpected_feature_ids = set()
        for row in ret_category_info:
            (feature_id, gene_uniquename) = (row[0], row[1])
            if feature_id not in gene_products.keys():
                unexpected_feature_ids.add(feature_id)
                continue
            # A gene product can be returned in many rows of the same query (one row per related gene).
            if feature_id not in category_members[label]:
                category_members[label].add(feature_id)
                gene_products[feature_id]['categories'].append(label)
            if gene_uniquename is not None:
                gene_products[feature_id]['genes'].add(gene_uniquename)
        log.info('Found {} gene products of category {}.'.format(len(category_members[label]), label))
        if unexpected_feature_ids:
            log.error('For category {}, found {} feature_ids absent from the set of current features having expression annotations.'
                      .format(label, len(unexpected_feature_ids)))
    # Gene products belonging to no category above get the catch-all label.
    category_members[UNKNOWN_LABEL] = set()
    for gene_product in gene_products.values():
        if not gene_product['categories']:
            gene_product['categories'].append(UNKNOWN_LABEL)
            category_members[UNKNOWN_LABEL].add(gene_product['feature_id'])
    log.info('Found {} gene products of category {}.'.format(len(category_members[UNKNOWN_LABEL]), UNKNOWN_LABEL))

    return category_members


def count_annotations(gene_products, feature_ids):
    """Sum up expression annotations for a set of gene products.

    Args:
        gene_products (dict): A dict of gene product dicts, keyed by feature_id.
        feature_ids (iterable): The feature_ids of the gene products to be counted.

    Returns:
        An int representing the number of expression annotations for these gene products.

    """
    return sum([gene_products[i]['n_annotations'] for i in feature_ids])


def report_group_counts(gene_products, group_name, groups):
    """Log gene product and annotation counts for named groups of gene products.

    Args:
        gene_products (dict): A dict of gene product dicts, keyed by feature_id.
        group_name (str): A description of the sort of group being counted: e.g., "category".
        groups (dict): A dict of sets of gene product feature_ids, keyed by group name.

    """
    log.info('{:<48}{:>12}{:>16}'.format(group_name.upper(), 'GENE_PRODUCTS', 'ANNOTATIONS'))
    for group_label, feature_ids in groups.items():
        log.info('{:<48}{:>12}{:>16}'.format(group_label, len(feature_ids), count_annotations(gene_products, feature_ids)))

    return


def report_feature_types(gene_products):
    """Log counts of the types of features having expression annotations.

    Args:
        gene_products (dict): A dict of gene product dicts, keyed by feature_id.

    """
    log.info('##### FEATURE TYPES HAVING EXPRESSION ANNOTATIONS #####')
    groupings = {
        'FlyBase ID type': 'id_type',
        'chado feature type': 'type',
        'organism': 'organism',
    }
    for grouping_name, grouping_key in groupings.items():
        groups = {}
        for gene_product in gene_products.values():
            if grouping_key == 'id_type':
                id_type_match = re.search(r'^(FB[a-z]{2})[0-9]+$', gene_product['uniquename'])
                group_label = id_type_match.group(1) if id_type_match else 'other'
            else:
                group_label = gene_product[grouping_key]
            try:
                groups[group_label].add(gene_product['feature_id'])
            except KeyError:
                groups[group_label] = {gene_product['feature_id']}
        sorted_groups = dict(sorted(groups.items(), key=lambda i: len(i[1]), reverse=True))
        report_group_counts(gene_products, 'gene products by {}'.format(grouping_name), sorted_groups)

    return


def report_category_counts(gene_products, category_members):
    """Log gene product and annotation counts per category.

    Args:
        gene_products (dict): A dict of gene product dicts, keyed by feature_id.
        category_members (dict): A dict of sets of gene product feature_ids, keyed by category label.

    """
    log.info('##### GENE PRODUCTS AND ANNOTATIONS PER CATEGORY #####')
    log.info('Categories are not mutually exclusive, so these counts sum to more than the totals reported above.')
    report_group_counts(gene_products, 'category', category_members)
    # Categories per gene product: how often is a gene product multiply categorized?
    counts_by_n_categories = {}
    for gene_product in gene_products.values():
        n_categories = len(gene_product['categories'])
        try:
            counts_by_n_categories[n_categories].add(gene_product['feature_id'])
        except KeyError:
            counts_by_n_categories[n_categories] = {gene_product['feature_id']}
    sorted_counts = {'{} categories'.format(k): v for k, v in sorted(counts_by_n_categories.items())}
    report_group_counts(gene_products, 'gene products by number of categories', sorted_counts)

    return


def report_category_overlaps(gene_products, category_members):
    """Log gene product and annotation counts for all pairwise category overlaps.

    Args:
        gene_products (dict): A dict of gene product dicts, keyed by feature_id.
        category_members (dict): A dict of sets of gene product feature_ids, keyed by category label.

    """
    log.info('##### PAIRWISE CATEGORY OVERLAPS #####')
    overlaps = {}
    for label_pair in combinations(category_members.keys(), 2):
        overlap_label = '{} + {}'.format(label_pair[0], label_pair[1])
        overlaps[overlap_label] = category_members[label_pair[0]].intersection(category_members[label_pair[1]])
    sorted_overlaps = dict(sorted(overlaps.items(), key=lambda i: len(i[1]), reverse=True))
    report_group_counts(gene_products, 'category overlap', sorted_overlaps)

    return


def report_gene_mapping(gene_products, category_members):
    """Log counts of gene products by the number of Dmel genes to which they map.

    Args:
        gene_products (dict): A dict of gene product dicts, keyed by feature_id.
        category_members (dict): A dict of sets of gene product feature_ids, keyed by category label.

    """
    log.info('##### GENE PRODUCT MAPPING TO CURRENT DMEL GENES #####')
    groups = {'no gene': set(), 'one gene': set(), 'many genes': set()}
    max_genes = 0
    for gene_product in gene_products.values():
        n_genes = len(gene_product['genes'])
        if n_genes == 0:
            group_label = 'no gene'
        elif n_genes == 1:
            group_label = 'one gene'
        else:
            group_label = 'many genes'
        groups[group_label].add(gene_product['feature_id'])
        max_genes = max(max_genes, n_genes)
    report_group_counts(gene_products, 'gene products by number of related genes', groups)
    log.info('The largest number of current Dmel genes related to a single gene product is {}.'.format(max_genes))
    log.info('{} gene products having {} annotations map to MANY current Dmel genes.'
             .format(len(groups['many genes']), count_annotations(gene_products, groups['many genes'])))
    aberrant_ids = category_members['ABERRANT']
    log.info('{} gene products having {} annotations represent aberrant expression of a current Dmel gene.'
             .format(len(aberrant_ids), count_annotations(gene_products, aberrant_ids)))

    return


def report_supplementary_counts():
    """Log the results of supplementary count queries."""
    log.info('##### SUPPLEMENTARY COUNTS #####')
    for count_query in COUNT_QUERIES:
        ret_count_info = connect(count_query['sql'], 'no_query', CONN)
        log.info('{} {}'.format(count_query['question'], ret_count_info))

    return


def report_invariants():
    """Check the curation invariants that the CURATED and ABERRANT categories depend upon."""
    log.info('##### CURATION INVARIANT CHECKS #####')
    for invariant_query in INVARIANT_QUERIES:
        ret_exceptions = connect(invariant_query['sql'], 'no_query', CONN)
        if ret_exceptions:
            log.error('INVARIANT BROKEN: {}. Found {} exception(s): {}'
                      .format(invariant_query['invariant'], len(ret_exceptions), ret_exceptions))
        else:
            log.info('INVARIANT HOLDS: {}.'.format(invariant_query['invariant']))

    return


def get_report_notes():
    """Build a list of notes describing the category labels for the report file header.

    Returns:
        A list of note strings.

    """
    notes = ['Category labels in the "categories" column are not mutually exclusive.']
    for category in CATEGORIES:
        notes.append('{}: {}'.format(category['label'], category['description']))
    notes.append('{}: The gene product fits none of the categories above.'.format(UNKNOWN_LABEL))

    return notes


def process_gene_products(gene_products):
    """Convert gene product dicts into a list of dicts for TSV output.

    Args:
        gene_products (dict): A dict of gene product dicts, keyed by feature_id.

    Returns:
        A list of dicts representing gene products and their categories.

    """
    log.info('Processing {} gene products for TSV export.'.format(len(gene_products)))
    data_list = []
    for gene_product in sorted(gene_products.values(), key=lambda i: i['uniquename']):
        data_list.append({
            'feature_id': gene_product['feature_id'],
            'uniquename': gene_product['uniquename'],
            'name': gene_product['name'],
            'categories': ','.join(gene_product['categories']),
            'n_categories': len(gene_product['categories']),
            'n_expression_annotations': gene_product['n_annotations'],
            'n_dmel_genes': len(gene_product['genes']),
            'dmel_genes': ','.join(sorted(gene_product['genes'])),
        })

    return data_list


if __name__ == "__main__":
    main()
