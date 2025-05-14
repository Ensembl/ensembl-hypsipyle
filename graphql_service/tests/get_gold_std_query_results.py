#!/usr/bin/env python3
"""
Script to generate gold standard query results for a list of test-case variants. 

Description:
    This script reads a list of variant IDs (one per line) from a text file and a genome ID string from the command line.
    For each variant, it executes the master GraphQL query against the in-memory schema and saves the query result 
    as a JSON file in the specified output directory.

Example usage:
    python -m graphql_service.tests.get_gold_std_query_results \
        /app/graphql_service/tests/gold_std_test_cases.txt \
        "a7335667-93e7-11ec-a39d-005056b38ce3" \
        /app/graphql_service/tests/gold_std_query_results
"""

import argparse
import asyncio
import json
from string import Template
from ariadne import graphql
from .test_utils import setup_test
import os

executable_schema, context = setup_test()

# -----------------------------------------------------------------------
# Define helper functions

# Define master query template
def master_query_template():
    """Query template"""
    
    query = """{
        variant(
            by_id: { genome_id: "$genome_id", variant_id: "$variant_id" }
        ) {
            name
            alternative_names { __typename }
            primary_source { __typename }
            type
            allele_type { __typename }
            slice  { __typename }
            alleles { 
            name
            allele_sequence
            reference_sequence
            alternative_names { __typename }
            type
            allele_type { __typename }
            slice { __typename }
            phenotype_assertions {
                feature
                feature_type { __typename }
                phenotype {
                name
                source { __typename }
                ontology_terms { __typename }
                }
                evidence {
                source { __typename }
                assertion { __typename }
                }
            }
            prediction_results {
                score
                result
                classification { __typename }
                analysis_method { __typename }
            }
            population_frequencies {
                population_name
                allele_frequency
                allele_count
                allele_number
                is_minor_allele
                is_hpmaf
            }
            predicted_molecular_consequences {
                allele_name
                stable_id
                feature_type { __typename }
                consequences { __typename }
                prediction_results {
                score
                result
                classification { __typename }
                analysis_method { __typename }
                }
                gene_stable_id
                gene_symbol
                protein_stable_id
                transcript_biotype
                cdna_location {
                relation { __typename }
                start
                end
                length
                percentage_overlap
                ref_sequence
                alt_sequence
                }
                cds_location {
                relation { __typename }
                start
                end
                length
                percentage_overlap
                ref_sequence
                alt_sequence
                }
                protein_location {
                relation { __typename }
                start
                end
                length
                percentage_overlap
                ref_sequence
                alt_sequence
                }
            }
            ensembl_website_display_data {
                count_transcript_consequences
                count_overlapped_genes
                count_regulatory_consequences
                count_variant_phenotypes
                count_gene_phenotypes
                representative_population_allele_frequency
            }
            }
            prediction_results {
            score
            result
            classification { __typename }
            analysis_method { __typename }
            }
            ensembl_website_display_data {
            count_citations
            }
        }
    }"""
    return Template(query)

async def execute_query(genome_id, variant_id):
    """
    Execute the query for a given variant_id and genome_id according to template
    query, returning tne query string, success status, and result.
    """
    template = master_query_template()
    query = template.substitute(genome_id=genome_id, variant_id=variant_id)
    query_data = {"query": query}
    success, result = await graphql(
        executable_schema, query_data, context_value=context(request={})
    )
    return query, success, result

# -----------------------------------------------------------------------
# Generate gold standard query results

async def main(variant_file, genome_id, output_dir):
    """
    Run master query against schema for all test cases.
    """
    
    with open(variant_file, "r") as f:
        variants = [line.strip() for line in f if line.strip()]

    os.makedirs(output_dir, exist_ok=True)

    # Run master query for each test case
    for variant_id in variants:
        query, success, result = await execute_query(genome_id, variant_id)
        if success:
            print(f"Query executed successfully for variant {variant_id}")
            output_file = f"{output_dir}/{variant_id}.json"
            with open(output_file, "w") as f_out:
                json.dump(result, f_out, indent=4)
            print(f"Result saved to {output_file}")
        else:
            print(f"Query execution failed for variant {variant_id}.\nQuery: {query}\nResult: {result}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate JSON files from GraphQL query results for variants")
    parser.add_argument("variant_file", help="Path to text file with variant IDs (one per line, format CHR:POS:ID)")
    parser.add_argument("genome_id", help="Genome ID string to use in the query")
    parser.add_argument("output_dir", help="Directory to save JSON result files")
    args = parser.parse_args()
    
    asyncio.run(main(args.variant_file, args.genome_id, args.output_dir))