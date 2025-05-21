"""
.. See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at
       http://www.apache.org/licenses/LICENSE-2.0
   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

from string import Template
from ariadne import graphql
import os

def get_test_case_ids(genome_id: str, gold_std_dir: str):
    """
    Dynamically generate a list of test cases from 
    gold standard query result JSON files in gold_std_dir.
    
    Expected file naming: <variant_id>.json
    Uses the same genome_id for all test cases.
    """
    cases = []
    for filename in os.listdir(gold_std_dir):
        if filename.endswith(".json"):
            variant_id = filename.replace(".json", "")
            cases.append((variant_id, genome_id))
    return cases

def master_query_template():
    """
    Template for the master query.
    """
    
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

async def execute_query(schema_and_context, genome_id: str, variant_id: str):
    """
    Execute the query for a given variant_id and genome_id according to template
    query, returning tne query string, success status, and result.
    """
    executable_schema, context = schema_and_context
    template = master_query_template()
    query = template.substitute(genome_id=genome_id, variant_id=variant_id)
    success, result = await graphql(
        executable_schema,
        {"query": query},
        context_value=context(request={})
    )
    return query, success, result