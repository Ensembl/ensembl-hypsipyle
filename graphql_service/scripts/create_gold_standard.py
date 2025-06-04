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
from graphql_service.tests.test_utils import execute_query, build_schema_context
import os

schema_and_context = build_schema_context()

async def main(variant_file, genome_id, output_dir):
    """
    Run master query against schema for all test cases.
    """
    
    with open(variant_file, "r") as f:
        variants = [line.strip() for line in f if line.strip()]

    os.makedirs(output_dir, exist_ok=True)

    # Run master query for each test case
    for variant_id in variants:
        query, success, result = await execute_query(schema_and_context, genome_id, variant_id)
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