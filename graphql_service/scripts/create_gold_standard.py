#!/usr/bin/env python3
"""Generate gold standard query results for a list of test-case variants.

Description:
    This script reads a list of variant IDs (one per line) from a text file and
    a genome ID string from the command line. For each variant, it executes the
    selected GraphQL query against the in-memory schema and saves the query
    result as a JSON file in the specified output directory.

Example usage:
    python -m graphql_service.scripts.create_gold_standard \
        /app/graphql_service/tests/test_cases/predefined.txt \
        "a7335667-93e7-11ec-a39d-005056b38ce3" \
        /app/graphql_service/tests/gold_standard/a7335667-93e7-11ec-a39d-005056b38ce3

    python -m graphql_service.scripts.create_gold_standard \
        /app/graphql_service/tests/test_cases/structural.txt \
        "a7335667-93e7-11ec-a39d-005056b38ce3" \
        /app/graphql_service/tests/gold_standard_structural/a7335667-93e7-11ec-a39d-005056b38ce3 \
        --variant-type structural_variant
"""

import argparse
import asyncio
import json
import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

VARIANT_TYPES = ("variant", "structural_variant")


async def main(variant_file, genome_id, output_dir, variant_type):
    """
    Run selected query against schema for all test cases.
    """
    from graphql_service.tests.test_utils import (
        build_schema_context,
        execute_query,
        execute_structural_variant_query,
    )

    schema_and_context = build_schema_context()
    query_executors = {
        "variant": execute_query,
        "structural_variant": execute_structural_variant_query,
    }
    execute_variant_query = query_executors[variant_type]

    with open(variant_file, "r") as f:
        variants = [line.strip() for line in f if line.strip()]

    os.makedirs(output_dir, exist_ok=True)

    # Run selected query for each test case
    for variant_id in variants:
        query, success, result = await execute_variant_query(
            schema_and_context, genome_id, variant_id
        )
        if success:
            print(f"Query executed successfully for {variant_type} {variant_id}")
            output_file = os.path.join(output_dir, f"{variant_id}.json")
            with open(output_file, "w") as f_out:
                json.dump(result, f_out, indent=4)
            print(f"Result saved to {output_file}")
        else:
            print(
                f"Query execution failed for {variant_type} {variant_id}."
                f"\nQuery: {query}\nResult: {result}"
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate JSON files from GraphQL query results for variants"
    )
    parser.add_argument(
        "variant_file",
        help="Path to text file with variant IDs (one per line, format CHR:POS:ID)",
    )
    parser.add_argument("genome_id", help="Genome ID string to use in the query")
    parser.add_argument("output_dir", help="Directory to save JSON result files")
    parser.add_argument(
        "--variant-type",
        choices=VARIANT_TYPES,
        default="variant",
        help="GraphQL query type to execute",
    )
    args = parser.parse_args()

    asyncio.run(
        main(args.variant_file, args.genome_id, args.output_dir, args.variant_type)
    )
