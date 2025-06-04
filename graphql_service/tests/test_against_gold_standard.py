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

import pytest
import os
import json
import difflib
from .test_utils import get_test_case_ids, execute_query

GOLD_STD_QUERY_RESULTS_PATH = "/app/graphql_service/tests/gold_standard"
GENOME_ID = "a7335667-93e7-11ec-a39d-005056b38ce3"
TEST_CASES = get_test_case_ids(GOLD_STD_QUERY_RESULTS_PATH, GENOME_ID)

# Run query for each test case, generating a new query result and comparing it against the stored gold standard
@pytest.mark.asyncio
@pytest.mark.parametrize("variant_id, genome_id", TEST_CASES)
async def test_query_results_against_gold_std(
    schema_context: tuple, 
    gold_standard_loader,
    variant_id: str, 
    genome_id: str
    ):
    """Test present query result for variant_id X, genome_id Y against 
    stored gold standard query result for variant_id X, genome_id Y."""

    # Run live query
    query, success, result = await execute_query(schema_context, genome_id, variant_id)
    
    # Load gold standard result
    gold_standard = gold_standard_loader(GOLD_STD_QUERY_RESULTS_PATH, genome_id, variant_id)
    
    # Compare
    result_str = json.dumps(result.get("data", {}), indent=4, sort_keys=True)
    gold_str = json.dumps(gold_standard.get("data", {}), indent=4, sort_keys=True)
    diff = ''.join(
        difflib.unified_diff(
            gold_str.splitlines(keepends=True),
            result_str.splitlines(keepends=True),
            fromfile="Gold Standard",
            tofile="Result"
        )
    )
    
    # Log the diff if there are any differences
    if result_str != gold_str:
        print(f"Diff for variant {variant_id}:\n{diff}")
    
    assert result_str == gold_str