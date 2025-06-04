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
from .test_utils import build_schema_context
import os
import json

@pytest.fixture(scope="session") 
def schema_context():
    """
    Session-scoped fixture reused by all tests.
    """
    return build_schema_context()

@pytest.fixture(scope="session")
def gold_standard_loader():
    """
    Returns a function that, given a variant_id, will read
    and parse the corresponding JSON gold-standard file.
    """
    def _loader(path: str, genome_id: str, variant_id: str) -> dict:
        path = os.path.join(path, genome_id, f"{variant_id}.json")
        with open(path, "r") as f:
            return json.load(f)
    return _loader