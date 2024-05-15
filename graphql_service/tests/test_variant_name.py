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

"""Testing of variant name"""

from string import Template
from ariadne import graphql
import pytest
from .test_utils import setup_test

executable_schema, context = setup_test()

def get_generic_query_template():
    query = """{
            variant(
            by_id: {genome_id: "", variant_id: "$variant_id"}
        ) {
            
            name
            
        }
        }"""
    template = Template(query)
    return template

async def test_variant_id():
    template = get_generic_query_template()
    query = template.substitute(variant_id="1:10007:rs1639538116")
    query_data = {"query": query}
    (success, result) = await graphql(
        executable_schema, query_data, context_value=context
    )
    assert success
