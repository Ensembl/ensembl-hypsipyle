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
from typing import Dict

from graphql import GraphQLError


class FieldNotFoundError(GraphQLError):
    """
    Custom error to be raised if a field cannot be found by id
    """

    def __init__(self, field_type: str, key_dict: Dict[str, str]):
        self.extensions = {"code": f"{field_type.upper()}_NOT_FOUND"}
        ids_string = ", ".join([f"{key}={val}" for key, val in key_dict.items()])
        message = f"Failed to find {field_type} with ids: {ids_string}"
        self.extensions.update(key_dict)
        super().__init__(message, extensions=self.extensions)



class VariantNotFoundError(FieldNotFoundError):
    """
    Custom error to be raised if variant is not found
    """
    def __init__(
        self, variant_id: str
    ):
        super().__init__("variant_id", {"variant_id": variant_id})

