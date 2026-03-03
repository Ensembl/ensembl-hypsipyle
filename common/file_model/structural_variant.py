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

from typing import Any, Mapping, List, Union
import re
import os
import json

from common.file_model.base_variant import BaseVariant


class StructuralVariant(BaseVariant):
    """StructuralVariant model – inherits shared behaviour from BaseVariant."""

    def __init__(self, record: Any, header: Any, genome_uuid: str) -> None:
        """Initialise SV-specific attributes and delegate shared setup.

        The `length` attribute is derived from the VCF INFO `SVLEN` where
        available; otherwise zero.
        """
        super().__init__(record, header, genome_uuid)
        self.type = "StructuralVariant"
        # SVLEN may be a list or integer depending on the generator; coerce to int
        svlen = self.info.get("SVLEN") if isinstance(self.info, dict) else None
        try:
            if isinstance(svlen, (list, tuple)) and svlen:
                self.length = int(svlen[0])
            elif svlen is not None:
                self.length = int(svlen)
            else:
                self.length = 0
        except Exception:
            self.length = 0
