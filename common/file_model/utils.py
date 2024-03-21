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


def minimise_allele(alt: str, ref:str) -> str:
        """
        VCF file has the representation without anchoring bases
        for prediction scores in INFO column. This function is useful
        in matching the SPDI format in VCF with the allele in memory
        """
        minimised_allele_string = alt
        if ref[0] == alt[0]:
            minimised_allele_string = alt[1:] if len(alt) > 1 else "-"
        return minimised_allele_string