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


def minimise_allele(alt: str, ref: str) -> str:
    """Converts a VCF allele string into a minimised SPDI format.

    This function converts a VCF allele representation that omits anchoring bases for prediction
    scores in the INFO column. It removes the anchoring base if the first base of the reference
    and alternate alleles is identical. If the resulting allele is empty, a hyphen ('-') is returned.

    Args:
        alt (str): The alternate allele from the VCF data.
        ref (str): The reference allele from the VCF data.

    Returns:
        str: The minimised allele in SPDI format.
    """
    minimised_allele_string = alt
    if ref[0] == alt[0]:
        minimised_allele_string = alt[1:] if len(alt) > 1 else "-"
    return minimised_allele_string
