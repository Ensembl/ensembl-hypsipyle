"""Provide schema, dummy variant data, and API mocks for resolver tests."""

from io import StringIO
from unittest.mock import Mock

import pytest
import vcfpy

from common.file_model.base_variant import BaseVariant
from common.file_model.structural_variant import StructuralVariant
from common.file_model.variant import Variant
from graphql_service.ariadne_app import prepare_executable_schema


@pytest.fixture
def genome_id():
    """Provide a genome represented in the population test metadata."""
    return "a7335667-93e7-11ec-a39d-005056b38ce3"


@pytest.fixture
def schema():
    """Load the repository's GraphQL schema and actual resolver bindings."""
    return prepare_executable_schema()


@pytest.fixture
def variants(monkeypatch, genome_id):
    """Construct real models from small VCF records without external data."""
    monkeypatch.setattr(BaseVariant, "variant_sources", {})
    header = (
        "##fileformat=VCFv4.2\n"
        "##VEP=v110\n"
        '##source="test" description="Test database" url="https://example.org/"\n'
        '##INFO=<ID=SOURCE,Number=1,Type=String,Description="Source">\n'
        "##INFO=<ID=CSQ,Number=.,Type=String,"
        'Description="Format: Allele|Consequence">\n'
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    models = {}
    for model_class, alt in ((Variant, "T"), (StructuralVariant, "<DEL>")):
        with vcfpy.Reader.from_stream(
            StringIO(
                header + f"1\t100\ttest1\tA\t{alt}\t.\tPASS\tSOURCE=test;SVLEN=1\n"
            )
        ) as reader:
            record = next(reader)
            record.INFO["CSQ"] = []
            models[model_class.__name__] = model_class(record, reader.header, genome_id)
    return models


@pytest.fixture
def region_response():
    """Provide complete region metadata for schema-contract tests."""
    return {
        "name": "13",
        "length": 114364328,
        "code": "chromosome",
        "topology": "linear",
        "assembly": {
            "id": "assembly-id",
            "name": "GRCh38.p14",
            "accession_id": "GCA_000001405.29",
            "accessioning_body": "INSDC",
            "default": True,
            "regions": [],
            "organism": {
                "id": "organism-id",
                "scientific_name": "Homo sapiens",
                "assemblies": [],
                "species": {
                    "scientific_name": "Homo sapiens",
                    "taxon_id": 9606,
                    "alternative_names": [],
                    "organisms": [],
                },
            },
        },
        "sequence": {"checksum": "test-checksum"},
        "metadata": {"ontology_terms": []},
    }


@pytest.fixture
def region_post(monkeypatch, region_response):
    """Replace HTTP requests with a controllable region response."""
    response = Mock()
    response.json.return_value = {"data": {"region": region_response}}
    post = Mock(return_value=response)
    monkeypatch.setattr("graphql_service.resolver.region_model.requests.post", post)
    return post
