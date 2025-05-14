import pytest
from common.file_client import FileClient
from graphql_service.ariadne_app import (
    prepare_executable_schema,
    prepare_context_provider,
)

@pytest.fixture(scope="module")
def schema_and_context():
    """
    Prepare the GraphQL executable schema and context provider for tests.
    """
    config = { "data_root": "/app/data" }
    schema = prepare_executable_schema()
    file_client = FileClient(config)
    context = prepare_context_provider({ "file_client": file_client })
    return schema, context