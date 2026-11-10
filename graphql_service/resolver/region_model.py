"""Resolve slice regions through the Ensembl core GraphQL API."""

import asyncio

import requests
from ariadne import ObjectType
from graphql import (
    GraphQLError,
    ast_from_value,
    get_named_type,
    is_abstract_type,
    print_ast,
)
from graphql.execution.collect_fields import collect_sub_fields
from graphql.execution.values import get_argument_values

SLICE_TYPE = ObjectType("Slice")
REGION_API_URL = "https://www.ensembl.org/data/graphql/api"


def _region_selection(info, parent_type, field_nodes):
    """Build nested upstream selections from GraphQL's collected fields.

    GraphQL handles fragment expansion and skip/include directives.
    Merge aliases by schema field name so normal local field resolution
    can read the returned dictionaries. Region metadata fields currently
    have no arguments; retain argument values when building selections.
    """
    if is_abstract_type(parent_type):
        return " ".join(
            f"... on {concrete.name} {{ "
            f"{_region_selection(info, concrete, field_nodes)} }}"
            for concrete in info.schema.get_possible_types(parent_type)
        )
    fields = collect_sub_fields(
        info.schema, info.fragments, info.variable_values, parent_type, field_nodes
    )
    by_name = {}
    for nodes in fields.values():
        by_name.setdefault(nodes[0].name.value, []).extend(nodes)
    selections = []
    for name, nodes in by_name.items():
        if name == "__typename":
            selections.append(name)
            continue
        field = parent_type.fields[name]
        selection = (
            "id: assembly_id"
            if parent_type.name == "Assembly" and name == "id"
            else name
        )
        arguments = get_argument_values(field, nodes[0], info.variable_values)
        if arguments:
            selection += (
                "("
                + ", ".join(
                    f"{key}: {print_ast(ast_from_value(value, field.args[key].type))}"
                    for key, value in arguments.items()
                )
                + ")"
            )
        if nodes[0].selection_set:
            children = _region_selection(info, get_named_type(field.type), nodes)
            selection += " { " + children + " }"
        selections.append(selection)
    return " ".join(selections) or "__typename"


def _fetch_region(query, variables):
    """Fetch a region or raise a GraphQL error for an upstream failure."""
    try:
        response = requests.post(
            REGION_API_URL,
            json={"query": query, "variables": variables},
            timeout=10,
        )
        response.raise_for_status()
        payload = response.json()
    except (requests.RequestException, ValueError) as error:
        raise GraphQLError("Unable to fetch region metadata from Ensembl.") from error
    if not isinstance(payload, dict) or payload.get("errors"):
        raise GraphQLError("Ensembl returned an error while resolving region metadata.")
    data = payload.get("data")
    region = data.get("region") if isinstance(data, dict) else None
    if not isinstance(region, dict):
        raise GraphQLError("Region metadata was not found in Ensembl.")
    return region


@SLICE_TYPE.field("region")
async def resolve_region(slice_value, info):
    """Fetch the requested region fields once per identical request query.

    The private genome key identifies the assembly without changing the
    public schema. Keep network I/O off the event loop and cache only
    within the current GraphQL request, including concurrent lookups.
    """
    genome_id = slice_value["_genome_id"]
    name = slice_value["region"]["name"]
    selection = _region_selection(
        info, info.schema.get_type("Region"), info.field_nodes
    )
    query = (
        "query RegionByName($genome_id: String!, $name: String!) { "
        "region(by_name: {genome_id: $genome_id, name: $name}) { " + selection + " } }"
    )
    variables = {"genome_id": genome_id, "name": name}
    cache = info.context.setdefault("region_queries", {})
    key = (genome_id, name, selection)
    if key not in cache:
        cache[key] = asyncio.create_task(
            asyncio.to_thread(_fetch_region, query, variables)
        )
    return await cache[key]
