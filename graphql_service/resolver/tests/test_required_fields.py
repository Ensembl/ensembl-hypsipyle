"""Check Python return values against GraphQL schema nullability."""

import inspect
from types import SimpleNamespace

import pytest
from graphql import (
    build_schema,
    default_field_resolver,
    get_named_type,
    is_leaf_type,
    is_list_type,
    is_non_null_type,
)
from graphql.language import FieldNode, NameNode, SelectionSetNode

from graphql_service.resolver.population_model import resolve_populations
from graphql_service.resolver.variant_model import resolve_api


def _selection_set(field_type, ancestors=()):
    """Build field selections for resolvers that delegate nested queries."""
    named_type = get_named_type(field_type)
    if is_leaf_type(named_type):
        return None
    if named_type.name in ancestors:
        return SelectionSetNode(
            selections=(FieldNode(name=NameNode(value="__typename")),)
        )
    return SelectionSetNode(
        selections=tuple(
            FieldNode(
                name=NameNode(value=name),
                selection_set=_selection_set(field.type, (*ancestors, named_type.name)),
            )
            for name, field in named_type.fields.items()
        )
    )


async def _nullability_errors(schema, field_type, value, path, ancestors=()):
    """Collect missing required values using the schema's nullability rules.

    Use registered resolvers or GraphQL's default key/attribute lookup.
    Inspect nullable objects when present, and distinguish a required list
    from required list items. Stop only actual object cycles. Empty lists
    and absent nullable parents do not exercise their child fields.
    """
    if is_non_null_type(field_type):
        if value is None:
            return [f"{path}: required {field_type}, got None or a missing field"]
        field_type = field_type.of_type
    elif value is None:
        return []

    if is_list_type(field_type):
        errors = []
        for index, item in enumerate(value):
            errors.extend(
                await _nullability_errors(
                    schema, field_type.of_type, item, f"{path}[{index}]", ancestors
                )
            )
        return errors
    if is_leaf_type(field_type):
        return []

    identity = (field_type.name, id(value))
    if identity in ancestors:
        return []
    ancestors = (*ancestors, identity)
    errors = []
    # Interfaces declare the same required fields as their implementations.
    for name, field in field_type.fields.items():
        field_path = f"{path}.{name}"
        info = SimpleNamespace(
            field_name=name,
            parent_type=field_type,
            return_type=field.type,
            schema=schema,
            context={},
            fragments={},
            variable_values={},
            field_nodes=[
                FieldNode(
                    name=NameNode(value=name), selection_set=_selection_set(field.type)
                )
            ],
        )
        resolver = field.resolve or default_field_resolver
        try:
            result = resolver(value, info)
            if inspect.isawaitable(result):
                result = await result
        except Exception as error:
            errors.append(f"{field_path}: {type(error).__name__}: {error}")
            continue
        errors.extend(
            await _nullability_errors(schema, field.type, result, field_path, ancestors)
        )
    return errors


@pytest.mark.parametrize(
    "type_name",
    ["Variant", "VariantAllele", "StructuralVariant", "StructuralVariantAllele"],
)
async def test_variant_required_fields(schema, variants, type_name, region_post):
    """Require every schema-marked value in model and resolver output."""
    if type_name.endswith("Allele"):
        values = variants[type_name.removesuffix("Allele")].get_alleles()
        assert values, "The fixture must exercise allele fields"
    else:
        values = [variants[type_name]]
    errors = []
    for index, value in enumerate(values):
        errors.extend(
            await _nullability_errors(
                schema, schema.get_type(type_name), value, f"{type_name}[{index}]"
            )
        )
    assert not errors, "\n".join(errors)


async def test_version_required_fields(schema):
    """Check the version dictionary using the schema's required fields."""
    value = resolve_api(None, None)
    assert value is not None
    errors = await _nullability_errors(
        schema, schema.get_type("Version"), value, "Version"
    )
    assert not errors, "\n".join(errors)


@pytest.mark.xfail(
    reason="Advisory: population metadata does not yet satisfy schema nullability",
    raises=AssertionError,
    strict=False,
)
async def test_population_required_fields(schema, genome_id):
    """Report population nullability gaps without failing the test run."""
    populations = resolve_populations(None, None, genome_id=genome_id)
    assert populations, "The fixture must exercise population fields"
    errors = await _nullability_errors(
        schema,
        schema.query_type.fields["populations"].type,
        populations,
        "populations",
    )
    assert not errors, "\n".join(errors)


@pytest.mark.parametrize(
    "payload,expected_paths",
    [
        ({"number": 0, "flag": False, "text": "", "items": []}, []),
        (SimpleNamespace(number=0, flag=False, text="", items=[]), []),
        ({"number": None, "flag": False, "text": "", "items": []}, ["Item.number"]),
        ({"flag": False, "text": "", "items": []}, ["Item.number"]),
        ({"number": 0, "flag": False, "text": "", "items": None}, ["Item.items"]),
        ({"number": 0, "flag": False, "text": "", "items": [None]}, ["Item.items[0]"]),
        (
            {"number": 0, "flag": False, "text": "", "items": [], "child": {}},
            ["Item.child.name"],
        ),
    ],
)
async def test_nullability_rules(payload, expected_paths):
    """Accept valid falsy values and detect missing nested and list values."""
    schema = build_schema(
        "type Query { item: Item } "
        "type Item { number: Int! flag: Boolean! text: String! "
        "items: [String!]! optional: String child: Child } "
        "type Child { name: String! }"
    )
    errors = await _nullability_errors(schema, schema.get_type("Item"), payload, "Item")
    assert [error.split(":", 1)[0] for error in errors] == expected_paths
