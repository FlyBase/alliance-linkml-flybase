"""Shared TSV writers for the AGR_data_retrieval_curation_* export scripts.

These helpers replace per-script copies of `generate_tsv_file` and
`generate_association_tsv_file` that previously lived in
AGR_data_retrieval_curation_cassette.py and
AGR_data_retrieval_curation_construct.py (and which appear in similar form
in the allele and transgenic_tool scripts).

DTO field names are derived from a `datatype` string (e.g. "cassette",
"construct"). All current scripts follow the regular convention
`{datatype}_full_name_dto`, `{datatype}_symbol_dto`,
`{datatype}_synonym_dtos`, etc.; if a future script breaks that
convention these helpers will need to be parameterized further.
"""

import os


NO_PUBS_SENTINEL = "NO PUBS"
EVIDENCE_DELIMITER = "|"

PRIMARY_TSV_HEADER = (
    "# Primary FBid\tValid symbol\tValid full name\tsecondary FBid(s)\tsynonyms\tinternal\n"
)
NOTES_TSV_HEADER = "# Primary FBid\ttype\tcomment\n"
COMPONENTS_TSV_HEADER = "# Primary FBid\tsymbol\trelation\ttaxon\tevidence\n"
# NB: existing tool_uses TSVs write rows as primary, tools, evidence; the
# header preserves that historic column order verbatim.
TOOL_USES_TSV_HEADER = "# Primary FBid\tevidence\ttool_uses\n"


def should_skip_obsolete():
    """Honor ADD_OBSOLETE=NO. Centralized so behavior is identical across scripts."""
    return os.environ.get('ADD_OBSOLETE') == 'NO'


def _is_excluded(entity_dict):
    return entity_dict.get('internal') or entity_dict.get('obsolete')


def write_primary_tsv(*, log, filename, entities, datatype):
    """Write the primary identifier TSV (used by every curation script)."""
    skip = should_skip_obsolete()
    if skip:
        log.info(f'ADD_OBSOLETE=NO: excluding obsolete/internal {datatype}s from TSV.')
    full_name_key = f'{datatype}_full_name_dto'
    symbol_key = f'{datatype}_symbol_dto'
    synonym_key = f'{datatype}_synonym_dtos'
    with open(filename, 'w') as outfile:
        outfile.write(PRIMARY_TSV_HEADER)
        for entity_dict in entities:
            if skip and _is_excluded(entity_dict):
                continue
            primary = entity_dict["primary_external_id"]
            symbol = ''
            name = ''
            secondary = []
            syns = []
            if full_name_key in entity_dict:
                name = entity_dict[full_name_key]["format_text"]
            if symbol_key in entity_dict:
                symbol = entity_dict[symbol_key]["format_text"]
            if synonym_key in entity_dict:
                for synonym in entity_dict[synonym_key]:
                    syns.append(synonym["format_text"])
            if "secondary_identifiers" in entity_dict:
                secondary = entity_dict["secondary_identifiers"]
            internal = entity_dict.get("internal", False)
            secondary_str = EVIDENCE_DELIMITER.join(secondary)
            syns_str = EVIDENCE_DELIMITER.join(syns)
            try:
                outfile.write(
                    f"{primary}\t{symbol}\t{name}\t{secondary_str}\t{syns_str}\t{internal}\n"
                )
            except TypeError:
                log.error(f"entity_dict: {entity_dict}")
                log.error(f"primary: {primary}")
                log.error(f"secondary {secondary}")
                log.error(f"symbol: {symbol}")
                log.error(f"name: {name}")
                log.error(f"syns: {syns}")
                log.error(f"internal: {internal}")
                raise


def write_notes_tsv(*, filename, entities):
    """Write the per-entity notes TSV (`note_dtos` column)."""
    skip = should_skip_obsolete()
    with open(filename, 'w') as outfile:
        outfile.write(NOTES_TSV_HEADER)
        for entity_dict in entities:
            if skip and _is_excluded(entity_dict):
                continue
            primary = entity_dict["primary_external_id"]
            for note in entity_dict.get("note_dtos", []):
                outfile.write(f"{primary}\t{note['note_type_name']}\t{note['free_text']}\n")


def write_components_tsv(*, filename, entities, datatype):
    """Write the component-slot TSV (`{datatype}_component_dtos`)."""
    skip = should_skip_obsolete()
    component_key = f'{datatype}_component_dtos'
    with open(filename, 'w') as outfile:
        outfile.write(COMPONENTS_TSV_HEADER)
        for entity_dict in entities:
            if skip and _is_excluded(entity_dict):
                continue
            primary = entity_dict["primary_external_id"]
            for comp in entity_dict.get(component_key, []):
                if 'evidence_curies' in comp:
                    evidence = EVIDENCE_DELIMITER.join(comp['evidence_curies'])
                else:
                    evidence = ""
                outfile.write(
                    f"{primary}\t{comp['component_symbol']}\t{comp['relation_name']}\t"
                    f"{comp['taxon_curie']}\t{evidence}\n"
                )


def write_tool_uses_tsv(*, filename, entities, datatype, no_pubs_sentinel=NO_PUBS_SENTINEL):
    """Write the tool-uses TSV (`{datatype}_use_dtos`)."""
    skip = should_skip_obsolete()
    use_key = f'{datatype}_use_dtos'
    with open(filename, 'w') as outfile:
        outfile.write(TOOL_USES_TSV_HEADER)
        for entity_dict in entities:
            if skip and _is_excluded(entity_dict):
                continue
            primary = entity_dict["primary_external_id"]
            for comp in entity_dict.get(use_key, []):
                if 'evidence_curies' in comp:
                    evidence = EVIDENCE_DELIMITER.join(comp['evidence_curies'])
                else:
                    evidence = no_pubs_sentinel
                tools = EVIDENCE_DELIMITER.join(comp["use_curies"])
                outfile.write(f"{primary}\t{tools}\t{evidence}\n")


def write_association_tsv(
    *,
    filename,
    rows,
    first_field,
    second_field,
    extra_fields=None,
    no_pubs_sentinel=NO_PUBS_SENTINEL,
):
    """Write an association TSV (subject, relation, object, evidence [, extras]).

    `extra_fields` is a list of `(header_label, dict_key, default_when_missing)`
    tuples; list/tuple values are joined with EVIDENCE_DELIMITER, missing keys
    fall back to `default_when_missing`.

    `no_pubs_sentinel` controls the cell written when `evidence_curies` is
    absent: cassette uses the default ("NO PUBS"); construct passes "" to
    preserve its historical empty-cell behavior.
    """
    skip = should_skip_obsolete()
    extras = extra_fields or []
    extra_headers = "".join(f"\t{label}" for label, _, _ in extras)
    with open(filename, 'w') as outfile:
        outfile.write(
            f"#{first_field}\tRelationship\t{second_field}\tEvidence{extra_headers}\n"
        )
        for entity_dict in rows:
            if skip and _is_excluded(entity_dict):
                continue
            sub = entity_dict[first_field]
            obj = entity_dict[second_field]
            rel_type = entity_dict['relation_name']
            if 'evidence_curies' in entity_dict:
                pubs = EVIDENCE_DELIMITER.join(entity_dict['evidence_curies'])
            else:
                pubs = no_pubs_sentinel
            extra_parts = []
            for _label, key, default in extras:
                if key in entity_dict and isinstance(entity_dict[key], (list, tuple)):
                    extra_parts.append(EVIDENCE_DELIMITER.join(entity_dict[key]))
                elif key in entity_dict:
                    extra_parts.append(str(entity_dict[key]))
                else:
                    extra_parts.append(default)
            extras_str = "".join(f"\t{p}" for p in extra_parts)
            outfile.write(f"{sub}\t{rel_type}\t{obj}\t{pubs}{extras_str}\n")
