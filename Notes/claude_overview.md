What This Repo Does

  This is an ETL (Extract-Transform-Load) pipeline that exports FlyBase chado database data into Alliance of Genome Resources
  LinkML-compliant JSON. FlyBase is a Drosophila Model Organism Database (MOD), and the Alliance is a consortium that aggregates data
  across MODs.

  Architecture

  Class Hierarchy

  DataHandler (handler.py) - base: DB connections, lookups, common logic
    └── PrimaryEntityHandler (entity_handler.py) - entities with FB curies
        └── FeatureHandler (feature_handler.py) - chado Feature table entities
            ├── GeneHandler (gene_handler.py)
            ├── ConstructHandler (construct_handler.py)
            ├── CassetteHandler (cassette_handler.py)
            └── ExperimentalToolHandler (transgenic_tool_handler.py)

  Separate handler files also exist for alleles, AGMs (strains/genotypes), disease annotations, and expression.

  Data Type Classes

  - fb_datatypes.py - FlyBase entity classes (FBGene, FBAllele, FBConstruct, etc.) wrapping chado objects
  - agr_datatypes.py - Alliance LinkML DTO classes (GeneDTO, AlleleDTO, ConstructDTO, etc.) matching the Alliance schema

  The 6 Stages of Data Flow

  Stage 1: Entry Point & Setup

  Each AGR_data_retrieval_curation_*.py script:
  - Parses CLI args (-l linkml version, -r reference DB for incremental updates)
  - Reads env vars (SERVER, DATABASE, RELEASE, LINKML_VERSION)
  - Creates SQLAlchemy sessions to PostgreSQL chado DB

  Stage 2: Build Lookups (get_general_data)

  Handler builds reference data structures:
  - Bibliography - pub_id → PMID/FBrf curies
  - CVterm lookup - cvterm_id → name, cv_name, curie
  - Organism lookup - organism_id → genus, species, taxon
  - Feature lookup - feature_id → uniquename, type, organism

  Stage 3: Query Chado (get_datatype_data)

  Each handler queries chado tables for its specific data type:
  - Primary entities (genes, alleles, constructs, etc.)
  - Relationships between entities
  - Publications, synonyms, cross-references
  - Properties and CV term annotations
  - Audit timestamps

  Stage 4: Synthesize (synthesize_info)

  Raw chado data is processed:
  - Timestamps consolidated
  - Synonyms organized by type
  - Secondary IDs collected
  - NCBITaxon IDs resolved
  - Handler-specific synthesis (e.g., construct components, gene types)

  Stage 5: Map to Alliance (map_fb_data_to_alliance)

  FlyBase objects are converted to Alliance LinkML DTOs:
  - Create DTO instances (GeneDTO, AlleleDTO, etc.)
  - Map identifiers, names, synonyms
  - Map cross-references and secondary IDs
  - Map timestamps and data provider info
  - Handler-specific mappings (gene type SO terms, allele mutation types, construct components)

  Stage 6: Export (generate_export_dict + generate_export_file)

  - Filter entities (exclude internal, non-exportable)
  - Convert DTOs to Python dicts
  - Write JSON with sorted keys, 2-space indent
  - Optionally generate TSV files for curator review (cassettes, tools)

