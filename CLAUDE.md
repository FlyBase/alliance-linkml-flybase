# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This Python codebase controls the export of FlyBase data to Alliance LinkML-compliant JSON format. It transforms FlyBase chado database data into Alliance-compatible objects for submission to the Alliance of Genome Resources.

## Architecture

The codebase follows a handler-based architecture:

- **Data Type Classes**: `fb_datatypes.py` (FlyBase) and `agr_datatypes.py` (Alliance) define the data structures
- **Handler Classes**: Type-specific handlers in `*_handler.py` files (e.g., `gene_handler.py`, `allele_handlers.py`) that transform FlyBase data to Alliance format
- **Retrieval Scripts**: `AGR_data_retrieval_curation_*.py` scripts that orchestrate the data transformation process
- **Base Handler**: `handler.py` provides the core `DataHandler` class that all specific handlers inherit from

### Key Components

- **FBExportEntity/FBDataEntity**: Base classes for FlyBase data entities with Alliance mapping capabilities
- **AuditedObjectDTO/SubmittedObjectDTO**: Base classes for Alliance LinkML objects
- **DataHandler**: Core handler with database connection, logging, and common transformation methods
- **harvdev_utils**: External dependency for database utilities and character conversions

## Development Commands

### Environment Setup
```bash
# Install dependencies
pip install -r requirements.txt

# The project uses Docker for production deployment
docker build -t alliance-linkml-flybase .
```

### Running Scripts
```bash
# Run data retrieval scripts from the src/ directory
python AGR_data_retrieval_curation_gene.py
python AGR_data_retrieval_curation_allele.py
python AGR_data_retrieval_curation_agm.py
python AGR_data_retrieval_curation_disease.py
python AGR_data_retrieval_curation_construct.py
```

### Data Validation
```bash
# Validate generated JSON files using Alliance schema
python validate_agr_schema.py <json_file>
```

## Database Configuration

Scripts expect these environment variables:
- `SERVER`: Database server (e.g., `flysql25`)
- `DATABASE`: Database name (e.g., `production_chado`)
- `RELEASE`: FlyBase release (e.g., `2022_05`)
- `LINKML_VERSION`: Alliance LinkML version (e.g., `v1.3.1`)
- `ALLIANCETOKEN`: API token for Alliance uploads

## File Structure

- `src/`: Main source code directory
  - `AGR_data_retrieval_*.py`: Data retrieval scripts
  - `*_handler.py`: Type-specific data handlers
  - `fb_datatypes.py`: FlyBase data type definitions
  - `agr_datatypes.py`: Alliance data type definitions
  - `handler.py`: Base handler class
  - `utils.py`: Utility functions
  - `retired/`: Deprecated code

## Production Workflow

The codebase is designed to run in GoCD pipelines using Docker containers. The typical workflow:

1. Scripts query FlyBase chado database
2. Data is transformed using handlers
3. JSON files are generated and validated
4. Files are uploaded to Alliance persistent store

## Development Notes

- Work on release branches corresponding to Alliance LinkML schema versions (e.g., `v1.3.1`)
- Validate files locally before uploading to Alliance
- The `testing` parameter in handlers enables test mode
- Use the Alliance Curation API to verify uploads
- Dependencies exist between data types (e.g., alleles must be loaded before disease annotations)