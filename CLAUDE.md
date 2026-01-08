# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

FastAPI web application for generating BED files from gene symbols. Users enter gene names (e.g., BRCA1, TP53), and the app fetches chromosomal positions from Ensembl REST API (GRCh38), then generates a downloadable BED file with optional region extension.

## Running the Application

```bash
# Install dependencies
pip install -r requirements.txt

# Start server
python3 app.py

# Development mode with auto-reload
uvicorn app:app --reload --host 127.0.0.1 --port 8000
```

Access at http://localhost:8000

## Architecture

Single-file FastAPI application (`app.py`) with embedded HTML/CSS/JS frontend. Uses Ensembl REST API for gene position lookup.

**Key components:**
- `fetch_gene_positions()`: Fetches gene coordinates from Ensembl API via POST /lookup/symbol/homo_sapiens
- `normalize_chromosome()`: Converts Ensembl format (17) to BED format (chr17)
- `is_standard_chromosome()`: Checks if chromosome is chr1-22, chrX, or chrY
- `process_bed_file()`: Merges, extends regions, and sorts by natural chromosome order

**API Endpoints:**
- `GET /`: Web interface
- `POST /fetch-positions`: Fetches gene positions from Ensembl (form field: `genes`)
- `POST /generate-bed`: Generates BED file (form fields: `genes_json`, `extend_bp`, `output_filename`)
- `GET /health`: Health check

## Application Flow

1. User enters gene symbols (comma or newline separated)
2. App calls Ensembl API to fetch positions
3. Genes on non-standard chromosomes (chrMT, patches, etc.) shown separately for user selection
4. Not found genes displayed as warning
5. User can extend regions by N base pairs
6. BED file generated with natural chromosome sorting

## Important Patterns

- Use `natsort.index_natsorted()` for chromosome sorting (chr1, chr2, ..., chr10)
- When extending regions, clip start to 0 to prevent negative coordinates
- Reverse strand positions are automatically normalized (start < end)
- Ensembl uses 1-based coordinates; convert to 0-based for BED format if needed
- Standard chromosomes: chr1-22, chrX, chrY (defined in `STANDARD_CHROMOSOMES` set)
