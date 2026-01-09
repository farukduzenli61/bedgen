# BED File Generator

A FastAPI web application for generating BED files from gene symbols. Enter gene names and automatically fetch chromosomal positions from Ensembl API (GRCh38).

## Features

- **Gene Symbol Input** - Enter genes as comma-separated or line-by-line (e.g., BRCA1, TP53, EGFR)
- **Ensembl API Integration** - Automatically fetch positions from Ensembl REST API (GRCh38)
- **Non-standard Chromosome Selection** - Choose to include/exclude genes on chrMT, patches, etc.
- **Extension Preview** - See before/after comparison for 3-5 random genes before downloading
- **Natural Chromosome Sorting** - Proper ordering (chr1, chr2, ..., chr10, chr11)
- **Reverse Strand Handling** - Automatically normalize positions (start < end)

## Quick Start

### Local Development

```bash
# Clone repository
git clone https://github.com/farukduzenli61/bedgen.git
cd bedgen

# Install dependencies
pip install -r requirements.txt

# Run server
python3 app.py
```

Open http://localhost:8000

### Deploy to Railway

1. Fork this repository
2. Go to [railway.app](https://railway.app)
3. "New Project" → "Deploy from GitHub repo"
4. Select your forked repo
5. Done! Get your URL from Settings → Generate Domain

## Usage

### Step 1: Enter Gene Symbols
```
BRCA1, TP53, EGFR, ACTB
```
or
```
BRCA1
TP53
EGFR
ACTB
```

### Step 2: Review Results
- Found genes displayed in a table with positions
- Warnings shown for genes not found
- Non-standard chromosomes listed separately for selection

### Step 3: Preview Extension
- Set extension amount (bp)
- Click "Preview" to see comparison table
- Random 3-5 genes shown with original vs extended positions

### Step 4: Download BED File
- Set output filename
- Click download

## API Endpoints

### `GET /`
Web interface

### `POST /fetch-positions`
Fetch gene positions from Ensembl API

**Form data:**
- `genes`: Gene symbols (comma or newline separated)

**Response:**
```json
{
  "genes": [
    {"gene": "BRCA1", "chr": "chr17", "start": 43044292, "end": 43170245, "strand": -1, "is_standard": true}
  ],
  "non_standard": [],
  "not_found": ["INVALIDGENE"]
}
```

### `POST /generate-bed`
Generate and download BED file

**Form data:**
- `genes_json`: JSON array of gene objects
- `extend_bp`: Base pairs to extend (default: 1000)
- `output_filename`: Output filename (default: genes.bed)

**Response:** BED file download

### `GET /health`
Health check - returns `{"status": "healthy"}`

## Output Format

Standard BED format (tab-separated):
```
chr7    55017820    55212628    EGFR
chr17   7660779     7688546     TP53
chr17   43043292    43171245    BRCA1
```

## Technical Details

- **Assembly:** GRCh38 (hg38)
- **API:** Ensembl REST API (rest.ensembl.org)
- **Coordinates:** 1-based (Ensembl native)
- **Standard Chromosomes:** chr1-22, chrX, chrY

## Requirements

- Python 3.8+
- FastAPI
- pandas
- requests
- natsort

## License

Open source - available for academic and research purposes.
