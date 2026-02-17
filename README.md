# BlastDocker

A FastAPI-based microservice for NCBI BLAST+ sequence similarity searches. Supports protein and nucleotide searches against yeast (SGD) and fungal genome databases.

## Features

- **NCBI BLAST+**: Full support for blastn, blastp, blastx, tblastn, and tblastx programs
- **Multiple Databases**: Search against SGD yeast datasets or multiple fungal genomes
- **Rich Output**: HTML-formatted results with SGD locus pages and JBrowse genome browser links
- **Graphical Hits**: Returns parsed hits for graphical display of alignment coverage
- **Sequence Retrieval**: Look up sequences by gene name or SGD ID

## Quick Start

### Using Docker

```bash
# Build the image
docker build -t blast-fastapi .

# Run with data volumes
docker run -d \
  -p 8000:8000 \
  -v /path/to/blast/data:/data/blast:ro \
  blast-fastapi
```

### Using Docker Compose

```bash
# Set data paths and run
DATA_DIR=/path/to/blast/data \
docker-compose up -d
```

### Local Development

```bash
# Install dependencies
pip install -r requirements.txt

# Run the application
cd www/app
uvicorn main:app --reload --host 0.0.0.0 --port 8000
```

## API Endpoints

### Health Check

```
GET /
```

Returns: `"Hello, we all love SGD!!"`

### BLAST Search

```
GET/POST /blast_search
```

**Search Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `seq` | string | required | Query sequence |
| `database` | string | required | Database(s) to search |
| `program` | string | required | BLAST program: blastn, blastp, blastx, tblastn, tblastx |
| `seqname` | string | `unknown` | Sequence name |
| `blastType` | string | `sgd` | Blast type: sgd or fungal |
| `outFormat` | string | null | Output format (ungapped for ungapped alignment) |
| `matrix` | string | `BLOSUM62` | Scoring matrix |
| `threshold` | string | null | Threshold value |
| `cutoffScore` | string | null | E-value cutoff |
| `alignToShow` | string | null | Number of alignments to show |
| `wordLength` | string | null | Word size |
| `filter` | string | null | Low complexity filter: on or off |

**Example:**

```bash
# Search for nucleotide sequence
curl "http://localhost:8000/blast_search?seq=ATGCATGC&database=YeastORF-Genomic&program=blastn"

# Search for protein sequence
curl "http://localhost:8000/blast_search?seq=MFVL&database=YeastORF&program=blastp"

# Search fungal databases
curl "http://localhost:8000/blast_search?seq=ATGCATGC&database=Afumigatus&program=blastn&blastType=fungal"
```

**Response:**

```json
{
  "cmd": "blastn -query ... -db '/data/blast/YeastORF-Genomic.fsa' ...",
  "result": "<pre>...</pre>",
  "hits": [
    {
      "query_length": 100,
      "name": "YAL001C: p=1e-50 s=200 Subunit...",
      "value": 95,
      "start": 1,
      "end": 95,
      "strand": 1,
      "exp": -50,
      "same_row": 0
    }
  ],
  "totalHits": 10,
  "showHits": 10
}
```

### Sequence Retrieval

```
GET /blast_search?name=<gene_name>&type=<sequence_type>
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `name` | string | Gene name (e.g., ACT1, YFL039C) or SGD ID |
| `type` | string | Sequence type: protein, pep, or dna |

**Example:**

```bash
# Get protein sequence
curl "http://localhost:8000/blast_search?name=ACT1&type=protein"

# Get DNA sequence
curl "http://localhost:8000/blast_search?name=YFL039C"
```

**Response:**

```json
{
  "seq": "MDSEVAALVIDNGSGMCKAGFAGDDAPRAVFPSIVGRP..."
}
```

### Configuration

```
GET /blast_search?conf=<config_name>
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `conf` | string | Configuration name: blast-sgd or blast-fungal |

**Example:**

```bash
curl "http://localhost:8000/blast_search?conf=blast-sgd"
```

## BLAST Programs

| Program | Query | Database | Description |
|---------|-------|----------|-------------|
| `blastn` | Nucleotide | Nucleotide | Nucleotide-nucleotide search |
| `blastp` | Protein | Protein | Protein-protein search |
| `blastx` | Nucleotide | Protein | Translated query vs protein database |
| `tblastn` | Protein | Nucleotide | Protein query vs translated database |
| `tblastx` | Nucleotide | Nucleotide | Translated query vs translated database |

## Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `DATA_DIR` | `/data/blast/` | Path to BLAST database files |
| `BIN_PATH` | `/tools/blast/bin/` | Path to NCBI BLAST+ binaries |
| `TMP_DIR` | `/var/tmp/` | Temporary file directory |
| `CONF_DIR` | `/var/www/conf/` | Configuration file directory |

## Data Files

### BLAST Databases (`/data/blast/`)

- `YeastORF.pep` - Yeast protein sequences
- `YeastORF-Genomic.fsa` - Yeast genomic sequences
- `Sc_nuclear.fsa` - S. cerevisiae nuclear chromosome sequences
- `Sc_mito_chr.fsa` - S. cerevisiae mitochondrial sequences
- `fungi/` - Fungal genome databases

### Configuration Files (`/var/www/conf/`)

- `blast-sgd.json` - SGD BLAST configuration
- `blast-fungal.json` - Fungal BLAST configuration

## Project Structure

```
BlastDocker/
├── Dockerfile              # Docker image definition
├── docker-compose.yml      # Docker Compose configuration
├── requirements.txt        # Python dependencies
├── README.md               # This file
└── www/
    ├── app/                # FastAPI application
    │   ├── main.py         # Application entry point
    │   ├── schemas.py      # Pydantic models
    │   ├── blast_service.py    # BLAST search service
    │   └── blast_markup.py     # Output formatting with links
    └── conf/
        ├── blast-sgd.json      # SGD configuration
        └── blast-fungal.json   # Fungal configuration
```

## Development

### Running Tests

```bash
# Install dev dependencies
pip install pytest httpx

# Run tests
pytest tests/
```

### Building for Production

```bash
# Build optimized image
docker build -t blast-fastapi:latest .

# Tag for registry
docker tag blast-fastapi:latest your-registry/blast-fastapi:latest

# Push to registry
docker push your-registry/blast-fastapi:latest
```

## API Compatibility

This service maintains backward compatibility with the legacy Flask-based API. All endpoints accept the same parameters and return the same response format.

## License

Copyright (c) Stanford University / SGD Project

## Related Projects

- [SGD (Saccharomyces Genome Database)](https://www.yeastgenome.org/)
- [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
