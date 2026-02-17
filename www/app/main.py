"""FastAPI application for BLAST searches."""

from contextlib import asynccontextmanager
from typing import Optional

from fastapi import FastAPI, Query, Form, Request
from fastapi.middleware.cors import CORSMiddleware

from blast_service import run_blast, get_seq, get_config


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application lifespan handler."""
    yield


app = FastAPI(
    title="BLAST Search API",
    description="FastAPI-based BLAST search service for SGD",
    version="2.0.0",
    lifespan=lifespan
)

# CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/")
def hello():
    """Health check endpoint."""
    return "Hello, we all love SGD!!"


@app.get("/blast_search")
@app.post("/blast_search")
async def blast_search(
    request: Request,
    # Query parameters
    name: Optional[str] = Query(default=None, description="Gene name for sequence lookup"),
    type: Optional[str] = Query(default=None, description="Sequence type: protein, pep, or dna"),
    conf: Optional[str] = Query(default=None, description="Configuration name: blast-sgd or blast-fungal"),
    seq: Optional[str] = Query(default=None, description="Query sequence"),
    database: Optional[str] = Query(default=None, description="Database(s) to search"),
    program: Optional[str] = Query(default=None, description="BLAST program"),
    seqname: Optional[str] = Query(default="unknown", description="Sequence name"),
    blastType: Optional[str] = Query(default=None, description="Blast type: sgd or fungal"),
    outFormat: Optional[str] = Query(default=None, description="Output format"),
    matrix: Optional[str] = Query(default=None, description="Scoring matrix"),
    threshold: Optional[str] = Query(default=None, description="Threshold value"),
    cutoffScore: Optional[str] = Query(default=None, description="E-value cutoff"),
    alignToShow: Optional[str] = Query(default=None, description="Number of alignments to show"),
    wordLength: Optional[str] = Query(default=None, description="Word size"),
    filter: Optional[str] = Query(default=None, description="Low complexity filter"),
):
    """
    BLAST search endpoint.

    Supports three modes:
    1. Sequence lookup: ?name=S000000001&type=protein
    2. Configuration lookup: ?conf=blast-sgd
    3. BLAST search: ?seq=ATGC...&database=YeastORF-Genomic&program=blastn
    """
    # Handle form data for POST requests
    form_data = {}
    if request.method == "POST":
        try:
            form_data = await request.form()
            form_data = dict(form_data)
        except Exception:
            pass

    # Override query params with form data if present
    seq = form_data.get('seq', seq)
    database = form_data.get('database', database)
    program = form_data.get('program', program)
    seqname = form_data.get('seqname', seqname) or "unknown"
    blastType = form_data.get('blastType', blastType)
    outFormat = form_data.get('outFormat', outFormat)
    matrix = form_data.get('matrix', matrix)
    threshold = form_data.get('threshold', threshold)
    cutoffScore = form_data.get('cutoffScore', cutoffScore)
    alignToShow = form_data.get('alignToShow', alignToShow)
    wordLength = form_data.get('wordLength', wordLength)
    filter = form_data.get('filter', filter)

    # Mode 1: Sequence lookup
    if name and program is None:
        data = get_seq(name, type)
        return data

    # Mode 2: Configuration lookup
    if conf and program is None:
        data = get_config(conf)
        return data

    # Mode 3: BLAST search
    if seq and database and program:
        data = run_blast(
            seq=seq,
            database=database,
            program=program,
            seqname=seqname,
            blast_type=blastType,
            outFormat=outFormat,
            matrix=matrix,
            threshold=threshold,
            cutoffScore=cutoffScore,
            alignToShow=alignToShow,
            wordLength=wordLength,
            filter=filter
        )
        return data

    return {"error": "Missing required parameters: seq, database, and program"}


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
