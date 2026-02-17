"""Pydantic models for BLAST API requests and responses."""

from typing import Optional, List, Dict, Any
from pydantic import BaseModel, Field


class BlastHit(BaseModel):
    """Individual BLAST hit for graphical display."""
    query_length: int
    name: str
    value: int
    start: int
    end: int
    strand: Optional[int] = None
    exp: float
    same_row: int = 0


class BlastRequest(BaseModel):
    """BLAST search request parameters."""
    seq: str = Field(..., description="Query sequence")
    database: str = Field(..., description="Database(s) to search")
    program: str = Field(..., description="BLAST program: blastn, blastp, blastx, tblastn, tblastx")
    seqname: Optional[str] = Field(default="unknown", description="Sequence name")
    blastType: Optional[str] = Field(default=None, description="Blast type: sgd or fungal")
    outFormat: Optional[str] = Field(default=None, description="Output format")
    matrix: Optional[str] = Field(default=None, description="Scoring matrix")
    threshold: Optional[str] = Field(default=None, description="Threshold value")
    cutoffScore: Optional[str] = Field(default=None, description="E-value cutoff")
    alignToShow: Optional[str] = Field(default=None, description="Number of alignments to show")
    wordLength: Optional[str] = Field(default=None, description="Word size")
    filter: Optional[str] = Field(default=None, description="Low complexity filter: on or off")


class BlastResponse(BaseModel):
    """BLAST search response."""
    cmd: str = Field(..., description="BLAST command executed")
    result: str = Field(..., description="HTML formatted BLAST output")
    hits: List[BlastHit] = Field(default_factory=list, description="Parsed hits for graphical display")
    totalHits: int = Field(default=0, description="Total number of hits found")
    showHits: int = Field(default=0, description="Number of hits displayed in graph")


class SeqResponse(BaseModel):
    """Sequence retrieval response."""
    seq: str = Field(..., description="Retrieved sequence")


class ConfigResponse(BaseModel):
    """Configuration response."""
    databasedef: Optional[Dict[str, str]] = None
    datagroup: Optional[Dict[str, str]] = None
    database: Optional[List[Dict[str, Any]]] = None
