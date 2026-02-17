"""BLAST service module for running NCBI BLAST+ searches."""

import json
import os
import subprocess
from typing import Dict, List, Tuple, Any, Optional

from blast_markup import markup_output

# Configuration
DATASET_PATH = os.environ.get('DATA_DIR', '/data/blast/')
BIN_PATH = os.environ.get('BIN_PATH', '/tools/blast/bin/')
TMP_PATH = os.environ.get('TMP_DIR', '/var/tmp/')
CONF_DIR = os.environ.get('CONF_DIR', '/var/www/conf/')

MAX_GRAPH_HITS = 50
LINE_WIDTH = 70
PVALUE_CUTOFF = 0.05

PROTEIN_LOOKUP_FILE = DATASET_PATH + 'YeastORF.pep'
DNA_LOOKUP_FILE = DATASET_PATH + 'YeastORF-Genomic.fsa'


def get_config(conf: str) -> Dict[str, Any]:
    """Load configuration from JSON file."""
    config_path = os.path.join(CONF_DIR, f"{conf}.json")
    with open(config_path, 'r') as f:
        return json.load(f)


def get_seq(name: str, seq_type: Optional[str] = None) -> Dict[str, str]:
    """Retrieve sequence by name from FASTA file."""
    lookup_file = DNA_LOOKUP_FILE
    if seq_type and seq_type in ['protein', 'pep']:
        lookup_file = PROTEIN_LOOKUP_FILE

    with open(lookup_file, 'r') as f:
        seq = ''
        found = False
        for line in f:
            if line.startswith('>'):
                if found:
                    break
                pieces = line.replace('>', '').split(' ')
                name_list = [pieces[0].upper()]
                if len(pieces) > 1 and pieces[1]:
                    name_list.append(pieces[1].upper())
                if len(pieces) > 2:
                    name_list.append(pieces[2].replace('SGDID:', '').replace(',', ''))
                if name.upper() in name_list:
                    found = True
            elif found:
                seq = seq + line.strip()

    return {'seq': seq}


def set_dataset_mapping(conf_file: str) -> Dict[str, Any]:
    """Build dataset mapping from configuration file."""
    dataset_mapping: Dict[str, Any] = {}
    datasetList4labelNm: Dict[str, List[str]] = {}
    groupLabel2groupNm: Dict[str, str] = {}

    conf_data = get_config(conf_file)
    dataset_mapping['databasedef'] = conf_data.get('databasedef', {})
    datagroup = conf_data.get('datagroup', {})

    for label in datagroup:
        dataset_list = datagroup[label].split(',')
        datasetList4labelNm[label] = dataset_list

    database = conf_data.get('database', [])
    for d in database:
        db = d['dataset']
        db_type = d['type']
        desc = d['label']
        mapping = dataset_mapping.get('dbList', [])
        mapping.append(db)
        dataset_mapping['dbList'] = mapping
        if 'dbType' not in dataset_mapping:
            dataset_mapping['dbType'] = {}
        dataset_mapping['dbType'][db] = db_type
        if db.startswith('label'):
            desc = '...' + desc
        if 'dbLabel' not in dataset_mapping:
            dataset_mapping['dbLabel'] = {}
        dataset_mapping['dbLabel'][db] = desc
        if db.startswith('label'):
            label = db.strip()
            groupNm = desc.lower()
            groupLabel2groupNm[label] = groupNm

    dataset_mapping['groupLabel2name'] = groupLabel2groupNm
    dataset_mapping['groupLabel2datasets'] = datasetList4labelNm

    return dataset_mapping


def combine_datasets(dataset_list: List[str], program: str, conf_file: str) -> Tuple[str, List[str]]:
    """Combine datasets into groups when possible."""
    dataset_mapping = set_dataset_mapping(conf_file)

    dataset_passed_in = {dataset: 1 for dataset in dataset_list}

    groupLabel2group = dataset_mapping['groupLabel2name']
    groupLabel2datasets = dataset_mapping['groupLabel2datasets']

    for label in groupLabel2datasets:
        dataset4label = groupLabel2datasets[label]
        not_found = False
        for dataset in dataset4label:
            if dataset not in dataset_list:
                not_found = True
                break
        if not_found:
            continue

        groupNm = groupLabel2group.get(label)
        if groupNm:
            dataset_passed_in[groupNm] = 1
            for dataset in dataset4label:
                if dataset in dataset_passed_in:
                    del dataset_passed_in[dataset]

    processed_datasets = ''
    processed_dataset_list = []

    for dataset in sorted(dataset_passed_in.keys()):
        if program in ['blastp', 'blastx']:
            dataset = dataset + ".pep"
        else:
            dataset = dataset + ".fsa"
        processed_dataset_list.append(dataset)
        processed_datasets = processed_datasets + ' ' + DATASET_PATH + "fungi/" + dataset

    return (processed_datasets.strip(), processed_dataset_list)


def prepare_datasets(datasets: str, program: str, conf_file: str) -> Tuple[str, List[str]]:
    """Prepare dataset paths for BLAST search."""
    datasets = datasets.replace(',', ' ').replace('+', ' ').replace('  ', ' ')
    dataset_list = datasets.split(' ')

    if "fungal" in conf_file and len(dataset_list) > 20:
        return combine_datasets(dataset_list, program, conf_file)

    processed_datasets = ''
    processed_dataset_list = []

    for dataset in dataset_list:
        if 'fungal' in conf_file:
            if program in ['blastp', 'blastx']:
                if '.pep' not in dataset:
                    dataset = dataset + '.pep'
            elif '.fsa' not in dataset:
                dataset = dataset + '.fsa'
        else:
            if program in ['blastp', 'blastx']:
                if '_cds' in dataset:
                    dataset = dataset.replace('_cds', '_pep')
                elif dataset.endswith('fsa'):
                    dataset = dataset.replace('.fsa', '.pep')
                elif '.pep' not in dataset:
                    dataset = dataset + ".pep"
            elif '.fsa' not in dataset:
                dataset = dataset + '.fsa'

        processed_dataset_list.append(dataset)
        if 'fungal' in conf_file:
            processed_datasets = processed_datasets + ' ' + DATASET_PATH + 'fungi/' + dataset
        else:
            processed_datasets = processed_datasets + ' ' + DATASET_PATH + dataset

    return (processed_datasets.strip(), processed_dataset_list)


def create_tmp_seq_file(query_file: str, seq: str, seqname: str) -> None:
    """Write sequence to temporary FASTA file."""
    with open(query_file, 'w') as fw:
        while len(seq) > LINE_WIDTH:
            if seq.startswith('>'):
                seq = ''
            else:
                fw.write(seq[0:LINE_WIDTH] + "\n")
                seq = seq[LINE_WIDTH:]
        fw.write(seq + "\n")


def _pvalue_to_exp(pvalue: float) -> float:
    """Convert p-value to exponent for sorting."""
    if pvalue == 0:
        return -500
    elif '-' in str(pvalue):
        return 0 - int(str(pvalue).split('-')[1])
    else:
        return pvalue


def _get_id_desc(line: str) -> Tuple[str, str]:
    """Extract ID and description from BLAST header line."""
    pieces = line.replace('>', '').strip().split(' ')
    id_str = pieces[0]
    desc = pieces[1] if len(pieces) > 1 else ''
    if len(pieces) > 2:
        desc = desc + ' ' + pieces[2].replace(',', '')
    return (id_str, desc)


def _get_coords(start: Optional[int], end: Optional[int]) -> Tuple[int, int, int]:
    """Get normalized coordinates and length."""
    if start is None or end is None:
        return (0, 0, 0)
    if start > end:
        (start, end) = (end, start)
    length = end - start + 1
    return (start, end, length)


def _set_name(id_str: str, pvalue: float, score: str, desc: str) -> str:
    """Format hit name with p-value and score."""
    return id_str + ': p=' + str(pvalue) + ' s=' + score + ' ' + desc


def _set_strand(line: str) -> int:
    """Determine strand from Sbjct line."""
    pieces = line.replace('Sbjct ', '').strip().split(' ')
    if int(pieces[0]) > int(pieces[-1]):
        return -1
    return 1


def parse_hits(blast_outfile: str) -> Tuple[int, int, List[Dict[str, Any]]]:
    """Parse BLAST output file for graphical display."""
    records = []
    total_hits = 0
    show_hits = 0

    start = None
    end = None
    id_str = None
    query_length = None
    strand = None
    same_row = 0
    pvalue: float = 0.0
    score = ''
    desc = ''

    with open(blast_outfile, 'r') as f:
        for line in f:
            if line.startswith('Length=') and query_length is None:
                query_length = int(line.strip().replace('Length=', ''))
            if line.startswith('>'):
                total_hits = total_hits + 1
                if id_str and start and end:
                    exp = _pvalue_to_exp(pvalue)
                    (start_norm, end_norm, length) = _get_coords(start, end)
                    if pvalue <= PVALUE_CUTOFF and show_hits < MAX_GRAPH_HITS:
                        show_hits = show_hits + 1
                        records.append({
                            'query_length': query_length,
                            'name': _set_name(id_str, pvalue, score, desc),
                            'value': length,
                            'start': start_norm,
                            'end': end_norm,
                            'strand': strand,
                            'exp': exp,
                            'same_row': same_row
                        })
                        start = None
                        end = None
                        id_str = None
                        strand = None
                        same_row = 0
                        pvalue = 0.0
                        score = ''
                        desc = ''
                (id_str, desc) = _get_id_desc(line)
            elif "Score = " in line:
                if start and end:
                    exp = _pvalue_to_exp(pvalue)
                    (start_norm, end_norm, length) = _get_coords(start, end)
                    if pvalue <= PVALUE_CUTOFF:
                        records.append({
                            'query_length': query_length,
                            'name': _set_name(id_str, pvalue, score, desc),
                            'value': length,
                            'start': start_norm,
                            'end': end_norm,
                            'strand': strand,
                            'exp': exp,
                            'same_row': same_row
                        })
                        same_row = 1
                    start = None
                    end = None
                    strand = None
                score = line.split('Score = ')[1].split(' ')[0]
                pvalue = float(line.split('Expect = ')[1].split(' ')[0].replace(',', ''))
            elif line.startswith("Query "):
                pieces = line.replace('Query ', '').strip().split(' ')
                if start is None:
                    start = int(pieces[0])
                end = int(pieces[-1])
            elif line.startswith("Sbjct ") and strand is None:
                strand = _set_strand(line)

    if id_str and start and end:
        exp = _pvalue_to_exp(pvalue)
        (start_norm, end_norm, length) = _get_coords(start, end)
        if pvalue <= PVALUE_CUTOFF and show_hits < MAX_GRAPH_HITS:
            show_hits = show_hits + 1
            records.append({
                'query_length': query_length,
                'name': _set_name(id_str, pvalue, score, desc),
                'value': length,
                'start': start_norm,
                'end': end_norm,
                'strand': strand,
                'exp': exp,
                'same_row': same_row
            })

    return (total_hits, show_hits, records)


def get_blast_options(
    program: str,
    database: Optional[str] = None,
    outFormat: Optional[str] = None,
    matrix: Optional[str] = None,
    threshold: Optional[str] = None,
    cutoffScore: Optional[str] = None,
    alignToShow: Optional[str] = None,
    wordLength: Optional[str] = None,
    filter: Optional[str] = None
) -> str:
    """Build BLAST command line options."""
    options = ''

    if cutoffScore:
        options = "-evalue " + cutoffScore
    if alignToShow:
        options = options + " -num_alignments " + alignToShow

    if database == 'Sc_mito_chr' and program in ['blastn', 'tblastx']:
        options = options + " -query_genetic_code 3"

    if program != 'blastn' and threshold and threshold != 'default':
        options = options + " -threshold " + threshold

    if program != 'blastn' and matrix and matrix != "BLOSUM62":
        options = options + " -matrix " + matrix

    if wordLength:
        if wordLength != 'default':
            options = options + " -word_size " + wordLength
        else:
            if program == 'blastn':
                options = options + " -word_size 11"
            else:
                options = options + " -word_size 3"

    options = options + " -outfmt 0"

    if outFormat and outFormat.startswith("ungapped"):
        options = options + " -ungapped"

    if filter:
        if filter == 'on':
            if program == 'blastp':
                options = options + " -seg yes"
        else:
            if program == 'blastn':
                options = options + " -dust 'no'"
            elif program != 'blastp':
                options = options + " -seg 'no'"

    return options


def run_blast(
    seq: str,
    database: str,
    program: str,
    seqname: str = "unknown",
    blast_type: Optional[str] = None,
    outFormat: Optional[str] = None,
    matrix: Optional[str] = None,
    threshold: Optional[str] = None,
    cutoffScore: Optional[str] = None,
    alignToShow: Optional[str] = None,
    wordLength: Optional[str] = None,
    filter: Optional[str] = None
) -> Dict[str, Any]:
    """Execute BLAST search and return results."""
    conf_file = 'blast-sgd'
    if blast_type == 'fungal':
        conf_file = 'blast-fungal'

    (processed_datasets, processed_dataset_list) = prepare_datasets(database, program, conf_file)

    program_path = BIN_PATH + program

    options = get_blast_options(
        program=program,
        database=database,
        outFormat=outFormat,
        matrix=matrix,
        threshold=threshold,
        cutoffScore=cutoffScore,
        alignToShow=alignToShow,
        wordLength=wordLength,
        filter=filter
    )

    pid = os.getpid()
    query_file = TMP_PATH + 'tmp.fsa.' + str(pid)
    create_tmp_seq_file(query_file, seq, seqname)

    blast_outfile = TMP_PATH + 'blast.out.' + str(pid)

    cmd = f"{program_path} -query {query_file} -db '{processed_datasets}' -out {blast_outfile} {options}"

    subprocess.run(cmd, shell=True, check=False)

    with open(blast_outfile, 'r') as f:
        output = f.read()

    # Clean up temp files
    try:
        os.remove(query_file)
        os.remove(blast_outfile)
    except OSError:
        pass

    if "** No hits found **" in output:
        return {
            "cmd": cmd,
            "result": "<font size=+1><pre>" + output + "</pre></font>",
            "hits": [],
            "totalHits": 0,
            "showHits": 0
        }

    output = markup_output(database, output)

    records = []
    total_hits = 0
    show_hits = 0

    try:
        # Re-read file for parsing (need to keep it temporarily)
        query_file_tmp = TMP_PATH + 'tmp.fsa.' + str(pid)
        blast_outfile_tmp = TMP_PATH + 'blast.out.' + str(pid)
        create_tmp_seq_file(query_file_tmp, seq, seqname)
        subprocess.run(cmd, shell=True, check=False)
        (total_hits, show_hits, records) = parse_hits(blast_outfile_tmp)
        try:
            os.remove(query_file_tmp)
            os.remove(blast_outfile_tmp)
        except OSError:
            pass
    except Exception:
        pass

    return {
        "cmd": cmd,
        "result": "<pre><p font-size='11px'>" + output + "</p></pre>",
        "hits": records,
        "totalHits": total_hits,
        "showHits": show_hits
    }
