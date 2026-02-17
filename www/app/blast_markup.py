"""BLAST output markup utilities for SGD/JBrowse links."""

import string

ROOT_URL = 'https://www.yeastgenome.org/'
LOCUS = ROOT_URL + 'locus/'
SEQ = ROOT_URL + 'seqTools?chr='
SEQAN = ROOT_URL + 'seqTools?seqname='
GBROWSE = 'https://jbrowse.yeastgenome.org/?loc='
NCBI_URL = 'https://www.ncbi.nlm.nih.gov/nuccore/'


def number2roman() -> dict:
    """Map chromosome numbers to Roman numerals."""
    num2rom = {}
    i = 0
    for roman in ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI',
                  'XII', 'XIII', 'XIV', 'XV', 'XVI', 'mt']:
        i = i + 1
        num2rom[str(i)] = roman
    num2rom['2-micron'] = '2-micron'
    return num2rom


def roman2number() -> dict:
    """Map Roman numerals to chromosome numbers."""
    num2rom = number2roman()
    rom2num = {}
    for num in num2rom:
        rom = num2rom[num]
        if num == '17':
            rom = 'Mito'
        rom2num[rom] = num
    return rom2num


def letter2number() -> dict:
    """Map letters to chromosome numbers."""
    letter2num = {}
    for letter in list(string.ascii_uppercase):
        if letter == 'Q':
            break
        letter2num[letter] = str(ord(letter) - 64)
    letter2num['Q'] = 'MT'
    letter2num['R'] = '2-micron'
    return letter2num


def record_line(line: str) -> str:
    """Add NCBI link to record line."""
    pieces = line.split("|")
    ncbi_id = pieces[1]
    url = NCBI_URL + ncbi_id
    return "<a href='" + url + "' target='_new'>" + pieces[0] + "|" + pieces[1] + "</a>|" + "|".join(pieces[2:])


def link_out_for_feature(line: str) -> str:
    """Generate SGD/JBrowse links for a feature."""
    piece = line.split(', Chr ')
    name = piece[0].split(' ')[0].replace('>', '')
    chr_coords = piece[1].split(' from ')
    chr_name = chr_coords[0]
    coords = chr_coords[1].split(' ')[0]
    coord_list = coords.split(',')
    beg = 0
    end = 0

    for coord_pair in coord_list:
        if '-' not in coord_pair:
            continue
        [tmp_beg, tmp_end] = coord_pair.split('-')
        if int(tmp_beg) > int(tmp_end):
            (tmp_beg, tmp_end) = (tmp_end, tmp_beg)
        if beg == 0:
            beg = int(tmp_beg)
        elif beg > int(tmp_beg):
            beg = int(tmp_beg)
        if end == 0:
            end = int(tmp_end)
        elif end < int(tmp_end):
            end = int(tmp_end)

    gsr_link = "<a href='" + SEQAN + name + "' target='_new'>Retrieve Sequence</a>"

    start = beg - 5000
    if start < 1:
        start = 1
    stop = end + 5000

    jbrowse_link = (
        "<a href='" + GBROWSE + "chr" + chr_name + "%3A" + str(start) + ".." + str(stop) +
        "&tracks=DNA%2CAll%20Annotated%20Sequence%20Features&highlight=chr" + chr_name +
        "%3A" + str(beg) + ".." + str(end) + "'  target='_new'>Genome Browser</a>"
    )

    lsp_link = "<a href='" + LOCUS + name + "' target='_blastout'>SGD Locus page</a> ]</b>"

    return "  <b>[ " + gsr_link + " / " + jbrowse_link + " / " + lsp_link + " ]</b>"


def link_out(chr_name: str, chrnum: str, beg: str, end: str) -> str:
    """Generate sequence retrieval and JBrowse links."""
    gsr_link = (
        "<a href='" + SEQ + chrnum + "&start=" + beg + "&end=" + end +
        "&submit2=Submit+Form" + "' target='_new'>Retrieve Sequence</a>"
    )

    beg_int = int(beg)
    end_int = int(end)
    if beg_int > end_int:
        (beg_int, end_int) = (end_int, beg_int)
    start = beg_int - 5000
    if start < 1:
        start = 1
    stop = end_int + 5000

    jbrowse_link = (
        "<a href='" + GBROWSE + "chr" + chr_name + "%3A" + str(start) + ".." + str(stop) +
        "&tracks=DNA%2CAll%20Annotated%20Sequence%20Features&highlight=chr" + chr_name +
        "%3A" + str(beg_int) + ".." + str(end_int) + "'  target='_new'>Genome Browser</a>"
    )

    return "  [ " + gsr_link + " / " + jbrowse_link + " ]"


def markup_chromosomal_coord(blast_output: str) -> str:
    """Add links to chromosomal coordinate records."""
    lines = blast_output.split('\n')

    rom2num = roman2number()
    newoutput = ''
    start_record = 0
    record_before_score = ''
    scoreline = ''
    record_after_score = ''
    chr_name = ''
    chrnum = ''
    beg = ''
    end = ''
    sub_hit = 0

    for line in lines:
        if line.startswith('>'):
            start_record = 1
        if start_record == 0:
            if "ref|NC_" in line or "gi|" in line:
                newline = record_line(line)
                newoutput = newoutput + newline + "\n"
            else:
                newoutput = newoutput + line + "\n"
        else:
            if line.startswith('>'):
                if record_before_score != '':
                    scoreline = scoreline + link_out(chr_name, chrnum, beg, end) + "\n"
                    newoutput = newoutput + record_before_score + "\n" + scoreline + "\n" + record_after_score + "\n"
                chr_name = ''
                chrnum = ''
                beg = ''
                end = ''
                sub_hit = 0
                newline = record_line(line)
                record_before_score = newline + "\n"
                record_after_score = ''
            elif " Score = " in line:
                scoreline = line
            elif "[chromosome=" in line:
                record_before_score = record_before_score + line + "\n"
                pieces = line.split("[chromosome=")
                chr_name = pieces[-1].replace(']', '')
                chrnum = rom2num.get(chr_name)
            elif scoreline == '':
                record_before_score = record_before_score + line + "\n"
            else:
                record_after_score = record_after_score + line + "\n"
                if beg != '' and end != '' and 'Identities =' in line:
                    sub_hit = 1
                if sub_hit == 1:
                    continue
                if line.startswith('Sbjct '):
                    pieces = line.split(' ')
                    tmp_beg = ''
                    tmp_end = ''
                    for item in pieces:
                        if item.isdigit():
                            if tmp_beg == '':
                                tmp_beg = item
                            else:
                                tmp_end = item
                    if beg == '':
                        beg = tmp_beg
                    end = tmp_end

    if record_before_score != '':
        scoreline = scoreline + link_out(chr_name, chrnum, beg, end) + "\n"
        newoutput = newoutput + record_before_score + "\n" + scoreline + "\n" + record_after_score + "\n"

    return newoutput


def markup_output(dataset: str, blast_output: str) -> str:
    """Add SGD/JBrowse links to BLAST output."""
    output_lines = blast_output.split('\n')

    blast_output = ''
    end = 0
    databases = ''
    finish_db_line = 0
    prev_line = ''
    summary_end = 0

    for line in output_lines:
        if line.startswith('BLAST'):
            continue
        if prev_line == '' and line == '':
            continue
        prev_line = line
        if "Sequences producing significant alignment" in line:
            end = 1
        if line.startswith('>'):
            summary_end = 1
        if end == 0:
            if '/data/blast/' in line:
                if finish_db_line == 1:
                    continue
                line = line.strip().replace("/data/blast/fungi/", "").replace("/data/blast/", "")
                if databases == '':
                    databases = line
                else:
                    databases = databases + " " + line
                    if databases.endswith(';'):
                        databases = databases + " etc"
                    finish_db_line = 1
                continue
            if " sequences; " in line and " total letters" in line:
                line = databases.strip() + "\n" + line
            if 'Query=' in line:
                continue
            if 'Length=' in line:
                line = 'Query ' + line
        if end == 1 and summary_end == 0 and "Sequences producing significant alignment" not in line:
            items = line.strip().split(' ')
            if len(items) > 3:
                name = items[0]
                if '|' in name:
                    name = name.split('|')[1]
                items[-1] = "<a href='#" + name + "'>" + items[-1] + "</a>"
                line = ' '.join(items)

        if ".fsa" in line and "Database:" not in line:
            continue

        if summary_end == 1:
            if line.startswith('>'):
                name = line.replace('>', '').strip().split(' ')[0]
                if '|' in name:
                    name = name.split('|')[1]
                blast_output = blast_output + "\n" + "<a name='" + name + "'>\n"

        blast_output = blast_output + "\n" + line

    if ("Sc_nuclear" in dataset or "Sc_mito" in dataset) and "[chromosome=2-micron]" not in blast_output:
        blast_output = markup_chromosomal_coord(blast_output)
    else:
        if "NotFeature" not in dataset:
            data = []
            lines = blast_output.split("\n")
            link_list = ''
            for line in lines:
                if line.startswith('>') and "SGDID:" in line:
                    data.append("<hr><p>")
                    data.append(line)
                    link_list = link_out_for_feature(line)
                elif line.startswith('>gi') or line.startswith('gi|'):
                    newline = record_line(line)
                    data.append(newline)
                elif line.startswith("Length="):
                    data.append("<br>" + link_list)
                    data.append("<br>" + line)
                    link_list = ''
                else:
                    data.append(line)
            blast_output = "\n".join(data)

            blast_output = blast_output.replace(
                " Verified ORF", " <font color='green'>Verified ORF</font>"
            )
            blast_output = blast_output.replace(
                " Uncharacterized ORF", " <font color='green'>Uncharacterized ORF</font>"
            )
            blast_output = blast_output.replace(
                " Dubious ORF", "<br> <font color='red'>** Dubious ORF **</font>"
            )

    return blast_output
