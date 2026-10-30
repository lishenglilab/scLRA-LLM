#!/usr/bin/env python3
"""Assemble observed, per-sample QC metrics and a concise offline baseline."""
import argparse
import csv
import gzip
import json
import math
import re
import statistics
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path


def read_cells(path):
    with gzip.open(str(path), 'rt', newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if not rows or len({r['cell_id'] for r in rows}) != len(rows):
        raise ValueError('Expected non-empty, unique filtered cell IDs')
    for row in rows:
        for key in list(row):
            if key not in ('cell_id', 'sample', 'doublet_label'):
                row[key] = None if row[key] in ('', 'NA') else float(row[key])
    return rows


def write_cells(rows, path):
    keys = list(rows[0])
    with gzip.open(str(path), 'wt', newline='') as handle:
        writer = csv.DictWriter(handle, keys, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: 'NA' if v is None else v for k, v in row.items()})


def finite_values(rows, key):
    return [r[key] for r in rows if isinstance(r.get(key), (float, int)) and math.isfinite(r[key])]


def correlation(x, y):
    if len(x) < 3:
        return None
    ax, ay = statistics.mean(x), statistics.mean(y)
    xx = math.fsum((v-ax)**2 for v in x)
    yy = math.fsum((v-ay)**2 for v in y)
    if not xx or not yy:
        return None
    return max(-1.0, min(1.0, math.fsum((a-ax)*(b-ay) for a, b in zip(x, y))/math.sqrt(xx*yy)))


def load_read_summary(path, cells):
    """An explicitly supplied BAM-derived table must cover the exact QC cohort."""
    selected = {r['cell_id'] for r in cells}
    result, seen = {}, set()
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(str(path), 'rt', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not {'CellBarcode', 'unique_reads'}.issubset(reader.fieldnames or []):
            raise ValueError('Read summary requires CellBarcode and unique_reads columns')
        for row in reader:
            name = row['CellBarcode']
            if name in seen:
                raise ValueError('Duplicate CellBarcode in the supplied read summary')
            seen.add(name)
            if name not in selected:
                continue
            values = {}
            for source, target in [('unique_reads', 'unique_read_ids'), ('unique_UMIs', 'unique_umi_strings')]:
                if source not in row:
                    continue
                value = float(row[source])
                if not math.isfinite(value) or value < 0 or value != int(value):
                    raise ValueError('Read/UMI summary values must be non-negative integer counts')
                values[target] = int(value)
            result[name] = values
    if set(result) != selected:
        raise ValueError('Supplied read summary does not cover every retained report cell')
    return result, {'source': 'BAM-derived summary', 'stage': 'BAM-derived read IDs',
                    'file': Path(path).name, 'cells': len(result)}


UMI_HEADER = re.compile(r'^(?P<barcode>[ACGTN]+)_(?P<umi>[ACGTN]+)#[^\s]+_[+-]$')


def summarize_read_groups(path, cells, bam_paths=()):
    """Count distinct ReadIDs, with disk sorting instead of an in-memory ID set.

    BAM query names, when supplied, are joined to the filtered read-group mapping
    on disk. Otherwise the source is filtered FASTQs before alignment. UMI-string
    counts are available only for validated BLAZE read names.
    """
    selected = {r['cell_id'] for r in cells}
    result = {}
    umi_valid = True
    with tempfile.TemporaryDirectory(prefix='qc-read-sort-', dir='.') as tmp:
        raw, ordered = Path(tmp)/'selected.tsv', Path(tmp)/'ordered.tsv'
        with Path(path).open() as handle, raw.open('w') as output:
            for line_number, line in enumerate(handle, 1):
                fields = line.rstrip('\r\n').split('\t')
                if len(fields) != 2 or not all(fields):
                    raise ValueError('Malformed read-groups row {}'.format(line_number))
                read_id, cell = fields
                if cell in selected:
                    output.write((read_id+'\t'+cell if bam_paths else cell+'\t'+read_id)+'\n')
        if bam_paths:
            import pysam  # Supplied by the sclong environment; not needed for offline rendering.
            bam_ids = Path(tmp)/'bam_ids.txt'
            with bam_ids.open('w') as handle:
                for bam in bam_paths:
                    with pysam.AlignmentFile(str(bam), 'rb', check_sq=False) as alignments:
                        for record in alignments.fetch(until_eof=True):
                            if record.query_name:
                                handle.write(record.query_name+'\n')
            sorted_map, sorted_ids = Path(tmp)/'mapping.tsv', Path(tmp)/'ids.txt'
            for source, destination in [(raw, sorted_map), (bam_ids, sorted_ids)]:
                subprocess.run(['sort', '-u', '-S', '256M', '-T', tmp, '-o', str(destination), str(source)], check=True)
            joined = Path(tmp)/'joined.tsv'
            with joined.open('w') as handle:
                subprocess.run(['join', '-t', '\t', str(sorted_map), str(sorted_ids)], stdout=handle, check=True)
            with joined.open() as handle, raw.open('w') as output:
                for line in handle:
                    read_id, cell = line.rstrip('\n').split('\t')
                    output.write(cell+'\t'+read_id+'\n')
        subprocess.run(['sort', '-u', '-S', '256M', '-T', tmp, '-o', str(ordered), str(raw)], check=True)
        current, n_reads, umis = None, 0, set()
        with ordered.open() as handle:
            for line in handle:
                cell, read_id = line.rstrip('\n').split('\t')
                if cell != current:
                    if current is not None:
                        result[current] = {'unique_read_ids': n_reads, 'unique_umi_strings': len(umis)}
                    current, n_reads, umis = cell, 0, set()
                n_reads += 1
                match = UMI_HEADER.fullmatch(read_id)
                if match and match.group('barcode') == cell.split('_', 1)[1]:
                    umis.add(match.group('umi'))
                else:
                    umi_valid = False
            if current is not None:
                result[current] = {'unique_read_ids': n_reads, 'unique_umi_strings': len(umis)}
    if set(result) != selected:
        raise ValueError('Filtered read groups do not cover every retained report cell')
    if not umi_valid:
        for row in result.values():
            row.pop('unique_umi_strings')
    return result, {'source': 'IsoQuant BAM query names joined to filtered read groups' if bam_paths else 'Filtered FASTQ read groups',
                    'stage': 'BAM-derived read IDs' if bam_paths else 'Filtered FASTQs before alignment',
                    'file': Path(path).name, 'bam_files': [Path(p).name for p in bam_paths],
                    'cells': len(result), 'validated_umi_headers': umi_valid,
                    'bam_policy': 'All stored alignment records, including unmapped; distinct query names per cell' if bam_paths else None}


def unique_extreme(records, key, high=False):
    candidates = [r for r in records if r.get(key) is not None]
    if len(candidates) < 2 or len({r[key] for r in candidates}) < 2:
        return None
    value = (max if high else min)(r[key] for r in candidates)
    hits = [r for r in candidates if r[key] == value]
    return hits[0] if len(hits) == 1 else None


def assemble(directory, read_groups=None, read_summary=None, project='scLRA-LLM', isoquant_dir=None, read_source='auto'):
    directory = Path(directory)
    cells = read_cells(directory/'cells.tsv.gz')
    records = json.loads((directory/'sample_metrics.json').read_text())
    manifest = json.loads((directory/'extraction.json').read_text())
    if read_summary:
        support, read_source = load_read_summary(read_summary, cells)
    elif read_groups:
        bams = sorted(Path(isoquant_dir).rglob('*.bam')) if isoquant_dir and read_source != 'fastq' else []
        if read_source == 'bam' and not bams:
            raise ValueError('BAM read support requested but no BAM exists in the IsoQuant output')
        support, read_source = summarize_read_groups(read_groups, cells, bams)
    else:
        support, read_source = {}, None
    for row in cells:
        row.update(support.get(row['cell_id'], {}))
    for record in records:
        d = [r for r in cells if r['sample'] == record['sample']]
        if len(d) != record['cells']:
            raise ValueError('Cell count differs between sample and per-cell tables')
        reads = finite_values(d, 'unique_read_ids')
        if len(reads) == len(d):
            record.update(read_ids_total=sum(reads), read_ids_mean=statistics.mean(reads),
                          read_ids_median=statistics.median(reads),
                          r_read_ids_gene_counts=correlation(reads, [r['gene_counts'] for r in d]))
        umis = finite_values(d, 'unique_umi_strings')
        if len(umis) == len(d):
            record['umi_strings_median'] = statistics.median(umis)
    low_key = 'read_ids_median' if all(r.get('read_ids_median') is not None for r in records) else 'gene_counts_median'
    low = unique_extreme(records, low_key)
    high_mt = unique_extreme(records, 'mt_median_pct', True)
    attention = []
    if low:
        attention.append({'sample': low['sample'], 'label': 'Lower depth',
                          'detail': 'Lowest median {} in this run.'.format('Read IDs' if low_key == 'read_ids_median' else 'gene counts')})
    if high_mt:
        attention.append({'sample': high_mt['sample'], 'label': 'Higher MT',
                          'detail': 'Highest median mitochondrial fraction in this run.'})
    reviews = {}
    for record in records:
        core = []
        if low and record['sample'] == low['sample']:
            core.append('Lowest median depth; review sample balance.')
        if high_mt and record['sample'] == high_mt['sample']:
            core.append('Highest median MT; review the distribution.')
        review = {'core': ' '.join(core) or ('Single-sample descriptive QC; no between-sample comparison.' if len(records) == 1 else
                                           'Descriptive sample summary; no unique review priority identified.'),
                  'attention': bool(core)}
        if record.get('read_ids_median') is not None:
            review['reads'] = 'Read-ID {} is {} the median.'.format('mean',
              'above' if record['read_ids_mean'] > record['read_ids_median'] else
              'below' if record['read_ids_mean'] < record['read_ids_median'] else 'equal to')
        else:
            review['reads'] = 'Detection ratio describes feature complexity.'
        if record.get('double_cell_available') is not None:
            review['reads'] += ' Existing doublet labels are descriptive.'
        review['correlations'] = 'Within-sample aggregate relationships; no accuracy score is assigned.'
        review['structure'] = 'Detected feature counts are descriptive.'
        if record.get('annotation_coverage_pct') is not None:
            review['structure'] = ('Complete' if record['annotation_coverage_pct'] == 100 else 'Partial') + ' structural annotation coverage; classes do not establish isoform accuracy.'
        reviews[record['sample']] = review
        record['assessment'] = review['core']
    n_cells = sum(r['cells'] for r in records)
    summary = 'Descriptive QC covers {:,} filtered cells across {} sample{}.'.format(n_cells, len(records), '' if len(records) == 1 else 's')
    summary += ' Review the within-run differences highlighted below.' if attention else ' No universal pass/fail score is assigned.'
    statistics_rows = [{'label': 'Retained cells', 'value': '{:,}'.format(n_cells), 'exact': n_cells}]
    if support:
        total = sum(r['read_ids_total'] for r in records)
        median = statistics.median(r['unique_read_ids'] for r in cells)
        statistics_rows.extend([{'label': 'Read IDs · total', 'value': '{:,}'.format(total), 'exact': total},
                                {'label': 'Read IDs / cell · median', 'value': '{:,.1f}'.format(median).rstrip('0').rstrip('.'), 'exact': median}])
    mt = finite_values(cells, 'mitochondrial_percent')
    if mt:
        median = statistics.median(mt)
        statistics_rows.append({'label': 'MT · median', 'value': '{:.2f}%'.format(median), 'exact': median})
    definitions = [
        ['Gene / transcript counts', 'Column totals in the corresponding linear RNA count matrices; not assumed to be UMI molecules.'],
        ['Detected genes / transcripts', 'Features with counts greater than zero in a cell.'],
        ['Detection ratio', 'Detected transcripts / detected genes in the same cell.'],
        ['Pearson r / Spearman ρ', 'Across cells within each sample, for the pairs stated in the headers; omitted when undefined.'],
        ['Detected features', 'Distinct features with positive counts in that sample.'],
    ]
    if support:
        definitions.insert(0, ['Read IDs', 'Distinct read IDs per retained cell from {}. Sample totals sum only retained cells.'.format(read_source['stage'])])
    if any(r.get('umi_strings_median') is not None for r in records):
        definitions.append(['UMI strings', 'Distinct UMI strings pooled across genes within a cell; not gene-aware molecule counts.'])
    if mt:
        definitions.append(['MT', 'Mitochondrial gene counts / total gene counts × 100, recomputed from the RNA counts.'])
    doublet_source = manifest.get('doublet_source')
    if doublet_source:
        definitions.append([doublet_source, 'Existing {} records in the filtered gene object. False is a label, not a measured residual doublet rate.'.format(doublet_source)])
    has_classes = any(r.get('annotation_coverage_pct') is not None for r in records)
    if has_classes:
        definitions.extend([
            ['Annotation coverage', 'Detected transcript features matched to a structural-classification record / all detected transcript features × 100.'],
            ['FSM / ISM', 'Reference-matching / incomplete reference-matching splice structures.'],
            ['NIC / NNC', 'New combinations of known splice sites / models containing new splice sites.'],
            ['Other structural classes', 'Antisense: opposite-strand overlap; intergenic: no gene overlap; fusion: multiple genes; genic: other gene overlap.'],
            ['Model length / exons', 'Length in bp and exon count among detected features with available annotation; equally weighted. Model length is not read length.'],
            ['IQR / range', '25th–75th percentiles / minimum–maximum.'],
        ])
    labels = {'double_cell': '{}|True / valid records'.format(doublet_source)} if doublet_source else {}
    export = [{'key': 'sample', 'label': 'Sample'}]
    for key in records[0]:
        if key not in ('sample', 'assessment') and any(r.get(key) is not None for r in records):
            export.append({'key': key, 'label': key})
    # A metric may be available only in a later sample.
    known = {x['key'] for x in export}
    for record in records[1:]:
        for key in record:
            if key not in known and key != 'assessment' and any(r.get(key) is not None for r in records):
                export.append({'key': key, 'label': key})
                known.add(key)
    export.append({'key': 'assessment', 'label': 'Assessment'})
    data = {'schema_version': '1.0', 'run': {'project': project, 'title': 'Single-cell long-read RNA',
            'date': datetime.now(timezone.utc).strftime('%d %B %Y'),
            'scope': 'Filtered cells · {} sample{}'.format(len(records), '' if len(records) == 1 else 's')},
            'summary': summary, 'attention': attention, 'statistics': statistics_rows,
            'records': records, 'sample_reviews': reviews, 'export_columns': export,
            'column_labels': labels,
            'group_notes': {'core': 'Relative sample comparisons are descriptive and depend on sequencing depth and cell composition.',
                            'reads': 'Read-ID source: {}.'.format(read_source['stage']) if read_source else '',
                            'correlations': 'These are within-platform aggregate relationships, not cross-platform validation.'},
            'structure_note': 'Class percentages use all detected transcript features as denominator; unmatched or other classes account for any remainder.' if has_classes else '',
            'definitions': definitions,
            'method_notes': [manifest['cohort_policy']+'. Gene and transcript columns use the same ordered cell IDs.',
                             'Per-cell values are sample medians unless a mean or total is stated.',
                             'Assessments are descriptive; no universal pass/fail score is assigned.'],
            'figure': 'qc_distributions.svg',
            'figure_caption': 'Retained-cell distributions by sample. Count and detection axes use log10(1 + value); MT is linear. All {:,} report cells are included.'.format(n_cells),
            'interpretation_label': 'Metric-based interpretation',
            'provenance': {'extraction': manifest, 'read_support': read_source,
                           'llm': {'status': 'disabled', 'numeric_metrics_modified': False}}}
    write_cells(cells, directory/'cells.tsv.gz')
    (directory/'qc_report.json').write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False)+'\n')
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data_dir', required=True, type=Path)
    parser.add_argument('--read_groups', type=Path)
    parser.add_argument('--read_summary', type=Path)
    parser.add_argument('--project', default='scLRA-LLM')
    parser.add_argument('--isoquant_dir', type=Path)
    parser.add_argument('--read_source', choices=['auto', 'bam', 'fastq'], default='auto')
    args = parser.parse_args()
    data = assemble(args.data_dir, args.read_groups, args.read_summary, args.project, args.isoquant_dir, args.read_source)
    print('Prepared {} sample records for {} retained cells'.format(len(data['records']), sum(r['cells'] for r in data['records'])))


if __name__ == '__main__':
    main()
