#!/usr/bin/env python3
"""Render a portable offline English QC report using only the Python stdlib.

Usage: python render_qc_report.py --data qc_report.json --output qc_report.html
Supply run metadata, observed metrics and concise assessments in the JSON input.
The renderer does not call an LLM or assign quality thresholds.
"""
import argparse
import base64
import csv
import html
import io
import json
import math
from pathlib import Path

ESC = lambda value: html.escape(str(value), quote=True)


def fmt(value, decimals=0):
    if value is None:
        return '—'
    text = f'{value:,.{decimals}f}'
    return text.rstrip('0').rstrip('.') if decimals else text


def cell_value(row, key):
    if key == 'double_cell':
        return f"{fmt(row.get('double_cell_true'))} / {fmt(row.get('double_cell_available'))}"
    if key in ('length_iqr', 'exons_iqr'):
        prefix = 'model_length' if key == 'length_iqr' else 'model_exons'
        return f"{fmt(row.get(prefix+'_median'))} ({fmt(row.get(prefix+'_q25'))}–{fmt(row.get(prefix+'_q75'))})"
    if key in ('length_range', 'exons_range'):
        prefix = 'model_length' if key == 'length_range' else 'model_exons'
        return f"{fmt(row.get(prefix+'_min'))}–{fmt(row.get(prefix+'_max'))}"
    if key.startswith(('r_', 'rho_')):
        return '—' if row.get(key) is None else f'{row[key]:.9f}'
    if key.endswith('_pct'):
        return fmt(row.get(key), 2)
    if key == 'detection_ratio_median':
        return fmt(row.get(key), 4)
    if key == 'sample':
        return row[key]
    return fmt(row.get(key), 1)


TABLES = [
    ('core', 'Overview', [
        ('sample', 'Sample'), ('cells', 'Cells'), ('read_ids_median', 'Read IDs|median / cell'),
        ('gene_counts_median', 'Gene counts|median / cell'), ('transcript_counts_median', 'Transcript counts|median / cell'),
        ('genes_detected_median', 'Detected genes|median / cell'), ('transcripts_detected_median', 'Detected transcripts|median / cell'),
        ('mt_median_pct', 'MT (%)|median / cell'),
    ]),
    ('reads', 'Read IDs & labels', [
        ('sample', 'Sample'), ('read_ids_total', 'Read IDs|sample total'), ('read_ids_mean', 'Read IDs|mean / cell'),
        ('umi_strings_median', 'UMI strings|median / cell'), ('detection_ratio_median', 'Transcript / gene|detection ratio · median'),
        ('double_cell', 'double.cell|True / valid records'),
    ]),
    ('correlations', 'Correlations', [
        ('sample', 'Sample'), ('r_gene_transcript_detection', 'Detected genes vs|detected transcripts|Pearson r'),
        ('rho_gene_transcript_detection', 'Detected genes vs|detected transcripts|Spearman ρ'),
        ('r_gene_transcript_counts', 'Gene total counts vs|transcript total counts|Pearson r'),
        ('r_read_ids_gene_counts', 'Read IDs vs|gene total counts|Pearson r'),
    ]),
    ('structure', 'Transcript models', [
        ('sample', 'Sample'), ('gene_features', 'Gene features|detected / sample'),
        ('transcript_features', 'Transcript features|detected / sample'),
        ('annotation_coverage_pct', 'Annotation|coverage (%)'), ('length_iqr', 'Model length (bp)|median (IQR)'),
        ('exons_iqr', 'Model exons|median (IQR)'), ('fsm_pct', 'FSM (%)'), ('ism_pct', 'ISM (%)'),
        ('nic_pct', 'NIC (%)'), ('nnc_pct', 'NNC (%)'),
    ]),
]
EXTRA_STRUCTURE = [('sample', 'Sample'), ('antisense_pct', 'Antisense (%)'), ('intergenic_pct', 'Intergenic (%)'),
                   ('fusion_pct', 'Fusion (%)'), ('genic_pct', 'Genic (%)'),
                   ('length_range', 'Model length (bp)|range'), ('exons_range', 'Model exons|range')]


def field_present(records, key):
    proxies = {'double_cell': 'double_cell_available', 'length_iqr': 'model_length_median',
               'exons_iqr': 'model_exons_median', 'length_range': 'model_length_min', 'exons_range': 'model_exons_min'}
    return any(r.get(proxies.get(key, key)) is not None for r in records)


def render_table(records, reviews, group, columns, extra=False):
    columns = [(key, label) for key, label in columns if field_present(records, key)]
    headers = ''.join('<th scope="col">'+''.join(f'<span>{ESC(part)}</span>' for part in label.split('|'))+'</th>' for key, label in columns)
    headers += '<th scope="col" class="assessment">Assessment</th>'
    rows = []
    for record in records:
        sample = record['sample']
        review = reviews.get(sample, {})
        cells = []
        for key, label in columns:
            if key == 'sample':
                cells.append(f'<th scope="row">{ESC(sample)}</th>')
            else:
                cells.append(f'<td>{ESC(cell_value(record, key))}</td>')
        assessment = review.get(group, record.get('assessment', ''))
        cells.append(f'<td class="assessment">{ESC(assessment)}</td>')
        rows.append(f'<tr data-sample="{ESC(sample)}" class="{("attention-row" if review.get("attention") and group=="core" else "")}">'+''.join(cells)+'</tr>')
    return '<div class="table-scroll" tabindex="0" role="region" aria-label="'+ESC(group)+' sample metrics"><table class="'+('wide ' if len(columns)>8 else '')+('extra ' if extra else '')+'"><thead><tr>'+headers+'</tr></thead><tbody>'+''.join(rows)+'</tbody></table></div>'


def validate(data):
    if data.get('schema_version') != '1.0' or not data.get('records'):
        raise ValueError('Expected schema_version 1.0 and at least one sample record')
    names = [r['sample'] for r in data['records']]
    if any(not isinstance(n, str) or not n.strip() for n in names) or len(names)!=len(set(names)):
        raise ValueError('Sample names must be non-empty, unique strings')
    for record in data['records']:
        if record.get('cells', -1)<0:
            raise ValueError('Each sample requires a non-negative cell count')
        for key, value in record.items():
            if isinstance(value, (int, float)) and (not math.isfinite(value) or value<0 and not key.startswith(('r_', 'rho_'))):
                raise ValueError(f'Invalid numeric value for {key}')
    if not data.get('summary'):
        raise ValueError('A concise summary is required')
    return data


def render(data, data_path):
    records = data['records']
    reviews = data.get('sample_reviews', {})
    tables = [(key, label, cols) for key, label, cols in TABLES
              if key=='core' or any(field_present(records, k) for k, _ in cols if k!='sample')]
    tabs, panels = [], []
    for key, label, columns in tables:
        tabs.append(f'<button type="button" role="tab" id="tab-{key}" data-panel="panel-{key}" aria-controls="panel-{key}" aria-selected="{str(key=="core").lower()}" tabindex="{0 if key=="core" else -1}">{ESC(label)}</button>')
        columns = [(k, data.get('column_labels', {}).get(k, label)) for k, label in columns]
        table = render_table(records, reviews, key, columns)
        note = data.get('group_notes', {}).get(key)
        if note:
            table += f'<p class="group-note">{ESC(note)}</p>'
        if key=='structure':
            if data.get('structure_note'):
                table = f'<p class="table-note">{ESC(data["structure_note"])}</p>'+table
            if any(field_present(records, k) for k, _ in EXTRA_STRUCTURE if k!='sample'):
                table += '<details class="extra-models"><summary>Other categories &amp; model ranges</summary>'+render_table(records, reviews, 'structure', EXTRA_STRUCTURE, True)+'</details>'
        panels.append(f'<div class="tab-panel" id="panel-{key}" role="tabpanel" aria-labelledby="tab-{key}"{(" hidden" if key!="core" else "")}>{table}</div>')

    statistics = ''.join(f'<div class="stat"><span>{ESC(s["label"])}</span><strong title="{ESC(s.get("exact",s["value"]))}">{ESC(s["value"])}</strong></div>' for s in data.get('statistics', []))
    attention = ''.join(f'<div class="flag"><b>{ESC(item["sample"])}</b><span class="flag-label">{ESC(item["label"])}</span><span>{ESC(item["detail"])}</span></div>' for item in data.get('attention', []))
    options = '<option value="all">All samples</option>'+''.join(f'<option value="{ESC(r["sample"])}">{ESC(r["sample"])}</option>' for r in records)
    definitions = ''.join(f'<div><dt>{ESC(term)}</dt><dd>{ESC(definition)}</dd></div>' for term,definition in data.get('definitions', []))
    notes = ''.join(f'<li>{ESC(note)}</li>' for note in data.get('method_notes', []))
    figure = ''
    if data.get('figure'):
        figure_path = (data_path.parent/data['figure']).resolve()
        suffix = figure_path.suffix.lower()
        mime = {'.svg':'image/svg+xml', '.png':'image/png', '.jpg':'image/jpeg'}.get(suffix)
        if not mime:
            raise ValueError('Figure must be SVG, PNG, or JPG')
        encoded = base64.b64encode(figure_path.read_bytes()).decode()
        figure = f'''<section class="figure-section" id="distributions"><div class="section-head"><h2>QC distributions</h2><button class="text-button" type="button" id="open-figure">Enlarge figure ↗</button></div><figure><img id="qc-figure" src="data:{mime};base64,{encoded}" alt="Retained-cell QC distributions by sample"><figcaption>{ESC(data.get('figure_caption',''))}</figcaption></figure></section>'''
    run = data.get('run', {})
    substitutions = {
        'PROJECT': ESC(run.get('project', 'Single-cell RNA')),
        'TITLE': ESC(run.get('title', 'Single-cell long-read RNA')),
        'DATE': ESC(run.get('date', '')), 'SCOPE': ESC(run.get('scope', 'Filtered cells')),
        'SUMMARY': ESC(data['summary']), 'STATISTICS': statistics, 'ATTENTION': attention,
        'SAMPLE_OPTIONS': options, 'TABS': ''.join(tabs), 'PANELS': ''.join(panels), 'FIGURE': figure,
        'DEFINITIONS': definitions, 'METHOD_NOTES': notes,
        'INTERPRETATION': ESC(data.get('interpretation_label', 'QC assessment')),
        'JSON_DATA': json.dumps(data, ensure_ascii=False, separators=(',',':'), allow_nan=False).replace('</','<\\/'),
    }
    template = Path(__file__).with_name('report_template.html').read_text()
    for key, value in substitutions.items():
        template = template.replace('@@'+key+'@@',value)
    return template


def write_tsv(data, path):
    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter='\t', lineterminator='\n')
    columns = data['export_columns']
    writer.writerow([c['label'] for c in columns])
    for row in data['records']:
        writer.writerow([row.get(c['key'], '') for c in columns])
    path.write_text(buffer.getvalue(), encoding='utf-8-sig')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    data_path = args.data.resolve()
    data = validate(json.loads(data_path.read_text()))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(render(data, data_path))
    write_tsv(data, args.output.with_suffix('.tsv'))
    print(json.dumps({'samples': len(data['records']), 'cells': sum(r['cells'] for r in data['records']), 'html_bytes': args.output.stat().st_size}))


if __name__ == '__main__':
    main()
