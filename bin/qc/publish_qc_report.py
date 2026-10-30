#!/usr/bin/env python3
"""Package the standalone HTML and its reproducible numerical inputs."""
import argparse
import json
import shutil
from pathlib import Path

from render_qc_report import render, validate, write_tsv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data_dir', required=True, type=Path)
    parser.add_argument('--interpreted', type=Path)
    parser.add_argument('--output_dir', required=True, type=Path)
    args = parser.parse_args()
    source = args.interpreted or args.data_dir/'qc_report.json'
    data = validate(json.loads(source.read_text()))
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    shutil.copytree(str(args.data_dir), str(output/'data'))
    (output/'data'/'qc_report.json').write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False)+'\n')
    if args.interpreted:
        audit = args.interpreted.parent/'llm'
        if audit.is_dir():
            shutil.copytree(str(audit), str(output/'llm'))
    # The legacy text output remains available even with interpretation disabled.
    lines = [data['run']['project']+' — filtered-cell QC', data['summary'], '']
    lines.extend('{}: {}'.format(r['sample'], data['sample_reviews'][r['sample']]['core']) for r in data['records'])
    lines.extend(['', data['interpretation_label']])
    (output/'sc_qc_report.txt').write_text('\n'.join(lines)+'\n')
    (output/'scLRA_LLM_QC_report.html').write_text(render(data, output/'data'/'qc_report.json'))
    write_tsv(data, output/'scLRA_LLM_QC_report.tsv')
    print('Packaged filtered-cell QC report: {} samples, {} cells'.format(len(data['records']), sum(r['cells'] for r in data['records'])))


if __name__ == '__main__':
    main()
