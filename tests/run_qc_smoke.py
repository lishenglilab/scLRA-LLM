#!/usr/bin/env python3
"""Run the actual Nextflow QC processes offline, with local mock LLM responses.

python tests/run_qc_smoke.py --output /tmp/new-qc-smoke --rscript /path/to/Rscript
"""
import argparse
import json
import os
import subprocess
import sys
import sys
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

from test_qc_report import interpretation

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--rscript', default=os.environ.get('RSCRIPT_BIN', 'Rscript'))
    parser.add_argument('--nextflow', default='nextflow')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    fixture = output/'fixture'
    subprocess.run([args.rscript, str(ROOT/'tests'/'make_qc_fixture.R'), str(fixture)], check=True)
    state, requests = {'mode': 'success'}, []
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *args):
            pass
        def do_POST(self):
            request = json.loads(self.rfile.read(int(self.headers['Content-Length'])))
            requests.append(request)
            payload = json.loads(request['messages'][1]['content'])
            samples = [r['sample'] for r in payload['records']]
            content = json.dumps(interpretation(samples)) if state['mode'] == 'success' else 'malformed JSON'
            body = json.dumps({'id': 'local-qc', 'object': 'chat.completion', 'created': 0,
                'model': 'mock-model', 'choices': [{'index': 0, 'message': {'role': 'assistant', 'content': content},
                'finish_reason': 'stop'}]}).encode()
            self.send_response(200)
            self.send_header('Content-Type', 'application/json')
            self.send_header('Content-Length', str(len(body)))
            self.end_headers()
            self.wfile.write(body)
    server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    results, baseline = [], None
    try:
        for mode in ('disabled', 'success', 'fallback'):
            state['mode'] = mode
            run = output/mode
            run.mkdir()
            env = dict(os.environ, NXF_OFFLINE='true', LLM_API_KEY='local-mock-key')
            env['PATH'] = str(Path(sys.executable).parent) + os.pathsep + env.get('PATH', '')
            env['PATH'] = str(Path(sys.executable).parent) + os.pathsep + env.get('PATH', '')
            command = [args.nextflow, '-c', str(ROOT/'tests'/'qc_smoke.config'), 'run',
                str(ROOT/'tests'/'qc_report_smoke.nf'), '-ansi-log', 'false', '-work-dir', str(run/'work'),
                '--fixture', str(fixture), '--pipeline_bin', str(ROOT/'bin'), '--outdir', str(run/'results'),
                '--rscript_bin', args.rscript, '--llm_enabled', 'false' if mode == 'disabled' else 'true',
                '--llm_base_url', 'http://127.0.0.1:{}/v1'.format(server.server_port),
                '--llm_model', 'mock-'+mode, '-with-trace', str(run/'trace.tsv')]
            # Isolate Nextflow session/cache state from production and other tests.
            proc = subprocess.run(command, cwd=str(run), env=env, stdout=subprocess.PIPE,
                                  stderr=subprocess.STDOUT, universal_newlines=True)
            (run/'nextflow.log').write_text(proc.stdout)
            if proc.returncode:
                raise RuntimeError('Nextflow QC smoke failed: '+str(run/'nextflow.log'))
            report = run/'results'/'PreliminarySingleCellAnalysis'/'qc_report'
            data = json.loads((report/'data'/'qc_report.json').read_text())
            numeric = [{k: v for k, v in row.items() if k != 'assessment'} for row in data['records']]
            if baseline is None:
                baseline = numeric
            assert numeric == baseline, 'LLM changed numerical metrics'
            assert sum(r['cells'] for r in data['records']) == 6
            status = data['provenance']['llm']['status']
            assert status == {'disabled': 'disabled', 'success': 'completed', 'fallback': 'failed'}[mode]
            assert (report/'scLRA_LLM_QC_report.html').stat().st_size > 10000
            assert (report/'scLRA_LLM_QC_report.tsv').exists()
            assert (report/'data'/'qc_distributions.svg').exists()
            if mode != 'disabled':
                assert (run/'results'/'PreliminarySingleCellAnalysis'/'deepseek_qc'/'sc_qc_report.txt').exists()
                assert (report/'llm'/'request.json').exists()
            results.append({'mode': mode, 'status': 'passed', 'llm_status': status,
                            'samples': len(data['records']), 'cells': 6, 'report': str(report)})
            print('Nextflow QC smoke: {} passed'.format(mode), flush=True)
        assert len(requests) == 2, 'Disabled QC made an unexpected LLM request'
        assert all('cell_id' not in json.dumps(r) and 'alpha_AAAA' not in json.dumps(r) for r in requests)
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)
    (output/'summary.json').write_text(json.dumps({'status': 'passed', 'runs': results,
        'aggregate_only_requests': True, 'numeric_metrics_unchanged': True}, indent=2)+'\n')


if __name__ == '__main__':
    main()
