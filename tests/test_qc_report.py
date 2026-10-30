"""Offline regression tests for scientific denominators and QC report integration.

Run: RSCRIPT_BIN=/path/to/Rscript python -m unittest discover -s tests -p test_qc_report.py -v
"""
import copy
import csv
import gzip
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import threading
import unittest
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'bin'/'qc'))
import prepare_qc_report as prepare
import render_qc_report as renderer
import interpret_qc as interpret


def interpretation(samples):
    return {'summary': 'Observed sample differences warrant descriptive review.',
            'sample_reviews': {s: {'core': 'Review the observed depth and MT distributions.',
                'reads': 'Read identifiers describe the stated processing stage.',
                'correlations': 'Aggregate relationships do not establish accuracy.',
                'structure': 'Interpret structural categories with their annotation coverage.'} for s in samples},
            'attention': [], 'group_notes': {'reads': 'Source stages are explicit.',
                'correlations': 'Within-platform correlations are descriptive.',
                'structure': 'Class percentages use detected features.'}}


class QCReportTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory(prefix='sclra-qc-test-')
        cls.path = Path(cls.tmp.name)
        cls.fixture = cls.path/'fixture'
        cls.rscript = os.environ.get('RSCRIPT_BIN', 'Rscript')
        subprocess.run([cls.rscript, str(ROOT/'tests'/'make_qc_fixture.R'), str(cls.fixture)],
                       check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        cls.metrics = cls.path/'metrics'
        cls.extract(cls.metrics)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    @classmethod
    def extract(cls, output, gene='gene.rds', transcript='transcript.rds', classification=None, check=True):
        cmd = [cls.rscript, str(ROOT/'bin'/'qc'/'extract_qc_metrics.R'),
               '--gene', str(cls.fixture/gene), '--transcript', str(cls.fixture/transcript),
               '--isoquant_dir', str(cls.fixture/'isoquant'), '--reference_gtf', str(cls.fixture/'reference.gtf'),
               '--output_dir', str(output)]
        if classification:
            cmd.extend(['--classification', str(cls.fixture/classification)])
        return subprocess.run(cmd, check=check, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

    def test_linear_counts_alignment_and_mt_recomputation(self):
        records = json.loads((self.metrics/'sample_metrics.json').read_text())
        self.assertEqual([r['sample'] for r in records], ['alpha', 'beta'])
        self.assertEqual(records[0]['gene_counts_median'], 7)
        self.assertEqual(records[0]['transcript_counts_median'], 9)
        self.assertAlmostEqual(records[0]['mt_median_pct'], 200/7)
        self.assertEqual(records[0]['double_cell_true'], 1)  # Existing TRUE labels are not silently removed.
        self.assertEqual(records[0]['double_cell_available'], 3)
        self.assertEqual(records[0]['annotation_coverage_pct'], 75)
        self.assertEqual(records[0]['fsm_pct'], 25)
        self.assertEqual(records[0]['model_length_median'], 200)

    def test_valid_double_cell_records_define_cohort(self):
        output = self.path/'missing-flag'
        self.extract(output, gene='gene_missing_flag.rds')
        manifest = json.loads((output/'extraction.json').read_text())
        self.assertEqual(manifest['report_cells'], 5)
        self.assertEqual(manifest['excluded_without_valid_double_cell'], 1)

    def test_skipped_doublet_detection_is_not_reported_as_negative(self):
        output = self.path/'skipped'
        self.extract(output, gene='gene_skipped_doublets.rds')
        self.assertIsNone(json.loads((output/'extraction.json').read_text())['doublet_source'])
        self.assertTrue(all('double_cell_available' not in r for r in json.loads((output/'sample_metrics.json').read_text())))
        self.extract(self.path/'performed', gene='gene_performed_doublets.rds')
        performed = json.loads((self.path/'performed'/'sample_metrics.json').read_text())
        self.assertEqual(performed[0]['double_cell_available'], 3)
        self.assertEqual(performed[0]['double_cell_true'], 0)

    def test_seurat_v5_split_counts_and_reference_mt_ids(self):
        output = self.path/'v5'
        self.extract(output, gene='gene_v5.rds')
        record = json.loads((output/'sample_metrics.json').read_text())[0]
        self.assertEqual(record['gene_counts_median'], 7)
        self.assertAlmostEqual(record['mt_median_pct'], 200/7)

    def test_no_supported_mt_metric_is_omitted(self):
        output = self.path/'no-mt'
        self.extract(output, gene='gene_no_mt.rds')
        self.assertTrue(all('mt_median_pct' not in r for r in json.loads((output/'sample_metrics.json').read_text())))
        data = prepare.assemble(output)
        self.assertNotIn('MT', data['summary'])
        self.assertTrue(all('MT' not in r['core'] for r in data['sample_reviews'].values()))

    def test_transcript_cohort_mismatch_fails(self):
        result = self.extract(self.path/'mismatch', transcript='transcript_missing_cell.rds', check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('does not cover', result.stderr)

    def test_headered_and_headerless_classification_agree(self):
        output = self.path/'headered'
        self.extract(output, classification='classification.tsv')
        self.assertEqual(json.loads((output/'sample_metrics.json').read_text()),
                         json.loads((self.metrics/'sample_metrics.json').read_text()))

    def test_read_group_distinct_ids_and_distinct_umi_strings(self):
        cells = prepare.read_cells(self.metrics/'cells.tsv.gz')
        values, source = prepare.summarize_read_groups(self.fixture/'read_groups.tsv', cells)
        self.assertEqual(len(values), 6)
        self.assertTrue(all(r == {'unique_read_ids': 3, 'unique_umi_strings': 2} for r in values.values()))
        self.assertEqual(source['stage'], 'Filtered FASTQs before alignment')
        self.assertTrue(source['validated_umi_headers'])

    def test_unrecognized_umi_headers_do_not_fabricate_umi_counts(self):
        cells = prepare.read_cells(self.metrics/'cells.tsv.gz')
        path = self.path/'unknown-headers.tsv'
        path.write_text(''.join('read{}\t{}\n'.format(i, r['cell_id']) for i, r in enumerate(cells)))
        values, source = prepare.summarize_read_groups(path, cells)
        self.assertFalse(source['validated_umi_headers'])
        self.assertTrue(all('unique_umi_strings' not in v for v in values.values()))

    @unittest.skipUnless(importlib.util.find_spec('pysam'), 'BAM test requires pysam from sclong')
    def test_bam_read_support_joins_ids_and_deduplicates_alignments(self):
        import pysam
        cells = prepare.read_cells(self.metrics/'cells.tsv.gz')
        bam = self.path/'fixture.bam'
        mapping = (self.fixture/'read_groups.tsv').read_text().splitlines()
        with pysam.AlignmentFile(str(bam), 'wb', header={'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chr1', 'LN': 100}]}) as handle:
            for row in cells:
                ids = [line.split('\t')[0] for line in mapping if line.endswith('\t'+row['cell_id'])][:2]
                for name in ids+[ids[0]]:
                    read = pysam.AlignedSegment()
                    read.query_name = name
                    read.query_sequence = 'AAAA'
                    read.flag = 4  # Stored unmapped records are counted by the stated BAM policy.
                    read.query_qualities = pysam.qualitystring_to_array('IIII')
                    handle.write(read)
        values, source = prepare.summarize_read_groups(self.fixture/'read_groups.tsv', cells, [bam, bam])
        self.assertTrue(all(v == {'unique_read_ids': 2, 'unique_umi_strings': 2} for v in values.values()))
        self.assertEqual(source['stage'], 'BAM-derived read IDs')

    def test_incomplete_or_duplicate_read_summary_fails(self):
        cells = prepare.read_cells(self.metrics/'cells.tsv.gz')
        path = self.path/'bad-summary.tsv'
        path.write_text('CellBarcode\tunique_reads\nalpha_AAAA\t2\n')
        with self.assertRaisesRegex(ValueError, 'every retained'):
            prepare.load_read_summary(path, cells)
        path.write_text(path.read_text()+'alpha_AAAA\t3\n')
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            prepare.load_read_summary(path, cells)

    def test_render_package_and_bam_metrics(self):
        import shutil
        output = self.path/'package-input'
        shutil.copytree(str(self.metrics), str(output))
        data = prepare.assemble(output, read_summary=self.fixture/'bam_summary.tsv', project='QC <fixture>')
        self.assertEqual(data['records'][0]['read_ids_total'], 60)
        self.assertEqual(data['records'][0]['read_ids_median'], 20)
        self.assertEqual(data['statistics'][2]['exact'], 35)
        self.assertEqual(data['provenance']['read_support']['stage'], 'BAM-derived read IDs')
        subprocess.run([self.rscript, str(ROOT/'bin'/'qc'/'plot_qc.R'), str(output/'cells.tsv.gz'),
                        str(output/'qc_distributions.svg')], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        package = self.path/'package'
        subprocess.run([sys.executable, str(ROOT/'bin'/'qc'/'publish_qc_report.py'), '--data_dir', str(output),
                        '--output_dir', str(package)], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        html = (package/'scLRA_LLM_QC_report.html').read_text()
        self.assertIn('QC &lt;fixture&gt;', html)
        self.assertIn('data:image/svg+xml;base64,', html)
        self.assertNotIn('@@', html)
        with (package/'scLRA_LLM_QC_report.tsv').open(encoding='utf-8-sig') as handle:
            rows = list(csv.reader(handle, delimiter='\t'))
        self.assertEqual(len(rows), 3)
        self.assertEqual(rows[0][-1], 'Assessment')

    def test_no_reference_ranges_and_undefined_correlations(self):
        self.assertIsNone(prepare.correlation([1,1,1], [2,3,4]))
        self.assertIsNone(prepare.correlation([1,2], [2,3]))
        self.assertEqual(prepare.correlation([1,2,3], [3,2,1]), -1)
        self.assertIsNone(prepare.unique_extreme([{'v': 1}, {'v': 1}], 'v'))
        self.assertIsNone(prepare.unique_extreme([{'v': 1}], 'v'))

    def test_llm_validation_and_numeric_immutability(self):
        import shutil
        output = self.path/'llm-base'
        shutil.copytree(str(self.metrics), str(output))
        data = prepare.assemble(output)
        value = interpretation(['alpha', 'beta'])
        parsed = interpret.parse_interpretation(json.dumps(value), ['alpha', 'beta'])
        updated = interpret.apply_interpretation(data, parsed)
        self.assertEqual([{k:v for k,v in r.items() if k != 'assessment'} for r in data['records']],
                         [{k:v for k,v in r.items() if k != 'assessment'} for r in updated['records']])
        value['sample_reviews']['invented'] = value['sample_reviews'].pop('alpha')
        with self.assertRaisesRegex(ValueError, 'exactly'):
            interpret.parse_interpretation(json.dumps(value), ['alpha', 'beta'])
        value = interpretation(['alpha', 'beta'])
        value['summary'] = 'word '*46
        with self.assertRaisesRegex(ValueError, 'overlong'):
            interpret.parse_interpretation(json.dumps(value), ['alpha', 'beta'])

    def test_mock_api_success_failure_and_strict_mode(self):
        import shutil
        output = self.path/'api-input'
        shutil.copytree(str(self.metrics), str(output))
        data = prepare.assemble(output)
        captured = []
        response = {'content': json.dumps(interpretation(['alpha', 'beta']))}
        class Handler(BaseHTTPRequestHandler):
            def log_message(self, *args):
                pass
            def do_POST(self):
                captured.append(json.loads(self.rfile.read(int(self.headers['Content-Length']))))
                body = json.dumps({'id': 'mock-qc', 'object': 'chat.completion', 'created': 0,
                    'model': 'mock-model', 'choices': [{'index': 0, 'message': {'role': 'assistant',
                    'content': response['content']}, 'finish_reason': 'stop'}]}).encode()
                self.send_response(200)
                self.send_header('Content-Type', 'application/json')
                self.send_header('Content-Length', str(len(body)))
                self.end_headers()
                self.wfile.write(body)
        server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        env = dict(os.environ, LLM_API_KEY='local-mock-key')
        cmd = [sys.executable, str(ROOT/'bin'/'qc'/'interpret_qc.py'), '--data', str(output/'qc_report.json'),
               '--base_url', 'http://127.0.0.1:{}/v1'.format(server.server_port), '--model', 'mock-model']
        try:
            destination = self.path/'api-success'
            subprocess.run(cmd+['--output_dir', str(destination)], env=env, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            result = json.loads((destination/'qc_report.json').read_text())
            self.assertEqual(result['provenance']['llm']['status'], 'completed')
            self.assertNotIn('cell_id', json.dumps(captured))
            self.assertNotIn('alpha_AAAA', json.dumps(captured))
            self.assertNotIn('local-mock-key', ''.join(p.read_text() for p in destination.rglob('*.json')))
            response['content'] = 'invalid json'
            destination = self.path/'api-fallback'
            subprocess.run(cmd+['--output_dir', str(destination)], env=env, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            result = json.loads((destination/'qc_report.json').read_text())
            self.assertEqual(result['records'], data['records'])
            self.assertEqual(result['provenance']['llm']['status'], 'failed')
            self.assertIn('LLM unavailable', result['interpretation_label'])
            strict = subprocess.run(cmd+['--output_dir', str(self.path/'api-strict'), '--strict'], env=env,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertNotEqual(strict.returncode, 0)
        finally:
            server.shutdown()
            server.server_close()
            thread.join(timeout=5)


if __name__ == '__main__':
    unittest.main()
