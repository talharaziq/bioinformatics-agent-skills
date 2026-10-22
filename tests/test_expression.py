import tempfile, unittest
from pathlib import Path
from bioskills.expression import expression_qc

ROOT=Path(__file__).resolve().parents[1]
class ExpressionQCTests(unittest.TestCase):
    def test_qc_writes_outputs_and_aligns_metadata(self):
        with tempfile.TemporaryDirectory() as d:
            result=expression_qc(ROOT/'examples/expression/counts.tsv',d,ROOT/'examples/expression/metadata.tsv')
            self.assertEqual(result['genes_or_features'],8)
            self.assertEqual(result['samples'],4)
            self.assertTrue((Path(d)/'sample_qc.tsv').is_file())
            self.assertTrue((Path(d)/'pairwise_similarity.tsv').is_file())
            self.assertTrue((Path(d)/'report.md').is_file())
            self.assertEqual(len(result['pairwise_pearson']),6)
            self.assertIn('condition',result['sample_qc'][0])

    def test_fractional_values_warn(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'m.tsv'; p.write_text('gene_id\ts1\ts2\ng1\t1.5\t2.0\ng2\t0\t4.5\n')
            result=expression_qc(p,Path(d)/'out')
            self.assertEqual(result['fractional_values'],2)
            self.assertTrue(any('Fractional' in w for w in result['warnings']))

    def test_metadata_mismatch_fails(self):
        with tempfile.TemporaryDirectory() as d:
            meta=Path(d)/'metadata.tsv'; meta.write_text('sample_id\tcondition\ns1\tA\ns3\tB\n')
            with self.assertRaisesRegex(ValueError,'sample mismatch'):
                expression_qc(ROOT/'examples/expression/counts.tsv',Path(d)/'out',meta)

if __name__=='__main__': unittest.main()
