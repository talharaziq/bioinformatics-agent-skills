import tempfile, unittest
from pathlib import Path
from bioskills.manifest import validate_manifest
ROOT=Path(__file__).resolve().parents[1]
class ManifestTests(unittest.TestCase):
    def test_example_manifest_validates_fastq(self):
        result=validate_manifest(ROOT/'examples/expression/manifest.tsv',check_fastq=True)
        self.assertEqual(result['n_samples'],4)
        self.assertEqual(result['conditions'],{'control':2,'treatment':2})
        self.assertTrue(all(s['layout']=='single-end' for s in result['samples']))

    def test_duplicate_sample_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x.fq'; p.write_text('@r\nACGT\n+\nIIII\n')
            m=Path(d)/'manifest.tsv'; m.write_text('sample_id\tcondition\tfastq_r1\tfastq_r2\ns1\tA\tx.fq\t\ns1\tB\tx.fq\t\n')
            with self.assertRaisesRegex(ValueError,'duplicate sample_id'):
                validate_manifest(m)

    def test_missing_fastq_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            m=Path(d)/'manifest.tsv'; m.write_text('sample_id\tcondition\tfastq_r1\tfastq_r2\ns1\tA\tno.fq\t\ns2\tB\tno2.fq\t\n')
            with self.assertRaises(FileNotFoundError): validate_manifest(m)

    def test_paired_read_ids_are_checked(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d)
            (d/'a_R1.fq').write_text('@readA/1\nACGT\n+\nIIII\n')
            (d/'a_R2.fq').write_text('@readB/2\nACGT\n+\nIIII\n')
            (d/'b.fq').write_text('@other\nACGT\n+\nIIII\n')
            m=d/'manifest.tsv'
            m.write_text('sample_id\tcondition\tfastq_r1\tfastq_r2\ns1\tA\ta_R1.fq\ta_R2.fq\ns2\tB\tb.fq\t\n')
            with self.assertRaisesRegex(ValueError,'identifiers differ'):
                validate_manifest(m,check_fastq=True)

if __name__=='__main__': unittest.main()
