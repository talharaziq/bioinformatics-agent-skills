import gzip, json, tempfile, unittest
from pathlib import Path
from bioskills.sequence import SequenceFormatError, sequence_qc

ROOT=Path(__file__).resolve().parents[1]

class SequenceQCTests(unittest.TestCase):
    def test_fasta_summary(self):
        result=sequence_qc(ROOT/'examples/sequences/toy.fasta')
        self.assertEqual(result['format'],'fasta')
        self.assertEqual(result['records'],3)
        self.assertEqual(result['total_bases'],29)
        self.assertEqual(result['min_length'],8)
        self.assertEqual(result['max_length'],12)
        self.assertGreater(result['ambiguous_bases'],0)
        self.assertTrue(0 <= result['gc_percent'] <= 100)

    def test_fastq_gzip_and_per_record_tsv(self):
        with tempfile.TemporaryDirectory() as d:
            gz=Path(d)/'reads.fastq.gz'
            with gzip.open(gz,'wt') as f: f.write('@a\nACGT\n+\nIIII\n@b\nGGGG\n+\nHHHH\n')
            rec=Path(d)/'records.tsv'
            result=sequence_qc(gz,records_tsv=rec)
            self.assertEqual(result['records'],2)
            self.assertEqual(result['mean_phred33_quality'],39.5)
            self.assertTrue(rec.read_text().startswith('id\tlength'))

    def test_rna_uracil_is_not_marked_ambiguous(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'rna.fasta'; p.write_text('>rna\nACGUU\n')
            result=sequence_qc(p)
            self.assertEqual(result['ambiguous_bases'],0)
            self.assertEqual(result['gc_percent'],40.0)

    def test_protein_alphabet_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'protein.fasta'; p.write_text('>protein\nMKWVTFISLL\n')
            with self.assertRaisesRegex(SequenceFormatError,'not protein FASTA'):
                sequence_qc(p)

    def test_bad_fastq_raises(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'bad.fastq'; p.write_text('@read\nACGT\n+\nIII\n')
            with self.assertRaises(SequenceFormatError): sequence_qc(p)

    def test_empty_file_raises(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'empty.fasta'; p.write_text('')
            with self.assertRaises(SequenceFormatError): sequence_qc(p)

if __name__=='__main__': unittest.main()
