import tempfile, unittest
from pathlib import Path
from bioskills.ora import bh_adjust, hypergeom_sf, run_ora
ROOT=Path(__file__).resolve().parents[1]
class ORATests(unittest.TestCase):
    def test_hypergeometric_known_value(self):
        self.assertAlmostEqual(hypergeom_sf(2,4,2,2),1/6,places=10)
        self.assertEqual(hypergeom_sf(3,10,2,4),0.0)

    def test_bh_adjust_monotonic_in_sorted_order(self):
        p=[0.04,0.001,0.02]
        q=bh_adjust(p)
        self.assertAlmostEqual(q[1],0.003)
        self.assertTrue(all(0<=x<=1 for x in q))

    def test_ora_output_uses_explicit_universe(self):
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)/'ora.tsv'
            rows=run_ora(ROOT/'examples/pathways/query_genes.txt',ROOT/'examples/pathways/toy.gmt',out,ROOT/'examples/pathways/universe.txt',min_size=2)
            self.assertEqual(len(rows),3)  # zero-overlap term is tested for FDR but omitted from the table
            self.assertAlmostEqual(rows[0]['fdr_bh'], 4 * (1/495), places=8)
            self.assertTrue(out.is_file())
            self.assertIn('fdr_bh',out.read_text().splitlines()[0])

    def test_query_outside_universe_fails(self):
        with tempfile.TemporaryDirectory() as d:
            genes=Path(d)/'genes.txt'; genes.write_text('not_in_universe\n')
            with self.assertRaisesRegex(ValueError,'not present'):
                run_ora(genes,ROOT/'examples/pathways/toy.gmt',Path(d)/'out.tsv',ROOT/'examples/pathways/universe.txt',min_size=2)

if __name__=='__main__': unittest.main()
