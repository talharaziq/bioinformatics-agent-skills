import tempfile, unittest
from pathlib import Path
from unittest.mock import patch
from bioskills.bulk import run_salmon

ROOT=Path(__file__).resolve().parents[1]

class BulkOrchestratorTests(unittest.TestCase):
    def test_builds_commands_and_run_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp=Path(tmp); index=tmp/'index'; index.mkdir(); out=tmp/'out'
            def fake_run(cmd, stdout=None, stderr=None, text=None):
                if cmd[0]=='salmon':
                    (Path(cmd[cmd.index('-o')+1])).mkdir(parents=True)
                if stdout: stdout.write('mock tool output\n')
                class Result: returncode=0
                return Result()
            with patch('bioskills.bulk.shutil.which',return_value='/mock/bin/tool'), patch('bioskills.bulk.subprocess.run',side_effect=fake_run) as mocked:
                result=run_salmon(ROOT/'examples/expression/manifest.tsv',index,out,threads=2)
            self.assertEqual(result['samples'],4)
            self.assertEqual(sum(call.args[0][0]=='salmon' for call in mocked.call_args_list),4)
            self.assertEqual(sum(call.args[0][0]=='fastp' for call in mocked.call_args_list),4)
            self.assertEqual(sum(call.args[0][0]=='fastqc' for call in mocked.call_args_list),4)
            self.assertTrue((out/'run_manifest.tsv').is_file())
            self.assertIn('batch', (out/'run_manifest.tsv').read_text().splitlines()[0])
            self.assertTrue((out/'logs'/'ctrl1.salmon.log').is_file())

    def test_existing_sample_quant_stops_before_any_commands(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp=Path(tmp); index=tmp/'index'; index.mkdir(); out=tmp/'out'; (out/'salmon'/'treat2').mkdir(parents=True)
            with patch('bioskills.bulk.shutil.which',return_value='/mock/bin/tool'), patch('bioskills.bulk.subprocess.run') as mocked:
                with self.assertRaisesRegex(FileExistsError,'treat2'):
                    run_salmon(ROOT/'examples/expression/manifest.tsv',index,out)
            mocked.assert_not_called()

if __name__=='__main__': unittest.main()
