import subprocess, sys, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
class SkillsTests(unittest.TestCase):
    def test_skill_validator(self):
        p=subprocess.run([sys.executable,str(ROOT/'scripts/validate_skills.py')],capture_output=True,text=True)
        self.assertEqual(p.returncode,0,p.stdout+p.stderr)
        self.assertIn('5 skill folders',p.stdout)
if __name__=='__main__': unittest.main()
