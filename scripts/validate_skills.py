#!/usr/bin/env python3
"""Validate Agent Skill folders' required YAML frontmatter and naming."""
from pathlib import Path
import re, sys
ROOT=Path(__file__).resolve().parents[1]
NAME=re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")
def main():
    skills_dir=ROOT/"skills"; errors=[]; checked=0
    for skill in sorted(p for p in skills_dir.iterdir() if p.is_dir()):
        f=skill/"SKILL.md"
        if not f.is_file(): errors.append(f"{skill.name}: missing SKILL.md"); continue
        text=f.read_text(encoding="utf-8")
        if not text.startswith("---\n"): errors.append(f"{skill.name}: file must start with YAML frontmatter"); continue
        end=text.find("\n---\n",4)
        if end<0: errors.append(f"{skill.name}: missing closing frontmatter delimiter"); continue
        front=text[4:end]
        metadata={}
        for line in front.splitlines():
            if ":" not in line: continue
            key,value=line.split(":",1); metadata[key.strip()]=value.strip().strip("\"'")
        name=metadata.get("name",""); desc=metadata.get("description","")
        if not name: errors.append(f"{skill.name}: missing name")
        elif not NAME.fullmatch(name): errors.append(f"{skill.name}: invalid skill name {name!r}")
        elif name!=skill.name: errors.append(f"{skill.name}: metadata name must match folder name (got {name!r})")
        if not desc: errors.append(f"{skill.name}: missing description")
        elif len(desc)>1024: errors.append(f"{skill.name}: description exceeds 1024 characters")
        if "# " not in text[end+5:]: errors.append(f"{skill.name}: missing Markdown title after frontmatter")
        if len(text.splitlines())<20: errors.append(f"{skill.name}: instructions look incomplete; provide actionable workflow details")
        checked+=1
    if errors:
        print("Skill validation FAILED:",file=sys.stderr)
        for err in errors: print(f"- {err}",file=sys.stderr)
        return 1
    print(f"Skill validation passed: {checked} skill folders")
    return 0
if __name__=="__main__": raise SystemExit(main())
