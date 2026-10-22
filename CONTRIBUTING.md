# Contributing

Contributions should add a reproducible, reviewable method—not only a prompt. Keep skills focused, include an actionable `SKILL.md`, document assumptions and failure modes, and prefer small deterministic fixtures. Do not commit human sequencing data, credentials, licensed reference bundles without redistribution rights, or large generated outputs.

## Development

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
python -m unittest discover -s tests -v
python scripts/validate_skills.py
bash -n scripts/run_salmon_pipeline.sh
```

## Pull request checklist

- [ ] Changed behavior has tests, including at least one failure case.
- [ ] CLI help, input schema, output columns, warnings, and dependency requirements are documented.
- [ ] Skill metadata uses a valid `SKILL.md` YAML frontmatter block and names match the directory.
- [ ] Outputs are deterministic or any random seed is recorded.
- [ ] Scientific assumptions, reference versions, and citations are documented.
- [ ] Toy examples are explicitly labelled synthetic.
- [ ] No user data, large generated data, or secrets are included.
