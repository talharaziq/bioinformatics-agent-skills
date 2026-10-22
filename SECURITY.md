# Security

The core tools operate on local files and do not evaluate input data as code. The optional `geo-fetch` command makes network requests to NCBI GEO for the exact accession supplied by the user. Never pass untrusted scripts to a shell; scripts here build external commands with argument arrays rather than `shell=True`.

Report suspected vulnerabilities privately to the repository maintainer rather than placing exploit details in a public issue. Include the affected version, reproduction steps, and expected impact. Do not submit biological or personal data in a report.
