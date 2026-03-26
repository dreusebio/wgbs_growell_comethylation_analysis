# Mini cytosine report dataset
## Generate sample_ids.txt

Run the following commands to create the list of 49 sample IDs from the Baseline cytosine reports.

```bash
mkdir -p data/demo/comethyl_minidata

ls data/processed/08_cytosine_reports/Baseline/*CpG_report.txt.gz \
| sort \
| sed 's#.*/##' \
| cut -d'_' -f1 \
| head -49 \
> data/demo/comethyl_minidata/sample_ids.txt
```
This directory was generated for lightweight local testing of comethyl.

- Number of requested samples: 43
- Number of matched samples: 43
- Mode: filtered
- Chromosomes kept: chr22
- BED used: none
- Max lines per file: none
