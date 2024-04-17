# APAeval chapter

## Dummy

### Commands to generate ground truth BED sample tables

Collects ground-truth BEDs downloaded from APAeval manuscript Zenodo repository into sample tables (mapping names and groups to filepaths & within-group pairwise comparison sample tables)

1. Make sure current directory is `<your_machine>/thesis/apa_benchmarking`
2. Run `scripts/beds_dir_to_samplesheet.py` with the following commands:

```bash
python scripts/beds_dir_to_samplesheet.py data/apaeval_ground_truths/chr_all_ground_truth_beds/ data/apaeval_ground_truths/chr_all
python scripts/beds_dir_to_samplesheet.py data/apaeval_ground_truths/chr_te_ground_truth_beds/ data/apaeval_ground_truths/chr_te
```
