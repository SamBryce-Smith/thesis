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


### Commands to calculate consistency metrics from sample tables

Calculate consistency metrics for ground truth data based on the generated sample tables. Currently just jaccard index of pairwise comparisons

1. Make sure current directory is `<your_machine>/thesis/apa_benchmarking`
2. Run the following commands. Note that the number of pair-wise comparisons is quite high (122 * 4 window sizes) and the code is not optimised, so can take a few mins to run

```bash
mkdir -p processed/gt_consistency
python scripts/gt_replicate_consistency_metrics.py data/apaeval_ground_truths/chr_all.sample2path.tsv data/apaeval_ground_truths/chr_all.pairwise_comparisons.tsv processed/gt_consistency/chr_all --window-sizes 0 10 25 50
python scripts/gt_replicate_consistency_metrics.py data/apaeval_ground_truths/chr_te.sample2path.tsv data/apaeval_ground_truths/chr_te.pairwise_comparisons.tsv processed/gt_consistency/chr_te --window-sizes 0 10 25 50
```
