python scripts/run_full_pipeline.py \
  --samples-csv examples/pipeline_samples_test.csv \
  --pipeline-root /mnt/d/Project/Clonetracker/103951/tests \
  --cellranger-bin /mnt/d/Project/Software/cellranger-10.0.0/cellranger \
  --transcriptome /mnt/d/Project/Clonetracker/103951/data/refdata-gex-GRCh38-2024-A \
  --bc14-file /mnt/d/Project/Clonetracker/103951/data/Cellecta-CloneTrackerXP-50M-BC14-LNGS-300-Library-Design.txt \
  --bc30-file /mnt/d/Project/Clonetracker/103951/data/Cellecta-CloneTrackerXP-50M-BC30-LNGS-300-Library-Design.txt \
  --barcode-search-umi-cutoff 0 \
  --assignment-min-total-umi 3 \
  --assignment-min-top-umi 3 
