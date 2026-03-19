python clonetracker_batch.py \
  --samples_csv barcode_sample.csv \
  --out_root /mnt/project/Cellecta_scCRISPR/103951_Revmed/CloneTracker_out_combined \
  --bc14_file /mnt/project/Cellecta_scCRISPR/103951_Revmed/Cellecta-CloneTrackerXP-50M-BC14-LNGS-300-Library-Design.txt \
  --bc30_file /mnt/project/Cellecta_scCRISPR/103951_Revmed/Cellecta-CloneTrackerXP-50M-BC30-LNGS-300-Library-Design.txt \
  --best_sequence_umi_py /mnt/project/Pipeline/CloneTracker_single_cell/best_sequence_umi.py \
  --barcode_process_py /mnt/project/Pipeline/CloneTracker_single_cell/barcode_process_umi_5.py \
  --final_assignment_py /mnt/project/Pipeline/CloneTracker_single_cell/barcode_final_assignment_10_Genomics_16.py
