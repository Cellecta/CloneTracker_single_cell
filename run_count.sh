# /mnt/project/Pipeline/Softeware/cellranger-8.0.0/cellranger count --id=Control_GEX \
#                   --fastqs=/mnt/project/Cellecta_scCRISPR/103951_Revmed/data/Control_GEX_trimmed \
#                   --transcriptome=/mnt/project/Pipeline/Reference/10XGenomics/refdata-gex-GRCh38-2024-A \
# 				  --create-bam true \
#                   --include-introns false

# /mnt/project/Pipeline/Softeware/cellranger-8.0.0/cellranger count --id=drugA_GEX \
#                     --fastqs=/mnt/project/Cellecta_scCRISPR/103951_Revmed/data/drugA_GEX_trimmed \
#                   --transcriptome=/mnt/project/Pipeline/Reference/10XGenomics/refdata-gex-GRCh38-2024-A \
# 				  --create-bam true \
#                   --include-introns false

/mnt/project/Pipeline/Softeware/cellranger-8.0.0/cellranger count --id=drugB_GEX \
                  --fastqs=/mnt/project/Cellecta_scCRISPR/103951_Revmed/data/drugB_GEX_trimmed \
                  --transcriptome=/mnt/project/Pipeline/Reference/10XGenomics/refdata-gex-GRCh38-2024-A \
				  --create-bam true \
                  --include-introns false                  