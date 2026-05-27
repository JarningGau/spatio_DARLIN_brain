pixi run Rscript scripts_mrna/01-mrna_qc.R L0927_Brain TRUE
pixi run Rscript scripts_mrna/01-mrna_qc.R L126_Brain_s1 FALSE
pixi run Rscript scripts_mrna/01-mrna_qc.R L126_Brain_s2 FALSE
pixi run Rscript scripts_mrna/01-mrna_qc.R L126_Brain_s3 FALSE

pixi run Rscript scripts_mrna/02-cluster.R L0927_Brain
pixi run Rscript scripts_mrna/02-cluster.R L126_Brain_s1
pixi run Rscript scripts_mrna/02-cluster.R L126_Brain_s2
pixi run Rscript scripts_mrna/02-cluster.R L126_Brain_s3

pixi run Rscript scripts_mrna/03-consensus_clusters_domains.R
pixi run Rscript scripts_mrna/04-latent_embeddings.R

pixi run Rscript scripts_lineage/01-preprocess_lineage.R
pixi run Rscript scripts_lineage/02-QC_metrics.R
pixi run Rscript scripts_lineage/03-QC_plot.R