pixi run Rscript scripts/01-mrna_qc.R L0927_Brain TRUE
pixi run Rscript scripts/01-mrna_qc.R L126_Brain_s1 FALSE
pixi run Rscript scripts/01-mrna_qc.R L126_Brain_s2 FALSE
pixi run Rscript scripts/01-mrna_qc.R L126_Brain_s3 FALSE

pixi run Rscript scripts/02-cluster.R L0927_Brain
pixi run Rscript scripts/02-cluster.R L126_Brain_s1
pixi run Rscript scripts/02-cluster.R L126_Brain_s2
pixi run Rscript scripts/02-cluster.R L126_Brain_s3

pixi run Rscript scripts/03-consensus_clusters_domains.R
