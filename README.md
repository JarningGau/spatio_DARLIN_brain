## Installation

```{bash}
git clone git@github.com:JarningGau/spatio_DARLIN_brain.git
cd spatio_DARLIN_brain
pixi install --all
pixi run init
```

## Run
```{bash}
bash run.sh
```

## Results

mRNA outputs
1. Latent embeddings -> Clone2vec
2. Cell type annotation
3. Spatial domain annotation -> Interpret the Clone2vec results

lineage outputs
1. Basis QC
2. CA/TA/RA barcodes -> clone ID
3. Clone2vec results

advanced analysis outputs
1. Cortex clonal pattern
2. Hippocampus clonal pattern
3. Hypothalamus clonal pattern
