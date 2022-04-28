# PanExplorer_workflow

## Install

1- Git clone

```
git clone https://github.com/SouthGreenPlatform/PanExplorer_workflow.git
```

2- Define the PANEX_PATH environnement variable

```
cp -rf PanExplorer_workflow /usr/local/bin
export PANEX_PATH=/usr/local/bin/PanExplorer_workflow
```

3- Retrieve COG database

```
wget https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/little_endian/Cog_LE.tar.gz
tar -xzvf Cog_LE.tar.gz
cp -rf Cog.* $PANEX_PATH/COG
```

## Run the workflow

```
snakemake --cores 1 -s $PANEX_PATH/Snakemake_files/Snakefile_wget_roary_heatmap_upset_COG
```
