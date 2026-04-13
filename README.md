# Claude-BIG40-Glioma-PheWAS


```bash

# Default: downloads to data/big40/
bash 01_download_big40.sh

# Custom directory:
bash 01_download_big40.sh /scratch/big40_data

# Kicked off the network? Just rerun:
bash 01_download_big40.sh   # picks up where it left off


bash 01_download_big40.sh /destination/path

```


```
data/big40/
├── metadata/
│   ├── IDPs.html              # Phenotype descriptions
│   ├── snp_stats33k.txt.gz    # Variant annotation
│   └── idp_numbers.txt        # Auto-generated IDP list
├── stats33k/
│   ├── 0.txt.gz               # Per-IDP GWAS summary stats
│   ├── 1.txt.gz
│   ├── ...
│   └── 3934.txt.gz
└── download.log               # Timestamped log
```


```bash
bash 02_download_ebi.sh /destination/path

```



```bash
03_validate_ebi.sh /destination/path

```



