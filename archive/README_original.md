# Claude-BIG40-Glioma-PheWAS



```


```bash/

mkdir -p /francislab/data1/refs/BIG40

nohup bash 01_download_big40.sh /destination/path 8 > 01_download_big40.log 2>&1 &

```



We've implemented the bandwidth limiter in a more responsive way, you should get better performance now. Additional tip, use https://open.oxcin.ox.ac.uk rather than https://open.win.ox.ac.uk to avoid a redirect step.





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

nohup 02_download_ebi.sh /destination/path > 02_download_ebi.log 2>&1 &

```



```bash

03_validate_ebi.sh /destination/path

```



```bash

bash 04_validate_stats33k.sh /destination/path

```


```bash

bash 05_install_prscs.sh ~/.local/

```


```bash

bash 06_format_for_prscs.sh /destination/path 16

```


Prep input directory to meet expectations ...

```bash
basedir=/francislab/data1/working/20250800-AGS-CIDR-ONCO-I370-TCGA
mkdir -p input
cd input

for ds in onco i370 tcga ; do
ln -s ${basedir}/20251218-survival_gwas/imputed-umich-${ds}
done

ds=cidr
ln -s /francislab/data1/working/20250813-CIDR/20260320f-impute_genotypes/imputed-umich-${ds}
cd ..
```



```bash

07_extract_bim.sh ${PWD}/input /destination/path

```





The hdf5 files are in ~/.local/ld_ref/ldblk_1kg_eur/



```bash

BIM_DIR=/francislab/data1/refs/BIG40/target_bim
SNPINFO=~/.local/ld_ref/ldblk_1kg_eur/snpinfo_1kg_hm3
REMAP=~/github/ucsffrancislab/Claude-BIG40-Glioma-PheWAS/07c_remap_bim_rsids.py

for cohort in cidr i370 onco tcga; do
    python3 $REMAP \
        $BIM_DIR/imputed-umich-${cohort}.bim \
        $BIM_DIR/${cohort}.bim \
        $SNPINFO
    echo
done

```


patch ~/.local/PRScs

```bash

bash ~/github/ucsffrancislab/Claude-BIG40-Glioma-PheWAS/patch_prscs_numpy2.sh

```




08b_prep_target.sh


```bash

bash 08b_run_ct.sh ${PWD}/input ${PWD}/ct_output

```





```bash

bash 08_run_prscs.sh cidr ${PWD}/ct_output

```

This as it would take about 60 days of uninterupted compute time on our node.
We have opted to filter the 3935 traits to the top by literature and the top by our own C+T run which we are doing next.

/francislab/data1/refs/sources/fileserve.mrcieu.ac.uk/ld/




09_ct_association.R ${PWD}/ct_output ${PWD}/input 64"







08_run_prscs.sh ${PWD}/prscs_output "$cohort" ~/github/ucsffrancislab/Claude-BIG40-Glioma-PheWAS/prscs_idp_list.txt












