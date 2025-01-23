# `cenplot`
Library for producing centromere figures.

> WIP

### Test
`rm`
* /project/logsdon_shared/projects/HGSVC3/CenMAP_verkko_run6/results/hor_stv/repeats/all_cens_chrX.annotation.fa.out
`stv`
* /project/logsdon_shared/projects/HGSVC3/CenMAP_verkko_run6/results/hor_stv/bed/chrX_AS-HOR_stv_row.all.bed
`cdr`
* /project/logsdon_shared/projects/HGSVC3/CenMAP_verkko_run6/results/hor_stv/bed/chrX_cdrs.bed
`self-ident`
* /project/logsdon_shared/projects/HGSVC3/CenMAP_verkko_run6/results/moddotplot/original/all_chrX.bed

```bash
python scriptsplot_cens_all_stv_mpl.py \
-i /project/logsdon_shared/projects/HGSVC3/CenMAP_verkko_run6/results/hor_stv/bed/chrX_AS-HOR_stv_row.all.bed \
-t test/tracks.csv
```

```bash
Rscript /project/logsdon_shared/projects/HGSVC3/CenMAP_verkko_run6/workflow/scripts/plot_cens_moddotplot.R         --bed chrY_ident.bed         --hor chrY_stv.bed         --hor_ort chrY_stv_ort.bed         --sat chrY_sat_annot.bed         --cdr chrY_cdrs.bed         --methyl chrY_methyl.bed
```