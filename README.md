# ABOUT

### GOAL(S)

Takes the Single cell matrix outputs from running scNanoGPS on long read transcript sequences and convert them into a format compatable with Seurat. Then analyze those files via the Seurat workflows included here.

### DATA

`/vf/users/CARDPB/data/snRNA_longread/scNanoGPS-neuro/scNanoGPS_res`

### OPERATIONS

- CONVERSIONS - of the scNanoGPS_res matrix outputs
  - matrix.tsv + matrix_isoform.tsv --> barcodes.tsv, genes.tsv, matrix.mtx (for each file)
  - matrix.tsv + matrix_isoform.tsv --> merged_matrix_and_matrix_isoform.tsv --> barcodes.tsv, genes.tsv, matrix.mtx
- SEURAT - analyses and plots
  - matrix.tsv
  - matrix_isoform.tsv
  - merged_matrix_and_matrix_isoform.tsv
  - short_reads
    - The rest of the file operations were intended for long reads so this short read script was for comparison for the same Sample ID # that had both long read and short read data available
  - Seurat QC, filtering, celltyping and plotting
- (WIP) HARMONY INTEGRATION
  - matrix.tsv OR matrix_isoform.tsv OR merged_matrix_and_matrix_isoform.tsv

# STEPS TO RUN

### PRE-REQUISITES

1) Run thru the scNanoGPS pipeline at:

```text
/vf/users/CARDPB/data/snRNA_longread/scNanoGPS-neuro
```

2) Use the output files that were created in the `/scNanoGPS_res` folder (the results after running scNanoGPS). The 2 specific files needed for this part are:

```text
/scNanoGPS_res/
   matrix.tsv
   matrix_isoform.tsv
```

### CONVERSIONS

3) Run thru the conversion script for your input data and desired output:


| SBATCH                                                                                  | CALLS THIS SCRIPT                                                                | INPUT(S)                            | OUTPUT(S)                                                                                                |
| ----------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------- | ------------------------------------- | :--------------------------------------------------------------------------------------------------------- |
| CONVERSION-scNanoGPS_matrix_and_matrix_isoform_to_barcodes_genes_mtx-argparse-sbatch.sh | CONVERSION-scNanoGPS_matrix_and_matrix_isoform_to_barcodes_genes_mtx-argparse.py | matrix.tsv,<br />matrix_isoform.tsv | (For each input file:)<br />barcodes.tsv, <br />genes.tsv, <br />matrix.mtx, <br />id.txt                |
| CONVERSION-scNanoGPS-merged_matrix_and_matrix_isoform-argparse-sbatch.sh                 | CONVERSION-scNanoGPS-merged_matrix_and_matrix_isoform-argparse.py                 | matrix.tsv,<br />matrix_isoform.tsv | merged_matrix_and_matrix_isoform.tsv,<br />barcodes.tsv, <br />genes.tsv, <br />matrix.mtx, <br />id.txt |

### SEURAT

4) Using matrix/matrix_isoform/merged_matrix_and_matrix_isoform files, run the associated script to peform Seurat QC, filtering, celltyping and plotting. Individual scripts can be run separately for testing (see SCRIPT column). Alternatively, the SBATCH script (see SBATCH column) can run the `-argpase.r` scripts for:

- matrix.tsv AND matrix_isoform.tsv --> uses the same inputs and can run BOTH (or separately), putting the `/isoform` outputs together with the same Sample ID #'s output
- OR
- merged_matrix_and_matrix_isoform.tsv

  Just make sure to comment/uncomment the relevant sections (matrix and matrix_isoform sections OR the merged section)


| SBATCH                                                        | CALLS THIS SCRIPT                                                                                                            | INPUT(S)                                                                                                                                                | OUTPUT(S)                                                                                             |
| --------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------- | :------------------------------------------------------------------------------------------------------ |
| SEURAT-for_scNanoGPS-matrix+isoform+merged-argparse-sbatch.sh | SEURAT-for_scNanoGPS-matrix-argparse.R,<br />SEURAT-for_scNanoGPS-matrix.R                                                   | matrix.tsv,<br />id.txt                                                                                                                                 | `/output/long_reads/...` <br />log-{date_time}.txt, <br />(PLOTS), <br />{Seurat_object}.rds          |
| SEURAT-for_scNanoGPS-matrix+isoform+merged-argparse-sbatch.sh | SEURAT-for_scNanoGPS-matrix_isoform-argparse.R,<br />SEURAT-for_scNanoGPS-matrix_isoform.R                                   | matrix_isoform.tsv,<br />id.txt                                                                                                                         | `/output/long_reads/...` <br />log-{date_time}.txt, <br />(PLOTS), <br />{Seurat_object}.rds          |
| SEURAT-for_scNanoGPS-matrix+isoform+merged-argparse-sbatch.sh | SEURAT-for_scNanoGPS-merged-argparse.R,<br />SEURAT-for_scNanoGPS-merged.R                                                   | merged_matrix_and_matrix_isoform.tsv,<br />id.txt                                                                                                       | `/output/merged/...` <br />log-{date_time}.txt, <br />(PLOTS), <br />{Seurat_object}.rds              |
| (N/A)                                                         | SEURAT-for_short_reads.R                                                                                                     | (Arc Cellranger output folder:)<br />filtered_feature_bc_matrix.h5 <br />(FULL celltype marker list) <br />(Folder of individual celltype marker lists) | `/output/short_reads/...` <br />(PLOTS)                                                               |
| (N/A)                                                         | HARMONY-Integrate_ngn2_0624_res2-WIP-matrix.R,<br />HARMONY-Integrate_ngn2_0624_res2-WIP-merged_matrix_and_matrix_isoforms.R | (multiple:)<br />matrix.tsv, <br />matrix_isoform.tsv, <br />merged_matrix_and_matrix_isoform.tsv                                                       | `/output/merged/Harmony_integration/...` <br />log-{date_time}.txt, <br />{Seurat_object}.rds, <br /> |

Plot outputs will be saved to the `/output` folder, then subdivided by which script and input data you used (`/long_reads`, `/merged`, `/short_reads`) and then by Sample ID #, date, plots.

General output folder structure:

```text
/output
   /{input_data_type: long_reads, merged, short_reads}
      /{sample_ID_#}
         /{date}
            /{plots}

```

### (WIP) HARMONY (optional)

5) The goal here was to try different ways of merging the data:

- all matrices together
- all matrix_isoforms together
- all matrices + all matrix_isoforms together

and then running thru Seurat. The results are slightly better than running the matrix_isoforms individually. However, it is probably more logically sound to merge each sample's matrix + matrix_isoform instead, which is why I created the `merged_matrix_and_matrix_isoform` files.

### ABOUT THE FILE NAMING

- Files ending in `-sbatch.sh` files are the sbatch script files to submit to Biowulf
- The sbatch file runs the file of the SAME name but without the `sbatch.sh` ending, it's replaced with `-argparse` and the extension (`.py` or `.R`) instead
- The same name but ending in `.ipynb` are the initial Jupyter Notebook files to test/run in `sinteractive` mode (keep in mind: these don't have the `argparse` commands)

# SOFTWARE VERSIONS


| R PACKAGES   | VERSION |
| -------------- | --------- |
| R            | 4.4.2   |
| ggplot       | 2_3.5.1 |
| patchwork    | 1.2.0   |
| Seurat       | 5.1.0   |
| SeuratObject | 5.0.2   |
| sp           | 2.1-4   |
| dplyr        | 1.1.4   |
| argparse     | 2.2.3   |

| PYTHON LIBRARIES |
| ------------------ |
| pandas           |
| os               |
| numpy            |
| re               |
| datetime         |
| argparse         |
| scipy.sparse     |
| scipy.io.mmwrite |
