# HKGI Flagship Paper 2025
This repository organized the key analysis scripts, and scripts for regenerating the figures of the manuscript. 

## Analysis

Each folder section may contain these subfolders:
- `script` - scripts and commands arranged in the order of execution, e.g. 01.sv_pipeline.sh, 02.summary.ipynb. Files without numbers as prefix are the reference data and helper scripts.
- `data` - data tables used by the scripts.
- `results` - per-sample outputs .
- `summary` - aggregated and consolidated tables.

### Key sections
1. `analysis/primary`
    1. `raw_data` - Raw VCF of small variants and SVs (manta and cnvkit), probands and family members of probands.
    1. `smallvar_vep` - Ensembl VEP-annotated VCF, probands and family members of probands.
    1. `01.sv` - Identifying pathogenic/likely pathogenic SVs for primary findings.
1. `analysis/secondary`
    1. `raw_data` - Raw VCF of small variants and SVs (manta and cnvkit), HKGP Chinese cohort.
    1. `smallvar_vep` - Ensembl VEP-annotated VCF, HKGP Chinese cohort.
    1. `dominant`
        1. `01.known_small_variants` - Identify known clinVar small variants from dominant disorder genes.
        1. `02.novel_small_variants` - Identify novel small variants from dominant disorder genes.
        1. `03.filtering` - Consolidating small variants from 01 and 02. Performing hard filtering on the variants.
        1. `04.sv` - Identifying pathogenic/likely pathogenic SVs in dominant disorder genes.
    1. `recessive`
        1. `01.known_small_variants` - Identify known clinVar small variants from recessive disorder genes.
        1. `02.novel_small_variants` - Identify novel small variants from recessive disorder genes.
        1. `03.filtering` - Consolidating small variants from 01 and 02. Performing hard filtering on the variants.
        1. `04a.SMN` - Genotype SMN1 copy number.
        1. `04b.HBA` - Genotype HBA1/2.
        1. `04c.sv` - Identifying pathogenic/likely pathogenic SVs in recessive disorder genes.
        1. `05.gene_level` - Summarizing carriers of recessive disorder variants in gene level.
    1. `PGx`
        1. `01a.genotype.aldy`, `01b.genotype.Cyrius`, `01c.genotype.HLA-HD`, `01d.genotype.manual` - PGx genotyping with different callers.
        1. `02.consolidate_genotypes` - Consolidating results from the above PGx callers.
        1. `03.eval_prescription` - Calculate statistics for secondary - PGx in unrelated Chinese.


## Figures

### Steps to reproduce figures

1. Build an apptainer image containing Jupyter notebook with all R and python dependencies installed: `./figures/build/build.sh`
2. Run the created apptainer image. This will start a Jupyter notebook: `apptainer run r_hkgi_flagship_figures-v0.0.1.sif`
3. Open the http://localhost link that is printed to the terminal to connect to the Jupyter notebook.
4. Within Jupyter notebook, run the notebooks under `figures/src`
