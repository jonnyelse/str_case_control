# str_case_control

**str_case_control** is a DNAnexus applet for the UK Biobank Research Analysis Platform (RAP) designed to perform a case-control Short Tandem Repeat (STR) Genome-Wide Association Study (GWAS) using joint analysis. It automates the process of extracting and modeling STR genotype data for binary phenotypes, leveraging both bash and R for efficient, large-scale analysis on biobank data.

---

## Features

- Efficient STR GWAS for case-control cohorts
- Automated data handling: VCF parsing, multi-allelic STR support, and regression modeling
- Designed for DNAnexus execution on UK Biobank RAP
- Covariate-adjusted logistic regression for association analysis

---

## Inputs

| Name             | Type   | Required | Description                                      |
|------------------|--------|----------|--------------------------------------------------|
| `vcf_file`       | file   | Yes      | Input VCF file containing joint STR genotypes    |
| `vcf_file_index` | file   | Yes      | Index for the VCF file (e.g., CSI file)          |
| `chr`            | string | Yes      | Chromosome to analyze (e.g., "chr1")             |
| `min_ac`         | int    | Yes      | Minimum allele count threshold for grouping       |
| `covar_pheno`    | file   | Yes      | Phenotype and covariate table (TSV, incl. IID, FID, disease_status, covariates) |
| `output_folder`  | string | No       | Output folder for results                        |
| `project`        | string | No       | DNAnexus project (optional override)             |

---

## Outputs

| Name            | Type | Description                                       |
|-----------------|------|---------------------------------------------------|
| `merged_results`| file | Final TSV with combined regression results         |

---

## Dependencies

This applet requires the following tools and R packages (managed within the DNAnexus environment):

- `bcftools`
- `tabix`
- `r-base`
- R packages: `dplyr`, `purrr`, `stringr`, `tidyr`, `data.table`, `reshape2`, `brglm2`

---

## Usage

1. **Prepare your inputs:**
    - Joint STR VCF and index
    - Covariate/phenotype table (TSV with columns: FID, IID, disease_status, covariates)
    - Specify chromosome and minimum allele count threshold

2. **Run the applet on DNAnexus:**
    - Upload your files to DNAnexus and launch the applet via the platform or CLI.

    Example DNAnexus CLI usage:
    ```bash
    dx run str_case_control \
      -i vcf_file=your_genotypes.vcf.gz \
      -i vcf_file_index=your_genotypes.vcf.gz.csi \
      -i chr=chr1 \
      -i min_ac=5 \
      -i covar_pheno=phenotypes.tsv \
      -i output_folder=/project/results/
    ```

3. **Output:**
    - `merged_results` (TSV): Combined regression results for all analyzed STRs, including odds ratios, p-values, and covariate effects.

---

## Output File Format

The output TSV (`merged_results`) contains:
- `locus` (variant identifier)
- `Alleles` (tested alleles)
- `OR` (odds ratios)
- `Std_Errors`
- `Z_Values`
- `P_Values`
- `Smallest_P_Value`
- `N_Samples`

---

## Notes

- Multi-allelic STRs are handled via allele grouping and merging where necessary.
- **Low allele count variants are grouped together, not filtered out.** Only monomorphic variants (variants with no observed polymorphism in the data) are filtered out of the association analysis.

---

## Citation

If you use this applet, please cite the UK Biobank and any relevant STR GWAS or software references.

---

## Author

- [jonnyelse](https://github.com/jonnyelse)

---

## License

Add your license here (e.g., MIT, GPL-3.0, etc.)

---

## Contributing

Contributions and suggestions are welcome; please open an issue or submit a pull request.
