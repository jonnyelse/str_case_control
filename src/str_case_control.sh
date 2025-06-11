#!/bin/bash
# str_case_control 0.0.1
# DNAnexus applet entry point script.
# This script runs a case-control STR GWAS on the specified chromosome using joint analysis.
# Inputs are loaded as environment variables.

main() {
    # Create temporary and output directories
    tmp_dir="tmp"
    mkdir -p $tmp_dir
    mkdir -p "QC"
    mkdir -p "results/chunks"

    # Download required input files from DNAnexus project
    dx download "$vcf_file" -o input_file
    dx download "$vcf_file_index" -o input_file.csi
    dx download "$covar_pheno" -o "covar_pheno.tsv"

    # Extract the specified chromosome from the VCF and trim alternate alleles
    bcftools view --trim-alt-alleles input_file -r $chr --threads 2 -Oz -o QC/${chr}_trimmed.vcf.gz

    # Generate a list of sample IDs from the VCF
    bcftools query -l input_file > samples_list.txt

    # Export genotypes for all variants in the QC'ed VCF into a TSV for processing
    bcftools query -f '%VARID\t%ALT[\t%GT]\n' QC/${chr}_trimmed.vcf.gz > QC/full_QC.genotypes.tsv

    # Process the genotype file:
    #   1. Remove symbolic STR markers from ALT allele column
    #   2. Replace missing genotypes (./.) with NA/NA
    awk -F'\t' -v OFS='\t' '{
        gsub(/<STR/, "", $2);  # Remove <STR prefix from ALT
        gsub(/>/, "", $2);     # Remove > suffix from ALT
        for (i = 3; i <= NF; i++) {
            if ($i == "./.") $i = "NA/NA";
        }
        print;
    }' QC/full_QC.genotypes.tsv > QC/processed_genotypes_corrected.tsv

    # For each variant, count alleles for each sample and output as TSV
    awk 'BEGIN {FS="\t"; OFS="\t"}
    {
      n = split($2, alleles, ",");
      printf $1 OFS $2;
      for (i = 3; i <= NF; i++) {
        if ($i == "NA/NA") {
            na_out = "NA";
            for (j = 2; j <= n; j++) {
                na_out = na_out ",NA";
            }
            printf OFS na_out;
            continue;
        } else {
            split($i, g, "/");
            delete counts;
            for (j = 1; j <= n; j++) counts[j] = 0;
            if (g[1] > 0) counts[g[1]]++;
            if (g[2] > 0) counts[g[2]]++;
            printf OFS counts[1];
            for (k = 2; k <= n; k++) printf "," counts[k];
        }
      }
      printf "\n";
    }' QC/processed_genotypes_corrected.tsv > QC/allele_counts_genotypes.tsv

    # Identify polymorphic variants and split them by allele count for parallel processing
    awk 'BEGIN { FS="\t"}
    {
      n = split($2, alleles, ",")
      isMonomorphic = 1
      if (n > 1) { isMonomorphic = 0 }
      else {
        for (i = 3; i <= NF && isMonomorphic; i++) {
          split($i, counts, ",")
          for (j in counts) {
            if (counts[j] != "0") {
              isMonomorphic = 0
              break
            }
          }
          if (isMonomorphic == 0) break
        }
      }
      if (!isMonomorphic) {
        fileName = "QC/" n "_alt_alleles.txt"
        print $0 > fileName
      }
    }' QC/allele_counts_genotypes.tsv

    # Run the association analysis in R (logistic regression & allele grouping)
    Rscript assoc.R $min_ac


    # Upload final results to DNAnexus and register as output
    final_file="results/final_combined_regression_results.tsv"
    filename_after_upload=${chr}_regression_results.tsv
    output="${output_folder}/${filename_after_upload}"
    results_file=$(dx upload $final_file --brief --path "$output" -p)
    dx-jobutil-add-output merged_results "$results_file" --class=file
}
