# Python-qiime-microbiome-analysis
16S bacterial amplicon sequencing analysis workflow on QIIME

source activate qiime1
cd ~/yourworkingdirectory

0. Raw data inspection

```{r eval = FALSE}
#Inspect read count and read length with QIIME
mkdir 0_raw_fastq
count_seqs.py -i 0_raw_fastq/"*.fastq" -o count_seqs_raw.txt
```
1. Stitching together the forward and reverse reads with QIIME 

```{r eval = FALSE}
multiple_join_paired_ends.py -i 0_raw_fastq -o 1_paired --read1_indicator '_R1' --read2_indicator '_R2'
count_seqs.py -i 1_paired_renamed/"*.fastq" -o count_seqs_paired.txt
```
2. Primer removal 

```{r eval = FALSE}
extract_barcodes.py -f 1_paired_renamed/Tw-22_paired.fastq -o 2_debarcoded_paired/Tw-22 -c barcode_paired_stitched -l 17 -L 20
```
3. QC check with Trimmomatic (Java)
```{r eval = FALSE}
java -jar Trimmomatic-0.36/trimmomatic-0.36.jar SE 2_debarcoded_paired/Tw-22/reads.fastq 3_trimmed_debarcoded_paired/Tw-22_trm_dbc_paired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
4. 2nd QC trimming and concatenating all reads into a single file 
```{r eval = FALSE}
split_libraries_fastq.py -i 3_trimmed_debarcoded_paired/Tw-22_trm_dbc_paired.fastq,3_trimmed_debarcoded_paired/Tw-23_trm_dbc_paired.fastq,3_trimmed_debarcoded_paired/Tw-24_trm_dbc_paired.fastq,3_trimmed_debarcoded_paired.fastq --sample_ids Tw-22,Tw-23,Tw-24 -o 4_slout_all -m map_mg3.txt -q 19 --barcode_type 'not-barcoded'
```

5. OUT picking, taxonomic assignment, and phylogenetic analysis
```{r eval = FALSE}
#OUT picking, taxonomic assignment (using uclust), and phylogenetic analysis (using PyNAST) using the workflow script pick_open_reference_otus.py
pick_open_reference_otus.py -i 4_slout_all/seqs.fna -o 5_uclust_all -s 0.1 -a -O 5

# Taxa assignment (using RDP classifier)
assign_taxonomy.py -i 5_uclust_all/rep_set.fna -o 5_uclust_all/rdp_assigned_taxonomy -m rdp

# Chimera-free phylogenetic analysis (using PyNAST)

# 1. Identification of chimeras using script identify_chimeric_seqs.py 
parallel_identify_chimeric_seqs.py -i 5_uclust_all/pynast_aligned_seqs/rep_set_aligned.fasta -o 5_uclust_all/pynast_aligned_seqs/chimeric_seqs.txt -m ChimeraSlayer -a /macqiime/greengenes/core_set_aligned.fasta.imputed -O 3

# 2. Remove chimeric sequences using script filter_fasta.py
filter_fasta.py -f 5_uclust_all/pynast_aligned_seqs/rep_set_aligned.fasta -o 5_uclust_all/pynast_aligned_seqs/rep_set_aligned_no_chmc.fasta -s 5_uclust_all/pynast_aligned_seqs/chimeric_seqs.txt -n

# 3. Filtering alignment using script filter_alignment.py
filter_alignment.py -i 5_uclust_all/pynast_aligned_seqs/rep_set_aligned_no_chmc.fasta -o 5_uclust_all/pynast_aligned_seqs/rep_set_aligned_no_chmc_pfiltered.fasta

# 4. Building a chimera-free phylogenetic tree![image](https://user-images.githubusercontent.com/7209275/172734581-838cb923-8560-4105-987b-fa1d2814ef3e.png)
make_phylogeny.py -i 5_uclust_all/pynast_aligned_seqs/rep_set_aligned_no_chmc_pfiltered.fasta/rep_set_aligned_no_chmc_pfiltered.fasta -o 5_uclust_all/rep_set_no_chmc.tre

# Re-generate an OTU table based on the RDP-taxa and chimera-free tree

# 1. Re-generate a RDP-taxa, chimera-free OTU table
make_otu_table.py -i 5_uclust_all/final_otu_map_mc2.txt -o 5_uclust_all/otu_table_w_rdp_no_chmc.biom -t 5_uclust_all/rdp_assigned_taxonomy/rep_set_tax_assignments.txt -e 5_uclust_all/pynast_aligned_seqs/chimeric_seqs.txt

# 2. Check OTU table stats (before filtering)
mkdir 6_biom_summary_all | biom summarize-table -i 5_uclust_all/otu_table_w_rdp_no_chmc.biom -o 6_biom_summary_all/biom_summary_w_rdp_no_chmc.txt

# 3. Filter the RDP-taxa, chimera-free OTU table
filter_otus_from_otu_table.py -i 5_uclust_all/otu_table_w_rdp_no_chmc.biom -o 5_uclust_all/otu_table_w_rdp_no_chmc_filtered.biom --min_count_fraction 0.00005

# 4. Check OUT table stats (after filtering)
##uclust_chmc 
biom summarize-table -i 5_uclust_all/otu_table_mc2_w_tax_no_pynast_failures.biom -o 6_biom_summary_all/biom_summary_mc2_w_tax_no_pynast_failures.txt 
##rdp_no_chmc
biom summarize-table -i 5_uclust_all/otu_table_w_rdp_no_chmc_filtered.biom -o 6_biom_summary_all/biom_summary_w_rdp_no_chmc_filtered.txt
```

6. Taxa summaries (L1 to L6)

```{r eval = FALSE}
#uclust_chmc 
summarize_taxa.py -i 5_uclust_all/otu_table_mc2_w_tax_no_pynast_failures.biom -o 7_taxa_summaries_all/mc2_w_tax_no_pynast_failures/

#rdp_no_chmc
summarize_taxa.py -i 5_uclust_all/otu_table_w_rdp_no_chmc_filtered.biom -o 7_taxa_summaries_all/w_rdp_no_chmc_filtered/

summarize_taxa_through_plots.py -i 5_uclust_all/otu_table_w_rdp_no_chmc_filtered.biom -m map_mg3.txt -o 9_taxa_plots_all/rdp_no_chmc/ -s

```

7. Diversity analyses
Alpha diversity
```{r eval = FALSE}
#Multiple subsamplings/rarefyings on the input OTU tables
multiple_rarefactions.py -i 5_uclust_all/otu_table_w_rdp_no_chmc_filtered.biom -m 10000 -x 340000 -s 5000 -n 10 -o 10_alpha_diversity_all/subsamplings

#Calculate alpha diversity on each rarefied OTU table
alpha_diversity.py -i 10_alpha_diversity_all/subsamplings/ -o 10_alpha_diversity_all/alpha_rare/ -t 5_uclust_all/rep_set_no_chmc.tre -m observed_species,chao1,PD_whole_tree,shannon,simpson_e,gini_index

#Summarize the alpha diversity data
collate_alpha.py -i 10_alpha_diversity_all/alpha_rare/ -o 10_alpha_diversity_all/alpha_collated/
```

Beta diversity
```{r eval = FALSE}
jackknifed_beta_diversity.py -i 5_uclust_all/otu_table_w_rdp_no_chmc_filtered.biom -o 11_beta_diversity_all -e 9000 -m map_mg3.txt -t 5_uclust_all/rep_set_no_chmc.tre --master_tree consensus
```

8. Network analysis
```{r eval = FALSE}
make_otu_network.py -i 5_uclust_all/otu_table_w_rdp_no_chmc_filtered.biom -m map_mg3_network_all.txt -o 13_otu_network_all/

collapse_samples.py -b 5_uclust_all/otu_table_w_rdp_no_chmc_filtered_CUTforCyto.biom -m map_mg3_collapse.txt --output_biom_fp 12_collapse/collapsed.biom --output_mapping_fp 12_collapse/mg_collapsed_map.txt --collapse_mode sum --collapse_fields TimePoint
```
































