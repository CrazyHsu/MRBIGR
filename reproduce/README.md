# unpack MRBIGR_data.tar.gz
```bash
tar -xf ../MRBIGR_data/MRBIGR_data.tar.gz -C ../MRBIGR_data
```

# clump to 381013 SNPs
```bash
python ../MRBIGR.py geno -clump -g ../MRBIGR_data/MRBIGR_data/chr_HAMP -o ../MRBIGR_data/MRBIGR_data/chr_HAMP_380k -r 0.7
```

# phylotree 2000k and 380k
## plot phylotree from bed file
```bash
python ../MRBIGR.py plot_new phylotree -bed ../MRBIGR_data/MRBIGR_data/chr_HAMP -o MRBIGR_phylo_2000k_from_bed -plot_fmt pdf -group ../MRBIGR_data/MRBIGR_data/491.group_info.csv -group_sep "," -drop ../MRBIGR_data/MRBIGR_data/drop_list
python ../MRBIGR.py plot_new phylotree -bed ../MRBIGR_data/MRBIGR_data/chr_HAMP_380k -o MRBIGR_phylo_510_380k_from_bed -plot_fmt pdf -group ../MRBIGR_data/MRBIGR_data/491.group_info.csv -group_sep "," -drop ../MRBIGR_data/MRBIGR_data/drop_list
```

## plot phylotree from nwk file
```bash
python ../MRBIGR.py plot_new phylotree -nwk ../MRBIGR_data/MRBIGR_data/AMP_510.2000k.tree.nwk -o MRBIGR_phylo_2000k -plot_fmt pdf -group ../MRBIGR_data/MRBIGR_data/491.group_info.csv -group_sep "," -drop ../MRBIGR_data/MRBIGR_data/drop_list
python ../MRBIGR.py plot_new phylotree -nwk ../MRBIGR_data/MRBIGR_data/AMP_510.380k.tree.nwk -o MRBIGR_phylo_380k -plot_fmt pdf -group ../MRBIGR_data/MRBIGR_data/491.group_info.csv -group_sep "," -drop ../MRBIGR_data/MRBIGR_data/drop_list
```

# PCA 2000k and 380k
## PCA for 2000k SNP
```bash
python ../MRBIGR.py geno -pca -g ../MRBIGR_data/MRBIGR_data/chr_HAMP -o AMP_510.2000k_pca
python ../MRBIGR.py plot_new scatterps -i AMP_510.2000k_pca_pca.csv -g ../MRBIGR_data/MRBIGR_data/AMP_537female_info.csv -ps_type pca -o AMP_510.2000k_pca_scatter
```
## PCA for 380k SNP
```bash
python ../MRBIGR.py geno -pca -g ../MRBIGR_data/MRBIGR_data/chr_HAMP_380k_clump -o AMP_510.380k_pca
python ../MRBIGR.py plot_new scatterps -i AMP_510.380k_pca_pca.csv -g ../MRBIGR_data/MRBIGR_data/AMP_537female_info.csv -ps_type pca -o AMP_510.380k_pca_scatter
```

# heatmap for 325 strains
```bash
csvtk cut -f 1 ../MRBIGR_data/MRBIGR_data/E3_log2.normalized_phe.csv | csvtk grep -P - ../MRBIGR_data/MRBIGR_data/AMP_537female_info.TEMP_TST.csv | csvtk rename -n "ID,Subpop" -f 1,2 > 325.group
python ../MRBIGR.py plot_new heatmap -i ../MRBIGR_data/MRBIGR_data/E3_log2.normalized_phe.csv -row_select 325.group -col_select ../MRBIGR_data/MRBIGR_data/102_metabolites.groups -anno_row 325.group -anno_col ../MRBIGR_data/MRBIGR_data/102_metabolites.groups -cluster_rows -cluster_cols -o MRBIGR_heatmap_metablite -height 6 -width 8
```

# PC1 distribution
```bash
python ../MRBIGR.py rd -i ../MRBIGR_data/MRBIGR_data/E3_log2.normalized_phe.csv -group_for_phe ../MRBIGR_data/MRBIGR_data/102_metabolites.groups -s_strains 325.group -input_sep "," -group_sep "," -m pca -group_min_n 5 -o MRBIGR_rd
csvtk gather -k metabolite -v PC1 -f -1 MRBIGR_rd.PCA.csv | csvtk join -f "ID" - 325.group > MRBIGR_rd.PCA.melt.csv
echo -e "Flavonoid\nAmino acid\nLysophosphatide\nPhenolamides\nTryptophan metabolite\nVitamin\nFatty Acids" > metabolites_order
python ../MRBIGR.py plot_new group -i MRBIGR_rd.PCA.melt.csv -group_field metabolite -value_field PC1 -fill_field Subpop -o MRBIGR_group.PC -height 6 -width 8 -group_order metabolites_order -add_pvalue
```

# GWAS
```
csvtk rename -f 1 -n "ID" ../MRBIGR_data/MRBIGR_data/AMP_kernel_transcriptome_v4.gz | csvtk cut -f "ID,Zm00001d028854" > P1_exp.csv
```
## run gwas for Zm00001d028854
```bash
python ../MRBIGR.py gwas -gwas -model lmm -thread 6 -g ../MRBIGR_data/MRBIGR_data/chr_HAMP -p P1_exp.csv
```
## draw manhattan and QQ plot for Zm00001d028854
```bash
python ../MRBIGR.py plot_new manhattan -i output/Zm00001d028854.assoc.txt -od MRBIGR_gwas -data_from file -plot_fmt jpg
python ../MRBIGR.py plot_new manhattan -i output/Zm00001d028854.assoc.txt -od MRBIGR_gwas_region -data_from file -chrom 1 -start 0 -end 100100000 -hl_pos chr1.s_48424403 -plot_fmt jpg
python ../MRBIGR.py plot_new qq -i output/Zm00001d028854.assoc.txt -od MRBIGR_gwas -t 1 -data_from file -plot_fmt jpg
```

## run gwas for 14 metabolites
```bash
echo -e "ID,n1109,n1370,n0146,n0511,n1111,n1144,n1201,n1240,n1372,n1391,n1555,n1562,n1569,n1570" > 14_metabolites.lst
csvtk cut -f "ID,n1109,n1370,n0146,n0511,n1111,n1144,n1201,n1240,n1372,n1391,n1555,n1562,n1569,n1570" ../MRBIGR_data/MRBIGR_data/E3_log2.normalized_phe.csv > 14_metabolites.abundance.csv
python ../MRBIGR.py gwas -gwas -model lmm -thread 1 -g ../MRBIGR_data/MRBIGR_data/chr_HAMP -p 14_metabolites.abundance.csv
```

## draw manhattan and QQ plot for n1109
```bash
python ../MRBIGR.py plot_new manhattan -i output/n1109.assoc.txt -od MRBIGR_gwas -data_from file -plot_fmt jpg
python ../MRBIGR.py plot_new manhattan -i output/n1109.assoc.txt -od MRBIGR_gwas_region -data_from file -plot_fmt jpg -chrom 1 -start 0 -end 100100000
python ../MRBIGR.py plot_new qq -i output/n1109.assoc.txt -od MRBIGR_gwas -data_from file -plot_fmt jpg
```
## draw manhattan and QQ plot for n1370
```bash
python ../MRBIGR.py plot_new manhattan -i output/n1370.assoc.txt -od MRBIGR_gwas -data_from file -plot_fmt jpg
python ../MRBIGR.py plot_new manhattan -i output/n1370.assoc.txt -od MRBIGR_gwas_region -data_from file -plot_fmt jpg -chrom 1 -start 0 -end 100100000
python ../MRBIGR.py plot_new qq -i output/n1370.assoc.txt -od MRBIGR_gwas -data_from file -plot_fmt jpg
```

# MR 
## run MODAS for qtl identification
```bash
MODAS.py genoidx -g ../MRBIGR_data/MRBIGR_data/chr_HAMP -genome_cluster -p 10 -o chr_HAMP
```
### run MODAS for 102 metabolites
```bash
MODAS.py prescreen -g ../test_files/chr_HAMP -genome_cluster chr_HAMP.genome_cluster.csv -phe ../MRBIGR_data/MRBIGR_data/102_metabolites.E3_log2.normalized_phe.csv -p 6 -o metabolite
MODAS.py regiongwas -g ../test_files/chr_HAMP -phe ../MRBIGR_data/MRBIGR_data/102_metabolites.E3_log2.normalized_phe.csv -phe_sig_qtl metabolite.phe_sig_qtl.csv -p 6 -o metabolite
```
### run MODAS for all 25,558 genes
```bash
MODAS.py prescreen -g ../test_files/chr_HAMP -genome_cluster chr_HAMP.genome_cluster.csv -phe ../MRBIGR_data/MRBIGR_data/AMP_kernel_transcriptome_v4 -p 6 -o genes -gwas_suggest_pvalue 1e-5
MODAS.py regiongwas -g ../test_files/chr_HAMP -phe ../MRBIGR_data/MRBIGR_data/AMP_kernel_transcriptome_v4 -phe_sig_qtl genes.phe_sig_qtl.csv -p 6 -o genes -p1 1e-7 -p2 1e-6 
```
### run Mendelian Randomization for genes to genes
```bash
MODAS.py mr -g ../test_files/chr_HAMP -exposure ../MRBIGR_data/MRBIGR_data/AMP_kernel_transcriptome_v4 -outcome ../MRBIGR_data/MRBIGR_data/AMP_kernel_transcriptome_v4 -qtl genes.region_gwas_qtl_res.csv -mlm -o g2g
```
### run Mendelian Randomization for genes to metabolites
MODAS.py mr -g ../test_files/chr_HAMP -exposure ../MRBIGR_data/MRBIGR_data/AMP_kernel_transcriptome_v4 -outcome ../MRBIGR_data/MRBIGR_data/102_metabolites.E3_log2.normalized_phe.csv -qtl genes.region_gwas_qtl_res.csv -mlm -o g2m

### draw forest and scatter plot for metabolites
```bash
ln -sf ../MRBIGR_data/MRBIGR_data/g2m.MR.csv
csvtk grep -f 2 -p "Zm00001d028854" g2m.MR.csv > g2m.p1.MR.csv
python ../MRBIGR.py plot_new forest -i g2m.p1.MR.csv -order_by_group -height 12 -group ../MRBIGR_data/MRBIGR_data/102_metabolites.groups -o MRBIGR_meta_forest
python ../MRBIGR.py plot_new scattermr -i g2m.p1.MR.csv -group_file ../MRBIGR_data/MRBIGR_data/102_metabolites.groups -order_by_group -sig_p 0.00049 -o MRBIGR_meta_scatter
```
### draw scatter plot for P1 related genes
```bash
ln -sf ../MRBIGR_data/MRBIGR_data/g2g.MR.csv
csvtk filter -f "pvalue<0.0001" g2g.MR.csv | csvtk grep -f 3 -p "Zm00001d028854" --quiet | csvtk cut -f 2 | csvtk grep -P - <(csvtk filter -f "pvalue<0.0001" g2g.MR.csv | csvtk grep -f 2 -p "Zm00001d028854" --quiet) -f 3 | csvtk rename -f 1,2,3,4,5,6 -n "snp,mTrait,pTrait,effect,TMR,pvalue" | grep -v "Inf" > g2g.p1.MR.csv
csvtk sort -k 6:n g2g.p1.MR.csv | csvtk cut -f 3 | sed '1d' > g2g.p1.MR.gene_order.csv
python ../MRBIGR.py plot_new scattermr -i g2g.p1.MR.csv -group_file ../MRBIGR_data/MRBIGR_data/genes.grp -order_by_pvalue -sig_p 0.0001 -o MRBIGR_gene_scatter
```

# GO and net
```bash
python ../MRBIGR.py plot_new go -i g2g.p1.MR.gene_order.csv -bg ../MRBIGR_data/MRBIGR_data/maize.genes2go.txt -height 9 -width 7 -o MRBIGR_P1_related_genes
cat <(echo "Zm00001d028854") <(csvtk cut -f 3 g2g.p1.MR.csv | sed '1d') <(cat ../MRBIGR_data/MRBIGR_data/maize_kernel_flavonoid_gene.lst) | sort -u > maize_flavonoid_related_62_genes.lst
selected_genes=$(cat maize_flavonoid_related_62_genes.lst | sed '1iID' | csvtk transpose)
csvtk rename -f 1 -n "ID" ../MRBIGR_data/MRBIGR_data/AMP_kernel_transcriptome_v4.gz | csvtk cut -f $selected_genes > maize_flavonoid_related_62_genes.exp.csv
python ../MRBIGR.py mr -g ../MRBIGR_data/MRBIGR_data/chr_HAMP -gene_exp maize_flavonoid_related_62_genes.exp.csv -pairwise -mlm -qtl genes.region_gwas_qtl_res.csv -thread 12 -o pairwise_mr_out
python ../MRBIGR.py plot_new net -i pairwise_mr_out.MR.csv -sig_p 0.05 -height 8 -width 8 -o MRBIGR_P1_related_genes
```

# nucleotide density
```bash
python ../MRBIGR.py plot_new nd -bed ../MRBIGR_data/MRBIGR_data/chr_HAMP -group ../MRBIGR_data/MRBIGR_data/AMP_537female_info.TEMP_TST.csv -gff_anno /home/xufeng/xufeng/Genomes/Zea_mays.AGPv4.50/Zea_mays.B73_RefGen_v4.50.gtf.gz -chrom 1 -chrom_start 48589178 -chrom_end 48597585 -left_offset 10000 -right_offset 10000 -window_size 1000 -window_step 100 -smooth -plot_left_expand 1000 -plot_right_expand 9000 -o MRBIGR_nd -od MRBIGR_nd
```

# allele frequency
```bash
csvtk grep -f 1 -p chr1.s_48424403 ../MRBIGR_data/MRBIGR_data/chr_HAMP_female.hmp.gz -t --quiet | cut -f 1,12- | csvtk transpose -t | tr '\t' ',' > snp_48424403.genotype
csvtk cut -f 1 ../MRBIGR_data/MRBIGR_data/AMP_kernel_transcriptome_v4.gz | sed '1d' | csvtk grep -P - ../MRBIGR_data/MRBIGR_data/AMP_537female_info.TEMP_TST.csv | csvtk grep -P <(csvtk cut -f 1 snp_48424403.genotype) - | tee 341_no_mixed.group | csvtk cut -f 1 | sed '1d' > 341_no_mixed.strains
python ../MRBIGR.py allele_freq -bed ../MRBIGR_data/MRBIGR_data/chr_HAMP -snps chr1.s_48424403 -group 341_no_mixed.group -group_sep "," -o MRBIGR_af_341_no_mixed -plot -plot_fmt pdf -width 3 -height 4
```

# gene expression boxplot
```bash
python ../MRBIGR.py plot_new phebox -i ../MRBIGR_data/MRBIGR_data/AMP_kernel_transcriptome_v4_FPKM.gz -snp snp_48424403.genotype -group ../MRBIGR_data/MRBIGR_data/AMP_537female_info.TEMP_TST.csv -s_strains 341_no_mixed.group -plot_type single_pop -s_genes Zm00001d028854,Zm00001d017077,Zm00001d037382,Zm00001d037383,Zm00001d047452,Zm00001d052673  -o MRBIGR_phebox_341_no_mixed
```
