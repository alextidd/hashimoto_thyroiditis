# TODO: remove germline SNPs (from caveman output)
# TODO: remove variants in >50% of cells
# TODO: remove RNA editing sites
# TODO: remove variants with global VAF > 0.25
# TODO: remove variants with cell VAF < 0.25 or mutant reads < 5
# TODO: rerun Sequoia filters - exclude cells with copy number changes from contributing to the germline filter
# TODO: generate a heterozygous SNP histogram per cell to look for outliers

# libraries
library(magrittr)
library(dplyr)
library(ggplot2)
library(tibble)
library(ape)
source("bin/build_phylogeny.R")

# dirs
seq_dir <- "out/resolveome/sequoia/"
dir.create(seq_dir, recursive = TRUE, showWarnings = FALSE)

# load matrices
mats <-
  c("NV", "NR") %>%
  purrr::set_names() %>%
  purrr::map(function(i) {
    file.path(Sys.getenv("LUSTRE_125"), "projects/hashimoto_thyroiditis",
              "out/resolveome/basejumper/bj-somatic-variantcalling/dna/PD63118/",
              "PD63118_run/SOMATIC_VARIANT_WORKFLOW_Heuristic_Filter_SEQUOIA",
              paste0("Sequoia_group_null_bino-10_rhosnp0.4_rhoindel0.4_mincov10_maxcov500_both_",
                     i, "_filtered_all.txt")) %>%
      read.table()
  })

# calculate VAFs and global VAFs
vafs <- mats$NV / (mats$NR + mats$NV)
global_vafs <- rowSums(mats$NV) / (rowSums(mats$NR) + rowSums(mats$NV))

# plot global VAF distribution
tibble::enframe(global_vafs, name = "mut_id", value = "global_vaf") %>%
  ggplot(aes(x = global_vaf, fill = grepl("chr1_", mut_id))) +
  geom_histogram(bins = 100) +
  lims(x = c(NA, 1)) +
  scale_y_continuous(trans = scales::pseudo_log_trans()) +
  theme_bw() +
  geom_vline(xintercept = 0.2, linetype = "dashed", color = "red")

# count variants in >50% of cells
n_cells_w_var <- rowSums(mats$NV > 0) / ncol(mats$NV)
table(n_cells_w_var > 0.5)

# count calls with <= 5 reads
table(mats$NR <= 5)

# count calls with mut_vaf < 0.3
table(vafs > 0.3)

# apply filters
filtered_mats <-
  mats %>%
  purrr::map(function(i) {
    # Set values to 0 where filters fail
    i[mats$NR <= 5] <- 0
    i[vafs < 0.3] <- 0
    # Keep only variants with global VAF >= 0.2 AND present in <50% of cells
    i[global_vafs < 0.2, ]
    # i[global_vafs < 0.2 & n_cells_w_var < 0.5, ]
  })

# remove variants with no remaining calls
filtered_mats <-
  filtered_mats %>%
  purrr::map(~ .x[rowSums(filtered_mats$NR) > 0, ])

# save nr/nv
purrr::walk2(names(filtered_mats), filtered_mats, function(i, mat) {
  write.table(mat, file = file.path(seq_dir, paste0(i, "_filtered.tsv")),
              quote = FALSE, sep = "\t")
})

# set sequoia defaults
opt <- list(
  donor_id = 'Patient',
  input_nv = NULL,
  input_nr = NULL,
  cgpvaf_output = NULL,
  output_dir = "",
  beta_binom_shared = TRUE,
  ncores = 1,
  normal_flt = 'PDv37is',
  snv_rho = 0.1,
  indel_rho = 0.15,
  min_cov = 10,
  max_cov = 500,
  only_snvs = TRUE,
  split_trees = TRUE,
  keep_ancestral = FALSE,
  exclude_samples = NULL,
  cnv_samples = NULL,
  vaf_absent = 0.1,
  vaf_present = 0.3,
  mixmodel = FALSE,
  min_clonal_mut = 35,
  tree_mut_pval = 0.01,
  genotype_conv_prob = FALSE,
  min_pval_for_true_somatic = 0.05,
  min_variant_reads_shared = 2,
  min_vaf_shared = 2,
  create_multi_tree = TRUE,
  mpboot_path = "",
  germline_cutoff = -5,
  genomeFile = "/nfs/cancer_ref01/Homo_sapiens/37/genome.fa",
  plot_spectra = FALSE,
  max_muts_plot = 5000,
  lowVAF_filter = 0,
  lowVAF_filter_positive_samples = 0,
  VAF_treshold_mixmodel = 0.3
)
dp_pos=opt$lowVAF_filter_positive_samples
ncores=opt$ncores
lowVAF_threshold=opt$lowVAF_filter
normal_flt=opt$normal_flt
snv_rho=opt$snv_rho
genomeFile=opt$genomeFile
plot_spectra=opt$plot_spectra
VAF_treshold=opt$VAF_treshold_mixmodel
indel_rho=opt$indel_rho
min_cov=opt$min_cov
max_cov=opt$max_cov
output_dir=opt$output_dir
only_snvs=opt$only_snvs
germline_cutoff=opt$germline_cutoff
if(is.null(opt$exclude_samples)) {samples_exclude=NULL} else {samples_exclude=unlist(strsplit(x=opt$exclude_samples,split = ","))}
if(is.null(opt$cnv_samples)) {samples_with_CNVs=NULL} else {samples_with_CNVs=unlist(strsplit(x=opt$cnv_samples,split = ","))}
if(is.null(opt$cgpvaf_output)) {cgpvaf_paths=NULL} else {cgpvaf_paths=unlist(strsplit(x=opt$cgpvaf_output,split = ","))}
keep_ancestral=opt$keep_ancestral
patient_ID=opt$donor_id
output_dir=opt$output_dir
nv_path=opt$input_nv
nr_path=opt$input_nr
max_muts_plot=opt$max_muts_plot
VAF_present=opt$vaf_present
VAF_absent=opt$vaf_absent
mixmodel=opt$mixmodel
split_trees=opt$split_trees
genotype_conv_prob=opt$genotype_conv_prob
min_pval_for_true_somatic_SHARED = opt$min_pval_for_true_somatic
min_variant_reads_SHARED=opt$min_variant_reads_shared
min_vaf_SHARED=opt$min_vaf_shared
tree_mut_pval=opt$tree_mut_pval
beta_binom_shared=opt$b
create_multi_tree=opt$create_multi_tree
path_to_mpboot=opt$mpboot_path
min_clonal_mut=opt$min_clonal_mut

# set sequoia arguments
NV_filtered <- filtered_mats$NV
NR_filtered <- filtered_mats$NR
NR_flt_nonzero <- filtered_mats$NR
NR_flt_nonzero[NR_flt_nonzero == 0] <- 1
genotype_bin <- as.matrix(NV_filtered / NR_flt_nonzero)
gender <- "male"
mut_id <- "both"
output_dir <- seq_dir

# create genotype bins
XY_chromosomal <- grepl("X|Y",rownames(NR_filtered))
autosomal <- !XY_chromosomal
genotype_bin <- as.matrix(NV_filtered / NR_flt_nonzero)
if (gender == "male") {
  genotype_bin[autosomal,][genotype_bin[autosomal,]<VAF_absent]=0
  genotype_bin[autosomal,][genotype_bin[autosomal,]>=VAF_present]=1
  genotype_bin[XY_chromosomal,][genotype_bin[XY_chromosomal,]<(2*VAF_absent)]=0
  genotype_bin[XY_chromosomal,][genotype_bin[XY_chromosomal,]>=(2*VAF_present)]=1
  genotype_bin[genotype_bin>0&genotype_bin<1]=0.5
} else if(gender=="female"){
  genotype_bin[genotype_bin<VAF_absent]=0
  genotype_bin[genotype_bin>=VAF_present]=1
  genotype_bin[genotype_bin>0&genotype_bin<1]=0.5
}
present_vars_full <- rowSums(genotype_bin > 0) > 0

# create dummy fasta consisting of As (WT) and Ts (Mutant)
Ref = rep("A",nrow(genotype_bin))
Alt = rep("T",nrow(genotype_bin))
dna_strings = list()
dna_strings[1]=paste(Ref,sep="",collapse="") #Ancestral sample
for (n in 1:ncol(genotype_bin)){
  Mutations = Ref
  Mutations[genotype_bin[,n]==0.5] = '?'
  Mutations[genotype_bin[,n]==1] = Alt[genotype_bin[,n]==1]
  dna_string = paste(Mutations,sep="",collapse="")
  dna_strings[n+1]=dna_string
}
names(dna_strings) <- c("Ancestral", colnames(genotype_bin))
require(seqinr)
write.fasta(dna_strings, names=names(dna_strings),paste0(output_dir,patient_ID,"_",mut_id,"_for_MPBoot.fa"))

# build tree
system(paste0(Sys.getenv("NFS_TEAM"), "/bin/mpboot -s ",
              output_dir, patient_ID, "_", mut_id, "_for_MPBoot.fa -bb 1000"),
       ignore.stdout = TRUE)

# mapping mutations onto the tree
tree=read.tree(paste0(output_dir,patient_ID,"_",mut_id,"_for_MPBoot.fa.treefile"))
tree=drop.tip(tree,"Ancestral")
if(!keep_ancestral){
  tree$edge.length=rep(1,nrow(tree$edge))
  NR_tree=NR_filtered[present_vars_full,]
  NV_tree=NV_filtered[present_vars_full,]
  res=assign_to_tree(tree,
                     mtr=as.matrix(NV_tree),
                     dep=as.matrix(NR_tree))
}else{
  tree <- add_ancestral_outgroup(tree) #Re add the ancestral outgroup after making tree dichotomous - avoids the random way that baseline polytomy is resolved
  tree$edge.length = rep(1, nrow(tree$edge)) 

  NR_tree=NR_filtered[present_vars_full,]
  NR_tree$Ancestral=30
  NV_tree=NV_filtered[present_vars_full,]
  NV_tree$Ancestral=0

  p.error = rep(0.01,ncol(NV_tree))
  p.error[colnames(NV_tree)=="Ancestral"]=1e-6
  res=assign_to_tree(tree,
                     mtr=as.matrix(NV_tree),
                     dep=as.matrix(NR_tree),
                     error_rate = p.error)
}

edge_length_nonzero = table(res$summary$edge_ml[res$summary$p_else_where<tree_mut_pval])
edge_length = rep(0,nrow(tree$edge))
names(edge_length)=1:nrow(tree$edge)
edge_length[names(edge_length_nonzero)]=edge_length_nonzero
tree$edge.length=as.numeric(edge_length)

if(create_multi_tree){
  print("Converting to a multi-furcating tree structure")
  if(keep_ancestral) {
    #Maintain the dichotomy with the ancestral branch
    ROOT=tree$edge[1,1]
    current_length<-tree$edge.length[tree$edge[,1]==ROOT & tree$edge[,2]!=which(tree$tip.label=="Ancestral")]
    new_length<-ifelse(current_length==0,1,current_length)
    tree$edge.length[tree$edge[,1]==ROOT & tree$edge[,2]!=which(tree$tip.label=="Ancestral")]<-new_length
  }
  tree<-di2multi(tree) #Now make tree multifurcating
  #Re-run the mutation assignment algorithm from the new tree
  res=assign_to_tree(tree,
                     mtr=as.matrix(NV_tree),
                     dep=as.matrix(NR_tree))
  edge_length_nonzero = table(res$summary$edge_ml[res$summary$p_else_where<tree_mut_pval])
  edge_length = rep(0,nrow(tree$edge))
  names(edge_length)=1:nrow(tree$edge)
  edge_length[names(edge_length_nonzero)]=edge_length_nonzero
  tree$edge.length=as.numeric(edge_length)
}

saveRDS(res,paste0(output_dir,patient_ID,"_",mut_id,"_assigned_to_tree.Rdata"))
write.tree(tree, paste0(output_dir,patient_ID,"_",mut_id,"_tree_with_branch_length.tree"))

if(split_trees & mut_id=="both"){
  Muts_coord=matrix(ncol=4,unlist(strsplit(rownames(NV_filtered)[present_vars_full],split="_")),byrow = T)
  is.indel=nchar(Muts_coord[,3])>1|nchar(Muts_coord[,4])>1
  
  edge_length_nonzero = table(res$summary$edge_ml[res$summary$p_else_where<tree_mut_pval&!is.indel])
  edge_length = rep(0,nrow(tree$edge))
  names(edge_length)=1:nrow(tree$edge)
  edge_length[names(edge_length_nonzero)]=edge_length_nonzero
  tree$edge.length=as.numeric(edge_length)
  pdf(paste0(output_dir,patient_ID,"_snv_tree_with_branch_length.pdf"), height = 20, width = 20)
  plot(tree)
  axisPhylo(side = 1,backward=F)
  dev.off()
  write.tree(tree, paste0(output_dir,patient_ID,"_snv_tree_with_branch_length.tree"))
  
  edge_length_nonzero = table(res$summary$edge_ml[res$summary$p_else_where<tree_mut_pval&is.indel])
  edge_length = rep(0,nrow(tree$edge))
  names(edge_length)=1:nrow(tree$edge)
  edge_length[names(edge_length_nonzero)]=edge_length_nonzero
  tree$edge.length=as.numeric(edge_length)
  pdf(paste0(output_dir,patient_ID,"_indel_tree_with_branch_length.pdf"), height = 20, width = 20)
  plot(tree)
  axisPhylo(side = 1,backward=F)
  dev.off()
  write.tree(tree, paste0(output_dir,patient_ID,"_indel_tree_with_branch_length.tree"))
  
}else{
  pdf(paste0(output_dir,patient_ID,"_",mut_id,"_tree_with_branch_length.pdf"))
  plot(tree)
  axisPhylo(side = 1,backward=F)
  dev.off()
  
  tree_collapsed=tree
  tree_collapsed$edge.length=rep(1,nrow(tree_collapsed$edge))
  pdf(paste0(output_dir,patient_ID,"_",mut_id,"_tree_with_equal_branch_length.pdf"))
  plot(tree_collapsed)
  dev.off()
}

Mutations_per_branch=as.data.frame(matrix(ncol=4,unlist(strsplit(rownames(NR_tree),split="_")),byrow = T))
colnames(Mutations_per_branch)=c("Chr","Pos","Ref","Alt")
Mutations_per_branch$Branch = tree$edge[res$summary$edge_ml,2]
Mutations_per_branch=Mutations_per_branch[res$summary$p_else_where<tree_mut_pval,]
Mutations_per_branch$Patient = patient_ID
Mutations_per_branch$SampleID = paste(patient_ID,Mutations_per_branch$Branch,sep="_")
write.table(Mutations_per_branch,paste0(output_dir,patient_ID,"_",mut_id,"_assigned_to_branches.txt"),quote=F,row.names=F,sep="\t")