# args
msa_attributions_file <- commandArgs(trailingOnly=TRUE)[1]
seq_dir <- commandArgs(trailingOnly=TRUE)[2]
out_dir <- commandArgs(trailingOnly=TRUE)[3]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# msa_attributions_file <- "out/resolveome/signatures/sigprofiler/mutational_signatures/sigprofiler_branches/output/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt"; seq_dir <- "out/resolveome/sequoia/20250918" ; out_dir <- "out/resolveome/signatures/sigprofiler/mutational_signatures/sigprofiler_branches/output/sequoia/"

# libraries
library(ape)
library(stringr)

# load msa attributions
exposures <- read.table(msa_attributions_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
row.names(exposures)=exposures$Samples

# initiate an exposures pile up: 
exposures_pileup <- data.frame()

## load tree file for the patient
tree_file=paste0(seq_dir, "/Patient_snv_tree_with_branch_length.tree")
tree = read.tree(tree_file)
tree_df=ggtree::fortify(tree)

## calculate the number of mutations that were not assigned to mutational signatures : (this relates to MSA output)
node_muts <- tree_df %>% dplyr::select(node, branch.length)
exposures$node=as.numeric(str_extract(exposures$Samples,'(?<=\\_)\\d+'))
exposures=left_join(x=exposures, y= node_muts,by='node')
exposures<-exposures%>%dplyr::select(-node)  
row.names(exposures)=exposures$Samples

## convert the truncal branches to unassigned:
# exposures$SBSTrunk=ifelse(exposures$Samples%in%trunk_branches,exposures$SBS1+exposures$SBS5+exposures$SBS18+exposures$SBS88+exposures$SBSUnassigned,0) # where branch_ids refers to pile up of truncal branch ids
# exposures$SBS1=ifelse(exposures$SBSTrunk>0,0,exposures$SBS1)
# exposures$SBS5=ifelse(exposures$SBSTrunk>0,0,exposures$SBS5)
# exposures$SBS18=ifelse(exposures$SBSTrunk>0,0,exposures$SBS18)
# exposures$SBS88=ifelse(exposures$SBSTrunk>0,0,exposures$SBS88)
# exposures$SBSUnassigned=ifelse(exposures$SBSTrunk>0,0,exposures$SBSUnassigned) 

## construct a pile up of exposures to ensure consistency across trees and downstream analysis :
exposures_pileup<-rbind(exposures_pileup,exposures)
## list the extract signatures :
sigs=colnames(exposures)[2:length(colnames(exposures))]
#colours for your signatures
cols= c("wheat2","powderblue","darkslateblue","maroon3","ivory3","grey30","#ffac3c","darkorange4","springgreen3","#0b6623","darkolivegreen2","mediumorchid3","pink","red","ivory4")
cols <- cols[1:length(sigs)]
names(cols)=sigs

samples=exposures$Samples
branches=sub('.*_', '', samples) #get the branch number, adjust to your sample names

pdf(paste0(out_dir, "signature_annotations.pdf"), height = 10)
br_lengths=NULL
plot(tree, cex = 0.7, label.offset = 0.01 * max(tree_df$x))

axisPhylo(side=1, las=1, backward = F)
for (k in 1:length(samples)){ #length(samples)

  n=as.numeric(branches[k])
  x_end=tree_df$x[n]
  x_start=tree_df$x[tree_df$parent[n]]
  x_intv=x_end-x_start
  y=node.height(tree)[n]
  l=tree_df$branch.length[n]
  t=samples[k]
  br_lengths <- rbind(br_lengths,cbind(t, l))
  tipnum=sum(tree_df$isTip)
  for (s in sigs){
    x_end=x_start+ exposures[samples[k],s]
    rect(ybottom=y-min(0.02*tipnum,0.4),ytop=y+min(0.02*tipnum,0.4),xleft=x_start,xright=x_end,col=cols[s], border ='black' , lwd = 0.5)
    x_start=x_end
    
  }
}

legend("topright", title = "signatures", legend = names(cols),
       fill = cols, bty = "n", cex = 0.8, ncol = 1, xjust = 0.5)
dev.off()
