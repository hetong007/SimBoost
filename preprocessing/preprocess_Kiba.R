## TODO add something here to download the raw Kiba file directly from the authors website
## (paper: "making sense of large scale kinase inhibitor...", file is in supplementary material)
require("openxlsx")
kiba_data = read.xlsx("ci400709d_si_002.xlsx",5,startRow = 1)

## save the names of drugs and targets
## the first column holds the CHEMBL IDs of the compounds
compound_IDs = kiba_data[,1]
target_IDs = names(kiba_data[,-1])
##

dt_mat = kiba_data[,-1]
dt_mat = as.matrix(dt_mat)
class(dt_mat) = "numeric"

## put into triplet format
kiba_triplet = matrix(0, nrow = length(which(!is.na(dt_mat))), ncol = 3)
kiba_triplet[,1] = which(!is.na(dt_mat), arr.ind=T)[,1]
kiba_triplet[,2] = which(!is.na(dt_mat), arr.ind=T)[,2]
kiba_triplet[,3] = dt_mat[which(!is.na(dt_mat))]
##

## remove all entries with less than 10 entries or more than 200
threshold=10
flag = TRUE
while (flag)
{
  # row stage
  tb = table(kiba_triplet[,1])
  ind = which(tb<=threshold | tb>200)
  cat('removing ',length(ind),' drugs with less than ',threshold,' entries\n')
  # nms = as.numeric(names(tb)[ind])
  nms = names(tb)[ind]
  ind = match(kiba_triplet[,1],nms,0)
  ind = which(ind==0)
  kiba_triplet = kiba_triplet[ind,]
  # col stage
  tb = table(kiba_triplet[,2])
  if (min(tb)>threshold) {
    flag = FALSE
  } else {
    ind = which(tb<=threshold)
    cat('removing ',length(ind),' targets with less than ',threshold,' entries\n')
    # nms = as.numeric(names(tb)[ind])
    nms = names(tb)[ind]
    ind = match(kiba_triplet[,2],nms,0)
    ind = which(ind==0)
    kiba_triplet = kiba_triplet[ind,]
  }
  # check row
  tb = table(kiba_triplet[,1])
  if (min(tb)>threshold)
    flag = FALSE
}

## number of drugs and targets and density now:
length(unique(kiba_triplet[,1]))
length(unique(kiba_triplet[,2]))
nrow(kiba_triplet)/(length(unique(kiba_triplet[,2]))*length(unique(kiba_triplet[,1])))
##

################################################################
## now we need to get the compounds similarities
## first, match CHEMBL IDs to pubchem CIDs for the web interface
require('httr')
temp_compound_IDs = compound_IDs[unique(kiba_triplet[,1])]
df = data.frame(chembl_id=temp_compound_IDs, cid = NA)
## the following loop may have to be restarted in case of a timeout..
for (i in which(is.na(df[,2]))){
  temp = paste("http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",temp_compound_IDs[i],sep ="")
  temp = paste(temp, "/cids/TXT",sep="")
  cat(temp,'\n')
  ans = GET(temp)
  cat(status_code(ans),'\n')
  if (status_code(ans)==200){
    cat(content(ans, encoding = "UTF-8"),'\n')
    df[i,2] = strsplit(gsub("[\n]", " ", content(ans, encoding = "UTF-8"))," ")[[1]][1]
  } else {
    cat('CID not found','\n')
    df[i,2] = "undefined"
  }
}
## now some of the compounds could not be assigned CIDs, so some remapping is necessary
which(df[,2]=="undefined")
length(which(df[,2]=="undefined"))
ids_of_undefined = df[which(df[,2]=="undefined"),1]
temp = compound_IDs[kiba_triplet[,1]] 
inds = which(temp %in% ids_of_undefined)
kiba_triplet = kiba_triplet[-inds,]

## this is how many drugs and targets we have now
n_drugs = length(unique(kiba_triplet[,1]))
n_targets = length(unique(kiba_triplet[,2]))

## now we get the actual similarity matrix
temp = compound_IDs[unique(kiba_triplet[,1])] 
write.table(df[match(temp, df[,1]),2], file="compound_cids.txt", sep="\n", col.names = F, row.names = F, quote=F)
## now copy paste the content of compound_cids.txt into the web interface https://pubchem.ncbi.nlm.nih.gov/assay/assay.cgi?p=clustering
## download the clustering result and copy the file into this folder. Then:
compound_sim = read.csv("tanimoto_cluster_reqid_3913388554273586095.csv")
compound_ids_sim = compound_sim[,1]
##first coloumn contains ids, last column contains NAs
compound_sim = compound_sim[,c(-1,-ncol(compound_sim))]
## we need to do some mapping
cids = df[match(compound_IDs[unique(kiba_triplet[,1])], df[,1]),2]
length(which(cids%in%compound_ids_sim))
temp = match(cids, compound_ids_sim)
## this is our drug sim mat:
drug_sim_mat = as.matrix(compound_sim[temp, temp])

save(drug_sim_mat, file="kiba_drug_sim.rda")
################################################################


################################################################
## now we need to get the target similarity mat
kiba_target_IDs = target_IDs[unique(kiba_triplet[,2])]
paste(kiba_target_IDs, collapse=' ')
## first we need the target sequences. copy paste the identifiers here https://www.ncbi.nlm.nih.gov/protein/
## download them as fasta file and copy the resulting file here
require("Biostrings")
require("seqinr")
require("Rcpp")
seqs = read.fasta(file = "sequence.fasta.txt", seqtype = "AA")
## now compte the normalized Smith Waterman similarity

self_aligned_score = rep(NA, length(seqs))

for (i in 1:length(seqs)){
  cat('aligning ',i,'/',length(seqs),' with itself..\n')
  seq = paste(toupper(getSequence(seqs[i][[1]])), sep="", collapse="")
  align_score = pairwiseAlignment(seq, seq, scoreOnly=TRUE, gapExtension=0.5, type="local", substitutionMatrix = "BLOSUM50")
  self_aligned_score[i] = align_score
}
# now compute the alignments for all pairs of targets.. note that I choose alignment parameters type="local" to perform SW alignment
# and gapExtension = 0.5, type="local", substitutionMatrix = "BLOSUM50" to get the same results as the target-target similarities here
# http://staff.cs.utu.fi/~aatapa/data/DrugTarget/

all_combinations = combn(length(seqs),2)

sim_mat = matrix(NA, nrow = length(seqs), ncol = length(seqs))

for (i in 1:ncol(all_combinations)){
  target_a = all_combinations[1,i]
  target_b = all_combinations[2,i]
  cat('computing similarity for ',target_a,', ',target_b,'(',i,'/',ncol(all_combinations),')\n')
  seq_1 = paste(toupper(getSequence(seqs[target_a][[1]])), sep="", collapse="")
  seq_2 = paste(toupper(getSequence(seqs[target_b][[1]])), sep="", collapse="")
  align_score_ab = pairwiseAlignment(seq_1, seq_2, scoreOnly=TRUE, gapExtension = 0.5, type="local", substitutionMatrix = "BLOSUM50")
  similarity = align_score_ab/(sqrt(self_aligned_score[target_a])*sqrt(self_aligned_score[target_b]))
  sim_mat[target_a,target_b] = similarity
  cat('similarity: ',similarity,'\n')
}

for (i in 1:nrow(sim_mat)){
  for (k in 1:ncol(sim_mat)){
    if (!is.na(sim_mat[i,k])){
      sim_mat[k,i] = sim_mat[i,k]
    }
    if (i==k){
      sim_mat[i,k] = 1
    }
  }
}

target_sim_mat = sim_mat
save(target_sim_mat, file='kiba_target_sim.rda')

##########################################


kiba_triplet[,1] = match(kiba_triplet[,1], unique(kiba_triplet[,1]))
kiba_triplet[,2] = match(kiba_triplet[,2], unique(kiba_triplet[,2]))

## save everything in one file
kiba_target_sim = target_sim_mat
kiba_drug_sim = drug_sim_mat

save(kiba_triplet, kiba_drug_sim, kiba_target_sim, file='../data/kiba_data.Rda')
