## davis data
## first get the interaction matrix and put into triplet format
download.file(url="http://staff.cs.utu.fi/~aatapa/data/DrugTarget/drug-target_interaction_affinities_Kd__Davis_et_al.2011.txt", destfile = "davis_raw.txt")
davis_mat = read.table('davis_raw.txt')
davis_mat = as.matrix(davis_mat)
davis_triplet = matrix(0, nrow = length(which(!is.na(davis_mat))), ncol = 3)
davis_triplet[,1] = which(!is.na(davis_mat), arr.ind=T)[,1]
davis_triplet[,2] = which(!is.na(davis_mat), arr.ind=T)[,2]
davis_triplet[,3] = davis_mat[which(!is.na(davis_mat))]
davis_triplet[,3] = -log(davis_triplet[,3]/1e9,10)
## download the drug similarity matrix
download.file(url="http://staff.cs.utu.fi/~aatapa/data/DrugTarget/drug-drug_similarities_2D.txt", destfile = "davis_drug_sim_mat.txt")
davis_drug_sim = read.table('davis_drug_sim_mat.txt')
davis_drug_sim = as.matrix(davis_drug_sim)
## download the target similarity matrix
download.file(url="http://staff.cs.utu.fi/~aatapa/data/DrugTarget/target-target_similarities_WS_normalized.txt", destfile = "davis_target_sim_mat.txt")
davis_target_sim = read.table('davis_target_sim_mat.txt')
davis_target_sim = as.matrix(davis_target_sim)
davis_target_sim = davis_target_sim/100
save(davis_triplet, davis_drug_sim, davis_target_sim, file = "../data/davis_data.Rda")

## same thing for the metz data
download.file(url="http://staff.cs.utu.fi/~aatapa/data/DrugTarget/known_drug-target_interaction_affinities_pKi__Metz_et_al.2011.txt", destfile = "metz_raw.txt")
metz_mat = read.table('metz_raw.txt')
metz_mat = as.matrix(metz_mat)
metz_triplet = matrix(0, nrow = length(which(!is.na(metz_mat))), ncol = 3)
metz_triplet[,1] = which(!is.na(metz_mat), arr.ind=T)[,1]
metz_triplet[,2] = which(!is.na(metz_mat), arr.ind=T)[,2]
metz_triplet[,3] = metz_mat[which(!is.na(metz_mat))]
## download the drug similarity matrix
download.file(url="http://staff.cs.utu.fi/~aatapa/data/DrugTarget/drug-drug_similarities_2D__Metz_et_al.2011.txt", destfile = "metz_drug_sim_mat.txt")
metz_drug_sim = read.table('metz_drug_sim_mat.txt')
metz_drug_sim = as.matrix(metz_drug_sim)
metz_drug_sim = metz_drug_sim/100
## download the target similarity matrix
download.file(url="http://staff.cs.utu.fi/~aatapa/data/DrugTarget/target-target_similarities_WS_normalized__Metz_et_al.2011.txt", destfile = "metz_target_sim_mat.txt")
metz_target_sim = read.table('metz_target_sim_mat.txt')
metz_target_sim = as.matrix(metz_target_sim)
metz_target_sim = metz_target_sim/100
save(metz_triplet, metz_drug_sim, metz_target_sim, file = "../data/metz_data.Rda")

