# ANCOM #

# Load example data
data(dietswap)
pseq = dietswap
n_taxa = ntaxa(pseq)
n_samp = nsamples(pseq)
# Metadata
meta_data = meta(pseq)
# Taxonomy table
taxonomy = tax_table(pseq)
# Absolute abundances
otu_absolute = abundances(pseq)

feature.table = otu_absolute; sample.var = "sample"; group.var = "nationality"; 
zero.cut = 0.90; lib.cut = 1000; neg.lb = TRUE
pre.process = feature_table_pre_process(feature.table, meta_data, sample.var, 
                                        group.var, zero.cut, lib.cut, neg.lb)
feature.table = pre.process$feature.table
group.name = pre.process$group.name
group.ind = pre.process$group.ind
struc.zero = pre.process$structure.zeros

# Paras for ANCOM-BC
grp.name = group.name; grp.ind = group.ind; adj.method = "bonferroni"
tol.EM = 1e-5; max.iterNum = 100; perNum = 1000; alpha = 0.05

out = ANCOM_BC(feature.table, grp.name, grp.ind, struc.zero,
               adj.method, tol.EM, max.iterNum, perNum, alpha)
res = cbind(taxon = rownames(out$feature.table), out$res)
write_csv(res, "demo_two_group.csv")