source(".init.R")
load("data/M.Rdata")

# perform SVD on M data ..                          <
M.svd <- svd(x=M)

# .. signatures
n <- 1000
sigs <- extract.sigs(M.svd, n)

rownames(sigs) <- rownames(M)
save(sigs, file="data/sigs.Rdata")
