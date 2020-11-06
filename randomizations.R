
require(wrswoR)
source(".init.R")

load("data/M.Rdata")

# Perform the ttc randomization method along with row/col shuffling

# generate "tendency to change" (randomization structure)
ttc <- cbind(apply(abs(M),1,sum)) %*% apply(abs(M),2,sum)

# ordered values for randomization
M.s <- M[order(abs(M), decreasing=T)]

M.p <- which(ttc>0, arr.ind=T)
M.p <- cbind(M.p, ttc[M.p]^alpha)
M.p <- M.p[order(M.p[,3], decreasing=T),]

ttc.d <- vector("list", 250)
row.m <- ttc.d
col.m <- ttc.d

n <- 250   # number of randomizations to perform
alpha <- 2 # weights

print("Building randomized matrices...")
for (i in 1:n) {

    # ttc randomization
    rand <- array(0, dim(M))
    rand[M.p[,1:2]] <- M.s[sample_int_crank(length(M.s), length(M.s), prob=M.p[,3])]
    colnames(rand) <- colnames(M)
    rownames(rand) <- rownames(M)
    ttc.d[[i]] <- rand

    # shuffling rows
    row.m[[i]] <- t(apply(M,1,sample))

    # shuffling cols
    col.m[[i]] <- apply(M,2,sample)

    print(i)
}

save(ttc.d, file="randomizations/ttc_d.Rdata")
save(col.m, file="randomizations/col_m.Rdata")
save(row.m, file="randomizations/row_m.Rdata")







