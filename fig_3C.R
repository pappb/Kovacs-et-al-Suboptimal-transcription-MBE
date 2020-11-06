
require(gridExtra)
require(Hmisc)
source(".init.R")

# .. load signatures
load('data/sigs.Rdata')

# .. environments data
load("data/p.ds.imputed.Rdata")

# .. make tem all have the same rownames
common.rows <- intersect(rownames(p.ds), rownames(sigs))

p.ds <- p.ds[match(common.rows, rownames(p.ds)),]
p.sigs <- sigs[match(common.rows, rownames(sigs)),]

# .. load random data
load("../data_processed/ttc_m.Rdata")

for (i in 1:length(ttc.m)) {
    sigs.rand <- ttc.m[[i]]
    rownames(sigs.rand) <- rownames(sigs)
    sigs.rand <- sigs.rand[match(common.rows, rownames(sigs)),]
    sigs.rand <- sigs.rand[,1:15]
    ttc.m[[i]] <- sigs.rand
}

# correlations
get_env_mode_cors <- function(env, modes, sig_genes, nmodes)
    {
        ind.cors <- matrix(0, ncol=nmodes, nrow=ncol(env))
        ind.cors.p <- ind.cors

        for (i in 1:nmodes) {
            for (j in 1:ncol(env)) {
                aux <- cor.test(modes[,i], env[,j], method='s')
                ind.cors[j,i] <- aux$estimate
                ind.cors.p[j,i] <- aux$p.value
            }
        }
        rownames(ind.cors) <- colnames(env)
        rownames(ind.cors.p) <- colnames(env)

        return(list(ind.cors, ind.cors.p))
    }

real.cors <- get_env_mode_cors(p.ds, p.sigs[,1:15], sig.genes, nmodes=15)

# correlations with randoimized datasets
rand.cors <- lapply(ttc.m[1:600], get_env_mode_cors,
                    env=p.ds, sig_genes = sig.genes, nmodes=15)
saveRDS(rand.cors, 'rand_cors.rds')

aux.rand <- sapply(rand.cors, function(x) apply(abs(x[[1]]), 2, max))
aux.rand <- melt(aux.rand)
names(aux.rand) <- c('mode', 'iter', 'rho')

# rearrange data and plot
ind.cors <- melt(real.cors[[1]])
ind.p <- melt(real.cors[[2]])
names(ind.p)[3] <- 'pval'
ind.cors <- merge(ind.cors, ind.p, by=c('Var1', 'Var2'), all=TRUE)
names(ind.cors) <- c('env', 'mode', 'rho', 'p')
setDT(ind.cors)

ind.cors <- ind.cors[abs(rho)>0.2 & p<0.01]

p.direct <- ggplot(ind.cors) +
    xlab('Mode') + ylab(expression("Spearman "*rho)) +
    scale_x_continuous(breaks=1:15) +
    theme_bw() + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
                       aspect.ratio=1, legend.position = c(0.85, 0.85),
                       panel.grid.minor.x = element_blank()) +
    geom_jitter(aes(mode, abs(rho)), size=1, width=.3,
                color='orange3', alpha=.5) +
    geom_boxplot(data=aux.rand, aes(x=mode, group=mode, y=rho))



pdf("results/env_mode_correlations.pdf", 4,4)
print(p.direct)
dev.off()


