require(viridis)
require(ggdendro)
rm(list=ls())
source(".init.R")

##_____________________________________________________________________________
##
## .. % explained by first x modes vs same % in individual deletions
##_____________________________________________________________________________

## .. computing the percentage of explained variance in each column
load("data/M.Rdata")
load('data/sigs.Rdata')

M.svd <- svd(M)

i <- 1 # mode

aux <- M.svd$u[,i] %*% t(M.svd$d[i]*M.svd$v[,i]) # projection

## .. percentage of variation explained by mode i overall
sum(aux^2)/sum(M^2)

## .. and the same as:
(M.svd$d[i]^2)/sum((M.svd$d^2))

## ..
## .. now compute the percentage of variation explained column k
k <- 2
aux <- M.svd$u[,i] * M.svd$d[i]* M.svd$v[k,i] # projection

sum(aux^2)/sum(M[,k]^2)

# .. same as:
(M.svd$d[i]^2 * M.svd$v[k,i]^2) / sum(M[,k]^2)


## .. compute for k modes and all deletions
k <- 1484

M.svd <- svd(M)

pvar.gnt <- matrix(0, ncol=k, nrow=ncol(M))
for (i in 1:k) {
    aux <- M.svd$u[,i] %*% t(M.svd$d[i]*M.svd$v[,i]) # projection
    pvar.gnt[,i] <- apply(aux^2,2,sum) / apply(M^2,2,sum)
}

# correlation, to see if it is positive or negative
pcor.gnt <- cor(M, sigs[,1:15])

# two methods, either using explained variance or correlation
pv <- pvar.gnt[,1:15]

get.signif.genes <- function(x, cutoff, names) names[which(x>cutoff)]

sgnf.pv <- apply(pv, 2, get.signif.genes, cutoff = .1, names=rownames(pv))

aux <- pv
aux[which(pcor.gnt<0)] <- 0
sgnf.pv.pos <- apply(aux, 2, get.signif.genes, cutoff = .1, names=rownames(pv))

aux <- pv
aux[which(pcor.gnt>0)] <- 0
sgnf.pv.neg <- apply(aux, 2, get.signif.genes, cutoff = .1, names=rownames(pv))


sgnf.pv.tab.p <- data.table(ORF = unlist(sgnf.pv.pos),
                            Mode = rep(1:15, sapply(sgnf.pv.pos, length)))
sgnf.pv.tab.n <- data.table(ORF = unlist(sgnf.pv.neg),
                            Mode = rep(1:15, sapply(sgnf.pv.neg, length)))

# Compute GO enrichment
slim <- fread('data/go_slim_mapping.tab', header=F)

# .. consider only process slim categories
slim <- slim[V4=='P']
slim <- slim[V5!='biological_process']

# compute hypergeometric dist. based enrichment
get.slim.enrichment <- function(genes, background, slim_mapping)
    {
        sl <- slim[V1 %in% background]
        scat <- unique(sl$V5)
        enr <- data.table("CAT"=scat, "f"=NA, "p"=NA, 'pres'=NA)

        # total genes in background
        tot <- length(background)

        for (j in 1:nrow(enr)) {

            # total genes with slim annotation j
            pres <- nrow(sl[V5==enr$CAT[j]])

            # genes with annotation in target gene set
            target <- nrow(sl[V5==enr$CAT[j] & V1 %in% genes])

            enr$f[j] <- target/(length(genes)*(pres/tot))
            enr$p[j] <- 1-phyper(target, pres, tot-pres, length(genes))
            enr$pres[j] <- pres

        }
        enr$p <- p.adjust(enr$p, 'fdr')
        return(enr)
    }

# meaningfully cluster go slim categories
slim.cat <- dcast(slim, V5~V1)
rn <- slim.cat$V5
slim.cat[, V5:= NULL]
slim.cat[!is.na(slim.cat)] <- 1
slim.cat[is.na(slim.cat)] <- 0
slim.cat <- as.matrix(slim.cat)
rownames(slim.cat) <- rn
slim.clust <- hclust(dist(slim.cat))

pdf("results/fig_S4_tree.pdf", 5, 9)
ggdendrogram(slim.clust, rotate = TRUE, theme_dendro = FALSE)
dev.off()

## compute go enrichments
pv.enr.pos <- lapply(sgnf.pv.pos, get.slim.enrichment,
                     background = colnames(M),
                     slim_mapping = slim)
pv.enr.neg <- lapply(sgnf.pv.neg, get.slim.enrichment,
                     background = colnames(M),
                     slim_mapping = slim)

aux <- nrow(pv.enr.pos[[1]])
pv.enr.pos <- rbindlist(pv.enr.pos)
pv.enr.pos$mode <- rep(1:15, each=aux)
pv.enr.pos$CAT <- ordered(pv.enr.pos$CAT, slim.clust$labels[slim.clust$order])

aux <- nrow(pv.enr.neg[[1]])
pv.enr.neg <- rbindlist(pv.enr.neg)
pv.enr.neg$mode <- rep(1:15, each=aux)
pv.enr.neg$CAT <- ordered(pv.enr.neg$CAT, slim.clust$labels[slim.clust$order])

#
pv.enr.pos$tail <- '+'
pv.enr.neg$tail <- '-'

aux <- rbind(pv.enr.pos, pv.enr.neg)

fwrite(aux, 'table_fig_S4.csv')

pdf("results/fig_S4_heatmap.pdf", 5.5, 9)
ggplot(aux[p<.05 & f>1 & pres>3]) +
    geom_tile(aes(x = tail, y = CAT, fill = log(f))) +
    scale_fill_viridis_c(na.value = NA) +
    scale_x_discrete(limits = c('-', '+')) +
    facet_grid(~ factor(as.numeric(mode)),
               space='free_x', scales='free_x', switch='x', drop = FALSE) +
    theme_bw() +
    theme(strip.placement = 'outside',
          strip.background = element_rect(fill=NA,colour='grey50'),
          panel.spacing=unit(0,'cm'),
          legend.position = 'top',
          axis.text.y = element_text(size=7)) +
    xlab('Mode')
dev.off()


# to what extent are deletions triggering any of the first 15 modes
# transcriptional regulators?

all.triggering <- unique(unlist(sgnf.pv))

aux <- slim[V1 %in% all.triggering]

transcr.related <- c("transcription from RNA polymerase II promoter",
                     "chromatin organization",
                     "histone modification",
                     "DNA-templated transcription, elongation",
                     "DNA-templated transcription, initiation",
                     "transcription from RNA polymerase III promoter",
                     "transcription from RNA polymerase I promoter")

aux <- aux[V5 %in% transcr.related]


length(all.triggering)
(length(all.triggering)-length(unique(aux$V1)))
(length(all.triggering)-length(unique(aux$V1)))/length(all.triggering)
