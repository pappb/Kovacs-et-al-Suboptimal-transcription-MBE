## ..
require(ggplot2)
require(ggdendro)
require(viridis)
source(".init.R")

load('data/sigs.Rdata')

tf.gene <- fread('data/TFmatrix/RegulationMatrix_Documented_2020219_DNAbindingANDexpressionChange_Yeastract.csv')
tf.gene <- tf.gene[, !'Uncharacterized']
tf.gene <- melt(tf.gene)
setDT(tf.gene)
tf.gene <- tf.gene[value==1]
names(tf.gene) <- c('TF', 'ORF', 'bind')
tf.gene[, bind:=NULL]

tf.bnd <- fread('data/TFmatrix/RegulationMatrix_Documented_2020219_DNAbindingEvidence.csv')
tf.bnd <- tf.bnd[, !'Uncharacterized']
tf.bnd <- melt(tf.bnd)
setDT(tf.bnd)
tf.bnd <- tf.bnd[value==1]
names(tf.bnd) <- c('TF', 'ORF', 'bind')
tf.bnd[, bind:=NULL]

tf <- unique(rbind(tf.gene, tf.bnd))

gene.orf <- fread('data/TFmatrix/orftogene.csv', header=F)[, .(V1, V3)]
names(gene.orf) <- c('ORF', 'Gene')
tf$ORF <- gene.orf$ORF[match(tf$ORF, gene.orf$Gene)]

tf.rel <- tf[ORF %in% rownames(sigs)]


# responding genes in each mode
sg <- fread('data/responding_genes_prjs_FC_1.3.csv')
only.sig <- sg[above_FC_1.3==1]
only.sig[, expression := NULL]
only.sig <- unique(only.sig)

md.enr <- vector('list', 15)
for (i in 1:15) {
    aux <- vector('list', 2)
    r.genes <- only.sig[mode==i]$responding_gene
    r.values <- sigs[match(r.genes, rownames(sigs)),i]

    aux[[1]] <- r.genes[which(r.values>0)]
    aux[[2]] <- r.genes[which(r.values<0)]

    md.enr[[i]] <- aux
}

# compute hypergeometric dist. based enrichment
get.tf.enrichment <- function(genes, background, tf_mapping)
    {
        sl <- tf_mapping[ORF %in% background]
        scat <- unique(sl$TF)
        enr <- data.table("CAT"=scat, "f"=NA, "p"=NA, 'pres' = NA)

        # total genes in background
        tot <- length(background)

        for (j in 1:nrow(enr)) {

            # total genes with slim annotation j
            pres <- nrow(sl[TF==enr$CAT[j]])

            # genes with annotation in target gene set
            target <- nrow(sl[TF==enr$CAT[j] & ORF %in% genes])

            enr$f[j] <- target/(length(genes)*(pres/tot))
            enr$p[j] <- 1-phyper(target, pres, tot-pres, length(genes))
            enr$pres[j] <- pres

        }
        enr$p <- p.adjust(enr$p, 'fdr')
        return(enr)
    }

enr <- lapply(md.enr, function(x) lapply(x, get.tf.enrichment,
                                         background = rownames(sigs),
                                         tf_mapping = tf.rel))

enr <- rbindlist(lapply(enr, function(x) rbindlist(x)))
enr$M <- rep(1:15, each = length(unique(unique(tf.rel$TF)))*2)
enr$tail <- rep(rep(c('-', '+'), each = length(unique(unique(tf.rel$TF)))),
             15)

enr.sgn <- enr[(p<0.001 & f>2)]
enr.sgn <- enr.sgn[pres>10]

#
tf.genes.sgnf <- unique(enr.sgn$CAT)
sink('data_processed/tf.genes.sgnf.txt')
cat(tf.genes.sgnf)
sink()
# this list was processed at http://www.yeastract.com/orftogene.php
# to give us the ORF names for the TFs:

orf_to_tf <- fread('data/orftogene_equivalence_signif_tfs.csv',
                   header=F)[,1:3]
names(orf_to_tf) <- c('CAT', 'ORF', 'Gene')
enr.sgn <- merge(enr.sgn, orf_to_tf, by='CAT', all=T)

# now retrieve the go mapping for the tfs from DAVID
sink('data/tf_orf_for_david.txt')
cat(orf_to_tf$ORF)
sink()

# cluster TFs according to their similarity in GO Slim BP belonging profile
GO.tf <- fread('data/signif_tf_GOBP.txt') # from DAVID
GO.cat <- as.list(GO.tf$GOTERM_BP_DIRECT)

process.go.string <- function(x)
{
    x <- unlist(strsplit(x, 'GO:'))
    x <- sapply(x[grep('~', x)], function(x) gsub(',', '', x))
    names(x) <- NULL
    x
}

GO.cat <- lapply(GO.cat, process.go.string)
GO.cat <- data.table(TF = rep(orf_to_tf[match(GO.tf$ID, orf_to_tf$ORF)]$CAT,
                              sapply(GO.cat, length)),
                     GO = unlist(GO.cat))
GO.cat$val <- 1
go.mat <- dcast(GO.cat, TF~GO)
go.mat <- as.matrix(go.mat[,-1], rownames = go.mat$TF)
go.mat[is.na(go.mat)] <- 0

tf.go.clust <- hclust(dist(go.mat))
hcd <- as.dendrogram(tf.go.clust)

pdf(paste0('results/TF_dendrogram.pdf'), 2, 10)
ggdendrogram(tf.go.clust, rotate = TRUE)
dev.off()

# order of CAT according to clustering
enr.sgn[, CAT := ordered(CAT, tf.go.clust$labels[tf.go.clust$order])]

# .. Matrix with fold changes of each TF in each mode, upregulated
p <- ggplot(enr.sgn[f>3]) +
    geom_tile(aes(x = tail, y = CAT, fill = f)) +
    scale_fill_viridis_c(na.value = NA) +
    scale_x_discrete(limits = c('-', '+')) +
    facet_grid(~ factor(as.numeric(M)),
               space='free_x', scales='free_x', switch='x', drop = FALSE) +
    theme_bw() +
    theme(strip.placement = 'outside',
          strip.background = element_rect(fill=NA,colour='grey50'),
          panel.spacing=unit(0,'cm'),
          legend.position = 'top',
          axis.text.y = element_text(size=8)) +
    xlab('Mode')

pdf(paste0('results/TF_enrichment.pdf'), 5, 9)
print(p)
dev.off()


