
source(".init.R")
load("data/ttc_pvar.Rdata")

load("~/work_data/holstege_data/M.Rdata")
pvar.obs <- get.pov(svd(M))


ttc.stats <- data.table('mode' = 1:ncol(ttc.pvar),
                        'ttc.mean' = apply(ttc.pvar,2,mean),
                        'ttc.sd' = apply(ttc.pvar,2,sd),
                        'pvar.obs' = pvar.obs[1:1000])


p <- "results_new/fig_2B.pdf"
pdf(p, 4, 4, family="CM Sans")
ggplot(ttc.stats[1:100], aes(x=mode, y=ttc.mean)) +
geom_point(colour = 'gray30', shape=95) +
geom_linerange(aes(ymin=ttc.mean-(2*ttc.sd), ymax=ttc.mean+(2*ttc.sd))) +
geom_point(aes(y=pvar.obs), colour='#008066', shape=1) +
geom_vline(xintercept=32) +
scale_y_log10() + annotation_logticks(sides = "l") +
theme_bw() + theme(panel.grid=element_blank())
dev.off()
embed_fonts(p, outfile=p)
