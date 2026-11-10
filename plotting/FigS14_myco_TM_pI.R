# Mycobacterium transmembrane protein isoelectric point:
# paired cave-vs-reference resampled pI distributions per topology.
# Same downsampling as the delta figure (90% cave count, 1000 iterations,
# WITH replacement), but records each iteration's cave mean and reference mean
# rather than only their difference. Three topology pairs of violins.

# ---- Paths (see README, "Running the scripts") ----
REPO     <- Sys.getenv("SULFURCAVE_REPO", unset = here::here())  # repository root
FIGOUT_S14 <- file.path(REPO, "figures/supplementary/FigS14")  # figure output folder (manuscript numbering)
dir.create(FIGOUT_S14, recursive = TRUE, showWarnings = FALSE)
EXTERNAL <- Sys.getenv("SULFURCAVE_EXTERNAL", unset = file.path(REPO, "external_data"))  # inputs not shipped in the repo

suppressMessages({library(data.table); library(ggplot2); library(ggbeeswarm)})
base <- file.path(EXTERNAL, "mycobacterium_pangenome")

tm <- fread(file.path(base,"merged_TM_protein_properties.tsv"), sep="\t",
            select=c("id","description","isoelectric_point"))
og <- fread(file.path(base,"MYCO_Pan_gene_clusters_summary.txt"), sep="\t",
            select=c("unique_id","genome_name")); setnames(og,"unique_id","id")
myco <- merge(tm, og, by="id", all.x=TRUE)
myco[, topology := fifelse(grepl("I_",description,fixed=TRUE),"Intracellular",
                   fifelse(grepl("O_",description,fixed=TRUE),"Extracellular",
                   fifelse(grepl("M_",description,fixed=TRUE),"Membrane-spanning","Unknown")))]
myco[, origin := fifelse(genome_name %in% c("M_Bin2","M_Bin3","M_methanotrophicum"),"cave","reference")]
myco[, isoelectric_point := as.numeric(isoelectric_point)]
myco <- myco[!is.na(isoelectric_point) & topology != "Unknown"]

set.seed(42); N_ITER <- 1000; FRAC <- 0.90; P_RES <- 1/N_ITER
topos <- c("Intracellular","Membrane-spanning","Extracellular")
resample_means <- function(dt) {
  cave <- dt[origin=="cave", isoelectric_point]; ref <- dt[origin=="reference", isoelectric_point]
  n <- floor(FRAC*length(cave))
  data.table(cave      = vapply(seq_len(N_ITER), function(i) mean(sample(cave,n,replace=TRUE)), numeric(1)),
             reference = vapply(seq_len(N_ITER), function(i) mean(sample(ref, n,replace=TRUE)), numeric(1)))
}
means_df <- rbindlist(lapply(topos, function(tp){
  m <- resample_means(myco[topology==tp]); m[, `:=`(topology=tp, iter=seq_len(N_ITER))]; m}))
long <- melt(means_df, id.vars=c("topology","iter"),
             measure.vars=c("cave","reference"), variable.name="origin", value.name="mean_pI")
long[, topology := factor(topology, levels=topos)]
long[, origin := factor(origin, levels=c("cave","reference"))]

# one-tailed p from per-iteration delta (cave - reference)
fmt_p <- function(p) ifelse(p < P_RES, "p < 0.001", sprintf("p = %.3f", p))
sig <- means_df[, {d<-cave-reference; obs<-mean(d); frac<- if(obs>0) mean(d<=0) else mean(d>=0)
                   .(p_lab=fmt_p(frac), ytop=max(c(cave,reference)))}, by=topology]
sig[, topology := factor(topology, levels=topos)]

BASE_FAMILY <- "Arial"; BASE_SIZE <- 20
origin_colors <- c(Cave="#e7298a", Reference="#0072B2")

# two-tier x axis (Cave/Reference on top, topology below) — Fig6 panel c format
long[, origin2 := factor(ifelse(origin=="cave","Cave","Reference"), levels=c("Cave","Reference"))]
lev <- as.vector(t(outer(topos, c("Cave","Reference"), function(t,o) paste(o,t,sep="."))))
long[, xg := factor(interaction(origin2, topology, sep="."), levels=lev)]
sig[, xc := c(1.5,3.5,5.5)[match(as.character(topology), topos)]]

p <- ggplot(long, aes(x=xg, y=mean_pI, fill=origin2)) +
  geom_violin(alpha=0.35, trim=FALSE, colour=NA, scale="width") +
  ggbeeswarm::geom_quasirandom(aes(colour=origin2), size=0.7, alpha=0.4, width=0.3) +
  geom_text(data=sig, aes(x=xc, y=ytop, label=p_lab), inherit.aes=FALSE,
            family=BASE_FAMILY, size=BASE_SIZE/ggplot2::.pt, vjust=-0.8, colour="black") +
  annotate("text", x=c(1.5,3.5,5.5), y=-Inf, label=topos,
           vjust=3.4, family=BASE_FAMILY, size=BASE_SIZE/ggplot2::.pt) +
  scale_fill_manual(values=origin_colors) +
  scale_colour_manual(values=origin_colors) +
  scale_x_discrete(labels=rep(c("Cave","Reference"),3)) +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.14))) +
  coord_cartesian(clip="off") +
  labs(x=NULL, y="Isoelectric point (pH units)") +
  theme_minimal(base_size=BASE_SIZE, base_family=BASE_FAMILY) +
  theme(text=element_text(family=BASE_FAMILY,size=BASE_SIZE,colour="black"),
        axis.text=element_text(family=BASE_FAMILY,size=BASE_SIZE,colour="black"),
        axis.title=element_text(family=BASE_FAMILY,size=BASE_SIZE,colour="black"),
        axis.line=element_line(colour="black",linewidth=0.6),
        panel.grid.major=element_line(color="grey92",linewidth=0.5),
        panel.grid.minor=element_line(color="grey92",linewidth=0.5),
        panel.border=element_rect(color="grey85",fill=NA,linewidth=0.8),
        legend.position="none", plot.margin=margin(t=14,r=10,b=48,l=14))

ggsave(file.path(FIGOUT_S14, "FigS14_myco_TM_pI_paired_violin.svg"), plot=p, width=12, height=7, bg="white")
