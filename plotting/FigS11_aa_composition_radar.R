# =====================================================================
# AA composition PAIRED radar panels (one per protein type), styled
# after Fig 6 panel b.
#
# Per panel: two polygons of ABSOLUTE mean amino-acid % -
#   Cave (magenta #e7298a) and IMG/VR (blue #0072B2) - on a shared 0..M%
#   scale (rings at 0 / M/2 / M %). Amino acids ordered by biochemical
#   group around the circle.
# Radial layout (centre -> out):
#   polygons/grid -> significance heatmap tiles -> property-group labels
#   -> thick property arcs -> amino-acid letters (outermost).
# Significance heatmap: 20 arc tiles coloured by one-tailed permutation p
#   (red/orange scale + grey n.s.). Property arcs use the Paul-Tol palette.
# Text: Arial, size 20, regular. One shared legend.
#
# INPUTS
#   ../viromics/cave_imgvr_features_raw.csv - per-protein features for both
#     groups (cave = all proteins; IMG/VR = seeded 3000/type subsample),
#     recomputed with Biopython ProtParam. *_Percent columns are stored as
#     (fraction x 10000); divide by 100 for true %.
#
# SIGNIFICANCE: computed inline with the EXACT permutation method of
#   FigS9_physiochemical_features_paired (R, set.seed(42), 1000 iterations,
#   n = floor(0.90 * cave_count) drawn WITH replacement from BOTH the cave
#   and reference pools, statistic = per-iteration difference of subsample
#   means, p_one_tailed = fraction of iterations on the side OPPOSITE the
#   observed mean difference). Written out to aa_permutation_p.csv.
# =====================================================================
# ---- Paths (see README, "Running the scripts") ----
REPO     <- Sys.getenv("SULFURCAVE_REPO", unset = here::here())  # repository root
FIGOUT_S11 <- file.path(REPO, "figures/supplementary/FigS11")  # figure output folder (manuscript numbering)
dir.create(FIGOUT_S11, recursive = TRUE, showWarnings = FALSE)

suppressWarnings(suppressMessages({
  library(tidyverse); library(ggradar); library(patchwork); library(cowplot)
}))

RAW <- file.path(REPO, "viromics/cave_imgvr_features_raw.csv")
BASE_FONT <- "Arial"
GEOM_SZ   <- 20/.pt   # ~20 pt expressed in geom_text/ggradar size units

aa_order <- c("A","V","L","I","M","F","W","Y","H","R","K","D","E","S","T","N","Q","G","P","C")
aa_groups <- list(
  Aliphatic=list(aas=c("A","V","L","I","M"),color="#117733"),
  Aromatic =list(aas=c("F","W","Y","H"),    color="#332288"),
  Positive =list(aas=c("R","K"),            color="#88CCEE"),
  Negative =list(aas=c("D","E"),            color="#AA4499"),
  Polar    =list(aas=c("S","T","N","Q"),    color="#44AA99"),
  Special  =list(aas=c("G","P","C"),        color="#DDCC77"))
sig_pal    <- c("p<0.001"="#99000D","p<0.01"="#EF6548","p<0.05"="#FDBB84","n.s."="#7F7F7F")
sig_labels <- c("p < 0.001","p < 0.01","p < 0.05","n.s.")

proteins <- c("major_head_protein","minor_head_protein","major_tail_protein","minor_tail_protein",
              "tail_fiber_protein","tail_spike_protein","dna_polymerase","endolysin",
              "head_maturation_protease","integrase","terminase_large_subunit","terminase_small_subunit")
nice <- function(x){s<-gsub("_"," ",x);s<-paste0(toupper(substr(s,1,1)),substr(s,2,nchar(s)));sub("^Dna ","DNA ",s)}

# ---- data ----------------------------------------------------------
feat <- read.csv(RAW)
absm <- feat %>% group_by(origin,protein) %>%
  summarise(across(ends_with("_Percent"), ~mean(.x,na.rm=TRUE)), .groups="drop")
getv <- function(pr,org){ r<-absm[absm$protein==pr & absm$origin==org,]
  setNames(as.numeric(r[paste0(aa_order,"_Percent")]),aa_order)/100 }   # /100 -> true %
Mabs <- ceiling(max(sapply(proteins,function(pr) max(getv(pr,"sulfurcave"),getv(pr,"imgvr")))))

# ---- significance: EXACT FigS9 permutation method ------------------
set.seed(42); N_ITER <- 1000; FRAC <- 0.90
resample_means <- function(cave, ref){                 # identical to FigS9
  n <- floor(FRAC * length(cave))
  data.frame(cave      = vapply(seq_len(N_ITER), function(i) mean(sample(cave, n, replace = TRUE)), numeric(1)),
             reference = vapply(seq_len(N_ITER), function(i) mean(sample(ref,  n, replace = TRUE)), numeric(1)))
}
sig_recs <- list()
for (aa in aa_order) for (pr in proteins){
  col <- paste0(aa, "_Percent")
  cv <- feat[feat$origin=="sulfurcave" & feat$protein==pr, col]; cv <- cv[is.finite(cv)]
  rf <- feat[feat$origin=="imgvr"      & feat$protein==pr, col]; rf <- rf[is.finite(rf)]
  if (length(cv) < 2 || length(rf) < 2) next
  m <- resample_means(cv, rf); d <- m$cave - m$reference; obs <- mean(d)
  p1 <- if (obs > 0) mean(d <= 0) else mean(d >= 0)
  sig_recs[[paste(aa, pr)]] <- data.frame(protein = pr, aa = aa, p_one_tailed = p1)
}
sig_all <- do.call(rbind, sig_recs)
getsig <- function(pr){ s<-sig_all[sig_all$protein==pr,]
  setNames(as.numeric(s$p_one_tailed), s$aa)[aa_order] }
siglev <- function(p) ifelse(is.na(p),"n.s.",ifelse(p<0.001,"p<0.001",ifelse(p<0.01,"p<0.01",ifelse(p<0.05,"p<0.05","n.s."))))

# ---- one radar panel ----------------------------------------------
build_radar2 <- function(pr, M, title, show_legend=FALSE){
  cv<-getv(pr,"sulfurcave"); iv<-getv(pr,"imgvr"); sg<-getsig(pr)
  rdf <- data.frame(Group=c("Cave","IMG/VR"), rbind(cv/M, iv/M), check.names=FALSE); colnames(rdf)[-1]<-aa_order
  p <- ggradar(rdf, values.radar=c("0%",paste0(M/2,"%"),paste0(M,"%")),
    grid.min=0, grid.mid=0.5, grid.max=1, group.colours=c("#e7298a","#0072B2"),
    group.line.width=0.8, group.point.size=1.3, background.circle.colour="#F9F9F9",
    gridline.min.colour="grey85", gridline.mid.colour="grey85", gridline.max.colour="grey85",
    gridline.mid.linetype="solid",
    axis.label.size=GEOM_SZ, grid.label.size=GEOM_SZ*0.75, axis.label.offset=1.60) +
    ggtitle(title) + labs(colour="Origin") +
    theme(text=element_text(family=BASE_FONT, size=20),
          plot.title=element_text(family=BASE_FONT, face="plain", size=20, hjust=0.5),
          legend.position=if(show_legend)"bottom" else "none")
  ali <- which(sapply(p$layers,function(l)!is.null(l$data)&&"axis.no"%in%colnames(l$data)))[1]
  ao <- p$layers[[ali]]$data %>% mutate(dist=sqrt(x^2+y^2)) %>% group_by(axis.no) %>%
        slice_max(dist,n=1) %>% ungroup() %>% arrange(axis.no) %>% mutate(aa=aa_order)
  outer_r<-max(ao$dist); half<-pi/20
  ti<-outer_r*1.09; to<-outer_r*1.24
  tiles<-bind_rows(lapply(seq_len(nrow(ao)),function(i){
    th<-atan2(ao$x[i],ao$y[i]); ts<-seq(th-half*0.86,th+half*0.86,length.out=14)
    data.frame(x=c(to*sin(ts),rev(ti*sin(ts))), y=c(to*cos(ts),rev(ti*cos(ts))),
               tile=i, Significance=siglev(sg[ao$aa[i]]))}))
  lab_r <- outer_r*1.31
  br_in <- outer_r*1.37; br_out <- outer_r*1.45
  p <- p + geom_polygon(data=tiles, aes(x=x,y=y,group=tile,fill=Significance), colour=NA, inherit.aes=FALSE) +
    scale_fill_manual(name="One-tailed significance", values=sig_pal,
      limits=c("p<0.001","p<0.01","p<0.05","n.s."), drop=FALSE, labels=sig_labels)
  for(g in names(aa_groups)){
    ga<-ao%>%filter(aa%in%aa_groups[[g]]$aas); th<-atan2(ga$x,ga$y)
    if(max(th)-min(th)>pi) th<-ifelse(th<0,th+2*pi,th)
    ts<-seq(min(th)-half*0.55,max(th)+half*0.55,length.out=80)
    band<-data.frame(x=c(br_out*sin(ts),rev(br_in*sin(ts))), y=c(br_out*cos(ts),rev(br_in*cos(ts))))
    p <- p + annotate("polygon", x=band$x, y=band$y, fill=aa_groups[[g]]$color, colour=NA)
    thc<-mean(th); phi<-(thc*180/pi)%%360
    rot<- -phi; if(phi>135 && phi<225) rot<-rot+180
    rot<-((rot+180)%%360)-180
    p <- p + annotate("text", x=lab_r*sin(thc), y=lab_r*cos(thc), label=g,
                      colour=aa_groups[[g]]$color, angle=rot, fontface="plain",
                      family=BASE_FONT, size=GEOM_SZ*0.8)
  }
  lim<-outer_r*1.60*1.14
  p + scale_x_continuous(limits=c(-lim,lim)) + scale_y_continuous(limits=c(-lim,lim)) + coord_fixed(clip="off")
}

# ---- assemble 4x3 grid + one shared legend -------------------------
leg_theme <- theme(legend.position="bottom", legend.box="horizontal",
                   legend.title=element_text(family=BASE_FONT, face="plain", size=20),
                   legend.text =element_text(family=BASE_FONT, face="plain", size=20))
donor <- build_radar2(proteins[3], Mabs, "", show_legend=TRUE) + leg_theme
leg   <- cowplot::get_legend(donor)
panels_nl <- lapply(proteins, function(pr) build_radar2(pr, Mabs, nice(pr), show_legend=FALSE))
body_grid <- wrap_plots(panels_nl, ncol=4) + plot_annotation(tag_levels="a") &
  theme(plot.tag=element_text(family=BASE_FONT, face="plain", size=20))
final <- cowplot::plot_grid(body_grid, leg, ncol=1, rel_heights=c(1, 0.05))

ggsave(file.path(FIGOUT_S11, "FigS11_aa_composition_radar.svg"),
       final, device=svglite::svglite, width=28, height=25, limitsize=FALSE)
