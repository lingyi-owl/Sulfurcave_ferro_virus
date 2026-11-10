# =============================================================================
# Fig. 6 — tail fiber protein physicochemistry (a, c, d violins; b amino-acid radar)
# Figure numbering follows the manuscript (Cave_Ferro_Virus, Nov 2026).
# =============================================================================

# =====================================================================
# Fig 6 (paper layout) with panel c swapped to aromatic content F+W+Y+H.
#   a = isoelectric point            (unchanged)
#   b = amino-acid usage radar       (unchanged; Fig6_tail_fiber_features.R)
#   c = Percentage of F+W+Y+H (%)    (was F+W+Y; now includes histidine)
#   d = secondary-structure fraction (unchanged)
#
# All violin panels: paired Cave (#e7298a) vs Reference / IMG-VR (#0072B2)
# distributions of 1000 subsample means; each iteration draws
# n = floor(0.90 * cave_count) WITHOUT replacement from the cave and the
# reference (IMG/VR 70% cluster-rep) tail-fiber-protein pools (set.seed(42)).
#
# Panels a, c, d read the precomputed per-iteration means table
# (fig6_paired_means_acd.csv); the FWYH series was computed with the identical
# scheme as the other features (see the "FWYH" feature). Panel b is built by
# Fig6_tail_fiber_features.R. Output: figures/main/Fig6/panels/Fig6_acd.svg
# =====================================================================
# ---- Paths (see README, "Running the scripts") ----
REPO     <- Sys.getenv("SULFURCAVE_REPO", unset = here::here())  # repository root

suppressWarnings(suppressMessages({library(data.table); library(ggplot2); library(ggbeeswarm); library(patchwork)}))

PLOT <- file.path(REPO, "plotting")
FIGS <- file.path(REPO, "figures/main/Fig6/panels")
dir.create(FIGS, recursive = TRUE, showWarnings = FALSE)
setwd(PLOT)

# ---- panel b: radar (unchanged) -> defines `radar_plot` -------------
# (inlined from the former Fig6_tail_fiber_features.R; reads fig6_aa_profile.csv + viromics diff CSV; builds radar_plot)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>> panel b >>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Panel d (paired): dual cave-vs-reference AA-usage radar, WITH the original
# panel-d outer layers — biochemical-class arcs + curved class labels + the
# one-tailed significance tier ring. Only the central series differs from the
# delta radar: two usage profiles (Cave, Reference) instead of one delta line.
# ---- Paths (see README, "Running the scripts") ----
REPO     <- Sys.getenv("SULFURCAVE_REPO", unset = here::here())  # repository root
FIGOUT_FIG6 <- file.path(REPO, "figures/main/Fig6/panels")  # figure output folder (manuscript numbering)
dir.create(FIGOUT_FIG6, recursive = TRUE, showWarnings = FALSE)

suppressMessages({library(dplyr); library(stringr); library(ggplot2); library(data.table)})
if (!requireNamespace("ggradar", quietly=TRUE))
  remotes::install_github("ricardo-bion/ggradar", upgrade="never", quiet=TRUE)
if (!requireNamespace("geomtextpath", quietly=TRUE))
  install.packages("geomtextpath", repos="https://cloud.r-project.org", quiet=TRUE)
suppressMessages(library(ggradar))

BASE_FAMILY <- "Arial"; BASE_SIZE <- 26; ABC_SIZE <- 36
ABC_GEOM <- ABC_SIZE/ggplot2::.pt; P_RES <- 1/1000
oc <- c(Cave="#e7298a", Reference="#808080")

# ---- significance (per-AA one-tailed p from the delta analysis) -------------
significance_df_all <- read.csv2(file.path(REPO, "viromics/sulfurcave_imgvr_comparisons/sulfurcave_imgvr_features_diff.csv"), sep=',')
significance_df <- significance_df_all[significance_df_all$protein=='tail_fiber_protein', ] %>%
  mutate(feature = str_remove(feature, "_Percent"))

# ---- AA usage profiles (cave vs reference, %) -------------------------------
aa_order <- c("A","V","L","I","M","F","W","Y","H","R","K","D","E","S","T","N","Q","G","P","C")
prof <- fread("fig6_aa_profile.csv"); setDF(prof)
prof <- prof[match(aa_order, prof$aa), ]
radar_raw <- as.data.frame(setNames(as.list(prof$cave), aa_order))   # colnames carrier

aa_groups <- list(
  Aliphatic=list(aas=c("A","V","L","I","M"),color="#117733"),
  Aromatic =list(aas=c("F","W","Y","H"),    color="#332288"),
  Positive =list(aas=c("R","K"),            color="#88CCEE"),
  Negative =list(aas=c("D","E"),            color="#AA4499"),
  Polar    =list(aas=c("S","T","N","Q"),    color="#44AA99"),
  Special  =list(aas=c("G","P","C"),        color="#DDCC77"))
group_colors <- unlist(lapply(names(aa_groups), function(g){
  present <- intersect(aa_groups[[g]]$aas, colnames(radar_raw)); rep(aa_groups[[g]]$color, length(present))}))

gmax <- ceiling(max(prof$cave, prof$reference))          # ~13
resc <- function(v) v/gmax
radar_rescaled <- data.frame(Group=c("Cave","Reference"),
  rbind(resc(prof$cave), resc(prof$reference)))
colnames(radar_rescaled)[-1] <- aa_order

p <- ggradar(radar_rescaled,
  values.radar=c("0%", paste0(gmax/2,"%"), paste0(gmax,"%")),
  grid.min=0, grid.mid=0.5, grid.max=1.0,
  background.circle.colour="#F9F9F9",
  gridline.min.colour="grey80", gridline.mid.colour="grey80", gridline.max.colour="grey80",
  gridline.min.linetype="solid", gridline.mid.linetype="solid", gridline.max.linetype="solid",
  grid.line.width=0.4,
  group.line.width=2.0, group.point.size=3.5, fill=TRUE, fill.alpha=0.10,
  group.colours=unname(oc[c("Cave","Reference")]),
  axis.label.size=ABC_GEOM, grid.label.size=ABC_GEOM, legend.text.size=ABC_SIZE,
  axis.label.offset=1.62, font.radar=BASE_FAMILY, legend.title="Origin") +
  theme(text=element_text(family=BASE_FAMILY, size=ABC_SIZE, colour="black"))

# AA letters -> black (class colour is carried by the arcs, not the letters)
color_lookup <- setNames(group_colors, colnames(radar_raw))
for (i in seq_along(p$layers)) {
  layer <- p$layers[[i]]; ldata <- layer$data
  if (inherits(layer$geom,"GeomText") && !is.null(ldata) && "text" %in% colnames(ldata) &&
      any(ldata$text %in% names(color_lookup))) {
    if (!any(is.na(color_lookup[ldata$text]))) {
      p$layers[[i]]$aes_params$colour <- "black"; p$layers[[i]]$aes_params$fontface <- "plain" }
  }
}

# ---- geometry of the outer layers (relative to ggradar grid radius) ---------
axis_layer_idx <- which(sapply(p$layers, function(l) !is.null(l$data) && "axis.no" %in% colnames(l$data)))[1]
axis_outer <- p$layers[[axis_layer_idx]]$data %>%
  mutate(dist=sqrt(x^2+y^2)) %>% group_by(axis.no) %>% slice_max(dist,n=1) %>% ungroup() %>%
  arrange(axis.no) %>% mutate(aa=colnames(radar_raw))
outer_r <- max(axis_outer$dist); bracket_r <- outer_r*1.52; label_r <- outer_r*1.70

group_membership <- bind_rows(lapply(names(aa_groups), function(g){
  present <- intersect(aa_groups[[g]]$aas, colnames(radar_raw)); if(length(present)==0) return(NULL)
  data.frame(aa=present, group_name=g, color=aa_groups[[g]]$color, stringsAsFactors=FALSE)}))
group_labels <- axis_outer %>% left_join(group_membership, by="aa", relationship="many-to-many") %>%
  group_by(group_name,color) %>% summarise(cx=mean(x),cy=mean(y),n_aa=n(),.groups="drop")

# class arcs
arc_data <- bind_rows(lapply(seq_len(nrow(group_labels)), function(i){
  g<-group_labels$group_name[i]; ga<-axis_outer %>% filter(aa %in% aa_groups[[g]]$aas)
  if(nrow(ga)<2) return(NULL); th<-atan2(ga$x,ga$y)
  if(max(th)-min(th)>pi) th<-ifelse(th<0,th+2*pi,th); ts<-seq(min(th),max(th),length.out=50)
  data.frame(x=bracket_r*sin(ts), y=bracket_r*cos(ts), group=g, color=group_labels$color[i])}))
p <- p + annotate("path", x=arc_data$x, y=arc_data$y, group=arc_data$group,
                  colour=arc_data$color, linewidth=2.6)
# curved class labels
label_arc_r <- outer_r*1.40
label_arc_data <- bind_rows(lapply(seq_len(nrow(group_labels)), function(i){
  g<-group_labels$group_name[i]; ga<-axis_outer %>% filter(aa %in% aa_groups[[g]]$aas)
  if(nrow(ga)<2) return(NULL); th<-atan2(ga$x,ga$y)
  if(max(th)-min(th)>pi) th<-ifelse(th<0,th+2*pi,th); ts<-seq(min(th),max(th),length.out=60)
  data.frame(x=label_arc_r*sin(ts), y=label_arc_r*cos(ts), group=g, label=g, color=group_labels$color[i])}))
p <- p + geomtextpath::geom_textpath(data=label_arc_data,
  aes(x=x,y=y,group=group,label=label), colour=label_arc_data$color,
  text_only=TRUE, upright=TRUE, size=ABC_GEOM, family=BASE_FAMILY, fontface="plain", inherit.aes=FALSE)

# ---- one-tailed significance tier ring -------------------------------------
tier_levels <- c("p < 0.001","p < 0.01","p < 0.05","n.s. (p \u2265 0.05)")
tier_colors <- c("p < 0.001"="#99000D","p < 0.01"="#EF6548","p < 0.05"="#FDBB84","n.s. (p \u2265 0.05)"="#D9D9D9")
sig_ring <- axis_outer %>% left_join(significance_df, by=c("aa"="feature")) %>%
  mutate(p_raw=as.numeric(p_one_tailed),
         tier=cut(p_raw, breaks=c(-Inf,0.001,0.01,0.05,Inf), labels=tier_levels, right=FALSE),
         tier=factor(tier, levels=tier_levels), theta=atan2(x,y))
n_aa<-nrow(sig_ring); gap_frac<-0.10; half_w<-(pi/n_aa)*(1-gap_frac)
bar_base<-outer_r*1.14; ring_thk<-outer_r*0.14
tile_poly <- bind_rows(lapply(seq_len(n_aa), function(i){
  th<-sig_ring$theta[i]; r_out<-bar_base+ring_thk; arc<-seq(th-half_w,th+half_w,length.out=20)
  data.frame(x=c(bar_base*sin(arc), r_out*sin(rev(arc))),
             y=c(bar_base*cos(arc), r_out*cos(rev(arc))),
             tile=sig_ring$aa[i], tier=sig_ring$tier[i])}))
base_ring <- data.frame(x=bar_base*sin(seq(0,2*pi,length.out=200)), y=bar_base*cos(seq(0,2*pi,length.out=200)))
plot_limit <- max(abs(c(tile_poly$x,tile_poly$y)), bracket_r)*1.12

radar_plot <- p +
  geom_path(data=base_ring, aes(x=x,y=y), colour="grey70", linewidth=0.4, inherit.aes=FALSE) +
  geom_polygon(data=tile_poly, aes(x=x,y=y,group=tile,fill=tier), colour="white", linewidth=0.25) +
  scale_fill_manual(values=tier_colors, breaks=tier_levels, limits=tier_levels, drop=FALSE,
    name="One-tailed\nsignificance",
    guide=guide_legend(override.aes=list(colour="white"), title.position="top", title.hjust=0,
      keywidth=grid::unit(1.4,"lines"), keyheight=grid::unit(1.4,"lines"), order=2)) +
  guides(colour=guide_legend(order=1, title="Origin")) +
  scale_x_continuous(limits=c(-plot_limit,plot_limit)) +
  scale_y_continuous(limits=c(-plot_limit,plot_limit)) +
  coord_fixed(clip="off") +
  theme(legend.position="right",
        legend.title=element_text(size=ABC_SIZE, family=BASE_FAMILY, colour="black"),
        legend.text =element_text(size=ABC_SIZE, family=BASE_FAMILY, colour="black"))

ggsave(file.path(FIGOUT_FIG6, "Fig6b_aa_radar.png"), plot=radar_plot, width=15, height=13, dpi=110, bg="white")
ggsave(file.path(FIGOUT_FIG6, "Fig6b_aa_radar.pdf"), plot=radar_plot, width=15, height=13, bg="white", device=cairo_pdf)
if (!requireNamespace("svglite", quietly=TRUE)) install.packages("svglite", repos="https://cloud.r-project.org", quiet=TRUE)
ggsave(file.path(FIGOUT_FIG6, "Fig6b_aa_radar.svg"), plot=radar_plot, width=15, height=13, bg="white")
cat("gmax", gmax, "done\n")

# <<<<<<<<<<<<<<<<<<<<<<<<<<<< panel b <<<<<<<<<<<<<<<<<<<<<<<<<<<<
setwd(PLOT)


# ---- shared violin style --------------------------------------------
BASE_FAMILY <- "Arial"; BASE_SIZE <- 20
oc  <- c(Cave="#e7298a", Reference="#0072B2")
cap <- function(x) factor(ifelse(x=="cave","Cave","Reference"), levels=c("Cave","Reference"))
base_theme <- theme_minimal(base_size=BASE_SIZE, base_family=BASE_FAMILY) +
  theme(text=element_text(family=BASE_FAMILY,size=BASE_SIZE,colour="black"),
        axis.text=element_text(family=BASE_FAMILY,size=BASE_SIZE,colour="black"),
        axis.title=element_text(family=BASE_FAMILY,size=BASE_SIZE,colour="black"),
        axis.line=element_line(colour="black",linewidth=0.6),
        panel.grid.major=element_line(color="grey92",linewidth=0.5),
        panel.grid.minor=element_line(color="grey92",linewidth=0.5),
        panel.border=element_rect(color="grey85",fill=NA,linewidth=0.8),
        legend.position="none", plot.margin=margin(t=14,r=10,b=8,l=14))

pm <- fread("fig6_paired_means_acd.csv")
one_tailed <- function(dt){ d <- dt$cave - dt$reference; if (mean(d) > 0) mean(d <= 0) else mean(d >= 0) }
fmt_p <- function(p) if (p < 1/1000) "p < 0.001" else sprintf("p = %.3f", p)

mk_single <- function(dt, ytitle){
  plab <- fmt_p(one_tailed(dt))
  long <- melt(dt, id.vars="iter", measure.vars=c("cave","reference"),
               variable.name="origin", value.name="m"); long[, origin:=cap(origin)]
  ggplot(long, aes(origin, m, fill=origin)) +
    geom_violin(alpha=0.35, trim=FALSE, colour=NA, scale="width") +
    ggbeeswarm::geom_quasirandom(aes(colour=origin), size=0.6, alpha=0.4, width=0.3) +
    annotate("text", x=1.5, y=max(long$m), label=plab, family=BASE_FAMILY,
             size=BASE_SIZE/ggplot2::.pt, vjust=-0.8) +
    scale_fill_manual(values=oc) + scale_colour_manual(values=oc) +
    scale_y_continuous(expand=expansion(mult=c(0.05,0.16))) +
    labs(x=NULL, y=ytitle) + base_theme
}

pa <- mk_single(pm[feature=="IsoelectricPoint"], "Isoelectric point (pH units)")
pc <- mk_single(pm[feature=="FWYH"],             "Percentage of F+W+Y+H (%)")   # swapped panel

# ---- panel d: secondary structure (unchanged) ----------------------
dc <- pm[feature %in% c("HelixFraction","SheetFraction","TurnFraction")]
dc[, `:=`(cave=cave*100, reference=reference*100)]
dc[, structure := factor(sub("Fraction$","",feature), levels=c("Helix","Sheet","Turn"))]
lc <- melt(dc, id.vars=c("structure","iter"), measure.vars=c("cave","reference"),
           variable.name="origin", value.name="m"); lc[, origin:=cap(origin)]
lev <- as.vector(t(outer(levels(dc$structure), c("Cave","Reference"),
                         function(s,o) paste(o,s,sep="."))))
lc[, xg := factor(interaction(origin, structure), levels=lev)]
sc <- data.table(structure=c("Helix","Sheet","Turn"), xc=c(1.5,3.5,5.5))
sc$ytop <- sapply(sc$structure, function(s) max(lc[structure==s]$m))
sc$plab <- sapply(sc$structure, function(s) fmt_p(one_tailed(dc[structure==s])))
pd <- ggplot(lc, aes(xg, m, fill=origin)) +
  geom_violin(alpha=0.35, trim=FALSE, colour=NA, scale="width") +
  ggbeeswarm::geom_quasirandom(aes(colour=origin), size=0.6, alpha=0.4, width=0.3) +
  geom_text(data=sc, aes(x=xc, y=ytop, label=plab), inherit.aes=FALSE,
            family=BASE_FAMILY, size=BASE_SIZE/ggplot2::.pt, vjust=-0.8) +
  annotate("text", x=c(1.5,3.5,5.5), y=-Inf, label=c("Helix","Sheet","Turn"),
           vjust=3.4, family=BASE_FAMILY, size=BASE_SIZE/ggplot2::.pt) +
  scale_fill_manual(values=oc) + scale_colour_manual(values=oc) +
  scale_x_discrete(labels=rep(c("Cave","Reference"),3)) +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.16))) +
  coord_cartesian(clip="off") +
  labs(x=NULL, y="Secondary structure fraction (%)") +
  base_theme + theme(plot.margin=margin(t=14,r=10,b=48,l=14))

# ---- assemble: a top-left, b(radar) top-right, c mid-left, d bottom -
fig6 <- pa + radar_plot + pc + pd +
  plot_layout(design="AB\nCB\nDD", widths=c(1,2), heights=c(1,1,1.2)) +
  plot_annotation(tag_levels="a") &
  theme(plot.tag=element_text(size=BASE_SIZE+4, family=BASE_FAMILY))

ggsave(file.path(FIGS,"Fig6_acd.svg"), plot=fig6, device=svglite::svglite,
       width=18, height=16, bg="white", limitsize=FALSE)
