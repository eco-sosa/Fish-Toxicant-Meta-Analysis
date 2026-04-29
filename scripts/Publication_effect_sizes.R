#======== Brandon Sosa =========== Updated 2025
#======== Publication-Quality Meta-Analysis Figures =========

library(ggplot2)
library(cowplot)

#===================== Hardcoded Data ===============================================

# -------- INCLUSIVE APPROACH (Borenstein-corrected for shared controls) --------
# Updated with Borenstein correction per Song et al. (2022)
# r calculated from sample sizes using Gleser & Olkin formula

# Size Inclusive (correlation-based - no shared control correction needed)
SIdata <- data.frame(
  Toxicant = c("Hg", "PCBs"),
  Mean = c(0.653, 0.438),
  CI_Lower = c(0.625, 0.345),
  CI_Upper = c(0.682, 0.531),
  n = c(17, 6)
)

# Muscle Inclusive (correlation-based - no shared control correction needed)
MuIdata <- data.frame(
  Toxicant = c("Hg", "PCBs"),
  Mean = c(0.913, 1.379),
  CI_Lower = c(0.892, 1.233),
  CI_Upper = c(0.934, 1.524),
  n = c(25, 6)
)

# Hatching Inclusive (Borenstein-corrected)
HIdata <- data.frame(
  Toxicant = c("Hg", "PCBs"),
  Mean = c(-1.0844, -0.7353),
  CI_Lower = c(-1.6427, -1.2926),
  CI_Upper = c(-0.5262, -0.1781),
  n = c(8, 5)
)

# Mortality Inclusive (Borenstein-corrected)
MIdata <- data.frame(
  Toxicant = c("Hg", "PCBs"),
  Mean = c(1.6711, 0.6901),
  CI_Lower = c(1.1306, 0.2816),
  CI_Upper = c(2.2116, 1.0985),
  n = c(12, 12)
)

# Growth Inclusive (Borenstein-corrected)
GIdata <- data.frame(
  Toxicant = c("Hg", "PCBs"),
  Mean = c(-0.9209, 0.0654),
  CI_Lower = c(-1.3381, -0.4464),
  CI_Upper = c(-0.5037, 0.5772),
  n = c(10, 5)
)

# -------- CONSERVATIVE APPROACH --------

# Size Conservative
SCdata <- data.frame(
  Toxicant = c("Hg", "PCBs"),
  Mean = c(0.638, 0.432),
  CI_Lower = c(0.565, 0.283),
  CI_Upper = c(0.711, 0.580),
  n = c(7, 4)
)

# Muscle Conservative
MuCdata <- data.frame(
  Toxicant = c("Hg", "PCBs"),
  Mean = c(0.958, 1.042),
  CI_Lower = c(0.926, 0.841),
  CI_Upper = c(0.989, 1.244),
  n = c(15, 5)
)

# Hatching Conservative
HCdata <- data.frame(
  Toxicant = c("Hg", "PCBs"),
  Mean = c(-0.655, -0.614),
  CI_Lower = c(-0.912, -0.999),
  CI_Upper = c(-0.398, -0.229),
  n = c(8, 5)
)

# Mortality Conservative
MCdata <- data.frame(
  Toxicant = c("Hg", "PCBs"),
  Mean = c(0.995, 0.475),
  CI_Lower = c(0.809, 0.324),
  CI_Upper = c(1.181, 0.626),
  n = c(12, 12)
)

# Growth Conservative
GCdata <- data.frame(
  Toxicant = c("Hg", "PCBs"),
  Mean = c(-0.653, 0.488),
  CI_Lower = c(-0.836, 0.116),
  CI_Upper = c(-0.470, 0.861),
  n = c(10, 5)
)

#===================== Colorblind-Friendly Palette ===============================
# Wong palette - distinguishable for most forms of colorblindness
hg_color <- "#E69F00"    # Orange
pcb_color <- "#0072B2"   # Blue

#============================= Publication Theme ===================================

pub_theme <- theme_classic(base_size = 11, base_family = "serif") +
  theme(
    # Axis formatting
    axis.line = element_line(linewidth = 0.6, color = "black"),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    
    # Title formatting
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    
    # Legend - will be added separately
    legend.position = "none",
    
    # Clean background
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA),
    
    # Margins
    plot.margin = margin(10, 10, 10, 10)
  )

#============================= Plot Function ===================================

make_plot <- function(data, title, show_ylab = FALSE, y_limits = c(-2, 2.25), ylab_text = "Effect Size") {
  
  # Create labels with n
  data$Label <- paste0(data$Toxicant, "\n(n=", data$n, ")")
  
  p <- ggplot(data, aes(x = Label, y = Mean)) +
    
    # Null effect line - prominent
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.7, color = "gray40") +
    
    # Error bars - thicker for visibility
    geom_errorbar(
      aes(ymin = CI_Lower, ymax = CI_Upper, color = Toxicant),
      width = 0.15, 
      linewidth = 1.0
    ) +
    
    # Point estimates - diamond shape, larger
    geom_point(
      aes(fill = Toxicant), 
      shape = 23,  # Diamond
      size = 4,
      color = "black",  # Black outline
      stroke = 0.8
    ) +
    
    # Colorblind-friendly palette
    scale_color_manual(values = c("Hg" = hg_color, "PCBs" = pcb_color)) +
    scale_fill_manual(values = c("Hg" = hg_color, "PCBs" = pcb_color)) +
    
    # Axis settings
    scale_y_continuous(
      limits = y_limits,
      breaks = seq(floor(y_limits[1]), ceiling(y_limits[2]), by = 0.5),
      expand = c(0.02, 0)
    ) +
    
    # Labels
    xlab("") +
    ylab(if(show_ylab) ylab_text else "") +
    ggtitle(title) +
    
    # Apply theme
    pub_theme
  
  return(p)
}

#============================= Create All Plots ===================================

# Larval Success Plots (Hedges' d)
HI_plot <- make_plot(HIdata, "Hatching Success", show_ylab = TRUE, ylab_text = "Effect Size (Hedges' d)")
HC_plot <- make_plot(HCdata, "Hatching Success", show_ylab = FALSE)

MI_plot <- make_plot(MIdata, "Mortality", show_ylab = TRUE, ylab_text = "Effect Size (Hedges' d)")
MC_plot <- make_plot(MCdata, "Mortality", show_ylab = FALSE)

GI_plot <- make_plot(GIdata, "Growth", show_ylab = TRUE, ylab_text = "Effect Size (Hedges' d)")
GC_plot <- make_plot(GCdata, "Growth", show_ylab = FALSE)

# Maternal Characteristics Plots (Fisher's Z)
SI_plot <- make_plot(SIdata, "Size", show_ylab = TRUE, ylab_text = "Effect Size (Fisher's Z)")
SC_plot <- make_plot(SCdata, "Size", show_ylab = FALSE)

MuI_plot <- make_plot(MuIdata, "Muscle", show_ylab = TRUE, ylab_text = "Effect Size (Fisher's Z)")
MuC_plot <- make_plot(MuCdata, "Muscle", show_ylab = FALSE)

#=======Titles==========
Inc_title <- ggdraw() + 
  draw_label("Inclusive Approach", fontface = 'bold', size = 14, x = 0.5, hjust = 0.5) +
  theme(plot.background = element_rect(fill = "white", color = NA))

Con_title <- ggdraw() + 
  draw_label("Conservative Approach", fontface = 'bold', size = 14, x = 0.5, hjust = 0.5) +
  theme(plot.background = element_rect(fill = "white", color = NA))

#==================== Combined Plots ====================

# Larval Success (Hatching, Mortality, Growth)
larval_plot <- cowplot::plot_grid(
  Inc_title, Con_title,
  HI_plot, HC_plot,
  MI_plot, MC_plot,
  GI_plot, GC_plot,
  ncol = 2, nrow = 4, 
  rel_heights = c(0.08, 1, 1, 1),
  labels = c("", "", "A", "B", "C", "D", "E", "F"),
  label_size = 14,
  label_fontface = "bold",
  label_fontfamily = "serif"
)

# Maternal Characteristics (Size, Muscle)
maternal_plot <- cowplot::plot_grid(
  Inc_title, Con_title,
  SI_plot, SC_plot,
  MuI_plot, MuC_plot,
  ncol = 2, nrow = 3,
  rel_heights = c(0.1, 1, 1),
  labels = c("", "", "A", "B", "C", "D"),
  label_size = 14,
  label_fontface = "bold",
  label_fontfamily = "serif"
)

# Display plots
print(larval_plot)
print(maternal_plot)

# Save plots - multiple formats for journal submission
# PNG for review
ggsave("figures/larval_success_plot.png", larval_plot, width = 8, height = 10, dpi = 900, bg = "white")
ggsave("figures/maternal_characteristics_plot.png", maternal_plot, width = 8, height = 7, dpi = 900, bg = "white")

# PDF for publication (vector, infinitely scalable)
ggsave("figures/larval_success_plot.pdf", larval_plot, width = 8, height = 10, bg = "white")
ggsave("figures/maternal_characteristics_plot.pdf", maternal_plot, width = 8, height = 7, bg = "white")

# TIFF for some journals (high-res raster)
ggsave("figures/larval_success_plot.tiff", larval_plot, width = 8, height = 10, dpi = 600, compression = "lzw", bg = "white")
ggsave("figures/maternal_characteristics_plot.tiff", maternal_plot, width = 8, height = 7, dpi = 600, compression = "lzw", bg = "white")