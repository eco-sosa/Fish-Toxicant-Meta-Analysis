# Load required libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Set working directory to the repo root before sourcing this script.
# In RStudio: Session > Set Working Directory > To Project Directory.
# Otherwise run: setwd("path/to/Hg_PCBs_Meta_Repo")

# Function to create forest plot for larval success endpoints using Borenstein-corrected data
create_larval_success_plot_borenstein <- function(file_path, xlim = c(-20, 20)) {

  # Read the Borenstein-corrected data (CSV format)
  df <- read.csv(file_path, check.names = FALSE)

  # Filter for larval success endpoints
  df_larval <- df[df$Relationship %in% c("Hatching", "Growth", "Mortality"), ]

  # Calculate 95% confidence intervals
  df_larval$CI_Lower <- df_larval$`Effect size` - 1.96 * df_larval$`Standard Error`
  df_larval$CI_Upper <- df_larval$`Effect size` + 1.96 * df_larval$`Standard Error`

  # Identify and exclude extreme outliers (outside xlim range)
  excluded <- df_larval[df_larval$`Effect size` < xlim[1] | df_larval$`Effect size` > xlim[2], ]
  if (nrow(excluded) > 0) {
    cat("\nStudies excluded from plot (effect size outside", xlim[1], "to", xlim[2], "):\n")
    for (i in 1:nrow(excluded)) {
      cat("  -", excluded$Study[i], "(", excluded$Relationship[i], "):",
          "d =", round(excluded$`Effect size`[i], 1), "\n")
    }
  }

  # Filter to only include studies within range
  df_larval <- df_larval[df_larval$`Effect size` >= xlim[1] & df_larval$`Effect size` <= xlim[2], ]

  # Create study labels
  df_larval$Study_Clean <- df_larval$Study
  df_larval$Study_Clean <- trimws(df_larval$Study_Clean)
  df_larval$Final_Label <- df_larval$Study_Clean

  # Calculate weights (inverse of variance)
  df_larval$Weight <- ifelse(df_larval$Variance > 0, 1 / df_larval$Variance, 1)

  # Sort by relationship type, then by effect size within each type
  df_larval$Relationship <- factor(df_larval$Relationship, levels = c("Mortality", "Growth", "Hatching"))
  df_larval <- df_larval[order(df_larval$Relationship, df_larval$`Effect size`), ]

  # Create row numbers with spacing for subgroup headers
  df_larval$Row <- NA
  header_positions <- c()
  header_labels <- c()
  current_row <- 0

  for(rel in c("Mortality", "Growth", "Hatching")) {
    rel_data <- df_larval[df_larval$Relationship == rel, ]
    if(nrow(rel_data) > 0) {
      # Add spacing and record header position
      current_row <- current_row + 1.5
      header_positions <- c(header_positions, current_row)
      header_labels <- c(header_labels, paste0(rel, " (", nrow(rel_data), ")"))

      # Add study rows for this relationship
      rel_indices <- which(df_larval$Relationship == rel)
      study_rows <- current_row + (1:length(rel_indices))
      df_larval$Row[rel_indices] <- study_rows
      current_row <- max(study_rows) + 0.5
    }
  }


  # Create the plot
  p1 <- ggplot(df_larval, aes(x = `Effect size`, y = Row)) +

    # Add confidence intervals with color mapping (clip to xlim)
    geom_segment(aes(x = pmax(CI_Lower, xlim[1]), xend = pmin(CI_Upper, xlim[2]),
                     y = Row, yend = Row, color = `Toxi Group`),
                 linewidth = 1.2, alpha = 0.6) +

    # Add point estimates (colored by toxicant)
    geom_point(aes(size = Weight, color = `Toxi Group`),
               alpha = 0.7) +

    # Add vertical line at null effect (0)
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.6) +

    # Add relationship headers on the left side
    annotate("text",
             x = xlim[1] - 3,
             y = header_positions,
             label = header_labels,
             hjust = 0, vjust = 0.5,
             fontface = "bold", size = 3, color = "black") +

    # Color scheme - gold for Hg, navy for PCBs with larger legend points
    scale_color_manual(values = c("Hg" = "gold3", "PCBs" = "navy"),
                       name = "",
                       guide = guide_legend(override.aes = list(size = 5))) +

    # Size legend - remove weight legend, increase toxicant legend size
    scale_size_continuous(range = c(2, 8), guide = "none") +

    # Y-axis labels (only for study rows, not headers)
    scale_y_continuous(breaks = df_larval$Row,
                       labels = df_larval$Final_Label,
                       limits = c(0.5, max(df_larval$Row) + 0.5)) +

    # Set x-axis limits
    coord_cartesian(clip = "off", xlim = xlim) +
    scale_x_continuous(expand = expansion(mult = c(0.25, 0.05))) +

    # Labels and theme
    labs(x = "Effect Size (Hedges d)",
         y = "",
         title = "Larval Success Endpoints (Borenstein-Corrected)") +

    theme_minimal() +
    theme(axis.text.y = element_text(size = 7),
          axis.text.x = element_text(size = 10),
          plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(size = 9, hjust = 0.5),
          legend.position = "right",
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(color = "grey80"))

  return(p1)
}

# Function to create forest plot for larval success endpoints (original - kept for reference)
create_larval_success_plot <- function(file_path) {

  # Read the data
  df <- read_excel(file_path)
  
  # Filter for larval success endpoints
  df_larval <- df[df$Relationship %in% c("Hatching", "Growth", "Mortality"), ]
  
  # Create study labels - use original study names as they are properly formatted
  df_larval$Study_Clean <- df_larval$Study
  df_larval$Study_Clean <- trimws(df_larval$Study_Clean)
  
  # Check for duplicate studies within each relationship type
  df_larval$Final_Label <- df_larval$Study_Clean
  
  for(rel in unique(df_larval$Relationship)) {
    rel_subset <- df_larval[df_larval$Relationship == rel, ]
    study_counts_in_rel <- table(rel_subset$Study_Clean)
    duplicated_studies <- names(study_counts_in_rel)[study_counts_in_rel > 1]
    
    # Add site info only for studies that appear multiple times within the same relationship
    mask <- df_larval$Relationship == rel & 
      df_larval$Study_Clean %in% duplicated_studies & 
      !is.na(df_larval$Site) & 
      df_larval$Site != ""
    
    df_larval$Final_Label[mask] <- paste0(df_larval$Study_Clean[mask], " (", df_larval$Site[mask], ")")
  }
  
  # Calculate 95% confidence intervals
  df_larval$CI_Lower <- df_larval$`Effect size` - 1.96 * df_larval$`Standard Error`
  df_larval$CI_Upper <- df_larval$`Effect size` + 1.96 * df_larval$`Standard Error`
  
  # Calculate weights (inverse of variance)
  df_larval$Weight <- ifelse(df_larval$Variance > 0, 1 / df_larval$Variance, 1)
  
  # Sort by relationship type, then by effect size within each type
  df_larval$Relationship <- factor(df_larval$Relationship, levels = c("Mortality", "Growth", "Hatching"))
  df_larval <- df_larval[order(df_larval$Relationship, df_larval$`Effect size`), ]
  
  # Create row numbers with spacing for subgroup headers
  df_larval$Row <- NA
  header_positions <- c()
  header_labels <- c()
  current_row <- 0
  
  for(rel in c("Mortality", "Growth", "Hatching")) {
    rel_data <- df_larval[df_larval$Relationship == rel, ]
    if(nrow(rel_data) > 0) {
      # Add spacing and record header position
      current_row <- current_row + 1.5
      header_positions <- c(header_positions, current_row)
      header_labels <- c(header_labels, paste0(rel, " (", nrow(rel_data), ")"))
      
      # Add study rows for this relationship
      rel_indices <- which(df_larval$Relationship == rel)
      study_rows <- current_row + (1:length(rel_indices))
      df_larval$Row[rel_indices] <- study_rows
      current_row <- max(study_rows) + 0.5
    }
  }
  
  # Create the plot
  p1 <- ggplot(df_larval, aes(x = `Effect size`, y = Row)) +
    
    # Add confidence intervals with color mapping
    geom_segment(aes(x = CI_Lower, xend = CI_Upper, y = Row, yend = Row, color = `Toxi Group`),
                 size = 1.2, alpha = 0.6) +
    
    # Add point estimates (colored by toxicant)
    geom_point(aes(size = Weight, color = `Toxi Group`),
               alpha = 0.7) +
    
    # Add vertical line at null effect (0)
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.6) +
    
    # Add relationship headers on the left side
    annotate("text", 
             x = min(df_larval$`Effect size`, na.rm = TRUE) - 5, 
             y = header_positions, 
             label = header_labels,
             hjust = 0, vjust = 0.5, 
             fontface = "bold", size = 3, color = "black") +
    
    # Color scheme - gold for Hg, navy for PCBs with larger legend points
    scale_color_manual(values = c("Hg" = "gold3", "PCBs" = "navy"),
                       name = "",
                       guide = guide_legend(override.aes = list(size = 5))) +
    
    # Size legend - remove weight legend, increase toxicant legend size
    scale_size_continuous(range = c(2, 8), guide = "none") +
    
    # Y-axis labels (only for study rows, not headers)
    scale_y_continuous(breaks = df_larval$Row,
                       labels = df_larval$Final_Label,
                       limits = c(0.5, max(df_larval$Row) + 0.5)) +
    
    # Extend x-axis to make room for headers  
    coord_cartesian(clip = "off") +
    scale_x_continuous(expand = expansion(mult = c(0.4, 0.1))) +
    
    # Labels and theme
    labs(x = "Effect Size (Hedges d)",
         y = "",
         title = "Larval Success Endpoints",
         subtitle = "Point size proportional to study weight; horizontal lines show 95% CI") +
    
    theme_minimal() +
    theme(axis.text.y = element_text(size = 7),
          axis.text.x = element_text(size = 10),
          plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(size = 10, hjust = 0.5),
          legend.position = "right",
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(color = "grey80"))
  
  return(p1)
}

# Function to create forest plot for size and muscle endpoints
create_size_muscle_plot <- function(file_path) {
  
  # Read the data
  df <- read_excel(file_path)
  
  # Filter for size and muscle endpoints
  df_size_muscle <- df[df$Relationship %in% c("Size", "Muscle"), ]
  
  # Create study labels - use original study names as they are properly formatted
  df_size_muscle$Study_Clean <- df_size_muscle$Study
  df_size_muscle$Study_Clean <- trimws(df_size_muscle$Study_Clean)
  
  # Check for duplicate studies within each relationship type
  df_size_muscle$Final_Label <- df_size_muscle$Study_Clean
  
  for(rel in unique(df_size_muscle$Relationship)) {
    rel_subset <- df_size_muscle[df_size_muscle$Relationship == rel, ]
    study_counts_in_rel <- table(rel_subset$Study_Clean)
    duplicated_studies <- names(study_counts_in_rel)[study_counts_in_rel > 1]
    
    # Add site info only for studies that appear multiple times within the same relationship
    mask <- df_size_muscle$Relationship == rel & 
      df_size_muscle$Study_Clean %in% duplicated_studies & 
      !is.na(df_size_muscle$Site) & 
      df_size_muscle$Site != ""
    
    df_size_muscle$Final_Label[mask] <- paste0(df_size_muscle$Study_Clean[mask], " (", df_size_muscle$Site[mask], ")")
  }
  
  # Calculate 95% confidence intervals
  df_size_muscle$CI_Lower <- df_size_muscle$`Effect size` - 1.96 * df_size_muscle$`Standard Error`
  df_size_muscle$CI_Upper <- df_size_muscle$`Effect size` + 1.96 * df_size_muscle$`Standard Error`
  
  # Calculate weights (inverse of variance)
  df_size_muscle$Weight <- ifelse(df_size_muscle$Variance > 0, 1 / df_size_muscle$Variance, 1)
  
  # Sort by relationship type, then by effect size within each type
  df_size_muscle$Relationship <- factor(df_size_muscle$Relationship, levels = c("Size", "Muscle"))
  df_size_muscle <- df_size_muscle[order(df_size_muscle$Relationship, df_size_muscle$`Effect size`), ]
  
  # Create row numbers with spacing for subgroup headers
  df_size_muscle$Row <- NA
  header_positions <- c()
  header_labels <- c()
  current_row <- 0
  
  for(rel in c("Size", "Muscle")) {
    rel_data <- df_size_muscle[df_size_muscle$Relationship == rel, ]
    if(nrow(rel_data) > 0) {
      # Add spacing and record header position
      current_row <- current_row + 1.5
      header_positions <- c(header_positions, current_row)
      header_labels <- c(header_labels, paste0(rel, " (", nrow(rel_data), ")"))
      
      # Add study rows for this relationship
      rel_indices <- which(df_size_muscle$Relationship == rel)
      study_rows <- current_row + (1:length(rel_indices))
      df_size_muscle$Row[rel_indices] <- study_rows
      current_row <- max(study_rows) + 0.5
    }
  }
  
  # Create the plot
  p2 <- ggplot(df_size_muscle, aes(x = `Effect size`, y = Row)) +
    
    # Add confidence intervals with color mapping
    geom_segment(aes(x = CI_Lower, xend = CI_Upper, y = Row, yend = Row, color = `Toxi Group`),
                 size = 1.2, alpha = 0.6) +
    
    # Add point estimates (colored by toxicant)
    geom_point(aes(size = Weight, color = `Toxi Group`),
               alpha = 0.7) +
    
    # Add vertical line at null effect (0)
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.6) +
    
    # Add relationship headers on the left side
    annotate("text", 
             x = min(df_size_muscle$`Effect size`, na.rm = TRUE) - 0.3, 
             y = header_positions, 
             label = header_labels,
             hjust = 0, vjust = 0.5, 
             fontface = "bold", size = 3, color = "black") +
    
    # Color scheme - gold for Hg, navy for PCBs with larger legend points
    scale_color_manual(values = c("Hg" = "gold3", "PCBs" = "navy"),
                       name = "",
                       guide = guide_legend(override.aes = list(size = 5))) +
    
    # Size legend - remove weight legend, increase toxicant legend size
    scale_size_continuous(range = c(2, 8), guide = "none") +
    
    # Y-axis labels (only for study rows, not headers)
    scale_y_continuous(breaks = df_size_muscle$Row,
                       labels = df_size_muscle$Final_Label,
                       limits = c(0.5, max(df_size_muscle$Row) + 0.5)) +
    
    # Extend x-axis to make room for headers
    coord_cartesian(clip = "off") +
    scale_x_continuous(expand = expansion(mult = c(0.4, 0.1))) +
    
    # Labels and theme
    labs(x = "Effect Size (Fisher Z)",
         y = "",
         title = "Size and Muscle Relationships",
         subtitle = "Point size proportional to study weight; horizontal lines show 95% CI") +
    
    theme_minimal() +
    theme(axis.text.y = element_text(size = 7),
          axis.text.x = element_text(size = 10),
          plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(size = 10, hjust = 0.5),
          legend.position = "right",
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(color = "grey80"))
  
  return(p2)
}

# Usage
# Use Borenstein-corrected data for larval success endpoints
file_path_borenstein <- "data/forest_plot_borenstein.csv"
file_path_original <- "data/forest plot.xlsx"

# Create both plots - larval uses Borenstein-corrected, size/muscle uses original
larval_plot <- create_larval_success_plot_borenstein(file_path_borenstein)
size_muscle_plot <- create_size_muscle_plot(file_path_original)

# Create combined stacked plot
combined_plot <- grid.arrange(
  larval_plot, 
  size_muscle_plot, 
  nrow = 2,
  heights = c(1.2, 1)  # Give more space to larval plot since it has more studies
)

# Display plots
print(larval_plot)
print(size_muscle_plot)
print(combined_plot)

# Save plots with highest resolution
ggsave("figures/larval_success_forest_plot.png", larval_plot, 
       width = 12, height = 10, dpi = 1200)

ggsave("figures/size_muscle_forest_plot.png", size_muscle_plot, 
       width = 12, height = 8, dpi = 1200)

# Save combined plot with highest resolution
ggsave("figures/combined_forest_plot.png", combined_plot, 
       width = 12, height = 16, dpi = 1200)

# Print summary statistics
df_larval <- read.csv(file_path_borenstein, check.names = FALSE)
df_larval <- df_larval[df_larval$Relationship %in% c("Hatching", "Growth", "Mortality"), ]

cat("Larval Success Forest Plot Summary (Borenstein-Corrected):\n")
cat("Total studies:", nrow(df_larval), "\n")
cat("Mercury studies:", sum(df_larval$`Toxi Group` == "Hg"), "\n")
cat("PCB studies:", sum(df_larval$`Toxi Group` == "PCBs"), "\n")
cat("Effect size range:", round(min(df_larval$`Effect size`), 3), "to", round(max(df_larval$`Effect size`), 3), "\n\n")

# Size/muscle summary
df_original <- read_excel(file_path_original)
df_size_muscle <- df_original[df_original$Relationship %in% c("Size", "Muscle"), ]
cat("Size/Muscle Forest Plot Summary:\n")
cat("Total studies:", nrow(df_size_muscle), "\n")
cat("Mercury studies:", sum(df_size_muscle$`Toxi Group` == "Hg"), "\n")
cat("PCB studies:", sum(df_size_muscle$`Toxi Group` == "PCBs"), "\n")
cat("Effect size range:", round(min(df_size_muscle$`Effect size`), 3), "to", round(max(df_size_muscle$`Effect size`), 3), "\n")