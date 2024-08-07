library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

######### load data ##########

AD427_metadata = read.csv("AD427_metadata.csv")
ADMR_metadata = read.csv("ADMR_metadata.csv")
metadata = read.csv("metadata.csv")


######### helper functions ##########

get_significance_code <- function(p_value) {
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else if (p_value < 0.1) {
    return(".")
  } else {
    return("")
  }
}

get_p_value <- function(x, y, alternative = 'two.sided') {
  num_categories <- length(unique(x))
  if (num_categories > 2) {
    # ANOVA test if there are more than 2 categories
    anova_model <- aov(y ~ x, data=metadata)

    # summarize ANOVA results
    anova_summary <- summary(anova_model)
    print(anova_summary)
    p_value_info <- anova_summary[[1]]["Pr(>F)"][1]
    p_value <- as.numeric(p_value_info[[1]][1])
  }
  else {
    # t-test if there are 2 categories
    t_test_result <- t.test(y ~ x, alternative = alternative)
    p_value <- t_test_result$p.value
    print(p_value)
  }
  return (p_value)
}

is_categorical <- function(data, col, num_categories) {
  unique_values = unique(data[[col]])
  return(length(unique_values) <= num_categories)
}

key_func <- function(col) is_categorical(metadata, col, num_categories = 10)

is_significant <- function (x, y, alterative = 'two.sided'){
  return (get_p_value(x, y, alterative) <= 0.05)
}


######### variable definitions ##########

pathologies = c('apoe_genotype', 'cogdx', 'age_bl', 'age_death', 'msex', 'race', 'braaksc', 'gpath', 'pmi', 'amyloid', 'plaq_d', 'plaq_n', 'nft', 'tangles', 'arteriol_scler', 'AD_states', 'hypertension_bl', 'chf_bl', 'heart_bl', 'stroke_bl', 'cogdx_stroke', 'diabetes_sr_rx_bl', 'headinjrloc_bl', 'thyroid_bl', 'ceradsc', 'cvda_4gp2')
categorical_vars = names(metadata)[sapply(names(metadata), key_func)]
categorical_pathologies <- pathologies[sapply(pathologies, key_func)]
continuous_pathologies <- pathologies[!sapply(pathologies, key_func)]

colors = brewer.pal(8, "Pastel2")

######### box plots ##########
par(mfrow=c(5, 4), mar=c(5.5, 5, 1, 1))
for (col in categorical_pathologies) {
    print(paste("Plotting", col))
    x = as.factor(metadata[[col]])
    y = metadata$bmi_lv

    # do a significance test
    p_value <- get_p_value(x, y, alternative = 'two.sided')
    sig_code = get_significance_code(p_value) # asterisks for display

    # create boxplot
    xlabel = paste(col, '\n', sig_code, 'p =', signif(p_value, 3))
    boxplot(y ~ x, xlab = xlabel, ylab = "bmi_lv", na.action = na.omit)
  }


######### scatter plots ##########
par(mfrow=c(4, 4), mar=c(5.5, 5, 1, 1))
for (col in pathologies) {
  if (!is_categorical(metadata, col, 20)) {
    plot(x=metadata$bmi_lv, y=metadata[[col]], xlab="bmi_lv", ylab=col)
  }
}

######### pie charts ##########
plot_pie_charts <- function(data, pathologies, title=NA) {
  # set margins
  if (!is.na(title)) {
    par(mfrow=c(5, 4), mar=c(1, 1, 4, 2), oma=c(0, 0, 2, 0))  # Add space for the main title
  } else {
    par(mfrow=c(5, 4), mar=c(1, 1, 4, 2))
  }

  for (col in pathologies) {
    if (is_categorical(data, col, 20)) {
      counts_table <- table(data[[col]])
      pie(counts_table, main = col, col = colors, cex = 0.8)
      legend("topright", legend = names(counts_table), fill = colors, cex = 0.8, title = "Categories", inset=c(-0.1, 0), xpd=TRUE)
    }
  }
  if (!is.na(title)) {
    mtext(title, side = 3, line = 0, outer = TRUE, cex = 1.5)
  }
}

png("pie_charts.png", width=2800, height=2400, res=300)
plot_pie_charts(metadata, pathologies)
dev.off()


######### pie charts by BMI groups ##########

bmi_groups <- c("bmi_<20", "bmi_20-25", "bmi_25-30", "bmi_30+")
# create a vector of pathology variables that correlate significantly with bmi (p < 0.05)
bmi_varied_vars <- categorical_pathologies[sapply(categorical_pathologies, function (var) is_significant(metadata[[var]], metadata$bmi_lv))]
filtered_bmi_data <- metadata[metadata$bmi_groups %in% bmi_groups, ] # filter data containing non-empty bmi values

# make a figure for each bmi group
for (group in bmi_groups) {
  filename <- sprintf("pie_charts_%s.png", group)

  png(filename = filename, width=2800, height=2400, res=300)
  plot_pie_charts(filtered_bmi_data[filtered_bmi_data$bmi_groups == group, ], pathologies, title = group)

  print(paste("Saving figure to:", filename))
  dev.off()
}

# make a figure comparing distributions across all bmi groups
png(filename = "pie_charts_by_bmi_groups.png", width=2800, height=2400, res=300)

par(mfrow=c(length(bmi_varied_vars), length(bmi_groups)), mar=c(1, 1, 4, 2))
for (col in bmi_varied_vars) {
  for (group in bmi_groups) {
    filtered_data <- metadata[metadata$bmi_groups == group, ]
    counts_table <- table(filtered_data[[col]])

    pie(counts_table, main = sprintf("%s, %s", col, group), col = colors)
    legend("topright", legend = names(counts_table), fill = colors, cex = 0.8, inset=c(-0.1, 0), xpd=TRUE)
  }
}

dev.off()


######### pie charts by AD groups ##########

AD_groups <- c('earlyAD', 'lateAD', 'nonAD')
AD_varied_vars <- c()

# significance test for two discrete variables done in discrete_significance_test.R

AD_varied_vars <- c("apoe_genotype", "cogdx", "msex", "braaksc", "diabetes_sr_rx_bl", "ceradsc")
