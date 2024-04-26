# GGE Biplot
# modified from '07_GGE_analysis.R' from Falcon et al.(2019)
# phenotypes & environmental covariates from Lopez-Cruz et al. (2023)
# Authors: Elliot Braun, Eleanor Carr
# Date: 26 April 2024

# packages and functions
library(ggplot2) # for making plots (scatterplots, barplots, etc)
library(tidyverse) # for editing & manipulating datasets
library(GGEBiplots) # for making the biplots

# set working directory
# directory should include:
  # phenotype file 
  # set of environmental covariates
  # list of desired genotypes to filter by
getwd()
setwd("/Users/eleanor/Documents/Spring2024/Frontiers/Module_2_G2F")

# read in files
phenos <- read.csv("PHENO_original.csv") # phenotype file
e_cov <- read.csv("ECOV.csv") # environmental covariates file
checks <- read.csv("checks.csv") # list of the checks

# filter and reformat the phenotype data
phenos <- phenos[phenos$genotype %in% checks$YS, ] # filter to only the check genotypes
# filter to only a particular year range
# phenos <- phenos %>% 
#   filter(year==2020 | year==2021) 

# save some variables as factors (optional)
# phenos$year_loc <- factor(phenos$year_loc)
# phenos$irrigated <- factor(phenos$irrigated)
# phenos$region <- factor(phenos$region)

# calculate the average yield for a genotype x year-location combination
average_yield_df <- phenos %>%
  group_by(year_loc, genotype) %>%
  summarise(average_yield = mean(yield, na.rm = TRUE), .groups = 'drop')

# use pivot_wider() to make the GGEBiplots input file
gge_input <- average_yield_df %>%
pivot_wider(names_from = year_loc, values_from = average_yield)
genotype <- gge_input$genotype
gge_input <- select(gge_input,-c(genotype)) # remove the genotype column
row.names(gge_input) <- genotype # set the row names as the genotypes

# a note: later, when you want to calculate the row and column averages for imputation,
  # all values must be numeric, so you have to deal with the genotype names somehow.
  # How you do that is up to you, we have chosen to get rid of them as a particular 
  # genotype performance was not relevant to our analysis.


# count the number of cells in the input file that are empty/missing
na_count <- 0
for (r in 1:nrow(gge_input)) {
  for (c in 1:ncol(gge_input)) {
    if (is.na(gge_input[r,c])) {
      na_count <- na_count + 1
    }
  }
}

# Remove rows with >90% missing data
# gge_v2 <- gge_input[,colMeans(is.na(gge_v1)) < 0.9] # remove environments
# gge_v3 <- gge_v2[rowMeans(is.na(gge_v2)) < 0.9,] # remove genotypes

# we can extrapolate provided code to set our own threshold, in this case >65%
# we also want to know how much imputation is happening, so we'll create some intermediate variables to calculate this

# a note from the editor: why are we using only the second column onward? Is this a remnant from when there was a genotype column?
col_na_proportions <- colMeans(is.na(gge_input[, 2:ncol(gge_input)])) # average number of NAs per column
env_pass_thres <- gge_input[, c(TRUE, col_na_proportions < 0.65)] # remove environments
row_na_proportions <- rowMeans(is.na(env_pass_thres[, -1])) # average number of NAs per remaining row
both_pass_thresh <- env_pass_thres[row_na_proportions < 0.65, ] # remove genotypes

# estimate missing data as row/column mean
row_means <- rowMeans(both_pass_thresh, na.rm = TRUE)
col_means <- colMeans(both_pass_thresh, na.rm = TRUE)
impute_count = 0 # keep track of how many times a value is imputed
for (r in 1:nrow(both_pass_thresh)) {
  for (c in 1:ncol(both_pass_thresh)) {
    if (is.na(both_pass_thresh[r,c])) {
      impute_count <- impute_count + 1
      both_pass_thresh[r,c] <- mean(row_means[r], col_means[c])
    }
  }
}

# now that we have the data prepped, we can input it to GGEBiplots
data <- both_pass_thresh

# GGEPlot() requires a GGEModel as an input
GGE1 <- GGEModel(data, centering = "tester", scaling = 'sd', SVP = "column") # Model for discriminability vs representativeness plot

# generate the GGE biplot using GGEPlot() 
discrimVsRepresent <- GGEPlot(GGE1, type = 7)

# if you want to see the physical plot, uncomment and run the below code
# GGEPlot(GGE1, type = 7, colSegment = "grey91", textGen = element_text(family = "", face = 1, color = "grey91", size = 1.5, hjust = 0, vjust = 0, angle = 0), textEnv = element_text(family = "", face = 1, color = "blue", size = 1, hjust = 0, vjust = 0, angle = 0))

# calculate the environment vector lengths
nenv <- ncol(data)
ngen <- nrow(data)
vectorLengths <- matrix(nrow = nenv, ncol = 1, dimnames = list(colnames(data), "vectorLength"))
for (v in (ngen + 1):(ngen + nenv)) {
  vectorLengths[(v - ngen),1] <- sqrt(discrimVsRepresent$data[v,1]^2 + (discrimVsRepresent$data[v,2])^2)
}

# if desired, write the vector lengths out to a file
write.csv(vectorLengths, "discriminability_ranks.csv")

# read in discriminability ranks
discrims <- read.csv("discriminability_ranks.csv")

# give the 'X' column a proper name
colnames(discrims)[1] <- "year_loc" 

# depending on how the file was read in, you may need to reformat the year-location names using the code below:
# discrim_names <- discrims$year_loc
# new_discrim_names <- gsub("^X(\\d{4})\\.(.*)$", "\\1-\\2", discrim_names)
# even_newer_discrim_names <- new_discrim_names <- gsub("\\.(?=Dry|Early|Late)", "-", new_discrim_names, perl = TRUE)
# even_newer_discrim_names <- even_newer_discrim_names
# discrims$year_loc <- even_newer_discrim_names
# discrims <- discrims %>% 
  # arrange(vectorLength) %>%
  # filter(year_loc != 'X')

# we chose to separate 'high' from 'low' based on the median discriminability value
hist(discrims$vectorLength)
med <- median(discrims$vectorLength)
high_discrims <- discrims %>% 
  filter(vectorLength > med)
low_discrims <- discrims %>% 
  filter(vectorLength < med)

colnames(e_cov)[1] <- "year_loc"

## perform t-test on year_loc groups

# set up high and low e_cov
# vector that is given by subsetting e_cov into covariate X and corresponding groups Y
group_covs1 <- e_cov %>% 
  filter(year_loc %in% high_discrims$year_loc) 

# vector that is given by subsetting e_cov into covariate X and corresponding groups Z
group_covs2 <- e_cov %>% 
  filter(year_loc %in% low_discrims$year_loc) 

# an empty df for placing the e_cov in
cov_df <- data.frame()

# for each of the environmental covariates:
for (i in 1:length(colnames(group_covs1 %>% select(-year_loc)))) {
  col_name <- colnames(group_covs1 %>% select(-year_loc))[i] # get the e_cov name
  covector1 <- group_covs1 %>% # filter the high e_covs for the current e_cov name
    select(-year_loc) %>% 
    select(all_of(col_name))
  covector2 <- group_covs2 %>% # filter the low e_covs for the current e_cov name
    select(-year_loc) %>% 
    select(all_of(col_name))
  t_test <- t.test(covector1, covector2) # perform a t-test on these two groups
  p_value <- t_test$p.value # extract the p-value
  new_row <- c(col_name, p_value) # create a df row with the name & p-value
  cov_df <- rbind(cov_df, new_row) # append the new row to the output df
}

# rename the output columns
colnames(cov_df) <- c("covariate", "p_value")

# we need to do some corrections on these p-values to correct for false positives
# see for more info: https://www.publichealth.columbia.edu/research/population-health-methods/false-discovery-rate 

# write a function to group the e_cov
# input: a list of the e_cov that belong to a group
# output: a dataframe containing only the columns containing the desired covariates
grouper <- function(string_list) {
  list_df <- list()
  
  for (starter in string_list) {
    cov_subsetted <- cov_df %>% 
      filter(str_starts(covariate, starter))
    
    list_df[[starter]] <- cov_subsetted
  }
  
  cov_full_group <- bind_rows(list_df)
  return(cov_full_group)
}

# use the defined function to group the e_covs
# groups are defined as:
  # related to heat
  # related to water
  # related to plant growth
temp_group <- grouper(c("HI30", "CumHI30", "TT"))
water_group <- grouper(c("Eo", "Eos", "Es", "ESW", "Flow", "Flux", "Infiltration", "Pot", "Runoff", "SW", "T_", "Water"))
plant_growth_group <- grouper(c( "Cover", "LAI", "biomass", "yield"))

# correct for false discovery based on the above grouping
p_adjusted_temp <- p.adjust(temp_group$p_value, method = "BH")
temp_group$p_adjusted_FDR <- p_adjusted_temp

p_adjusted_water <- p.adjust(water_group$p_value, method = "BH")
water_group$p_adjusted_FDR <- p_adjusted_water

p_adjusted_growth <- p.adjust(plant_growth_group$p_value, method = "BH")
plant_growth_group$p_adjusted_FDR <- p_adjusted_growth

water_df_low_p_vals <- water_group %>% 
  filter(p_adjusted_FDR <.05)
# if desired, write the lowest p-values out to .csv
# write.csv(water_df_low_p_vals, "fdr_adjusted_water_covariates.csv", row.names = F)

# generate a barplot to compare soil water between high and low discriminability groups
sw_highs <- group_covs1 %>% 
  select(year_loc, starts_with("SW"))
sw_lows <- group_covs2 %>%
  select(year_loc, starts_with("SW"))
row_means_highs <- rowMeans(sw_highs[,-1])
row_means_lows <- rowMeans(sw_lows[,-1])
row_sds_highs <- apply(sw_highs[,-1], 1, sd)
row_sds_lows <- apply(sw_lows[,-1], 1, sd)

sw_highs$mean <- row_means_highs
sw_lows$mean <- row_means_lows
sw_highs$sd <- row_sds_highs
sw_lows$sd <- row_sds_lows
#write.csv(sw_highs, "sw_highs.csv")
#write.csv(sw_lows, "sw_lows.csv")

sw_all <- read.csv("csv_SW_covariate_vals_per_year_loc.csv")
sw_all <- sw_all %>%
  arrange(Discriminibility.Group, year_loc) %>%
  mutate(year_loc = factor(year_loc, levels = unique(year_loc)))

plot <- ggplot(sw_all, aes(x = year_loc, y = mean, fill = Discriminibility.Group)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, position = position_dodge(0.9)) +
  #scale_fill_manual(values = c("High" = "red", "Low" = "blue")) +
  labs(x = "Year Location", y = "Mean Soil Water Covariate Value", fill = "Discriminability Group") +
  ggtitle("Soil Water Values at Each Year-Location") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1),  # Rotate and shrink x-axis labels
    axis.text.y = element_text(size = rel(.5))  # Set y-axis text to two-thirds the default size
  )

# Print the plot
print(plot)

print(colnames(sw_all))

# generate a barplot to compare mean yield between high and low discriminability groups
# input: GGEBiplots input file (above)
# output: ggplot barplot

# Step 1: get the mean of each column
col_means <- colMeans(data)
# col_means <- data.frame(col_means)

# Step 2: get the standard deviation of each column
stds <- sapply(data, function(x) if(is.numeric(x)) sd(x, na.rm = TRUE))
# stds <- data.frame(stds)

# Step 3: merge the 2 vectors together to make a table
out <- cbind(col_means, stds)
# out <- data.frame(out)

# Step 4: assign each year-loc its designated 'high' or 'low' designation
groups <- read.csv('output.csv', row.names = 1)
final <- merge(out, groups, by = 0)
colnames(final)[1] <- "year_loc"

# Step 5: make the barplot
final$year_loc <- as.factor(final$year_loc)
# Create an ordering index based on group within each year_loc
final$year_loc <- factor(final$year_loc, levels = unique(final$year_loc[order(final$group)]))
ggplot(final, aes(x = year_loc, y = col_means, fill = group)) +
  geom_col() +
  geom_errorbar(aes(ymin = col_means - stds, ymax = col_means + stds), width = 0.2) +
  labs(title = "Mean Yield by Year-Location", x = "Year-Location", y = "Mean Yield", fill = "Discriminability Score") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5) , legend.position = "bottom" # Slant labels at 45 degrees
  )

