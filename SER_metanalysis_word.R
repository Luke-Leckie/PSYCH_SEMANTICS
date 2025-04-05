rm(list=ls())
gc()
library(gam)
library(interp)
library(cowplot)
library(ggplot2)
library(HelpersMG)
library(glmm)
library(lme4)
library(lmerTest)
library(MASS)
library(vegan)
library(multcomp)
library(segmented)
library(splines)
library(lspline)
library(caret)
library(dplyr)
library(ggbreak)
library(report)
library(rsq)
library(car)
library(emmeans)
library(segmented)
library(lme4)
library(dplyr)
library(ggbeeswarm)
library(ggplot2)
library(multcompView)


z_norm<-function(df, var){
  df_placebo<-subset(df, df$Drug=='PL')
  mean_placebo<-mean(df_placebo[[var]])
  std_placebo<-sd(df_placebo[[var]])
  df$centered_var<-df[[var]]-mean_placebo
  df$z_score<-df$centered_var/std_placebo 
  return(df)
}
#var




working_dir='/home/ll16598/Documents/POSTDOC/'
data_save_dir=paste0(working_dir, 'semantic_distance_output/')
stats_save_dir='/home/ll16598/Documents/POSTDOC/TDA/stats_output/'
image_save_dir='/home/ll16598/Documents/POSTDOC/TDA/significant_plots/'





#'PEM_df', 
z_normalise=TRUE



###############################################################################
# Example R snippet for computing Hedges' g and running a random-effects meta-analysis
###############################################################################
library(dplyr)
library(metafor)

#install.packages("clubSandwich") # if not installed
library(clubSandwich)
# Replace with your list of data frames. 

# Directory where you keep your CSVs:



analysis='word_freq'
analysis='trajectory'

if (analysis=='word_freq'){
  variable_list<- c(#'noun_entropy', 'adj_entropy', 'verb_entropy', #'noun_frequency',
                    #'adj_frequency', 
                   # 'verb_frequency',
                   # 'avg_freq_and_gap',
                   'avg_frequency',
                    'avg_gap'
                   )
  working_dir <- "/home/ll16598/Documents/POSTDOC/PSYCH_SEMANTICS/word_count_results/"
} else if (analysis=='trajectory'){
  variable_list<- c('trajectory_length',
                    'semantic_speed', 'semantic_acceleration',
                    'mean_curvature', 'recurrence_rate',#'determinism', 'laminarity','permutation_entropy',
                    'fractal_dimension',
                    'RR', 'DET', 'L',
                    'Lmax', 'DIV', 'Lentr', 'LAM', 'V', 'Vmax', 'Ventr', 'W', 'Wentr'#,#,
                  #  'DET/RR', 'LAM/DET',
                   # 'diagonal_freq', 
                  # 'vertical_freq',#,
                 #   'white_vertical_freq'
                    )
  working_dir <- "/home/ll16598/Documents/POSTDOC/trajectory_output/"
}


#

#variable_list<-'modal_participation_ratio'
{small_talk=TRUE
JUST_SER=FALSE
RANDOM_TASK=FALSE
window=80
overlap=0.1
POOL_SER=FALSE ##TO POOL THE SER EXPERIMENTAL CONDITIONS
z_normalise=TRUE
dim_reduction='50'
just_small_talk=FALSE
}
initial_var
if(just_small_talk==TRUE){
  df_names=c('MASM','cleaned_DEI')#'SER_monologs',
  
}else if(JUST_SER==TRUE){
df_names = c("SER_IPSP", "PEM_df")}else if (small_talk==TRUE){
  df_names=c(  'SER_IPSP','PEM_df','MASM','cleaned_DEI')#,'SER_monologs')
  }else{
    df_names=c('SER_IPSP','SER_monologs','PEM_df')#,
    
}

z_norm<-function(df, var){
  df[[var]]<-df[[var]]+10
  df_placebo<-subset(df, df$Drug=='PL')
  df_placebo <- df_placebo %>% filter(!is.na(.[[var]]))
  
  mean_placebo<-mean(df_placebo[[var]])
  std_placebo<-sd(df_placebo[[var]])
  df$centered_var<-df[[var]]-mean_placebo
  df$z_score<-df$centered_var/std_placebo 
  return(df)
}

remove_bracs<-function(df, var){
  if (!is.numeric(df[[var]])) {
    df[[var]] <- as.numeric(gsub("\\[|\\]", "", df[[var]]))
  }
  return(df)
}
fix_ser<-function(df_with_tda){
df_drug=read.csv('/home/ll16598/Documents/POSTDOC/Context-DATM/atom_assigned_dfs/df_monolog_244.csv')
df_with_tda$Participant<-df_with_tda$participant
df_with_tda$Session<-df_with_tda$session
df_with_tda$Participant <- gsub("[^0-9]", "", df_with_tda$Participant)
df_with_tda$Participant <- as.numeric(gsub("[^0-9]", "", df_with_tda$Participant))
df_with_tda$Drug<-df_with_tda$Drug.x
df_with_tda <- df_with_tda %>%
  left_join(df_drug, by = c("Participant", "Session"))
return(df_with_tda)}

#ser_df_with_tda$Drug
#df_SER2$Drug
all_meta=list()
#f_with_tda$D
#df_with_tda[[initial_var]]
filter_length=5
initial_var
#df_with_tda[[initial_var]]
df_with_tda$Session
meta_save_dir='/home/ll16598/Documents/POSTDOC/PSYCH_SEMANTICS/metanalysis_results/'
for (span in c(80,70,60,50)){
for(initial_var in variable_list){
  # We assume these variables are already defined somewhere above:
  # working_dir, df_names, window, overlap, z_normalise, initial_var, etc.
  # We'll create meta_data to accumulate effect sizes from each df_name
  meta_data <- data.frame()
 # print(initial_var)
  for (df_name in df_names) {
  #  print(df_name)
    ## Decide on the CSV filename as in your original code
    
    wind <- as.character(window)
    if (analysis=='word_freq'){
      csv_file <- paste0(working_dir, df_name, "_word_freq.csv")                  ##MEAN if TDA
    }else if (analysis=='trajectory'){
        csv_file <- paste0(working_dir, df_name,'_',as.character(span),'_',dim_reduction, "_utterance_trajectory_results.csv")                  ##MEAN if TDA
    }
    df_with_tda <- read.csv(csv_file)
    if(df_name=='MASM'){
      df_with_tda$Drug<-df_with_tda$condition
      df_with_tda$Drug[df_with_tda$Drug == "['PLC']"] <- "PL"
      df_with_tda$Drug[df_with_tda$Drug == "['MDMA']"] <- "MDMA"
      df_with_tda$Drug[df_with_tda$Drug == "['MA']"] <- "MA"
    }else if(df_name=='cleaned_DEI'){
      df_with_tda$Drug<-df_with_tda$condition
      df_with_tda$Drug[df_with_tda$Drug == "['PLC.txt']"] <- "PL"
      df_with_tda$Drug[df_with_tda$Drug == "['MDMA.txt']"] <- "MDMA"
      df_with_tda$Drug[df_with_tda$Drug == "['MA.txt']"] <- "MA"
      df_with_tda$participant_numbers_only <- as.numeric(gsub("\\D+", "", df_with_tda$participant))
    #  df_with_tda<-subset(df_with_tda,df_with_tda$participant_numbers_only<521)
    }else if(df_name=='SER_IPSP'){
      df_with_tda<-fix_ser(df_with_tda)
    }
    nrow(df_with_tda)
    df_with_tda<-subset(df_with_tda,df_with_tda$length>=filter_length)
    df_with_tda$Drug[df_with_tda$Drug == 0.00] <- "PL"
    
    if (POOL_SER==TRUE){
      df_with_tda$Drug[df_with_tda$Drug == 0.75] <- "1.5"}
    df_with_tda$Drug <- as.factor(df_with_tda$Drug)
    
    # Remove bracket symbols if needed, then z-score based on placebo
    df_with_tda <- remove_bracs(df_with_tda, initial_var)
    df_with_tda <- z_norm(df_with_tda, initial_var)
    nrow(df_with_tda)
    # We'll analyze the z-scored variable if z_normalise == TRUE
    if (z_normalise == TRUE) {
      var <- "z_score"
    } else {
      var <- initial_var
    }
    
    # Subset to the relevant experimental group
    if (df_name == "PEM_df" ) {
      exp_group <- "MDMA"
      df_sub <- df_with_tda %>% filter(Drug %in% c("PL", "MDMA"))
    } else if (df_name == "MASM" ) {
      exp_group <- "MDMA"
      df_sub <- df_with_tda %>% filter(Drug %in% c("PL", "MDMA"))
    } else if (df_name == "cleaned_DEI" ) {
      exp_group <- "MDMA"
      df_sub <- df_with_tda %>% filter(Drug %in% c("PL", "MDMA"))}else {
      exp_group <- "1.5"
      df_sub <- df_with_tda %>% filter(Drug %in% c("PL", "1.5"))
    }
    nrow(df_sub)
    # Compute group-level descriptive stats (PL vs. experimental)
    PL_mean <- mean(df_sub[[var]][df_sub$Drug == "PL"], na.rm = TRUE)
    PL_sd   <- sd(df_sub[[var]][df_sub$Drug == "PL"], na.rm = TRUE)
    PL_n    <- sum(!is.na(df_sub[[var]][df_sub$Drug == "PL"]))
    
    EX_mean <- mean(df_sub[[var]][df_sub$Drug == exp_group], na.rm = TRUE)
    EX_sd   <- sd(df_sub[[var]][df_sub$Drug == exp_group], na.rm = TRUE)
    EX_n    <- sum(!is.na(df_sub[[var]][df_sub$Drug == exp_group]))
    
    # Skip if missing data in either group
    if (PL_n == 0 || EX_n == 0) next
    
    # Calculate Hedges' g (unbiased standardized mean difference)
    tmp_es <- escalc(measure = "SMD",
                     m1i = EX_mean, sd1i = EX_sd, n1i = EX_n,
                     m2i = PL_mean, sd2i = PL_sd, n2i = PL_n,
                     vtype = "UB",
                     data = data.frame(id = 1))  # dummy data frame for escalc
    
    # Label this effect size row with the dataset name
    tmp_es$Study <- df_name
    
    # Assign a "EXP_ID" (cluster) label so that 
    # the first two studies share the same label,
    # and the third is different. (Adjust if needed!)
    if (df_name %in% c("SER_IPSP", "SER_monologs")) {
      tmp_es$EXP_ID <- "SER"
    } else if (df_name == "PEM_df") {
      tmp_es$EXP_ID <- "PEM"
    } else {
      tmp_es$EXP_ID <- df_name
    }
    if (df_name %in% c("SER_IPSP", "PEM_df")) {
      tmp_es$task <- "semi"
    } else if (df_name == "SER_monologs") {
      tmp_es$task <- "unstructured"
    } else if (df_name %in% c("cleaned_DEI", "MASM")) {
      tmp_es$task <- 'small-t'
    } 
    
    # Append to the master data frame
    meta_data <- rbind(meta_data, tmp_es)
  }
  
  # Now meta_data has columns: yi, vi, Study, and Participant
  # Convert Participant to a factor
  meta_data$EXP_ID <- as.factor(meta_data$EXP_ID)
  meta_data$task <- as.factor(meta_data$task)
  
 # print(meta_data)
  
  # Fit a multilevel random-effects model, accounting for the shared participant cluster
  if (RANDOM_TASK == TRUE) {
    # Try the two-random-effects model
    res <- tryCatch(
      {
        rma.mv(
          yi,
          V = vi,
          random = list(~ 1 | EXP_ID, ~ 1 | task),
          data = meta_data,
          method = "REML"
        )
      },
      error = function(e) {
        message("Error with the two-random-effects model: ", e$message)
        message("Fitting simpler model with only random intercept by study instead.")
        # Fallback: single random effect
        rma.mv(
          yi,
          V = vi,
          random = list(~ 1 | EXP_ID),
          data = meta_data,
          method = "REML"
        )
      }
    )
  } else {
    # If RANDOM_TASK is FALSE, just use the simpler model
    res <- rma.mv(
      yi,
      V = vi,
      random = list(~ 1 | EXP_ID),
      data = meta_data,
      method = "REML"
    
    )
  }
  
  # Then proceed with `res` normally, e.g. summary(res), etc.
  
  
  # Quick forest plot
  res_sum <- summary(res)
  
  # If there's only one coefficient (intercept), its p-value is usually at index [1]
  p_val <- res_sum$pval[1]
  if (p_val<0.10){
    if (analysis=='trajectory'){
  print(res_sum)
  print(initial_var)
  # Now insert that p-value into the forest plot title:
  forest(
    res,
    slab = meta_data$Study,
    main = paste0(
      initial_var,
      " (p=", round(p_val, 4), ")", span
    )
  )}else{
    print(res_sum)
    print(initial_var)
    # Now insert that p-value into the forest plot title:
    forest(
      res,
      slab = meta_data$Study,
      main = paste0(
        initial_var,
        " (p=", round(p_val, 4), ")"
      )
    )
  }
  
  }
  # Export forest-plot data if desired
  study_data <- data.frame(
    Study  = meta_data$Study,
    yi     = meta_data$yi,
    ci.lb  = meta_data$yi - 1.96 * sqrt(meta_data$vi),
    ci.ub  = meta_data$yi + 1.96 * sqrt(meta_data$vi),
    weight = weights(res)
  )
  
  overall_row <- data.frame(
    Study  = "Overall (Random-Effects)",
    yi     = res$b[1],
    ci.lb  = res$ci.lb,
    ci.ub  = res$ci.ub,
    weight = sum(weights(res))  # sum of study weights (or NA)
  )
  
  forest_data <- rbind(study_data, overall_row)
 # forest_data$span<-span
  forest_data$p_val <- rep(p_val, nrow(forest_data))
  
 # all_meta=c(all_meta,forest_data)
  name_to_save= paste0(meta_save_dir, initial_var,'_',as.character(span), "_forest_plot_data.csv")

  write.csv(
    forest_data, name_to_save, 
    row.names = FALSE)
}}

name_to_save


###


# Initialize a data frame to store the p-values, along with span and variable name
p_vals_df <- data.frame(span = numeric(),
                        variable = character(),
                        p_value = numeric(),
                        stringsAsFactors = FALSE)

library(dplyr)
library(lme4)
library(car)
library(multcomp)

# Initialize a data frame to store the p-values


z_norm<-function(df, var){
  df[[var]]<-df[[var]]+1
  df_placebo<-subset(df, df$Drug=='PL')
  df_placebo <- df_placebo %>% filter(!is.na(.[[var]]))
  
  mean_placebo<-mean(df_placebo[[var]])
  std_placebo<-sd(df_placebo[[var]])
  df$centered_var<-df[[var]]-mean_placebo
  df$z_score<-df$centered_var/std_placebo 
  return(df)
}

extract_participant_info<-function(df){
  df$Participant <- sub(".*[Ss][Ee][Rr]([0-9]{3}).*", "\\1", df$classification)
  df$Session <- sub(".*\\.s([0-9]+).*", "\\1", df$classification)
  return(df)
}

df_names=c('SER_IPSP','SER_monologs','PEM_df', 'SER_combined', 'Semi_combined')#'SER_monologs',

df_name<-df_names[[4]]
image_save_dir='/home/ll16598/Documents/POSTDOC/TDA/significant_plots/'

###
combine_SER<-function(){
  csv_file1 <- paste0(working_dir, 'SER_monologs', "_", span, "_", dim_reduction, "_utterance_trajectory_results.csv")
  df_with_tda1 <- read.csv(csv_file1)
  df_with_tda1<-extract_participant_info(df_with_tda1)
  df_with_tda1$Drug[df_with_tda1$Drug == 0.00] <- "PL"
  df_with_tda1$Drug <- factor(df_with_tda1$Drug, levels = c("PL", "0.75","1.5"))
  
  csv_file2 <- paste0(working_dir, 'SER_IPSP', "_", span, "_", dim_reduction, "_utterance_trajectory_results.csv")

  df_with_tda2 <- read.csv(csv_file2)
  df_with_tda2$Drug[df_with_tda2$Drug == 0.00] <- "PL"
  df_with_tda2$Drug <- factor(df_with_tda2$Drug, levels = c("PL", "0.75","1.5"))
  df_with_tda1$task<-'monolog'
  df_with_tda2$task<-'semi'
  common_cols <- intersect(names(df_with_tda1), names(df_with_tda2))
  df_with_tda1 <- df_with_tda1[common_cols]
  df_with_tda2 <- df_with_tda2[common_cols]
  df_with_tda1<-z_norm(df_with_tda1,var)#I normalise separately if it is combined df
  df_with_tda2<-z_norm(df_with_tda2,var)
  df_with_tda <- rbind(df_with_tda1, df_with_tda2)
  return(df_with_tda)}

combine_semi<-function(){
  csv_file1 <- paste0(working_dir, 'PEM_df', "_", span, "_", dim_reduction, "_utterance_trajectory_results.csv")
  df_with_tda1 <- read.csv(csv_file1)
  df_with_tda1$Drug <- factor(df_with_tda1$Drug, levels = c("PL", "MDMA"))
  
  csv_file2 <- paste0(working_dir, 'SER_IPSP', "_", span, "_", dim_reduction, "_utterance_trajectory_results.csv")
  df_with_tda2 <- read.csv(csv_file2)
  df_with_tda2$Drug[df_with_tda2$Drug == 0.00] <- "PL"
  df_with_tda2$Drug[df_with_tda2$Drug == 1.5] <- "MDMA"
  df_with_tda2$Drug <- factor(df_with_tda2$Drug, levels = c("PL", "MDMA"))
  df_with_tda1$task<-'monolog'
  df_with_tda2$task<-'semi'
  
  common_cols <- intersect(names(df_with_tda1), names(df_with_tda2))
  df_with_tda1 <- df_with_tda1[common_cols]
  df_with_tda2 <- df_with_tda2[common_cols]
  df_with_tda1<-z_norm(df_with_tda1,var)#I normalise separately if it is combined df
  df_with_tda2<-z_norm(df_with_tda2,var)
  df_with_tda <- rbind(df_with_tda1, df_with_tda2)
  return(df_with_tda)}


p_vals_df <- data.frame(span = numeric(),
                        variable = character(),
                        anova_p_value = numeric(),
                        tukey_PL_vs_0.75 = numeric(),
                        tukey_PL_vs_1.5  = numeric(),
                        tukey_0.75_vs_1.5 = numeric(),
                        stringsAsFactors = FALSE)
#df_with_tda$DET
for (span in c(80,70,60,50,40)) {
  if (!grepl("combined", df_name)) {
    csv_file <- paste0(working_dir, df_name, "_", span, "_", dim_reduction, "_utterance_trajectory_results.csv")
    df_with_tda <- read.csv(csv_file)
    nrow(df_with_tda)}
  
  if (df_name=='SER_monologs'){
    df_with_tda<-extract_participant_info(df_with_tda)
    df_with_tda$Drug[df_with_tda$Drug == 0.00] <- "PL"
    df_with_tda$Drug <- factor(df_with_tda$Drug, levels = c("PL", "0.75","1.5"))}else if (df_name=='SER_IPSP'){
      df_with_tda$Drug[df_with_tda$Drug == 0.00] <- "PL"
      df_with_tda$Drug <- factor(df_with_tda$Drug, levels = c("PL", "0.75","1.5"))}
  
  for (var in variable_list) {
    if (df_name=='SER_combined') {
      df_with_tda<-combine_SER()
    }else if (df_name=='Semi_combined') {
      df_with_tda<-combine_semi()
    }
  #  df_with_tda<-subset(df_with_tda,df_with_tda$nodes>3)
    # Recode and order Drug: change 0.00 to "PL" and set levels as PL, 0.75, 1.5
    
    df_with_tda <- df_with_tda %>% filter(!is.na(.[['Drug']]))
    
    df_with_tda$Participant <- as.factor(df_with_tda$Participant)
    # if (!grepl("combined", df_name)) {
    #   df_with_tda$edge_density      <- df_with_tda$edges / choose(df_with_tda$nodes, 2)
    #   df_with_tda$triangle_density  <- df_with_tda$tris  / choose(df_with_tda$nodes, 3)
    #   df_with_tda$tetra_density     <- df_with_tda$tetra / choose(df_with_tda$nodes, 4)}
    if (!grepl("combined", df_name)) {
      df_with_tda <- z_norm(df_with_tda, var)}
    
    df_with_tda <- df_with_tda %>% filter(!is.na(.[[var]]))
    
    df_with_tda$variable <- df_with_tda$z_score
    print(length(df_with_tda$variable))
    # df_with_tda$variable <- df_with_tda[[var]]
    if (!grepl("combined", df_name)) {
      model <- lmer(variable ~ Drug + (1 | Participant), data = df_with_tda)}else{
        model <- lmer(variable ~ Drug + task+(1 | Participant), data = df_with_tda)}
    anova_model <- car::Anova(model, type = "II")
    anova_p_value <- anova_model[1, 3]
    
    # Perform Tukey post hoc comparisons with BH correction
    tukey_comp <- glht(model, linfct = mcp(Drug = "Tukey"), correction = "BH")
    tukey_summary <- summary(tukey_comp)
    tukey_pvalues <- tukey_summary$test$pvalues
    
    # Extract individual p-values (names should be like "0.75 - PL", "1.5 - PL", "1.5 - 0.75")
    p_tukey_PL_vs_0.75 <- tukey_pvalues[1]
    p_tukey_PL_vs_1.5  <- tukey_pvalues[2]
    p_tukey_0.75_vs_1.5 <- tukey_pvalues[3]
    # Save the results to the data frame
    p_vals_df <- rbind(p_vals_df, data.frame(span = span,
                                             variable = var,
                                             anova_p_value = anova_p_value,
                                             tukey_PL_vs_0.75 = p_tukey_PL_vs_0.75,
                                             tukey_PL_vs_1.5  = p_tukey_PL_vs_1.5,
                                             tukey_0.75_vs_1.5 = p_tukey_0.75_vs_1.5,
                                             stringsAsFactors = FALSE))
    
    # If overall anova is significant, print results and generate plot
    if (anova_p_value < 0.05) {
      print(var)
      print(anova_model)
      print(tukey_summary)
      test_plot(model, var, df_with_tda, save_path = paste0(image_save_dir, df_name, '_', var, '_', as.character(span), '.png'))
    }
  }
}

dir_save_p='/home/ll16598/Documents/POSTDOC/PSYCH_SEMANTICS/stats_results/'
# Save the collected p-values to a CSV file
write.csv(p_vals_df, file = paste0(dir_save_p, df_name, "_p_values.csv"), row.names = FALSE)
dim_reduction=50



p_vals_df <- data.frame(span = numeric(),
                        variable = character(),
                        anova_p_value = numeric(),
                        tukey_MDMAvsMA = numeric(),
                        tukey_PL_vs_MA  = numeric(),
                        tukey_PL_vs_MDMA = numeric(),
                        e_MDMAvsMA = numeric(),
                        e_PL_vs_MA  = numeric(),
                        e_PL_vs_MDMA = numeric(),
                        stringsAsFactors = FALSE)

analysis='TDA'
if (analysis=='TDA'){
  variable_list<- c('rt',
                    #'rt_centroid', 
                    'density',
                    'edges',
                    'nodes',
                    #   'diameter',
                    # 'wighted_diameter',
                    'shortest_path_weighted','shortest_path_unweighted',
                    #  'num_triangles', 'num_tetrahedra',
                    'modularity_louvain', 'modularity_unw','spectral_gap','num_comms_unw','num_comms',
                    'clustering_coefficient', 'max_degree',
                    'mean_degree', 'max_betweenness', 'mean_betweenness', 'max_strength',
                    'mean_strength', 
                    'fiedler_value', 'largest_laplacian_eigenvalue',
                    'death_rate_dim0', 'mean_persistence_dim0',
                    'max_persistence_dim0', 'std_persistence_dim0', 'skewness_dim0',
                    'kurtosis_dim0', 'entropy_dim0', 'number_dim0',
                    
                    'birth_rate_dim1',
                    'death_rate_dim1',
                    'mean_persistence_dim1', 'max_persistence_dim1',
                    'std_persistence_dim1', 'skewness_dim1', 'kurtosis_dim1',
                    'entropy_dim1', 'number_dim1')
  #  working_dir <- "/media/ll16598/One Touch/TDA_output_filtered/sparse_0.5/"
  working_dir<-'/mnt/onetouch/TDA/TDA_output/'
  
}
test_plot<-function(model, var,sp,df_with_tda,save_path=NULL){
  print(var)
  
  # Perform Tukey post hoc comparisons on the Drug factor
  tukey_comp <- glht(model, linfct = mcp(Drug = "Tukey"))
  tukey_summary <- summary(tukey_comp)
  print(tukey_summary)
  
  # Create the combined plot
  p <- ggplot(df_with_tda, aes_string(x = "Drug", y = var)) +
    # Boxplot without default outliers
    geom_boxplot(outlier.shape = NA) +
    # Beeswarm points, colored by Participant
    geom_beeswarm(aes(color = Participant), size = 3, cex = 3) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = paste("Boxplot and Beeswarm Plot for", var, ' ', sp),
         x = "Drug",
         y = var)
  
  print(p)
  if (!is.null(save_path)) {
    # Define the file name; here we use the variable name as part of the file name
    ggsave(filename = save_path, plot = p, width = 8, height = 6)
  }
}
df_with_tda
for (span in c(80,60,40,20)){
csv_file1 <- paste0(working_dir, 'MASM', "_", span, "_", dim_reduction, "_utterance_trajectory_results.csv")
csv_file2 <- paste0(working_dir, 'cleaned_DEI', "_", span, "_", dim_reduction, "_utterance_trajectory_results.csv")
if (analysis=='TDA'){
csv_file1 <- paste0(working_dir, 'MASM',"_", span, "_mean_", dim_reduction, "_window_", analysis, "_results.csv")
csv_file2 <- paste0(working_dir, 'cleaned_DEI', "_", span, "_mean_", dim_reduction, "_window_", analysis, "_results.csv")}

df_with_tda1 <- read.csv(csv_file1)
df_with_tda1$Drug<-df_with_tda1$condition
df_with_tda1$Drug[df_with_tda1$Drug == "['PLC']"] <- "PL"
df_with_tda1$Drug[df_with_tda1$Drug == "['MDMA']"] <- "MDMA"
df_with_tda1$Drug[df_with_tda1$Drug == "['MA']"] <- "MA"



df_with_tda2 <- read.csv(csv_file2)
df_with_tda2$Drug<-df_with_tda2$condition
df_with_tda2$Drug[df_with_tda2$Drug == "['PLC.txt']"] <- "PL"
df_with_tda2$Drug[df_with_tda2$Drug == "['MDMA.txt']"] <- "MDMA"
df_with_tda2$Drug[df_with_tda2$Drug == "['MA.txt']"] <- "MA"
df_with_tda1$study<-'MASM'
df_with_tda2$study<-'DEI'
common_cols <- intersect(names(df_with_tda1), names(df_with_tda2))
df_with_tda1 <- df_with_tda1[common_cols]
df_with_tda2 <- df_with_tda2[common_cols]


for (var in variable_list) {
  df_with_tda1<-z_norm(df_with_tda1,var)#I normalise separately if it is combined df
  df_with_tda2<-z_norm(df_with_tda2,var)
  df_with_tda <- rbind(df_with_tda1, df_with_tda2)
  df_with_tda$Drug <- factor(df_with_tda$Drug, levels = c("PL", "MDMA","MA"))
  df_with_tda <- df_with_tda %>% filter(!is.na(.[['Drug']]))
  df_with_tda$Participant <- as.factor(df_with_tda$participant)
  
  df_with_tda <- df_with_tda %>% filter(!is.na(.[[var]]))
  df_with_tda$variable <- df_with_tda$z_score
  #print(length(df_with_tda$variable))
  # df_with
  model <- lmer(variable ~ Drug + (1 |study/ Participant), data = df_with_tda)
  anova_model <- car::Anova(model, type = "II")
  anova_p_value <- anova_model[1, 3]

  tukey_comp <- glht(model, linfct = mcp(Drug = "Tukey"), correction = "BH")
  tukey_summary <- summary(tukey_comp)
  
  tukey_pvalues <- tukey_summary$test$pvalues
  tukey_estimates<-tukey_summary$test$coefficients
  # Extract individual p-values (names should be like "0.75 - PL", "1.5 - PL", "1.5 - 0.75")
  p_tukey_MDMAvsMA <- tukey_pvalues[3]
  p_tukey_PL_vs_MA  <- tukey_pvalues[2]
  p_tukey_PL_vs_MDMA <- tukey_pvalues[1]
  e_tukey_MDMAvsMA <- tukey_estimates[3]
  e_tukey_PL_vs_MA  <- tukey_estimates[2]
  e_tukey_PL_vs_MDMA <- tukey_estimates[1]
  # Save the results to the data frame
  p_vals_df <- rbind(p_vals_df, data.frame(span = span,
                                           variable = var,
                                           anova_p_value = anova_p_value,
                                           tukey_MDMAvsMA = p_tukey_MDMAvsMA,
                                           tukey_PL_vs_MA  = p_tukey_PL_vs_MA,
                                           tukey_PL_vs_MDMA = p_tukey_PL_vs_MDMA,
                                           e_MDMAvsMA = e_tukey_MDMAvsMA,
                                           e_PL_vs_MA  = e_tukey_PL_vs_MA,
                                           e_PL_vs_MDMA = e_tukey_PL_vs_MDMA,
                                           stringsAsFactors = FALSE))
  
  # If overall anova is significant, print results and generate plot
  if (anova_p_value < 0.05) {
    print(var)
    print(anova_model)
    print(tukey_summary)
    test_plot(model, var,as.character(span), df_with_tda, save_path = paste0(image_save_dir, df_name, '_', var, '_', as.character(span), '.png'))
  }
}}


dir_save_p='/home/ll16598/Documents/POSTDOC/PSYCH_SEMANTICS/stats_results/'
# Save the collected p-values to a CSV file
write.csv(p_vals_df, file = paste0(dir_save_p, "small_talk_p_values.csv"), row.names = FALSE)
p_vals_df



p_vals_df <- data.frame(span = numeric(),
                        variable = character(),
                        anova_p_value = numeric(),
                        tukey_MDMAvsMA = numeric(),
                        tukey_PL_vs_MA  = numeric(),
                        tukey_PL_vs_MDMA = numeric(),
                        e_MDMAvsMA = numeric(),
                        e_PL_vs_MA  = numeric(),
                        e_PL_vs_MDMA = numeric(),
                        stringsAsFactors = FALSE)



mat <- rbind('MDMA-PL'=c(0,1,0,0,0,0),
  'MA-PL'=c(0,0,1,0,0,0),
  'MA-MDMA'=c(0,-1,1,0,0,0)
)


for (var in variable_list) {
  all_dfs <- list()
for (span in c(50,40,30,20)){
  csv_file1 <- paste0(working_dir, 'MASM', "_", span, "_", dim_reduction, "_utterance_trajectory_results.csv")
  df_with_tda1 <- read.csv(csv_file1)
  df_with_tda1$Drug<-df_with_tda1$condition
  df_with_tda1$Drug[df_with_tda1$Drug == "['PLC']"] <- "PL"
  df_with_tda1$Drug[df_with_tda1$Drug == "['MDMA']"] <- "MDMA"
  df_with_tda1$Drug[df_with_tda1$Drug == "['MA']"] <- "MA"
  
  csv_file2 <- paste0(working_dir, 'cleaned_DEI', "_", span, "_", dim_reduction, "_utterance_trajectory_results.csv")
  
  df_with_tda2 <- read.csv(csv_file2)
  df_with_tda2$Drug<-df_with_tda2$condition
  df_with_tda2$Drug[df_with_tda2$Drug == "['PLC.txt']"] <- "PL"
  df_with_tda2$Drug[df_with_tda2$Drug == "['MDMA.txt']"] <- "MDMA"
  df_with_tda2$Drug[df_with_tda2$Drug == "['MA.txt']"] <- "MA"
  df_with_tda1$study<-'MASM'
  df_with_tda2$study<-'DEI'
  common_cols <- intersect(names(df_with_tda1), names(df_with_tda2))
  df_with_tda1 <- df_with_tda1[common_cols]
  df_with_tda2 <- df_with_tda2[common_cols]

    df_with_tda1<-z_norm(df_with_tda1,var)#I normalise separately if it is combined df
    df_with_tda2<-z_norm(df_with_tda2,var)
    df_with_tda <- rbind(df_with_tda1, df_with_tda2)
    df_with_tda$Drug <- factor(df_with_tda$Drug, levels = c("PL", "MDMA","MA"))
    df_with_tda <- df_with_tda %>% filter(!is.na(.[['Drug']]))
    df_with_tda$Participant <- as.factor(df_with_tda$participant)
    df_with_tda <- df_with_tda %>% filter(!is.na(.[[var]]))
    
    df_with_tda$variable <- df_with_tda$z_score
    df_with_tda$narrative <- factor(seq_len(nrow(df_with_tda)))
    df_with_tda$sp<-span
    all_dfs <- append(all_dfs, list(df_with_tda))
}
  combined_df <- dplyr::bind_rows(all_dfs)
  model <- lmer(variable ~ Drug*sp + (1 |study/Participant/narrative), data = combined_df)
  summary(model)
  anova_model <- car::Anova(model, type = "III")
  anova_p_value <- anova_model[4, 3]
  anova_drug <- anova_model[2, 3]
  if (anova_p_value>0.05){
    model <- lmer(variable ~ Drug+sp + (1 |study/Participant/narrative), data = combined_df)
    anova_model <- car::Anova(model, type = "II")
    anova_p_value <- anova_model[1, 3]
    tukey_comp <- glht(model, linfct = mcp(Drug = "Tukey"), correction = "BH")
    tukey_summary <- summary(tukey_comp)
    
    tukey_pvalues <- tukey_summary$test$pvalues
    tukey_estimates<-tukey_summary$test$coefficients
    p_tukey_MDMAvsMA <- tukey_pvalues[3]
    p_tukey_PL_vs_MA  <- tukey_pvalues[2]
    p_tukey_PL_vs_MDMA <- tukey_pvalues[1]
    e_tukey_MDMAvsMA <- tukey_estimates[3]
    e_tukey_PL_vs_MA  <- tukey_estimates[2]
    e_tukey_PL_vs_MDMA <- tukey_estimates[1]
  }else{
    
    X <- model.matrix(~ Drug * sp, data = combined_df)
    
    # Find average row for each Drug group
    drug_levels <- levels(combined_df$Drug)
    avg_X_by_drug <- sapply(drug_levels, function(d) {
      rows <- combined_df$Drug == d
      colMeans(X[rows, , drop = FALSE])
    })
    
    avg_X_by_drug <- t(avg_X_by_drug)  
    mat <- rbind(
      "MDMA - PL" = avg_X_by_drug["MDMA", ] - avg_X_by_drug["PL", ],
      "MA - PL"   = avg_X_by_drug["MA",   ] - avg_X_by_drug["PL", ],
      "MA - MDMA" = avg_X_by_drug["MA",   ] - avg_X_by_drug["MDMA", ]
    )
    differences <-glht(model, linfct = mat, correction='BH')
    sum_dif=summary(differences)
    e_tukey_MDMAvsMA <- sum_dif$test$coefficients[[3]]
    e_tukey_PL_vs_MA  <- sum_dif$test$coefficients[[2]]
    e_tukey_PL_vs_MDMA <- sum_dif$test$coefficients[[1]]
    p_tukey_MDMAvsMA <-  sum_dif$test$pvalues[[3]]
    p_tukey_PL_vs_MA  <-  sum_dif$test$pvalues[[2]]
    p_tukey_PL_vs_MDMA <-  sum_dif$test$pvalues[[1]]
  }
  plot(model)  # residuals vs fitted
  qqnorm(resid(model)); qqline(resid(model))  # normality of residuals
  
  

  # Save the results to the data frame
  p_vals_df <- rbind(p_vals_df, data.frame(span = span,
                                           variable = var,
                                           anova_p_value = anova_p_value,
                                           tukey_MDMAvsMA = p_tukey_MDMAvsMA,
                                           tukey_PL_vs_MA  = p_tukey_PL_vs_MA,
                                           tukey_PL_vs_MDMA = p_tukey_PL_vs_MDMA,
                                           e_MDMAvsMA = e_tukey_MDMAvsMA,
                                           e_PL_vs_MA  = e_tukey_PL_vs_MA,
                                           e_PL_vs_MDMA = e_tukey_PL_vs_MDMA,
                                           stringsAsFactors = FALSE))
  
  if (anova_p_value<0.05){
    print(var)
    print(anova_model)
    print(tukey_summary)

   # test_plot(model, var,as.character(span), combined_df, save_path = paste0(image_save_dir, df_name, '_', var, '_', as.character(span), '.png'))
  
  }
   } 
  
  dir_save_p='/home/ll16598/Documents/POSTDOC/PSYCH_SEMANTICS/stats_results/'
  # Save the collected p-values to a CSV file
  write.csv(p_vals_df, file = paste0(dir_save_p, "small_talk_p_values2.csv"), row.names = FALSE)
  p_vals_df
  

if(df_name=='MASM'){
  df_with_tda$Drug<-df_with_tda$condition
  df_with_tda$Drug[df_with_tda$Drug == "['PLC']"] <- "PL"
  df_with_tda$Drug[df_with_tda$Drug == "['MDMA']"] <- "MDMA"
  df_with_tda$Drug[df_with_tda$Drug == "['MA']"] <- "MA"
}else if(df_name=='cleaned_DEI'){

  df_with_tda$participant_numbers_only <- as.numeric(gsub("\\D+", "", df_with_tda$participant))}
  
  
  
####NOW MA
JUST_SER=FALSE
df_names=c('MASM', 'cleaned_DEI','SER1')#'SER_monologs',

analysis
window=80
window='utterances'
overlap=0.1
filter_length=5
POOL_SER=FALSE ##TO POOL THE SER EXPERIMENTAL CONDITIONS
z_normalise=TRUE
#var
#df_with_tda$rt
meta_save_dir='/home/ll16598/Documents/POSTDOC/PSYCH_SEMANTICS/metanalysis_results/'
for (span in c(80,60,40,20)){
  for(initial_var in variable_list){
    # We assume these variables are already defined somewhere above:
    # working_dir, df_names, window, overlap, z_normalise, initial_var, etc.
    # We'll create meta_data to accumulate effect sizes from each df_name
    meta_data <- data.frame()
    print(initial_var)
    for (df_name in df_names) {
      print(df_name)
      ## Decide on the CSV filename as in your original code
      
      wind <- as.character(window)
      if (analysis=='word_freq'){
        csv_file <- paste0(working_dir, df_name, "_word_freq.csv")                  ##MEAN if TDA
      }else if (analysis=='trajectory'){
        csv_file <- paste0(working_dir, df_name,'_',as.character(span),'_',dim_reduction, "_mahattan_utterance_trajectory_results.csv")                  ##MEAN if TDA
      }
    df_with_tda <- read.csv(csv_file)
    if(df_name=='MASM'){
      df_with_tda$Drug<-df_with_tda$condition
      df_with_tda$Drug[df_with_tda$Drug == "['PLC']"] <- "PL"
      df_with_tda$Drug[df_with_tda$Drug == "['MDMA']"] <- "MDMA"
      df_with_tda$Drug[df_with_tda$Drug == "['MA']"] <- "MA"
    }else if(df_name=='cleaned_DEI'){
      df_with_tda$Drug<-df_with_tda$condition
      df_with_tda$Drug[df_with_tda$Drug == "['PLC.txt']"] <- "PL"
      df_with_tda$Drug[df_with_tda$Drug == "['MDMA.txt']"] <- "MDMA"
      df_with_tda$Drug[df_with_tda$Drug == "['MA.txt']"] <- "MA"
      df_with_tda$participant_numbers_only <- as.numeric(gsub("\\D+", "", df_with_tda$participant))
   #   df_with_tda<-subset(df_with_tda,df_with_tda$participant_numbers_only>520)
      
    }
    if (analysis=='TDA'){
      df_with_tda$edge_density      <- df_with_tda$edges / choose(df_with_tda$nodes, 2)
      df_with_tda$triangle_density  <- df_with_tda$tris  / choose(df_with_tda$nodes, 3)
      df_with_tda$tetra_density     <- df_with_tda$tetra / choose(df_with_tda$nodes, 4)
      
    }
    nrow(df_with_tda)
    df_with_tda<-subset(df_with_tda,df_with_tda$length>=filter_length)
    print('NOT FILTERING LENGTH')
    print(nrow(df_with_tda))
    # Convert 0 to "PL" (placebo), ensure "Drug" is a factor
    df_with_tda$Drug[df_with_tda$Drug == 0.00] <- "PL"
    
    if (POOL_SER==TRUE){
      df_with_tda$Drug[df_with_tda$Drug == 0.75] <- "1.5"}
    df_with_tda$Drug <- as.factor(df_with_tda$Drug)
    
    # Remove bracket symbols if needed, then z-score based on placebo
    df_with_tda <- remove_bracs(df_with_tda, initial_var)
    df_with_tda <- z_norm(df_with_tda, initial_var)
    nrow(df_with_tda)
    # We'll analyze the z-scored variable if z_normalise == TRUE
    if (z_normalise == TRUE) {
      var <- "z_score"
    } else {
      var <- initial_var
    }
    
    # Subset to the relevant experimental group
    if (df_name == "MASM" ) {
      exp_group <- "MA"
      df_sub <- df_with_tda %>% filter(Drug %in% c("PL", "MA"))
    } else if (df_name == "cleaned_DEI" ) {
      exp_group <- "MA"
      df_sub <- df_with_tda %>% filter(Drug %in% c("PL", "MA"))
    }  else {
      exp_group <- "20"
      df_sub <- df_with_tda %>% filter(Drug %in% c("PL", "20"))
    }
    nrow(df_sub)
    # Compute group-level descriptive stats (PL vs. experimental)
    PL_mean <- mean(df_sub[[var]][df_sub$Drug == "PL"], na.rm = TRUE)
    PL_sd   <- sd(df_sub[[var]][df_sub$Drug == "PL"], na.rm = TRUE)
    PL_n    <- sum(!is.na(df_sub[[var]][df_sub$Drug == "PL"]))
    
    EX_mean <- mean(df_sub[[var]][df_sub$Drug == exp_group], na.rm = TRUE)
    EX_sd   <- sd(df_sub[[var]][df_sub$Drug == exp_group], na.rm = TRUE)
    EX_n    <- sum(!is.na(df_sub[[var]][df_sub$Drug == exp_group]))
    
    # Skip if missing data in either group
    if (PL_n == 0 || EX_n == 0) next
    
    # Calculate Hedges' g (unbiased standardized mean difference)
    tmp_es <- escalc(measure = "SMD",
                     m1i = EX_mean, sd1i = EX_sd, n1i = EX_n,
                     m2i = PL_mean, sd2i = PL_sd, n2i = PL_n,
                     vtype = "UB",
                     data = data.frame(id = 1))  # dummy data frame for escalc
    
    # Label this effect size row with the dataset name
    tmp_es$Study <- df_name
    
    # Assign a "EXP_ID" (cluster) label so that 
    # the first two studies share the same label,
    # and the third is different. (Adjust if needed!)
    tmp_es$EXP_ID <- df_name
    
    # Append to the master data frame
    meta_data <- rbind(meta_data, tmp_es)
  }
  
  # Now meta_data has columns: yi, vi, Study, and Participant
  # Convert Participant to a factor
  meta_data$EXP_ID <- as.factor(meta_data$EXP_ID)
  
  print(meta_data)
  
  # Fit a multilevel random-effects model, accounting for the shared participant cluster
  if (RANDOM_TASK == TRUE) {
    # Try the two-random-effects model
    res <- tryCatch(
      {
        rma.mv(
          yi,
          V = vi,
          random = list(~ 1 | EXP_ID, ~ 1 | task),
          data = meta_data,
          method = "REML"
        )
      },
      error = function(e) {
        message("Error with the two-random-effects model: ", e$message)
        message("Fitting simpler model with only random intercept by study instead.")
        # Fallback: single random effect
        rma.mv(
          yi,
          V = vi,
          random = list(~ 1 | EXP_ID),
          data = meta_data,
          method = "REML"
        )
      }
    )
  } else {
    # If RANDOM_TASK is FALSE, just use the simpler model
    res <- rma.mv(
      yi,
      V = vi,
      random = list(~ 1 | EXP_ID),
      data = meta_data,
      method = "REML"
    )
  }
  
  # Then proceed with `res` normally, e.g. summary(res), etc.
  
  
  # Quick forest plot
  res_sum <- summary(res)
  
  # If there's only one coefficient (intercept), its p-value is usually at index [1]
  p_val <- res_sum$pval[1]
  if (p_val<0.10){
    print(res_sum)
    
    # Now insert that p-value into the forest plot title:
    forest(
      res,
      slab = meta_data$Study,
      main = paste0(
        initial_var, " ", as.character(span),
        " (p=", round(p_val, 4), ")"
      )
    )}
  # Export forest-plot data if desired
  study_data <- data.frame(
    Study  = meta_data$Study,
    yi     = meta_data$yi,
    ci.lb  = meta_data$yi - 1.96 * sqrt(meta_data$vi),
    ci.ub  = meta_data$yi + 1.96 * sqrt(meta_data$vi),
    weight = weights(res)
  )
  
  overall_row <- data.frame(
    Study  = "Overall (Random-Effects)",
    yi     = res$b[1],
    ci.lb  = res$ci.lb,
    ci.ub  = res$ci.ub,
    weight = sum(weights(res))  # sum of study weights (or NA)
  )
  
  forest_data <- rbind(study_data, overall_row)

  name_to_save= paste0(meta_save_dir, initial_var, "MA_forest_plot_data.csv")

  write.csv(
    forest_data, name_to_save, 
    row.names = FALSE)
}}

# 3) Combine into a single table



mat <- rbind(
  "PL vs (10,20)"      = c(0, -0.5, -0.5,  0,    0),
  "PL vs (0.75,1.5)"   = c(0,  0,    0,   -0.5, -0.5),
  #"PL vs 0.75" = c(0,  0,  0, -1,  0),
  #"PL vs 1.5"  = c(0,  0,  0,  0, -1),
  #"PL vs 10"   = c(0, -1,  0,  0,  0),
  #"PL vs 20"   = c(0,  0, -1,  0,  0),
  "10 vs 20"           = c(0,  1,   -1,    0,    0),
  "0.75 vs 1.5"        = c(0,  0,    0,    1,   -1),
  # 'hiMa vsHi MDMA' =c(0,0,-1,0,1),
  #  'loMa vslo MDMA' =c(0,-1,0,1,0)
  "(10,20) vs (0.75,1.5)" = c(0,  0.5,  0.5,  -0.5, -0.5)
)
var
# Keeps digits and decimal points

z_norm<-function(df, var){
  df[[var]]<-df[[var]]+1
  df_placebo<-subset(df, df$Drug=='PL')
  df_placebo <- df_placebo %>% filter(!is.na(.[[var]]))
  
  mean_placebo<-mean(df_placebo[[var]])
  std_placebo<-sd(df_placebo[[var]])
  df$centered_var<-df[[var]]-mean_placebo
  df$z_score<-df$centered_var/std_placebo 
  return(df)
}
var

df_names=c('SER_combined')#,'SER_monologs','PEM_df','SER1')#'SER_monologs',
#'PEM_df', 
#

colnames(df_with_tda)
z_normalise=TRUE
continuous_drug=FALSE
colnames((df_with_tda))
for (df_name in df_names){
  results_df <- data.frame(data_name=character(),Variable = character(),Step=character(),Window=character(), P_value = character(), Estimate = numeric(), Sig_Code = character(), stringsAsFactors = FALSE)
  for (overlap in c(0.1)){
    for (window in c(80)){
      step=as.character(window*overlap)
      wind=as.character(window)
      if(df_name=='SER_combined'){
        if (analysis=='syntax'){
          df_with_tda1=read.csv(paste0(working_dir, 'SER_IPSP_', analysis, '_results.csv'))
          df_with_tda2=read.csv(paste0(working_dir, 'SER1_', analysis,'_results.csv'))
        }else{
        df_with_tda1=read.csv(paste0(working_dir, 'SER_IPSP_',wind,'_',step,'_', analysis, '_results.csv'))
        df_with_tda2=read.csv(paste0(working_dir, 'SER1_',wind,'_',step,'_', analysis,'_results.csv'))
        }
        nrow(df_with_tda1)
        df_with_tda1$study<-'1'
        df_with_tda2$study<-'2'
        df_with_tda1$Drug[df_with_tda1$Drug == 0.00] <- "PL"

        
      }else{
        if (analysis=='syntax'){
          df_with_tda=read.csv(paste0(working_dir, df_name,'_', analysis,'_results.csv'))
          
        }else{
          df_with_tda=read.csv(paste0(working_dir, df_name,'_',wind,'_',step,'_', analysis,'_results.csv'))
        }
        
        print(nrow(df_with_tda))

        # if(continuous_drug==FALSE){
        #   df_with_tda$Drug <- as.factor(df_with_tda$Drug)}else{df_with_tda$Drug <- as.numeric(df_with_tda$Drug)}
      }

            # Suppose your data frame is called df, and the column is df$var

      colnames(df_with_tda)
      #df_with_tda <- df_with_tda[df_with_tda$nodes >= 3, ]
      for (var in variable_list){
        
        
        if(df_name=='SER_combined'){
          df_with_tda1<-remove_bracs(df_with_tda1,var)
          df_with_tda2<-remove_bracs(df_with_tda2,var)
          df_with_tda1<-z_norm(df_with_tda1,var)#I normalise separately if it is combined df
          df_with_tda2<-z_norm(df_with_tda2,var)
          common_cols <- intersect(names(df_with_tda1), names(df_with_tda2))
          
          df_with_tda1 <- df_with_tda1[common_cols]
          df_with_tda2 <- df_with_tda2[common_cols]
          df_with_tda <- rbind(df_with_tda1, df_with_tda2)
          nrow(df_with_tda)
          df_with_tda$Drug <- factor(df_with_tda$Drug, levels=c('PL', '10','20','0.75','1.5'))
          nrow(df_with_tda)
          
        }else{
          df_with_tda<-remove_bracs(df_with_tda,var)
          df_with_tda<-z_norm(df_with_tda,var)}
        
        if(z_normalise==TRUE){
          print('Z-score normalised based on placebo')
          df_with_tda$variable<-df_with_tda$z_score
        }else{
          print('NO Z-normalisation')
          df_with_tda$variable<-df_with_tda[[var]]}
        
        
        df_with_tda <- df_with_tda %>% filter(!is.na(.[['Drug']]))
        nrow(df_with_tda)
        df_with_tda$Participant <- as.factor(df_with_tda$Participant)
        
        
        if(df_name=='SER_combined'){
          
          model <- lmer(variable~Drug+(1|Participant)+(1|session),
                        data = df_with_tda)

        }else{
          model <- lmer(variable~Drug+(1|Participant)+(1|session),
                        data = df_with_tda)}
        anova_model <- car::Anova(model, type = "II")
        
        p_value <- anova_model[1,3]
        print(p_value)
        if (p_value<0.05){
          print(var)
          if(df_name=='SER_combined'){
            differences <-glht(model, linfct = mat, correction='BH')
            test_plot(model, 'z_score',df_with_tda,save_path=paste0(image_save_dir,df_name,'_', var,'_',step,'_',wind,'.png'))
            print(summary(differences))}else if(continuous_drug==TRUE){
              test_plot_continuous(model, var,df_with_tda)
            }  else{test_plot(model, var,df_with_tda,save_path=paste0(image_save_dir,df_name,'_', var,'_',step,'_',wind,'.png'))}
        }
        
        # Extract estimate from the model summary
        estimate <- summary(model)$coefficients[2, "Estimate"]  # Adjust if necessary
        sig_code <- ifelse(p_value < 0.001, '***',
                           ifelse(p_value < 0.01, '**',
                                  ifelse(p_value < 0.05, '*',
                                         ifelse(p_value < 0.1, '.', ' '))))
        
        
        new_row <- data.frame(data_name=df_name, Variable = var, Step=step,Window=wind,P_value = p_value, Estimate = estimate, Sig_Code = sig_code)    
        # Bind this row to the results data frame
        results_df <- rbind(results_df, new_row)
      }
    }}
  if(analysis=='sytnax'){
    stats_name=paste0(stats_save_dir,df_name,'_', analysis,'.csv')
  }else{
    stats_name=paste0(stats_save_dir,df_name,'_', analysis,'_',step,'_',wind,'.csv')
  }
  write.csv(results_df,stats_name)
  
}


z_norm_rt<-function(df, var){
  df[[var]]<-df[[var]]+1
  df_placebo<-subset(df, df$Drug=='PL')
  df_placebo <- df_placebo %>% filter(!is.na(.[[var]]))
  
  mean_placebo<-mean(df_placebo[[var]])
  std_placebo<-sd(df_placebo[[var]])
  df$centered_var<-df[[var]]-mean_placebo
  df$z_score<-df$centered_var/std_placebo 
  return(df)}

dimension_to_test<-0
if (dimension_to_test==0){
  dim_var<-'alive_dim0'
  x_var<-'scales_dim0'
} else if (dimension_to_test==1){
  dim_var<-'alive_dim1'
  x_var<-'scales_dim1'
}else if (dimension_to_test==2){
  dim_var<-'alive_dim2'
  x_var<-'scales_dim2'
}
data_direc<-'/media/ll16598/One Touch/TDA_output_filtered/TDA_output/'
SER_MA_simplices_over_time=read.csv(paste0(data_direc, 'SER1_', as.character(100), '_', as.character(10), '_',dimension_to_test,'_simplices_over_time.csv'))
SER_MDMA_simplices_over_time=read.csv(paste0(data_direc, 'SER_IPSP_', as.character(100), '_', as.character(10), '_',dimension_to_test,'_simplices_over_time.csv'))


SER_MA_simplices_over_time$study<-'1'
SER_MDMA_simplices_over_time$study<-'2'
SER_MDMA_simplices_over_time$Drug[SER_MDMA_simplices_over_time$Drug == 0.00] <- "PL"


SER_MA_simplices_over_time<-z_norm_rt(SER_MA_simplices_over_time, dim_var)
SER_MDMA_simplices_over_time<-z_norm_rt(SER_MDMA_simplices_over_time, dim_var)


common_cols <- intersect(names(SER_MA_simplices_over_time), names(SER_MDMA_simplices_over_time))

SER_MA_simplices_over_time <- SER_MA_simplices_over_time[common_cols]
SER_MDMA_simplices_over_time <- SER_MDMA_simplices_over_time[common_cols]
df_with_tda <- rbind(SER_MA_simplices_over_time, SER_MDMA_simplices_over_time)
nrow(df_with_tda)
df_with_tda$Drug <- factor(df_with_tda$Drug, levels=c('PL', '10','20','0.75','1.5'))
df_with_tda$processed_text <- factor(df_with_tda$processed_text)

nrow(df_with_tda)

colnames(df_with_tda)

df_with_tda$Participant<-as.factor(df_with_tda$Participant)
df_with_tda$Session<-as.factor(df_with_tda$Session)

df_with_tda<-df_with_tda %>% filter(!is.na(.[[x_var]]))
df_with_tda<-df_with_tda %>% filter(!is.na(.[['z_score']]))





##DIM1
model <- gamm(
  z_score ~ s(scales_dim1, by = Drug)  + Drug,
  random = list(Participant = ~1, Session = ~1, processed_text~1),
  #random = ~ 1 | Participant/Session,
  corr = corARMA(form = ~ 1 | Participant/Session, p=1,q=1),
  data = df_with_tda, method = 'ML'
)

#acf(resid(model$lme), lag.max = 36, main = "ACF of Model Residuals")
acf(residuals(model$lme,type="normalized"),main="standardized residual ACF")
# Optionally, plot the partial autocorrelation function (PACF)
pacf(resid(model$lme), lag.max = 36, main = "Partial ACF of Model Residuals")

null_model_GAMM <- gamm(
  z_score ~ s(scales_dim1)+ Drug,
  random = list(Participant = ~1, Session = ~1),
  #random = ~ 1 | Participant/Session,
  corr = corARMA(form = ~ 1 | Participant/Session, p=1,q=1),
  data = df_with_tda, method = 'ML'
)

# Extract lme objects
lme_model <- model$lme
lme_null <- null_model_GAMM$lme
# Perform likelihood ratio test
likelihood_test<-anova(lme_null, lme_model)
likelihood_test
anova(model$gam, null_model_GAMM$gam, test = "Chisq")


poly_num=5
model <- lmer(z_score~Drug*ns(scales_dim1,poly_num)+(1|Participant/Session),
              data = df_with_tda)
car::Anova(model, type=3)
hist(resid(model))
plot(model)
newdata <- expand.grid(
  Drug = unique(df_with_tda$Drug),
  scales_dim1 = seq(min(df_with_tda[[x_var]], na.rm = TRUE),
                    max(df_with_tda[[x_var]], na.rm = TRUE),
                    length.out = 100)
)

# Generate predictions using only fixed effects (re.form = NA)
newdata$pred <- predict(model, newdata = newdata, re.form = NA)
X <- model.matrix(~ Drug * ns(scales_dim1, poly_num), data = newdata)
pred_var <- diag(X %*% vcov(model) %*% t(X))
newdata$se <- sqrt(pred_var)
INTERVALS='SE'
if (INTERVALS=='CI'){
newdata$upper <- newdata$pred + 1.96 * newdata$se
newdata$lower <- newdata$pred - 1.96 * newdata$se
}else if (INTERVALS=='SE'){
newdata$upper <- newdata$pred + newdata$se
newdata$lower <- newdata$pred - newdata$se}
# Plot the predictions along with confidence bands.
library(ggplot2)
ggplot(newdata, aes(x = scales_dim1, y = pred, color = Drug, group = Drug)) +
  geom_line(size = 1) +
  # geom_jitter(data = df_with_tda, aes(x = scales_dim1, y = z_score), 
  #             alpha = 0.3, width = 0.1, height = 0, inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Drug), alpha = 0.1, color = NA) +
  # Add points for the observed mean at each timepoint
  stat_summary(data = df_with_tda, 
               aes(x = scales_dim1, y = z_score, color = Drug, group = Drug),
               fun = mean, geom = "point", size = 2) +
  labs(
    title = "Model Predicted z_score by Drug",
    x = "scales_dim0 (Time)",
    y = "Predicted z_score"
  ) +
  theme_minimal()




##DIM1

library(mgcv)

# Ensure you have a variable that represents time/order. 
# Here we assume 'Time' is available; adjust as needed.
df_with_tda$narr <- paste(df_with_tda$Participant, df_with_tda$Session, sep = "_")

df_with_tda$diff_z <- c(NA, diff(df_with_tda$z_score))
model <- gamm(
  z_score ~ s(scales_dim0, by = Drug)  + Drug,
  random = list(Participant = ~1, Session = ~1, narr=~1),
  #random = ~ 1 | Participant/Session,
  corr = corARMA(form = ~ 1 | narr, p=1,q=1),
  data = df_with_tda, method = 'ML'
)

acf(resid(model$lme), lag.max = 36, main = "ACF of Model Residuals")
acf(residuals(model$lme,type="normalized"),main="standardized residual ACF")
# Optionally, plot the partial autocorrelation function (PACF)
pacf(resid(model$lme), lag.max = 36, main = "Partial ACF of Model Residuals")

null_model_GAMM <- gamm(
  z_score ~ s(scales_dim0),
  random = list(Participant = ~1, Session = ~1, narr=~1),
  #random = ~ 1 | Participant/Session,
  corr = corARMA(form = ~ 1 | narr, p=1,q=1),
  data = df_with_tda, method = 'ML'
)

# Extract lme objects
lme_model <- model$lme
lme_null <- null_model_GAMM$lme
# Perform likelihood ratio test
likelihood_test<-anova(lme_null, lme_model)
likelihood_test
anova(model$gam, null_model_GAMM$gam, test = "Chisq")


nrow(df_with_tda)

df_with_tda<-subset(df_with_tda,df_with_tda$scales_dim0<=0.75)
poly_num=4
model <- lmer(z_score~Drug*ns(scales_dim0,poly_num)+(1|Participant/Session),
              data = df_with_tda)
res <- resid(model)
acf(res, main = "ACF of Model Residuals")

car::Anova(model, type=3)
hist(resid(model))
plot(model)

library(e1071)
residuals <- resid(model)
skew_value <- skewness(residuals)
kurt_value <- kurtosis(residuals)

# Print the results
print(paste("Skewness:", skew_value))
print(paste("Kurtosis:", kurt_value))


newdata <- expand.grid(
  Drug = unique(df_with_tda$Drug),
  scales_dim0 = seq(min(df_with_tda[[x_var]], na.rm = TRUE),
                    max(df_with_tda[[x_var]], na.rm = TRUE),
                    length.out = 100)
)

# Generate predictions using only fixed effects (re.form = NA)
newdata$pred <- predict(model, newdata = newdata, re.form = NA)
X <- model.matrix(~ Drug * ns(scales_dim0, poly_num), data = newdata)
pred_var <- diag(X %*% vcov(model) %*% t(X))
newdata$se <- sqrt(pred_var)
INTERVALS='SE'
if (INTERVALS=='CI'){
  newdata$upper <- newdata$pred + 1.96 * newdata$se
  newdata$lower <- newdata$pred - 1.96 * newdata$se
}else if (INTERVALS=='SE'){
  newdata$upper <- newdata$pred + newdata$se
  newdata$lower <- newdata$pred - newdata$se}
# Plot the predictions along with confidence bands.
library(ggplot2)
ggplot(newdata, aes(x = scales_dim0, y = pred, color = Drug, group = Drug)) +
  geom_line(size = 1) +
  # geom_jitter(data = df_with_tda, aes(x = scales_dim1, y = z_score), 
  #             alpha = 0.3, width = 0.1, height = 0, inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Drug), alpha = 0.1, color = NA) +
  # Add points for the observed mean at each timepoint
  stat_summary(data = df_with_tda, 
               aes(x = scales_dim0, y = z_score, color = Drug, group = Drug),
               fun = mean, geom = "point", size = 2) +
  labs(
    title = "Model Predicted z_score by Drug",
    x = "RT",
    y = "Predicted z_score"
  ) +
  theme_minimal()


# Load ggplot2
library(ggplot2)

# Basic line plot with points, grouped by Drug
ggplot(df_with_tda, aes(x = df_with_tda[[x_var]], y = z_score, color = Drug, group = Drug)) +
 # geom_line() +
#  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  labs(
    title = "Time Series of alive_dim0 by Drug with Smooth Mean Lines",
    x = "scales_dim0 (Time)",
    y = "alive_dim0"
  ) +
  theme_minimal()


####
nrow(df_with_tda)
poly_num=3
model <- lmer(z_score~Drug*ns(scales_dim2,poly_num)+(1|Participant/Session),
              data = df_with_tda)
car::Anova(model, type=3)
hist(resid(model))
plot(model)
newdata <- expand.grid(
  Drug = unique(df_with_tda$Drug),
  scales_dim2 = seq(min(df_with_tda[[x_var]], na.rm = TRUE),
                    max(df_with_tda[[x_var]], na.rm = TRUE),
                    length.out = 100)
)

# Generate predictions using only fixed effects (re.form = NA)
newdata$pred <- predict(model, newdata = newdata, re.form = NA)
X <- model.matrix(~ Drug * ns(scales_dim2, poly_num), data = newdata)
pred_var <- diag(X %*% vcov(model) %*% t(X))
newdata$se <- sqrt(pred_var)
INTERVALS='SE'
if (INTERVALS=='CI'){
  newdata$upper <- newdata$pred + 1.96 * newdata$se
  newdata$lower <- newdata$pred - 1.96 * newdata$se
}else if (INTERVALS=='SE'){
  newdata$upper <- newdata$pred + newdata$se
  newdata$lower <- newdata$pred - newdata$se}
# Plot the predictions along with confidence bands.
library(ggplot2)
ggplot(newdata, aes(x = scales_dim2, y = pred, color = Drug, group = Drug)) +
  geom_line(size = 1) +
  # geom_jitter(data = df_with_tda, aes(x = scales_dim1, y = z_score), 
  #             alpha = 0.3, width = 0.1, height = 0, inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Drug), alpha = 0.1, color = NA) +
  # Add points for the observed mean at each timepoint
  stat_summary(data = df_with_tda, 
               aes(x = scales_dim2, y = z_score, color = Drug, group = Drug),
               fun = mean, geom = "point", size = 2) +
  labs(
    title = "Model Predicted z_score by Drug",
    x = "scales_dim0 (Time)",
    y = "Predicted z_score"
  ) +
  theme_minimal()



####log-log
##DIM1
df_with_tda<-subset(df_with_tda,df_with_tda$scales_dim0<=0.75)
df_with_tda$log_z_score<-log(df_with_tda$z_score+2)
df_with_tda$log_scales_dim0<-log(df_with_tda$scales_dim0+1)

model <- lmer(log_z_score~Drug*log_scales_dim0+(1|Participant/Session),
              data = df_with_tda)
car::Anova(model, type=3)
hist(resid(model))
plot(model)

library(e1071)
residuals <- resid(model)
skew_value <- skewness(residuals)
kurt_value <- kurtosis(residuals)

# Print the results
print(paste("Skewness:", skew_value))
print(paste("Kurtosis:", kurt_value))


newdata <- expand.grid(
  Drug = unique(df_with_tda$Drug),
  log_scales_dim0 = seq(min(df_with_tda$log_scales_dim0, na.rm = TRUE),
                        max(df_with_tda$log_scales_dim0, na.rm = TRUE),
                        length.out = 100)
)

# Generate predictions using only fixed effects (re.form = NA)
newdata$pred <- predict(model, newdata = newdata, re.form = NA)
X <- model.matrix(~ Drug * log_scales_dim0, data = newdata)
pred_var <- diag(X %*% vcov(model) %*% t(X))
newdata$se <- sqrt(pred_var)
INTERVALS='SE'
if (INTERVALS=='CI'){
  newdata$upper <- newdata$pred + 1.96 * newdata$se
  newdata$lower <- newdata$pred - 1.96 * newdata$se
}else if (INTERVALS=='SE'){
  newdata$upper <- newdata$pred + newdata$se
  newdata$lower <- newdata$pred - newdata$se}
# Plot the predictions along with confidence bands.
library(ggplot2)
ggplot(newdata, aes(x = exp(log_scales_dim0)-1, y = exp(pred)-2, color = Drug, group = Drug)) +
  geom_line(size = 1) +
  # geom_jitter(data = df_with_tda, aes(x = scales_dim1, y = z_score), 
  #             alpha = 0.3, width = 0.1, height = 0, inherit.aes = FALSE) +
  # geom_ribbon(aes(ymin = lower, ymax = upper, fill = Drug), alpha = 0.1, color = NA) +
  # Add points for the observed mean at each timepoint
  stat_summary(data = df_with_tda, 
               aes(x = scales_dim0, y = z_score, color = Drug, group = Drug),
               fun = mean, geom = "point", size = 2) +
  labs(
    title = "Model Predicted z_score by Drug",
    x = "RT",
    y = "Predicted z_score"
  ) +
  theme_minimal()



####

colnames((df_with_tda))
for (df_name in df_names){
  results_df <- data.frame(data_name=character(),Variable = character(),Step=character(),Window=character(), P_value = character(), Estimate = numeric(), Sig_Code = character(), stringsAsFactors = FALSE)

  step=as.character(window*overlap)
  wind=as.character(window)

  df_with_tda=read.csv(paste0(data_save_dir, df_name,'_',wind,'_',step,'_distance_results.csv'))
  print(nrow(df_with_tda))
  df_with_tda$Drug[df_with_tda$Drug == 0.00] <- "PL"
  
  
  df_with_tda <- df_with_tda %>% filter(!is.na(.[['Drug']]))
  
  df_with_tda$Participant <- as.factor(df_with_tda$Participant)
  colnames(df_with_tda)
  for (var in variable_list){

      df_with_tda<-remove_bracs(df_with_tda,var)
      df_with_tda<-z_norm(df_with_tda,var)
    
    if(z_normalise==TRUE){
      print('Z-score normalised based on placebo')
      df_with_tda$variable<-df_with_tda$z_score
    }else{
      print('NO Z-normalisation')
      df_with_tda$variable<-df_with_tda[[var]]}
    
    if(df_name=='SER_combined'){
      
      model <- lmer(variable~Drug+(1|Participant),
                    data = df_with_tda)
      
      
    }else{
      model <- lmer(variable~Drug+(1|Participant),
                    data = df_with_tda)}
    anova_model <- car::Anova(model, type = "II")
    
    p_value <- anova_model[1,3]
    print(p_value)
    if (p_value<0.05){
      print(var)
      if(df_name=='SER_combined'){
        differences <-glht(model, linfct = mat, correction='BH')
        test_plot(model, var,df_with_tda,save_path=paste0(image_save_dir,df_name,'_', var,'_',step,'_',wind,'.png'))
        print(summary(differences))}else if(continuous_drug==TRUE){
          test_plot_continuous(model, var,df_with_tda)
        }  else{test_plot(model, var,df_with_tda,save_path=paste0(image_save_dir,df_name,'_', var,'_',step,'_',wind,'.png'))}
    }
    
    # Extract estimate from the model summary
    estimate <- summary(model)$coefficients[2, "Estimate"]  # Adjust if necessary
    sig_code <- ifelse(p_value < 0.001, '***',
                       ifelse(p_value < 0.01, '**',
                              ifelse(p_value < 0.05, '*',
                                     ifelse(p_value < 0.1, '.', ' '))))
    
    
    new_row <- data.frame(data_name=df_name, Variable = var, Step=step,Window=wind,P_value = p_value, Estimate = estimate, Sig_Code = sig_code)    
    # Bind this row to the results data frame
    results_df <- rbind(results_df, new_row)
  }


write.csv(results_df,paste0(stats_save_dir,df_name,'_distances_',step,'_',wind,'.csv'))

}
