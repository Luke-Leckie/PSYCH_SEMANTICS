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
test_plot<-function(model, var,df_with_tda,save_path=NULL){
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
    # Connect the same Participant's data across conditions with a faint line,
    # with the same color as the participant.
    # geom_line(aes(group = Participant, color = Participant), alpha = 0.3) +
    # Add the significance letters above each boxplot
    # geom_text(data = pos_df, aes(x = Drug, y = label_y, label = sig),
    #           color = "black", size = 5) +
    # Remove the legend and use a clean theme
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = paste("Boxplot and Beeswarm Plot for", var),
         x = "Drug",
         y = var)
  
  print(p)
  if (!is.null(save_path)) {
    # Define the file name; here we use the variable name as part of the file name
    ggsave(filename = save_path, plot = p, width = 8, height = 6)
  }
}

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


variable_list<- c('semantic_density_mean_sentence_embeddings',
                  'conversational_distance_sentence_embeddings',
                  'av_distance_mean_sentence_embeddings',
                  'av_distance_entropy_sentence_embeddings',
                  'conversational_breadth_mean_sentence_embeddings')

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

remove_bracs<-function(df, var){
  if (!is.numeric(df[[var]])) {
    df[[var]] <- as.numeric(gsub("\\[|\\]", "", df[[var]]))
  }
  return(df)
}



#'PEM_df', 
z_normalise=TRUE

window=100
window='utterances'
overlap=0.1

###############################################################################
# Example R snippet for computing Hedges' g and running a random-effects meta-analysis
###############################################################################
library(dplyr)
library(metafor)
#install.packages("clubSandwich") # if not installed
library(clubSandwich)
# Replace with your list of data frames. 

# Directory where you keep your CSVs:

# This data frame will accumulate the effect sizes from each study.
variable_list<- c('semantic_density_mean_sentence_embeddings',
                  'conversational_distance_sentence_embeddings',
                  'av_distance_mean_sentence_embeddings',
                  'av_distance_entropy_sentence_embeddings',
                  'conversational_breadth_mean_sentence_embeddings')


analysis='distance'

analysis='DMD'

if (analysis=='distance'){
  variable_list<- c('length','semantic_density_mean_sentence_embeddings',
                  'conversational_distance_sentence_embeddings',
                  'av_distance_mean_sentence_embeddings',
                  #'av_distance_entropy_sentence_embeddings',
                  'conversational_breadth_mean_sentence_embeddings')
  working_dir <- "/home/ll16598/Documents/POSTDOC/semantic_distance_output/"
} else if (analysis=='DMD'){
                    variable_list<- c('energy_top2', 'energy_top3', 'largest_abs_eig',
                                      #'num_unstable_eigs',
                                      'num_stable_eigs','modal_participation_ratio'
                                      ,'spectral_gap','decay_spectral_index','mode_amplitude_entropy',
                                      'koopman_condition_number')
                    working_dir <- "/home/ll16598/Documents/POSTDOC/DMD_output/"
                    
                    }
#variable_list<-'modal_participation_ratio'
#variable_list<- c('energy_top2')

JUST_SER=FALSE
if(JUST_SER==TRUE){
df_names = c("SER_IPSP", "PEM_df")}else{
  df_names=c('SER_IPSP','SER_monologs','PEM_df')#'SER_monologs',
  
}

POOL_SER=TRUE ##TO POOL THE SER EXPERIMENTAL CONDITIONS
z_normalise=TRUE
meta_save_dir='/home/ll16598/Documents/POSTDOC/PSYCH_SEMANTICS/metanalysis_results/'
for(initial_var in variable_list){
  # We assume these variables are already defined somewhere above:
  # working_dir, df_names, window, overlap, z_normalise, initial_var, etc.
  # We'll create meta_data to accumulate effect sizes from each df_name
  meta_data <- data.frame()
  print(initial_var)
  for (df_name in df_names) {
    
    # Decide on the CSV filename as in your original code
    
    wind <- as.character(window)
    if(wind=='utterances'){csv_file <- paste0(working_dir, df_name, "_", wind, "_", analysis, "_results.csv")}else{
      step <- as.character(window * overlap)
    csv_file <- paste0(working_dir, df_name, "_", wind, "_", step, "_", analysis, "_results.csv")}
    df_with_tda <- read.csv(csv_file)
    print(nrow(df_with_tda))
    # Convert 0 to "PL" (placebo), ensure "Drug" is a factor
    df_with_tda$Drug[df_with_tda$Drug == 0.00] <- "PL"
    
    if (POOL_SER==TRUE){
      df_with_tda$Drug[df_with_tda$Drug == 0.75] <- "1.5"}
    df_with_tda$Drug <- as.factor(df_with_tda$Drug)
    
    # Remove bracket symbols if needed, then z-score based on placebo
    df_with_tda <- remove_bracs(df_with_tda, initial_var)
    df_with_tda <- z_norm(df_with_tda, initial_var)
    
    # We'll analyze the z-scored variable if z_normalise == TRUE
    if (z_normalise == TRUE) {
      var <- "z_score"
    } else {
      var <- initial_var
    }
    
    # Subset to the relevant experimental group
    if (df_name == "PEM_df") {
      exp_group <- "MDMA"
      df_sub <- df_with_tda %>% filter(Drug %in% c("PL", "MDMA"))
    } else {
      exp_group <- "1.5"
      df_sub <- df_with_tda %>% filter(Drug %in% c("PL", "1.5"))
    }
    
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
      # fallback in case there's some other dataset
      tmp_es$EXP_ID <- df_name
    }
    
    # Append to the master data frame
    meta_data <- rbind(meta_data, tmp_es)
  }
  
  # Now meta_data has columns: yi, vi, Study, and Participant
  # Convert Participant to a factor
  meta_data$EXP_ID <- as.factor(meta_data$EXP_ID)
  
  print(meta_data)
  
  # Fit a multilevel random-effects model, accounting for the shared participant cluster
  res <- rma.mv(
    yi, 
    V = vi,
    random = ~ 1 | EXP_ID,   # random intercept by participant cluster
    data = meta_data,
    
    method = "REML"
  )

  print(summary(res))
  # summary(res_robust)
  # Quick forest plot
  forest(res, slab = meta_data$Study)
  
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
  write.csv(
    forest_data, 
    paste0(meta_save_dir, initial_var, "_", as.character(window), "_forest_plot_data.csv"), 
    row.names = FALSE
  )
}
# 3) Combine into a single table





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
