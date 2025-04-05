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

test_plot_continuous <- function(model, var, df_with_tda) {
  # Print the variable name for reference
  print(var)
  
  # Create a new data frame with a grid of Drug values over the observed range.
  # (If your model includes additional predictors, you'll need to supply fixed values for them.)
  new_data <- data.frame(Drug = seq(
    min(df_with_tda$Drug, na.rm = TRUE),
    max(df_with_tda$Drug, na.rm = TRUE),
    length.out = 100
  ))
  
  # Use the model to predict responses (and standard errors) for these new Drug values.
  preds <- predict(model, newdata = new_data, se.fit = TRUE,re.form=NA)
  new_data$fit   <- preds$fit
  new_data$upper <- preds$fit + preds$se.fit
  new_data$lower <- preds$fit - preds$se.fit
  
  # Build the plot
  p <- ggplot(df_with_tda, aes_string(x = "Drug", y = var)) +
    # Plot the raw data points, coloring by Participant
    geom_point(aes(color = 'Participant')) +
    # Connect points from the same Participant with faint lines
    geom_line(aes(group = 'Participant', color = 'Participant'), alpha = 0.3) +
    # Add the fitted model's standard error ribbon
    geom_ribbon(data = new_data,
                aes(x = Drug, y = fit, ymin = lower, ymax = upper),
                inherit.aes = FALSE,
                fill = "grey70",
                alpha = 0.3) +
    # Overlay the fitted model line
    geom_line(data = new_data,
              aes(x = Drug, y = fit),
              color = "black",
              size = 1) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = paste("Plot with Model Fit for", var),
         x = "Drug",
         y = var)
  
  # Display the plot
  print(p)
}




df_names=c('SER_IPSP','SER1', 'SER_combined')#'SER_monologs',
working_dir='/media/ll16598/One Touch/TDA/'
data_save_dir=paste0(working_dir, 'TDA_output/')
stats_save_dir='/home/ll16598/Documents/POSTDOC/TDA/stats_output/'
image_save_dir='/home/ll16598/Documents/POSTDOC/TDA/significant_plots/'



variable_list<- c('rt','rt_centroid', 'density', 'edges',
                  'tris', 'tetra', 'penta', 'euler', 'nodes',
                  #'shortest_path_weighted','shortest_path_unweighted',
                  'num_triangles', 'num_tetrahedra',
                  'modularity_louvain', 'clustering_coefficient', 'max_degree',
                  'mean_degree', 'max_betweenness', 'mean_betweenness', 'max_strength',
                  'mean_strength', 'fiedler_value', 'largest_laplacian_eigenvalue',
                  'death_rate_dim0', 'mean_persistence_dim0',
                  'max_persistence_dim0', 'std_persistence_dim0', 'skewness_dim0',
                  'kurtosis_dim0', 'entropy_dim0', 'number_dim0')
                  #'birth_rate_dim1',
                  #'death_rate_dim1', 
                  #'

                  #'mean_persistence_dim1', 'max_persistence_dim1',
                  #'std_persistence_dim1', 'skewness_dim1', 'kurtosis_dim1',
                 # 'entropy_dim1', 'number_dim1')

continuous_drug=FALSE

mat <- rbind(
  "PL vs (10,20)"      = c(0, -0.5, -0.5,  0,    0),
  "PL vs (0.75,1.5)"   = c(0,  0,    0,   -0.5, -0.5),
  "10 vs 20"           = c(0,  1,   -1,    0,    0),
  "0.75 vs 1.5"        = c(0,  0,    0,    1,   -1),
  "(10,20) vs (0.75,1.5)" = c(0,  0.5,  0.5,  -0.5, -0.5)
)


z_norm<-function(df, var){
  df_placebo<-subset(df, df$Drug=='PL')
  mean_placebo<-mean(df_placebo[[var]])
  std_placebo<-sd(df_placebo[[var]])
  df$centered_var<-df[[var]]-mean_placebo
  df$z_score<-df$centered_var/std_placebo 
  return(df)
}
#var

df_with_tda[[var]]
z_normalise=TRUE
#df_names=c( 'SER_monologs')#,'SER_IPSP','SER1','SER_combined')#'SER_monologs',
colnames(df_with_tda1)
for (df_name in df_names){
  results_df <- data.frame(data_name=character(),Variable = character(),Step=character(),Window=character(), P_value = character(), Estimate = numeric(), Sig_Code = character(), stringsAsFactors = FALSE)
  for (overlap in c(0.1)){
    for (window in c(100)){
      step=as.character(window*overlap)
      wind=as.character(window)
      if(df_name=='SER_combined'){
        df_with_tda1=read.csv(paste0(data_save_dir, 'SER_IPSP_',wind,'_',step,'_TDA_results.csv'))
        df_with_tda2=read.csv(paste0(data_save_dir, 'SER1_',wind,'_',step,'_TDA_results.csv'))
        nrow(df_with_tda1)
        df_with_tda1$study<-'1'
        df_with_tda2$study<-'2'
        df_with_tda1$Drug[df_with_tda1$Drug == 0.00] <- "PL"
        
        
        
        #df_with_tda1 <- df_with_tda1[df_with_tda1$Drug != '0.75', ]
        # df_with_tda2 <- df_with_tda2[df_with_tda2$Drug != '10', ]
        common_cols <- intersect(names(df_with_tda1), names(df_with_tda2))
        
        
        # Subset both dataframes to keep only common columns
        df_with_tda1 <- df_with_tda1[common_cols]
        df_with_tda2 <- df_with_tda2[common_cols]
        # Bind the data frames together
        
        
      }else{
        df_with_tda=read.csv(paste0(data_save_dir, df_name,'_',wind,'_',step,'_TDA_results.csv'))
        print(nrow(df_with_tda))
        df_with_tda$Drug[df_with_tda$Drug == 0.00] <- "PL"
        
        
        
        if(continuous_drug==FALSE){
          df_with_tda$Drug <- as.factor(df_with_tda$Drug)}else{df_with_tda$Drug <- as.numeric(df_with_tda$Drug)}
      }
      continuous_drug=FALSE#whether drug should be treatd as a factor or continuous
      # Suppose your data frame is called df, and the column is df$var
      df_with_tda <- df_with_tda %>% filter(!is.na(.[['Drug']]))
      
      df_with_tda$Participant <- as.factor(df_with_tda$Participant)
      
      
      colnames(df_with_tda)
      #df_with_tda <- df_with_tda[df_with_tda$nodes >= 3, ]
      for (var in variable_list){
        if(df_name=='SER_combined'){
        df_with_tda1<-z_norm(df_with_tda1,var)#I normalise separately if it is combined df
        df_with_tda2<-z_norm(df_with_tda2,var)
        df_with_tda <- rbind(df_with_tda1, df_with_tda2)
        df_with_tda$Drug <- factor(df_with_tda$Drug, levels=c('PL', '10','20','0.75','1.5'))
        }else{
        df_with_tda<-z_norm(df_with_tda,var)}
        
        if(z_normalise==TRUE){
          print('Z-score normalised based on placebo')
          df_with_tda$variable<-df_with_tda$z_score
        }else{
        print('NO Z-normalisation')
        df_with_tda$variable<-df_with_tda[[var]]}

         model <- lmer(variable~Drug+(1|Participant),
                        data = df_with_tda)
        anova_model <- car::Anova(model, type = "II")
        
        p_value <- anova_model[1,3]
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
    }}

write.csv(results_df,paste0(stats_save_dir,df_name,'_TDA_',step,'_',wind,'.csv'))

}



df_with_tda_2<-subset(df_with_tda, df_with_tda$Drug=='PL')
model <- lmer(nodes~Drug+(1|Participant),
              data = df_with_tda2)

car::Anova(model, type=2)


working_dir='/home/ll16598/Documents/POSTDOC/'
data_save_dir=paste0(working_dir, 'semantic_distance_output/')
stats_save_dir='/home/ll16598/Documents/POSTDOC/TDA/stats_output/'
image_save_dir='/home/ll16598/Documents/POSTDOC/TDA/significant_plots/'


variable_list<- c('semantic_density_mean_sentence_embeddings',
                  'conversational_distance_sentence_embeddings',
                  'av_distance_mean_sentence_embeddings',
                 'av_distance_entropy_sentence_embeddings',
                  'conversational_breadth_mean_sentence_embeddings')

mat<-rbind('PlmMDMA'=c(0,0,0,-1,-1),
'PlmMA'=c(0,-1,-1,0,0),
'hlMDMA'=c(0,0,0,1,-1),
'hlMA'=c(0,1,-1,0,0),
'MDMAmMA'=c(0,-1,-1,1,1))

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
  #df$z_score<-df$centered_var
  return(df)
}

remove_bracs<-function(df, var){
  if (!is.numeric(df[[var]])) {
  df[[var]] <- as.numeric(gsub("\\[|\\]", "", df[[var]]))
  }
  return(df)
}

df_names=c('SER_IPSP','SER_monologs','PEM_df','SER1')#'SER_monologs',

df_names=c('SER_IPSP')#'SER_monologs',

#'PEM_df', 
z_normalise=TRUE
colnames((df_with_tda))
for (df_name in df_names){
  results_df <- data.frame(data_name=character(),Variable = character(),Step=character(),Window=character(), P_value = character(), Estimate = numeric(), Sig_Code = character(), stringsAsFactors = FALSE)
  for (overlap in c(0.1)){
    for (window in c(100)){
      step=as.character(window*overlap)
      wind=as.character(window)
      if(df_name=='SER_combined'){
        df_with_tda1=read.csv(paste0(data_save_dir, 'SER_IPSP_',wind,'_',step,'_distance_results.csv'))
        df_with_tda2=read.csv(paste0(data_save_dir, 'SER1_',wind,'_',step,'_distance_results.csv'))
        nrow(df_with_tda1)
        df_with_tda1$study<-'1'
        df_with_tda2$study<-'2'
        df_with_tda1$Drug[df_with_tda1$Drug == 0.00] <- "PL"
        
        
        
      }else{
        df_with_tda=read.csv(paste0(data_save_dir, df_name,'_',wind,'_',step,'_distance_results.csv'))
        print(nrow(df_with_tda))
        df_with_tda$Drug[df_with_tda$Drug == 0.00] <- "PL"
        
        
        
        if(continuous_drug==FALSE){
          df_with_tda$Drug <- as.factor(df_with_tda$Drug)}else{df_with_tda$Drug <- as.numeric(df_with_tda$Drug)}
      }
      continuous_drug=FALSE#whether drug should be treatd as a factor or continuous
      # Suppose your data frame is called df, and the column is df$var
      df_with_tda <- df_with_tda %>% filter(!is.na(.[['Drug']]))
      
      df_with_tda$Participant <- as.factor(df_with_tda$Participant)
      
      
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
          df_with_tda$Drug <- factor(df_with_tda$Drug, levels=c('PL', '10','20','0.75','1.5'))
        
          }else{
          df_with_tda<-remove_bracs(df_with_tda,var)
          df_with_tda<-z_norm(df_with_tda,var)}
        
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
    }}
  
  write.csv(results_df,paste0(stats_save_dir,df_name,'_distances_',step,'_',wind,'.csv'))
  
}
####



df_names=c('SER_combined', 'SER_monologs','SER_IPSP','PEM_df','SER1')#'SER_monologs',

working_dir='/home/ll16598/Documents/POSTDOC/'
data_save_dir=paste0(working_dir, 'DMD_output/')
stats_save_dir='/home/ll16598/Documents/POSTDOC/DMD_output/stats_output/'
image_save_dir='/home/ll16598/Documents/POSTDOC/DMD_output/significant_plots/'
continuous_drug=FALSE

######DISTANCE
variable_list<- c('energy_top2', 'energy_top3', 'largest_abs_eig',
                  'dominant_frequency','num_unstable_eigs','num_stable_eigs','modal_participation_ratio'
                  ,'spectral_gap','decay_spectral_index','mode_amplitude_entropy',
                  'koopman_condition_number')
variable_list<- c('energy_top2', 'energy_top3','mode_amplitude_entropy','num_stable_eigs')
df_names=c('SER_IPSP','PEM_df','SER1', 'SER_combined','PEM_df')#, 'SER_monologs','PEM_df')#'SER_monologs',
#df_names=c('SER_IPSP')#, 'SER_monologs','PEM_df')#'SER_monologs',
df_names=c('SER1','SER_IPSP','PEM_df')#'SER_monologs',
z_normalise=TRUE

df_with_tda[[var]]
for (df_name in df_names){
  results_df <- data.frame(data_name=character(),Variable = character(),Step=character(),Window=character(), P_value = character(), Estimate = numeric(), Sig_Code = character(), stringsAsFactors = FALSE)
  for (overlap in c(0.1)){
    for (window in c(100)){
      step=as.character(window*overlap)
      wind=as.character(window)
      if(df_name=='SER_combined'){
        df_with_tda1=read.csv(paste0(data_save_dir, 'SER_IPSP_',wind,'_',step,'_DMD_results.csv'))
        df_with_tda2=read.csv(paste0(data_save_dir, 'SER1_',wind,'_',step,'_DMD_results.csv'))
        nrow(df_with_tda1)
        df_with_tda1$study<-'1'
        df_with_tda2$study<-'2'
        df_with_tda1$Drug[df_with_tda1$Drug == 0.00] <- "PL"
        
        
        
        #df_with_tda1 <- df_with_tda1[df_with_tda1$Drug != '0.75', ]
        # df_with_tda2 <- df_with_tda2[df_with_tda2$Drug != '10', ]
        common_cols <- intersect(names(df_with_tda1), names(df_with_tda2))
        
        
        # Subset both dataframes to keep only common columns
        df_with_tda1 <- df_with_tda1[common_cols]
        df_with_tda2 <- df_with_tda2[common_cols]
        # Bind the data frames together
        
        
      }else{
        df_with_tda=read.csv(paste0(data_save_dir, df_name,'_',wind,'_',step,'_DMD_results.csv'))
        print(nrow(df_with_tda))
        df_with_tda$Drug[df_with_tda$Drug == 0.00] <- "PL"
        
        
        
        if(continuous_drug==FALSE){
          df_with_tda$Drug <- as.factor(df_with_tda$Drug)}else{df_with_tda$Drug <- as.numeric(df_with_tda$Drug)}
      }
      continuous_drug=FALSE#whether drug should be treatd as a factor or continuous
      # Suppose your data frame is called df, and the column is df$var
      df_with_tda <- df_with_tda %>% filter(!is.na(.[['Drug']]))
    #  df_with_tda <- df_with_tda %>% filter(!is.na(.[[var]]))
      
      df_with_tda$Participant <- as.factor(df_with_tda$Participant)
      
      
      colnames(df_with_tda)
      #df_with_tda <- df_with_tda[df_with_tda$nodes >= 3, ]
      for (var in variable_list){
        if(df_name=='SER_combined'){
          df_with_tda1<-z_norm(df_with_tda1,var)#I normalise separately if it is combined df
          df_with_tda2<-z_norm(df_with_tda2,var)
          df_with_tda <- rbind(df_with_tda1, df_with_tda2)
          df_with_tda$Drug <- factor(df_with_tda$Drug, levels=c('PL', '10','20','0.75','1.5'))
        }else{
          df_with_tda<-z_norm(df_with_tda,var)}
        
        if(z_normalise==TRUE){
          print('Z-score normalised based on placebo')
          df_with_tda$variable<-df_with_tda$z_score
        }else{
          print('NO Z-normalisation')
          df_with_tda$variable<-df_with_tda[[var]]}
          
        
        model <- lmer(variable~Drug+(1|Participant),
                        data = df_with_tda)
        anova_model <- car::Anova(model, type = "II")
        
        p_value <- anova_model[1,3]
        #print(p_value)
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
    }}
  
  write.csv(results_df,paste0(stats_save_dir,df_name,'_DMD_',step,'_',wind,'.csv'))
  
}


####encoder




###########
summary(model)

var
df_names=c('monolog', 'SER_IPSP')
working_dir='/home/ll16598/Documents/POSTDOC/Context-DATM/'
data_save_dir=paste0(working_dir, 'TDA_output/')

window<-'100'
step<-'20'

df_name=df_names[2]
df_betti=read.csv(paste0(data_save_dir, df_name,'_',window,'_',step,'_0_simplices_over_time.csv'))
df_betti$scales_dim0<-as.numeric(df_betti$scales_dim0)
df_betti$alive_dim0<-as.numeric(df_betti$alive_dim0)
df_betti <- df_betti %>% filter(!is.na(.[['scales_dim0']]))
df_betti <- df_betti %>% filter(!is.na(.[['alive_dim0']]))

continuous_drug=FALSE#whether drug should be treatd as a factor or continuous
if(continuous_drug==FALSE){
  df_betti$Drug <- as.factor(df_betti$Drug)}else{df_betti$Drug <- as.numeric(df_betti$Drug)}
df_betti <- df_betti %>% filter(!is.na(.[['Drug']]))
df_betti$Participant <- as.factor(df_betti$Participant)

colnames(df_betti)
model <- lmer(alive_dim0~Drug*poly(scales_dim0,2)+(1|Participant),
              data = df_betti)
anova_model <- car::Anova(model, type = "III")
anova_model
summary(model)
hist(resid(model))
plot(model)

df_betti=read.csv(paste0(data_save_dir, df_name,'_',window,'_',step,'_1_simplices_over_time.csv'))

continuous_drug=FALSE#whether drug should be treatd as a factor or continuous
if(continuous_drug==FALSE){
  df_betti$Drug <- as.factor(df_betti$Drug)}else{df_betti$Drug <- as.numeric(df_betti$Drug)}
df_betti <- df_betti %>% filter(!is.na(.[['Drug']]))
df_betti$Participant <- as.factor(df_betti$Participant)

colnames(df_betti)
model <- lmer(log(alive_dim1+1)~Drug+log(scales_dim1+1)+(1|Participant),
              data = df_betti)
anova_model <- car::Anova(model, type = "II")
anova_model
summary(model)

# window<-'100'
# step<-'20'
df_betti=read.csv(paste0(data_save_dir, df_name,'_',window,'_',step,'_centroid_simplices_over_time.csv'))
continuous_drug=FALSE#whether drug should be treatd as a factor or continuous
if(continuous_drug==FALSE){
  df_betti$Drug <- as.factor(df_betti$Drug)}else{df_betti$Drug <- as.numeric(df_betti$Drug)}
df_betti <- df_betti %>% filter(!is.na(.[['Drug']]))
df_betti$Participant <- as.factor(df_betti$Participant)
df_betti$centroid_alive_dim0
model <- lmer(log(centroid_alive_dim0+1)~Drug+log(centroid_scales_dim0+1)+(1|Participant),
              data = df_betti)
anova_model <- car::Anova(model, type = "II")
anova_model
summary(model)
hist(resid(model))
plot(model)
###


dream_dat=read.csv('/home/ll16598/Documents/altered_states/2024/final_network_results/corpus_results_130_18.csv')
#dat=read.csv('/home/ll16598/Documents/altered_states/2024/stats/df_stats.csv')

head(dream_dat)
na.omit(dream_dat)
dat_noblind<-subset(dream_dat, blind=='False')
dat<-subset(dream_dat, veterans=='False')
dat<-subset(dat, blind=='False')

variable_list<-c('Efficiency','Diameter','Clustering', 'Feedback.loops', 'DH','Transitivity','AP','Density', 'gini', 'Modularity')
dat_noblind[[var]]
for (var in variable_list){
  dat_noblind$variable<-dat_noblind[[var]]
  if(var=='Feedback.loops'){
    dat$variable<-log(dat$variable+1)}
  print(var)
  model <- lmer(median_Valence~variable+(1|dreamer)+(1|Nodes),
                        data = dat_noblind)
  print(car::Anova(model, type=3))
  print(summary(model))
}

dat=read.csv('/home/ll16598/Documents/altered_states/2024/final_network_results/corpus_results_130_DI_cosine.csv')
#dat$s_to_e_path
results_df <- data.frame(Response = character(), Variable = character(), P_value = character(), Estimate = numeric(), Sig_Code = character(), stringsAsFactors = FALSE)
response='median_Valence'
for (var in variable_list){
  dat$variable<-dat[[var]]
  # if(var=='s_to_e_path'){
  #   dat$variable<-log(dat$variable+1)}
  print(var)
  model <- lmer(median_Valence~variable+(1|dreamer)+(1|Nodes),
                data = dat)
  print(hist(resid(model)))
  anova_model <- car::Anova(model, type = "II")
  
  p_value <- anova_model[1,3]
  print(p_value)
  
  # Extract estimate from the model summary
  estimate <- summary(model)$coefficients[2, "Estimate"]  # Adjust if necessary
  sig_code <- ifelse(p_value < 0.001, '***',
                     ifelse(p_value < 0.01, '**',
                            ifelse(p_value < 0.05, '*',
                                   ifelse(p_value < 0.1, '.', ' '))))
  
  # Combine results into a data frame row
  new_row <- data.frame(Response = response, Variable = var, P_value = p_value, Estimate = estimate, Sig_Code = sig_code)    
  # Bind this row to the results data frame
  results_df <- rbind(results_df, new_row)
}

results_df

write.csv(results_df, '/home/ll16598/Documents/altered_states/2024/network_results_130_valence.csv')

dat_deg_val <- read.csv('/home/ll16598/Documents/altered_states/2024/valence_quant.csv', header=T, stringsAsFactors = F)
dat_deg_val$log_degree=log(dat_deg_val$degree)

model <- lmer(log_degree~prevalence*quantile+(1|node_id),
              data = dat_deg_val)

hist(resid(model))
car::Anova(model, type='III')
summary(model)
dat<-dat_deg_val
degree_range <- max(dat$degree) - min(dat$degree)
prevalence_range <- max(dat$prevalence) - min(dat$prevalence)

newdat <- expand.grid(
  prevalence = seq(min(dat$prevalence) - 0.05 * prevalence_range, max(dat$prevalence) + 0.05 * prevalence_range, length.out = 100),
  quantile = unique(dat$quantile)
)

# Predict using the new data frame
newdat$log_degree <- predict(model, newdat, re.form=NA)
newdat$degree<-exp(newdat$log_degree)
  
dat_aro_val <- read.csv('/home/ll16598/Documents/altered_states/2024/arousal_quant.csv', header=T, stringsAsFactors = F)
dat_aro_val$log_degree=log(dat_aro_val$degree)
model_aro <- lmer(log_degree~prevalence*quantile+(1|node_id),
              data = dat_aro_val)
hist(resid(model_aro))
car::Anova(model_aro, type='III')
summary(model_aro)

for (var in variable_list){
  dat$variable<-dat[[var]]
  print(var)
  model <- lmer(variable~blind+ (1|dreamer)+(1|gender),
                data = dat)
  print(car::Anova(model, type=3))
  print(summary(model))
}
for (var in variable_list){
  dat$variable<-dat[[var]]
  print(var)
  model <- lmer(variable~median_Valence*Arousal+ (1|dreamer)+(1|gender),
                data = dat)
  print(car::Anova(model, type=3))
  print(summary(model))
}

for (var in variable_list){##Femalesmore nodes and average path
  dat$variable<-dat[[var]]
  print(var)
  model <- lmer(variable~gender+ (1|dreamer),
                data = dat)
  print(car::Anova(model, type=3))
}
for (var in variable_list){
  dat$variable<-dat[[var]]
  print(var)
  model <- lmer(variable~Arousal+ (1|dreamer)+(1|gender),
                data = dat)
  print(car::Anova(model, type=3))
  print(summary(model))
}


for (var in variable_list){
  dream_dat$variable<-dream_dat[[var]]
  print(var)
  model <- lmer(variable~blind+ (1|dreamer)+(1|gender),
                data = dream_dat)
  print(car::Anova(model, type=3))
}
for (var in variable_list){
  dream_dat$variable<-dream_dat[[var]]
  print(var)
  model <- lmer(variable~dreamer+(1|gender),
                data = dream_dat)
  print(car::Anova(model, type=3))
}

model <- lmer(median_Valence~gender+ (1|dreamer),
              data = dat)
print(car::Anova(model, type=3))
model <- lmer(Arousal~gender+ (1|dreamer),
              data = dat)
print(car::Anova(model, type=3))
model <- lmer(Valence~dreamer+ (1|dreamer),
              data = dat)
print(car::Anova(model, type=3))
model <- lmer(Valence~blind+ (1|dreamer),
              data = dream_dat)
print(car::Anova(model, type=3))

dream_dat$dreamer
model_sqrt_v <- lmer( Nodes~ blind+ (1|dreamer),
                     data = dream_dat)
model_sqrt_v <- lmer( Nodes~ gender+ (1|dreamer),
                      data = dream_dat)
car::Anova(model_sqrt_v, type=3)
summary(model_sqrt_v)
