dataset6$Treatment <- gsub("\\.", " ", dataset6$Treatment)
dataset6$Treatment <- factor(dataset6$Treatment, 
                             levels = c("Forest interior", "Ecotone", "Retention clearcut"))
levels(dataset6$Treatment)

library(brms)
library(emmeans)
library(dplyr)
#control = list(adapt_delta = 0.95, max_treedepth = 15)
#save_pars = save_pars(all = TRUE)
model7 <- brm(
  formula = Number ~ Treatment * Movement.pattern + Season + Functional.group+(1 | Trap),
  data = dataset6,
  family = negbinomial(link = "log"),
  chains = 4,
  cores = 4,
  iter = 6000, control = list(adapt_delta = 0.99),
  seed = 123)

library(bayesplot)
library(ggplot2)
#Diagnostics 1: Posterior predictive checks plot
# How well model captures the distribution of the response variable?
dq<-pp_check(model7, type = "rootogram", ndraws = 200)
dq
# Diagnostics 2: Residuals vs. fitted and histogram of PPR
y_obs <- model7$data$Number 
y_pred <- posterior_predict(model7)
# Compute posterior predictive residuals
residuals <- sweep(y_pred, 2, y_obs)  
mean_resid <- colMeans(residuals)     
sd_resid <- apply(residuals, 2, sd)   
# Prepare data for plotting
fitted_vals <- fitted(model7)[, "Estimate"] 
df_plot <- data.frame(Fitted = fitted_vals, Residual = mean_resid)
# Plot Bayesian posterior predictive means: it describes the average deviation between observed and predicted values for each observation plotted against its predicted value.
dz <- ggplot(df_plot, aes(x = Fitted, y = Residual)) +
  geom_point(color = "black") +
  geom_hline(yintercept = 0, color = "red") +
  labs(
    x = "Fitted values",
    y = "Mean residuals (posterior predictive)",
    title = "Mean residuals compared to fitted values") +
  theme(
    plot.title = element_text(hjust = 0.5))

print(dz)

# Residuals were symmetrically distributed and centered around zero, indicating no systematic bias.
de<-hist(mean_resid,
         breaks = 50,
         main = "Histogram of posterior predictive residuals",
         xlab = "Mean residuals")
plot(de)
# R²
bayes_R2(model6)

tiff('Histogram.tiff',units="in",width=6,height=5,bg="white",res=300)
dz
dev.off()

dataset6$Treatment <- gsub("\\.", " ", dataset6$Treatment)
dataset6$Treatment <- factor(dataset6$Treatment, 
                             levels = c("Forest interior", "Ecotone", "Retention clearcut"))
levels(dataset6$Treatment)

emm <- emmeans(model7, ~ Movement.pattern * Treatment,
               re_formula = NA,      
               type = "response")

emm_df <- as.data.frame(emm)
d <- ggplot(emm_df, aes(x = Treatment, y = prob,
                        color = Movement.pattern,
                        group = Movement.pattern)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower.HPD, ymax = upper.HPD),  
                width = 0.2,
                linewidth = 0.8,
                position = position_dodge(width = 0.5)) +
  geom_line(aes(linetype = Movement.pattern),
            position = position_dodge(width = 0.5),
            linewidth = 0.8) +
  labs(
    x = "Treatment",
    y = "Total predicted N° of individuals"
  ) +
  scale_y_continuous(
    limits = c(0, 6),
    labels = function(x) ifelse(x == 0, "0", x)
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    strip.text = element_text(size = 15, family = "Arial"),
    strip.background = element_rect(fill = "grey90"),
    axis.text.x = element_text(size = 12, family = "Arial", angle = -45, hjust = 0, vjust = 1),
    axis.text.y = element_text(size = 12, family = "Arial"),
    axis.title.x = element_text(size = 15, family = "Arial"),
    axis.title.y = element_text(size = 15, family = "Arial", margin = margin(r = 10)),
    legend.text = element_text(size = 12, family = "Arial"),
    legend.title = element_blank(),
    legend.position = "right",
    legend.direction = "vertical",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black")
  ) +
  scale_color_manual(values = c("Across" = "grey60", "Along" = "black")) +
  scale_linetype_manual(values = c("Across" = "dashed", "Along" = "solid"))

print(d)
tiff("Total_abundance.tiff", units = "in", width = 6, height = 5, res = 300, bg = "white")
print(d)
dev.off()

# Generate emmeans for the interaction term (Movement.pattern * Treatment)
interaction_emm <- emmeans(model7, 
                           ~ Movement.pattern * Treatment)

# Perform pairwise comparisons only between `Along` and `Across` within each treatment
pairwise_interaction <- contrast(interaction_emm, 
                                 method = "pairwise", 
                                 by = "Treatment")

# Summarize the results
summary(pairwise_interaction)

#Model species richness using brms
model_full10 <- brm(formula=SpeciesRichness ~ Treatment * Movement.pattern +Season+(1 | Trap), 
                    data = species_richness_data , 
                    family = negbinomial(link = "log"),
                    chains = 4,
                    cores = 4,
                    iter = 6000, control = list(adapt_delta = 0.99),
                    seed = 123)
interaction_emm <- emmeans(model_full10, 
                           ~ Movement.pattern * Treatment)

# Perform pairwise comparisons only between `Along` and `Across` within each treatment
pairwise_interaction <- contrast(interaction_emm, 
                                 method = "pairwise", 
                                 by = "Treatment")

# Summarize the results
summary(pairwise_interaction)

emm <- emmeans(model_full10, ~ Treatment * Movement.pattern, 
               re_formula = NA, type = "response")
emm_df <- as.data.frame(emm)
ggplot(emm_df, aes(x = Treatment, y = prob, 
                   color = Movement.pattern, 
                   group = Movement.pattern)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower.HPD, ymax = upper.HPD),
                position = position_dodge(width = 0.5),
                width = 0.2, linewidth = 0.8) +
  geom_line(aes(linetype = Movement.pattern),
            position = position_dodge(width = 0.5),
            linewidth = 0.8) +
  labs(
    x = "Treatment",
    y = "Predicted Species Richness"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1, size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )
#####
# Perform abundance in each season for each functional group
library(bayestestR)
library(emmeans)
library(dplyr)
library(coda)

models_brm <- list()
functional_groups <- unique(dataset6$Functional.group)

for (group in functional_groups) {
  group_data <- dataset6 %>% filter(Functional.group == group)
  models_brm[[group]] <- list()
  
  for (season in unique(group_data$Season)) {
    season_data <- group_data %>% filter(Season == season)
    
    model <- brm(
      Number ~ Treatment * Movement.pattern + (1 | Trap),
      data = season_data,
      family = negbinomial(link = "log"),
      chains = 4,
      cores = 4,
      iter = 4000,
      control = list(adapt_delta = 0.95),
      seed = 123
    )
    
    # Save the model to disk
    filename <- paste0("model_brm_", group, "_", season, ".rds")
    saveRDS(model, file = filename)
    
    # Optional: also store in list if working interactively
    models_brm[[group]][[season]] <- model
    
    cat("Model saved for Functional Group:", group, ", and Season:", season, "\n")
  }
}

summary(readRDS("model_brm_Omnivore_Autumn.rds"))
omni_autumn <- readRDS("model_brm_Omnivore_Autumn.rds")
emm <- emmeans(omni_autumn, ~ Movement.pattern * Treatment, re_formula = NA, type = "link",posterior = TRUE)
contrasts <- contrast(emm, method = "revpairwise", by = "Treatment")
contrasts

# Step 1: Extract posterior samples from emmGrid
emm <- emmeans(omni_summer, ~ Movement.pattern * Treatment,
               re_formula = NA, type = "link", posterior = TRUE)

# Step 2: Convert to MCMC list before contrast
emm_mcmc <- as.mcmc(emm)  # ensures posterior samples are preserved

# Step 3: Apply contrast WITHIN the posterior draw object
contrasts <- contrast(emm, method = "revpairwise", by = "Treatment")

# Step 4: Manually extract posterior draws from contrast object
# `contrast()` still returns an `emmGrid`, so we force extract draws
contrast_draws <- as.mcmc(contrasts)

# Step 5: If you want per-treatment pd, split manually
treatment_levels <- unique(contrasts@grid$Treatment)
pd_summaries <- lapply(treatment_levels, function(treatment) {
  sub_contrast <- contrast(emm, method = "revpairwise", by = "Treatment") %>%
    subset(Treatment == treatment)
  
  sub_draws <- do.call(rbind, as.mcmc(sub_contrast))
  param_name <- colnames(sub_draws)[1]  # typically one column per contrast
  
  if (nrow(sub_draws) > 0 && !is.null(param_name)) {
    summary <- describe_posterior(sub_draws[, param_name], ci = 0.95, rope_range = NULL)
    summary$Parameter <- param_name
    summary$Treatment <- treatment
    return(summary)
  } else {
    NULL
  }
})

# Step 6: Combine all
pd_df <- bind_rows(pd_summaries)
print(pd_df)

#####
# Species richness
# Aggregate species richness data
species_richness_data <- dataset6 %>%
  group_by(Trap, Treatment, Season, Movement.pattern, Functional.group) %>%
  summarize(
    SpeciesRichness = n_distinct(Species),          # Number of unique species
    Abundance = sum(Number, na.rm = TRUE),          # Total abundance
    Number = sum(Number, na.rm = TRUE),             # Redundant, but retained
    .groups = "drop"
  )

# Initialize a list to store models
models_brm_species_richness <- list()

# Get unique functional groups
functional_groups2 <- unique(species_richness_data$Functional.group)

# Loop over functional groups
for (group2 in functional_groups2) {
  group_data2 <- species_richness_data %>%
    filter(Functional.group == group2)
  
  if (nrow(group_data2) == 0) {
    cat("Skipping Functional Group:", group2, "due to no data\n")
    next
  }
  
  seasons <- unique(group_data2$Season)
  models_brm_species_richness[[group2]] <- list()
  
  for (season in seasons) {
    season_data2 <- group_data2 %>%
      filter(Season == season)
    
    if (nrow(season_data2) < 5) {
      cat("Skipping Functional Group:", group2, "Season:", season, "due to insufficient data\n")
      next
    }
    
    tryCatch({
      # Step 5: Fit Bayesian model
      model <- brm(
        SpeciesRichness ~ Treatment * Movement.pattern + (1 | Trap),
        data = season_data2,
        family = negbinomial(link = "log"),
        chains = 4,
        cores = 4,
        iter = 4000,
        control = list(adapt_delta = 0.95),
        seed = 123
      )
      
      # Step 6: Save model to disk
      filename <- paste0("brm_species_richness_", group2, "_", season, ".rds")
      saveRDS(model, file = filename)
      
      # Step 7: Store model in list (optional)
      models_brm_species_richness[[group2]][[season]] <- model
      
      cat("brms model fitted and saved for Functional Group:", group2, ", Season:", season, "\n")
      
    }, error = function(e) {
      cat("Error fitting model for Functional Group:", group2, ", Season:", season, ": ", e$message, "\n")
    })
  }
}

summary(readRDS("brm_species_richness_Predator_Spring.rds"))
sp_pred_autumn <- readRDS("brm_species_richness_Predator_Autumn.rds")
emm <- emmeans(sp_pred_autumn, ~ Movement.pattern * Treatment, re_formula = NA, type = "link",posterior = TRUE)
contrasts <- contrast(emm, method = "revpairwise", by = "Treatment")
contrasts

# Step 1: Extract posterior samples from emmGrid
emm <- emmeans(sp_omni_summer, ~ Movement.pattern * Treatment,
               re_formula = NA, type = "link", posterior = TRUE)

# Step 2: Convert to MCMC list before contrast
emm_mcmc <- as.mcmc(emm)  # ensures posterior samples are preserved

# Step 3: Apply contrast WITHIN the posterior draw object
contrasts <- contrast(emm, method = "revpairwise", by = "Treatment")

# Step 4: Manually extract posterior draws from contrast object
# `contrast()` still returns an `emmGrid`, so we force extract draws
contrast_draws <- as.mcmc(contrasts)

# Step 5: If you want per-treatment pd, split manually
treatment_levels <- unique(contrasts@grid$Treatment)
pd_summaries <- lapply(treatment_levels, function(treatment) {
  sub_contrast <- contrast(emm, method = "revpairwise", by = "Treatment") %>%
    subset(Treatment == treatment)
  
  sub_draws <- do.call(rbind, as.mcmc(sub_contrast))
  param_name <- colnames(sub_draws)[1]  # typically one column per contrast
  
  if (nrow(sub_draws) > 0 && !is.null(param_name)) {
    summary <- describe_posterior(sub_draws[, param_name], ci = 0.95, rope_range = NULL)
    summary$Parameter <- param_name
    summary$Treatment <- treatment
    return(summary)
  } else {
    NULL
  }
})

# Step 6: Combine all
pd_df <- bind_rows(pd_summaries)
print(pd_df)