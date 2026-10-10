########Random forest######################
library(terra)
library(caret)
library(randomForest)


######Data
setwd("D:/DSM/tiff")

#####Stack
covStack <- rast(c(
  'Aspect.tif',
  #'cai.tif',
  #'Clayi.tif',
  'Dem.tif',
  #'GI.tif',
  'LST.tif',
  'MRRTF.tif',
  'MRVBF.tif',
  'MSI.tif',
  #'NDSI.tif',
  #'NDVI.tif',
  #'NDWI.tif',
  #'plan_crv.tif',
  'Profile_crv.tif',
  'Rain.tif',
  'River_dis.tif',
  'Road_dis.tif',
  'Slope.tif',
  #'SPI.tif',
  #'Tmax.tif',
  'Tmin.tif',
  #'TPI.tif',
  #'TRI.tif',
  #'TWI.tif',
  #'Valley_depth.tif',
  'Sand.tif',
  #'Clay.tif',
  'Silt.tif'
))
###########DEPENDENT Variable###############
pH_data <- read.table("soil.txt", header = TRUE, sep = "\t", na.strings = "NA", dec = ".", strip.white = TRUE)
str(pH_data)
######################Test and Train 
#set.seed(123)
training <- sample(nrow(pH_data), 0.8 * nrow(pH_data))

#################RF
ctrl <- trainControl(method = "boot", number = 10)

# Fit the model using caret (cross-validation is handled automatically).
rf_model <- randomForest(
  pH ~  Aspect+ Dem + LST + MRVBF +MRRTF+MSI+ 
    Profile_crv + River_dis + Road_dis + Sand + Silt +
    Slope + Rain + Tmin,
  data = pH_data,
  trControl = ctrl,
  importance = TRUE,
  ntree = 1000 
)
print(rf_model)
varImpPlot(rf_model, main = "Variable Importance - Random Forest pH Model")
importance_df <- data.frame(importance(rf_model))
importance_df$Variable <- rownames(importance_df)
importance_df <- importance_df[order(importance_df$IncNodePurity, decreasing = TRUE), ]

#####Validati

goof2 <- function(observed, predicted, title = "Model Validation") {
  observed <- as.numeric(observed)
  predicted <- as.numeric(predicted)
  
  # Remove missing (NA) values.
  valid <- complete.cases(observed, predicted)
  observed <- observed[valid]
  predicted <- predicted[valid]
  
  # Calculate evaluation metrics.
  RMSE <- sqrt(mean((observed - predicted)^2))
  MAE  <- mean(abs(observed - predicted))
  R2   <- cor(observed, predicted)^2
  BIAS <- mean(predicted - observed)
  NSE  <- 1 - (sum((observed - predicted)^2) / sum((observed - mean(observed))^2))
  RPD  <- sd(observed) / RMSE
  
 
  results <- data.frame(
    RMSE = RMSE,
    MAE  = MAE,
    R2   = R2,
    NSE  = NSE,
    BIAS = BIAS,
    RPD  = RPD
  )
  print(results)
  
  # 
  plot(observed, predicted,
       pch = 19, col = "steelblue",
       xlab = "Observed pH",
       ylab = "Predicted pH",
       main = title)
  abline(0, 1, col = "red", lwd = 2)  
  legend("topleft",
         legend = c(
           paste0("R² = ", round(R2, 3)),
           paste0("RMSE = ", round(RMSE, 3)),
           paste0("MAE = ", round(MAE, 3))
         ),
         bty = "n")
  
  return(invisible(results))
}

# =====================================================================
#  Validation (Training)
# =====================================================================
pred_train <- predict(rf_model, newdata = pH_data[training, ])
cat(" Internal Validation (Training Data)\n")
goof2(observed = pH_data$pH[training], predicted = pred_train,
      title = "Random Forest - Training (Internal Validation)")


# =====================================================================
# Validation (Testing)
# =====================================================================
pred_test <- predict(rf_model, newdata = pH_data[-training, ])
cat(" External Validation (Testing Data)\n")
goof2(observed = pH_data$pH[-training], predicted = pred_test,
      title = "Random Forest - Testing (External Validation)")

val_train <- goof2(pH_data$pH[training], pred_train)
val_test  <- goof2(pH_data$pH[-training], pred_test)
# =====================================================================
# =====================================================================
# Install the iml package if it is not already installed, then load the required package.
if (!requireNamespace("iml", quietly = TRUE)) install.packages("iml")
library(lime)
library(randomForest)  


model_type.randomForest <- function(x, ...) {
  return("regression")  # 
}

predict_model.randomForest <- function(x, newdata, ...) {
  data.frame(Response = predict(x, newdata = newdata))
}

# =====================================================================
# 
# Remove the target variable from the predictor data (pH in this example).
train_x <- pH_data[training, -which(names(pH_data) == "pH")]

# =====================================================================
# Explainer
explainer <- lime(
  x = train_x,
  model = rf_model
)

# =====================================================================
#  Sample Selection
# MAX 
case_max <- pH_data[which.max(pH_data$pH), !(names(pH_data) %in% c("EC","pH", "x", "y","SPI"))]
case_max

# MIN 
case_min <- pH_data[which.min(pH_data$pH), !(names(pH_data) %in% c("EC","pH", "x", "y","SPI"))]
case_min

# =====================================================================
# LIME
explanation_max <- explain(
  case_max,
  explainer,
  n_labels = 1,
  n_features = 15,   # Number of the most important features to display.
  kernel_width = 0.5
)

explanation_max
explanation_min <- explain(
  case_min,
  explainer,
  n_labels = 1,
  n_features = 15,
  kernel_width = 0.5
)


# =====================================================================
# Result LIME
par(mfrow = c(1, 2)) 
plot_features(explanation_max)
title(" Highest pH case (LIME explanation)")
plot_features(explanation_min)
title(" Lowest pH case (LIME explanation)")
par(mfrow = c(1, 1))
##################Prediction
pH_map <- predict(covStack, rf_model, filename = "pH_RF_Map.tif",
                  format = "GTiff", datatype = "FLT4S", overwrite = TRUE)
plot(pH_map, main = "Predicted Soil pH (Random Forest)")




X_train <- pH_data[training, !(names(pH_data) %in% c("pH","EC","x","y","SPI"))]
predictor_rf <- Predictor$new(
  model = rf_model,
  data = X_train,
  y = pH_data$pH[training]
)

library (iml)

# Create a predictor object for iml.
predictor_rf <- Predictor$new(
  model = rf_model,
  data = X_train,
  y = pH_data$pH[training]
)

# Select the case with the highest pH value.
case_max <- X_train[which.max(pH_data$pH[training]), , drop = FALSE]
# Select the case with the lowest pH value.
case_min <- X_train[which.min(pH_data$pH[training]), , drop = FALSE]

### Calculate SHAP values for the maximum and minimum pH cases.
shap_max <- Shapley$new(predictor_rf, x.interest = case_max)
shap_min <- Shapley$new(predictor_rf, x.interest = case_min)

# Plot feature contributions.
plot(shap_max) + ggtitle("Highest pH Case - SHAP Explanation")
plot(shap_min) + ggtitle("Lowest pH Case - SHAP Explanation")

#########################################################
###  SHAP Summary (Global Importance)
#########################################################
library (iml)
# Calculate SHAP values for a sample of observations (e.g., 100 randomly selected samples).
set.seed(1)
sample_idx <- sample(1:nrow(X_train), 50)
X_sample <- X_train[sample_idx, ]
shap_global <- Shapley$new(predictor_rf, x.interest = X_sample)

# Calculate the mean absolute SHAP value for each feature.
shap_values <- shap_global$results
shap_values
avg_importance <- aggregate(abs(phi) ~ feature, data = shap_values, FUN = mean)
avg_importance <- avg_importance[order(avg_importance$`abs(phi)`, decreasing = TRUE), ]
avg_importance
# Display the results as a plot.
ggplot(avg_importance, aes(x = reorder(feature, `abs(phi)`), y = `abs(phi)`)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(x = "Features", y = "Mean(|SHAP value|)", title = "Global Feature Importance (SHAP)") +
  theme_minimal()


#########################################################
##########################################################
### SHAP Summary Plot (Beeswarm – Correct way in iml)
#########################################################

library(iml)
library(ggplot2)
library(dplyr)

# Set the number of samples used for SHAP analysis.
set.seed(123)
n_shap <- min(50, nrow(X_train))
sample_idx <- sample(1:nrow(X_train), n_shap)
X_shap <- X_train[sample_idx, , drop = FALSE]

# ===============================
# Calculate SHAP values for each sample.
# ===============================
shap_list <- lapply(1:nrow(X_shap), function(i) {

  shap_i <- Shapley$new(
    predictor_rf,
    x.interest = X_shap[i, , drop = FALSE]
  )

  res <- shap_i$results
  res$sample_id <- i
  res$feature_value <- as.numeric(X_shap[i, res$feature])

  return(res)
})

# Combine the results.
shap_results <- bind_rows(shap_list)

# ===============================
# Sort the features.
# ===============================
feature_order <- shap_results %>%
  group_by(feature) %>%
  summarise(mean_abs_shap = mean(abs(phi))) %>%
  arrange(mean_abs_shap) %>%
  pull(feature)

shap_results$feature <- factor(shap_results$feature,
                               levels = feature_order)

# ===============================
# SHAP Summary (Beeswarm) Plot
# ===============================
ggplot(shap_results,
       aes(x = phi,
           y = feature,
           color = feature_value)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_viridis_c(option = "D") +
  labs(
    x = "SHAP value (impact on model output)",
    y = "Input features",
    color = "Feature value",
    title = "SHAP Summary Plot (Random Forest – Training samples)"
  ) +
  theme_minimal(base_size = 13)


#########################################################
### SHAP Dependence Plot (Correct R implementation)
#########################################################

library(ggplot2)
library(dplyr)

# Select the main feature and the interacting feature.
feat_main <- "LST"
feat_interact <- "Rain"

# Filter SHAP results for the main feature.
shap_dep <- shap_results %>%
  filter(feature == feat_main) %>%
  mutate(interaction_value = X_train[sample_idx, feat_interact])

# Create the SHAP dependence plot.
ggplot(shap_dep,
       aes(x = feature_value,
           y = phi,
           color = interaction_value)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_viridis_c(option = "D") +
  labs(
    x = paste(feat_main, "value"),
    y = "SHAP value",
    color = paste(feat_interact, "value"),
    title = paste("SHAP Dependence Plot:", feat_main, "vs", feat_interact)
  ) +
  theme_minimal(base_size = 13)


library(openxlsx)
output_path <- "D:/DSM/output_RF_pH_results.xlsx"
wb <- createWorkbook()

# Create output worksheets.
addWorksheet(wb, "Variable_Importance")
writeData(wb, "Variable_Importance", importance_df)

addWorksheet(wb, "Training_Validation")
writeData(wb, "Training_Validation", val_train)

addWorksheet(wb, "Testing_Validation")
writeData(wb, "Testing_Validation", val_test)

pred_all <- data.frame(
  Observed = c(pH_data$pH[training], pH_data$pH[-training]),
  Predicted = c(pred_train, pred_test),
  Type = c(rep("Training", length(pred_train)), rep("Testing", length(pred_test)))
)
addWorksheet(wb, "Observed_vs_Predicted")
writeData(wb, "Observed_vs_Predicted", pred_all)

addWorksheet(wb, "SHAP_Values")
writeData(wb, "SHAP_Values", shap_values)

addWorksheet(wb, "SHAP_Global_Importance")
writeData(wb, "SHAP_Global_Importance", avg_importance)

addWorksheet(wb, "LIME_Max_pH")
writeData(wb, "LIME_Max_pH", explanation_max)

addWorksheet(wb, "LIME_Min_pH")
writeData(wb, "LIME_Min_pH", explanation_min)

addWorksheet(wb, "SHAP_Max_pH")
writeData(wb, "SHAP_Max_pH", shap_max$results)

addWorksheet(wb, "SHAP_Min_pH")
writeData(wb, "SHAP_Min_pH", shap_min$results)

saveWorkbook(wb, output_path, overwrite = TRUE)

cat("\ All results saved successfully in:", output_path, "\n")
############################################################