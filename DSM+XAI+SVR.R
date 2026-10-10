######## Support Vector Regression (SVR) ######################
library(terra)
library(caret)
library(e1071)
library(ggplot2)
library(iml)
library(lime)
library(openxlsx)

###### Data
setwd("D:/DSM/tiff")

##### Stack
##### Stack (same as before)
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

########### DEPENDENT Variable ###############
pH_data <- read.table("soil.txt", header = TRUE, sep = "\t",
                      na.strings = "NA", dec = ".", strip.white = TRUE)
str(pH_data)

###################### Train / Test
#set.seed(123)
training <- sample(nrow(pH_data), 0.8 * nrow(pH_data))

###################### SVR Model
ctrl <- trainControl(method = "boot", number = 10)

svr_model <- train(
  pH ~ Aspect + Dem + LST + MRVBF + MRRTF + MSI + 
    Profile_crv + River_dis + Road_dis + Sand + Silt +
    Slope + Rain + Tmin,
  data = pH_data,
  method = "svmRadial",
  trControl = ctrl,
  preProcess = c("center", "scale"),
  tuneLength = 10
)

print(svr_model)

###################### Variable Importance (SVR)
importance_df <- varImp(svr_model, scale = TRUE)$importance
importance_df$Variable <- rownames(importance_df)
importance_df <- importance_df[order(importance_df$Overall, decreasing = TRUE), ]

ggplot(importance_df,
       aes(x = reorder(Variable, Overall), y = Overall)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Variable Importance - SVR", x = "", y = "Importance")

##### Validation Function
goof2 <- function(observed, predicted, title = "Model Validation") {

  observed <- as.numeric(observed)
  predicted <- as.numeric(predicted)
  valid <- complete.cases(observed, predicted)
  observed <- observed[valid]
  predicted <- predicted[valid]

  RMSE <- sqrt(mean((observed - predicted)^2))
  MAE  <- mean(abs(observed - predicted))
  R2   <- cor(observed, predicted)^2
  BIAS <- mean(predicted - observed)
  NSE  <- 1 - (sum((observed - predicted)^2) /
                 sum((observed - mean(observed))^2))
  RPD  <- sd(observed) / RMSE

  results <- data.frame(RMSE, MAE, R2, NSE, BIAS, RPD)
  print(results)

  plot(observed, predicted, pch = 19, col = "steelblue",
       xlab = "Observed pH", ylab = "Predicted pH", main = title)
  abline(0, 1, col = "red", lwd = 2)

  return(results)
}

# ===================== Validation =====================
pred_train <- predict(svr_model, newdata = pH_data[training, ])
cat(" Internal Validation (Training Data)\n")
val_train <- goof2(pH_data$pH[training], pred_train,
                   "SVR - Training")

pred_test <- predict(svr_model, newdata = pH_data[-training, ])
cat(" External Validation (Testing Data)\n")
val_test <- goof2(pH_data$pH[-training], pred_test,
                  "SVR - Testing")

###################### LIME
model_type.train <- function(x, ...) "regression"
predict_model.train <- function(x, newdata, ...)
  data.frame(Response = predict(x, newdata))

train_x <- pH_data[training, -which(names(pH_data) == "pH")]

explainer <- lime(train_x, svr_model)

case_max <- pH_data[which.max(pH_data$pH),
                    !(names(pH_data) %in% c("EC","pH","x","y","SPI"))]
case_min <- pH_data[which.min(pH_data$pH),
                    !(names(pH_data) %in% c("EC","pH","x","y","SPI"))]

explanation_max <- explain(case_max, explainer,
                           n_labels = 1, n_features = 15)
explanation_min <- explain(case_min, explainer,
                           n_labels = 1, n_features = 15)

par(mfrow = c(1,2))
plot_features(explanation_max)
plot_features(explanation_min)
par(mfrow = c(1,1))

###################### Prediction Map
svr_predict_fun <- function(model, data, ...) {
  pred <- predict(model, newdata = data)
  return(as.numeric(pred))
}
pH_map <- predict(
  covStack,
  svr_model,
  fun = svr_predict_fun,
  filename = "pH_SVR_Map.tif",
  datatype = "FLT4S",
  overwrite = TRUE,
  na.rm = TRUE
)

plot(pH_map, main = "Predicted Soil pH (SVR)")

###################### SHAP (iml)
X_train <- pH_data[training, !(names(pH_data) %in% c("pH","EC","x","y","SPI"))]

predictor_svr <- Predictor$new(
  model = svr_model,
  data = X_train,
  y = pH_data$pH[training]
)

case_max <- X_train[which.max(pH_data$pH[training]), , drop = FALSE]
case_min <- X_train[which.min(pH_data$pH[training]), , drop = FALSE]

shap_max <- Shapley$new(predictor_svr, x.interest = case_max)
shap_min <- Shapley$new(predictor_svr, x.interest = case_min)

plot(shap_max)
plot(shap_min)

###################### Global SHAP
set.seed(1)
sample_idx <- sample(1:nrow(X_train), 50)
X_sample <- X_train[sample_idx, ]

shap_global <- Shapley$new(predictor_svr, x.interest = X_sample)
shap_values <- shap_global$results

avg_importance <- aggregate(abs(phi) ~ feature,
                            data = shap_values, mean)
avg_importance <- avg_importance[
  order(avg_importance$`abs(phi)`, decreasing = TRUE), ]

ggplot(avg_importance,
       aes(x = reorder(feature, `abs(phi)`),
           y = `abs(phi)`)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Global Feature Importance (SHAP - SVR)",
       x = "Features", y = "Mean(|SHAP|)")
#########################################################
### SHAP Summary Plot (Beeswarm – Correct way in iml)
#########################################################

library(iml)
library(ggplot2)
library(dplyr)

# Number of samples for SHAP analysis
set.seed(123)
n_shap <- min(50, nrow(X_train))
sample_idx <- sample(1:nrow(X_train), n_shap)
X_shap <- X_train[sample_idx, , drop = FALSE]

# ===============================
# Calculate SHAP values for each sample
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

# Combine the results
shap_results <- bind_rows(shap_list)

# ===============================
# Sort the features
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
###################################################3
#########################################################
### SHAP Summary Plot (Beeswarm – Publication ready)
#########################################################

library(iml)
library(ggplot2)
library(dplyr)

# ===============================
# Plot
# ===============================
library(iml)
library(ggplot2)
library(dplyr)

p_shap <- ggplot(
  shap_results,
  aes(
    x = phi,
    y = feature,
    color = feature_value
  )
) +
  geom_point(alpha = 0.7, size = 2) +

  scale_color_viridis_c(
    option = "D",
    name = "Feature value",
    breaks = c(
      min(shap_results$feature_value, na.rm = TRUE),
      max(shap_results$feature_value, na.rm = TRUE)
    ),
    labels = c("Low", "High")
  ) +

  labs(
    x = "SHAP value (impact on model output)",
    y = "Input features",
    title = "SHAP Summary Plot (SVR)"
  ) +

  theme_minimal(base_size = 14) +

  theme(
    # ===== Set all text to black =====
    text = element_text(color = "black"),

    axis.title = element_text(color = "black"),
    axis.text  = element_text(color = "black"),

    plot.title = element_text(
      color = "black",
      face = "bold",
      hjust = 0.5
    ),

    legend.title = element_text(color = "black"),
    legend.text  = element_text(color = "black"),

    legend.position = "right",
    legend.direction = "vertical",
    legend.key.height = unit(1, "npc"),
    legend.key.width  = unit(0.4, "cm")
  ) +

  guides(
    color = guide_colorbar(
      barheight = unit(0.8, "npc"),
      barwidth  = unit(0.6, "cm"),
      ticks = FALSE
    )
  )

# Display the plot
print(p_shap)

# Display the map
print(p_shap)
ggsave(
  filename = "D:/DSM/SHAP_Summary_SVR_pH.png",
  plot = p_shap,
  width = 10,
  height = 10,
  dpi = 1000
)

###################### Save Outputs
library(openxlsx)

output_path <- "D:/DSM/output_SVR_pH_results.xlsx"
wb <- createWorkbook()

# =====================================================
# Variable Importance
# =====================================================
addWorksheet(wb, "Variable_Importance")
writeData(wb, "Variable_Importance", importance_df)

# =====================================================
# Training Validation
# =====================================================
addWorksheet(wb, "Training_Validation")
writeData(wb, "Training_Validation", val_train)

# =====================================================
# Testing Validation
# =====================================================
addWorksheet(wb, "Testing_Validation")
writeData(wb, "Testing_Validation", val_test)

# =====================================================
# Observed vs Predicted (Training + Testing)
# =====================================================
pred_all <- data.frame(
  Observed  = c(pH_data$pH[training], pH_data$pH[-training]),
  Predicted = c(pred_train, pred_test),
  Type = c(
    rep("Training", length(pred_train)),
    rep("Testing", length(pred_test))
  )
)

addWorksheet(wb, "Observed_vs_Predicted")
writeData(wb, "Observed_vs_Predicted", pred_all)

# =====================================================
# SHAP Values (Local – all calculated)
# =====================================================
addWorksheet(wb, "SHAP_Values")
writeData(wb, "SHAP_Values", shap_values)

# =====================================================
# SHAP Global Importance
# =====================================================
addWorksheet(wb, "SHAP_Global_Importance")
writeData(wb, "SHAP_Global_Importance", avg_importance)

# =====================================================
# LIME – Max pH
# =====================================================
addWorksheet(wb, "LIME_Max_pH")
writeData(wb, "LIME_Max_pH", explanation_max)

# =====================================================
# LIME – Min pH
# =====================================================
addWorksheet(wb, "LIME_Min_pH")
writeData(wb, "LIME_Min_pH", explanation_min)

# =====================================================
# SHAP – Max pH Case
# =====================================================
addWorksheet(wb, "SHAP_Max_pH")
writeData(wb, "SHAP_Max_pH", shap_max$results)

# =====================================================
# SHAP – Min pH Case
# =====================================================
addWorksheet(wb, "SHAP_Min_pH")
writeData(wb, "SHAP_Min_pH", shap_min$results)

# =====================================================
# Save Workbook
# =====================================================
saveWorkbook(wb, output_path, overwrite = TRUE)

cat("All SVR results saved successfully in:", output_path, "\n")

