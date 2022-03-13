# cross-biome soil respiration analysis
library(ggplot2)
library(RColorBrewer)
library(MASS)
library(deming)
library(sfsmisc)
library(dplyr)
library(ggsci)
library("scales")
library(wesanderson)
library(ggh4x)#set scales for facet_wrap

safe_colorblind_palette <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
                             "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
                             "#920000","#924900","#db6d00","#24ff24","#ffff6d")


CW = read.csv("C_W_1994_2020_updated.csv", stringsAsFactors = T)
table(CW$EcosystemType2)

CW <- filter(CW, !EcosystemType %in% c("WetSedgeTundra")) %>% droplevels() %>%
  select("RecordID", "DOY", "EcosystemType2", "Percent_change" ,"log_change" ,
         "change_in_T","cohensd","Stemp_Avg_C", "n_Replicates.1", "Source" ) %>% na.omit()
table(CW$Source)

CW$EcosystemType2 = factor(CW$EcosystemType2, levels = c("BorealForest", "TemperateGrassland", 
                                                       "TemperateForest", "TemperateAgriculture",
                                                       "Meadow","Shrubland","Desert","TropicalForest")) 

CW.clean <- lapply(split(CW, CW$EcosystemType2), function(dataframe){
  mod <- lm(dataframe[, 'Percent_change'] ~  DOY+ I(DOY^2), data = dataframe)
  cooksd <- cooks.distance(mod)
  influential <- as.numeric(names(cooksd)[(cooksd >  4 * mean(cooksd, na.rm = TRUE))])
  final <- filter(dataframe, !rownames(dataframe) %in% influential)
  final}) %>% do.call(rbind, .)

CW.clean$EcosystemType = CW$EcosystemType2[match(CW.clean$RecordID,CW$RecordID)]
table(CW.clean$EcosystemType2)

# with DOY ----
fit = lm(Percent_change ~ DOY + I(DOY^2) + change_in_T, CW.clean)
summary(fit)

##
fit = lm(Percent_change ~ DOY+ I(DOY^2)+change_in_T, weights = n_Replicates.1, filter( CW.clean, EcosystemType == "BorealForest"))
summary(fit)
plot(fit)
fit = lm(Percent_change ~ DOY+ I(DOY^2), filter( CW.clean, EcosystemType == "BorealForest"))
summary(fit)
plot(fit)
##
fit = lm(Percent_change ~ DOY+ I(DOY^2)+change_in_T, weights = n_Replicates.1, filter( CW.clean, EcosystemType == "TemperateForest"))
summary(fit)
plot(fit)
fit = lm(Percent_change ~ DOY+ I(DOY^2), weights = n_Replicates.1, filter( CW.clean, EcosystemType == "TemperateForest"))
summary(fit)
fit = lm(Percent_change ~ DOY+ I(DOY^2) +change_in_T, filter( CW.clean, EcosystemType == "TemperateForest"))
summary(fit)
plot(fit)
##
fit = lm(Percent_change ~ DOY+ I(DOY^2)+change_in_T, weights = n_Replicates.1, filter( CW.clean, EcosystemType == "TemperateGrassland"))
summary(fit)
plot(fit)
fit = lm(Percent_change ~ DOY+ I(DOY^2), weights = n_Replicates.1, filter( CW.clean, EcosystemType == "TemperateGrassland"))
summary(fit)
plot(fit)
##
fit = lm(Percent_change ~ DOY+ I(DOY^2)+change_in_T, weights = n_Replicates.1, filter( CW.clean, EcosystemType == "TemperateAgriculture"))
summary(fit)
plot(fit)
fit = lm(Percent_change ~ DOY+ I(DOY^2), weights = n_Replicates.1, filter( CW.clean, EcosystemType == "TemperateAgriculture"))
summary(fit)
fit = lm(Percent_change ~ DOY+ I(DOY^2)+change_in_T, filter( CW.clean, EcosystemType == "TemperateAgriculture"))
summary(fit)
plot(fit)
##
fit = lm(Percent_change ~ DOY+ I(DOY^2)+change_in_T, weights = n_Replicates.1, filter( CW.clean, EcosystemType == "Meadow"))
summary(fit)
plot(fit)
fit = lm(Percent_change ~ DOY+ I(DOY^2), weights = n_Replicates.1, filter( CW.clean, EcosystemType == "Meadow"))
summary(fit)
fit = lm(Percent_change ~ DOY+ I(DOY^2)+change_in_T, filter( CW.clean, EcosystemType == "Meadow"))
summary(fit)
plot(fit)
##
fit = lm(Percent_change ~ DOY+ I(DOY^2)+change_in_T, weights = n_Replicates.1, filter( CW.clean, EcosystemType == "Shrubland"))
summary(fit)
plot(fit)
fit = lm(Percent_change ~ DOY+ I(DOY^2), weights = n_Replicates.1, filter( CW.clean, EcosystemType == "Shrubland"))
summary(fit)
fit = lm(Percent_change ~ DOY+ I(DOY^2)+change_in_T, filter( CW.clean, EcosystemType == "Shrubland"))
summary(fit)
plot(fit)
##
fit = lm(Percent_change ~ DOY+ I(DOY^2)+change_in_T, weights = n_Replicates.1, filter( CW.clean, EcosystemType == "Desert"))
summary(fit)
plot(fit)
fit = lm(Percent_change ~ DOY+ I(DOY^2), weights = n_Replicates.1, filter( CW.clean, EcosystemType == "Desert"))
summary(fit)
fit = lm(Percent_change ~ DOY+ I(DOY^2)+change_in_T, filter( CW.clean, EcosystemType == "Desert"))
summary(fit)
plot(fit)
##
fit = lm(Percent_change ~ DOY+ I(DOY^2)+change_in_T, weights = n_Replicates.1, filter( CW.clean, EcosystemType == "TropicalForest"))
summary(fit)
plot(fit)
fit = lm(Percent_change ~ DOY+ I(DOY^2), weights = n_Replicates.1, filter( CW.clean, EcosystemType == "TropicalForest"))
summary(fit)
fit = lm(Percent_change ~ DOY+ I(DOY^2)+change_in_T, filter( CW.clean, EcosystemType == "TropicalForest"))
summary(fit)
plot(fit)

ggplot(CW.clean) +
  aes(x = DOY, y = Percent_change, 
      size = n_Replicates.1, weight = n_Replicates.1) +
  scale_size_continuous(name="Number of replicates",  range = c(0.5,4)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = "gray50", size = 0.5)+
  geom_point(shape = "circle open", alpha = 0.7, color = "gray50") +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se = T,aes(weight = n_Replicates.1, color = EcosystemType2,fill = EcosystemType2)) +
  scale_y_continuous(limits = c(-100, 150), breaks = c(-100, -50, 0, 50, 100, 150)) +
  labs(x = "DOY", y = "Percent change in Rs (%)") +
  ggthemes::theme_few() +
  facet_wrap(vars(EcosystemType2), ncol = 4) +
  scale_fill_manual(values = safe_colorblind_palette[c(3,4,7,6,8,10,13,15)]) +
  scale_color_manual(values = safe_colorblind_palette[c(3,4,7,6,8,10,13,15)]) +
  theme(strip.text = element_text(size = 14),
        legend.spacing.y = unit(0.5, 'cm'),
        legend.title = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 14),
        axis.text = element_text(color = "black", size = 14),
        axis.title.y = element_text(color = "black", size = 14, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        #legend.position = "none",
        axis.title.x = element_text(color = "black", size = 14, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1))


ggsave("PC_DOY_AllEco_bw_withtropics.pdf", width = 13.5, height = 6)

# with temperature ----
library(wCorr)

weights::wtd.cor(filter(CW.clean,EcosystemType == "BorealForest")$Percent_change,
        filter(CW.clean,EcosystemType == "BorealForest")$Stemp_Avg_C, 
        weight=filter(CW.clean,EcosystemType == "BorealForest")$n_Replicates.1)
cor(filter(CW.clean,EcosystemType == "BorealForest")$Percent_change, filter(CW.clean,EcosystemType == "BorealForest")$Stemp_Avg_C)

weights::wtd.cor(filter(CW.clean,EcosystemType == "TemperateGrassland")$Percent_change, 
        filter(CW.clean,EcosystemType == "TemperateGrassland")$Stemp_Avg_C, 
        weight=filter(CW.clean,EcosystemType == "TemperateGrassland")$n_Replicates.1)
cor(filter(CW.clean,EcosystemType == "TemperateGrassland")$Percent_change, filter(CW.clean,EcosystemType == "TemperateGrassland")$Stemp_Avg_C)

weights::wtd.cor(filter(CW.clean,EcosystemType == "TemperateForest")$Percent_change, 
        filter(CW.clean,EcosystemType == "TemperateForest")$Stemp_Avg_C, 
        weight=filter(CW.clean,EcosystemType == "TemperateForest")$n_Replicates.1)
cor(filter(CW.clean,EcosystemType == "TemperateForest")$Percent_change, filter(CW.clean,EcosystemType == "TemperateForest")$Stemp_Avg_C)

weights::wtd.cor(filter(CW.clean,EcosystemType == "TemperateAgriculture")$Percent_change, 
        filter(CW.clean,EcosystemType == "TemperateAgriculture")$Stemp_Avg_C, 
        weight=filter(CW.clean,EcosystemType == "TemperateAgriculture")$n_Replicates.1)
cor(filter(CW.clean,EcosystemType == "TemperateAgriculture")$Percent_change, filter(CW.clean,EcosystemType == "TemperateAgriculture")$Stemp_Avg_C)

weights::wtd.cor(filter(CW.clean,EcosystemType == "Meadow")$Percent_change, 
        filter(CW.clean,EcosystemType == "Meadow")$Stemp_Avg_C, 
        weight=filter(CW.clean,EcosystemType == "Meadow")$n_Replicates.1)
cor(filter(CW.clean,EcosystemType == "Meadow")$Percent_change, filter(CW.clean,EcosystemType == "Meadow")$Stemp_Avg_C)

weights::wtd.cor(filter(CW.clean,EcosystemType2 == "Shrubland")$Percent_change, 
        filter(CW.clean,EcosystemType2 == "Shrubland")$Stemp_Avg_C, 
        weight=filter(CW.clean,EcosystemType2 == "Shrubland")$n_Replicates.1)
cor(filter(CW.clean,EcosystemType == "N_Shrubland")$Percent_change, filter(CW.clean2,EcosystemType == "N_Shrubland")$Stemp_Avg_C)

weights::wtd.cor(filter(CW.clean,EcosystemType == "S_Shrubland")$Percent_change, 
        filter(CW.clean,EcosystemType == "S_Shrubland")$Stemp_Avg_C,
        weight=filter(CW.clean,EcosystemType == "S_Shrubland")$n_Replicates.1)
cor(filter(CW.clean,EcosystemType == "S_Shrubland")$Percent_change, filter(CW.clean,EcosystemType == "S_Shrubland")$Stemp_Avg_C)

weights::wtd.cor(filter(CW.clean,EcosystemType == "Desert")$Percent_change,
        filter(CW.clean,EcosystemType == "Desert")$Stemp_Avg_C, 
        weight=filter(CW.clean,EcosystemType == "Desert")$n_Replicates.1)
cor(filter(CW.clean,EcosystemType == "Desert")$Percent_change, filter(CW.clean,EcosystemType == "Desert")$Stemp_Avg_C)

weights::wtd.cor(filter(CW.clean,EcosystemType == "TropicalForest")$Percent_change,
                 filter(CW.clean,EcosystemType == "TropicalForest")$Stemp_Avg_C, 
                 weight=filter(CW.clean,EcosystemType == "TropicalForest")$n_Replicates.1)
cor(filter(CW.clean,EcosystemType == "TropicalForest")$Percent_change, filter(CW.clean,EcosystemType == "TropicalForest")$Stemp_Avg_C)

fit = lm(Percent_change ~ Stemp_Avg_C + change_in_T, weights = n_Replicates.1, filter( CW.clean, EcosystemType == "BorealForest"))
summary(fit)
fit = lm(Percent_change ~ Stemp_Avg_C + change_in_T, weights = n_Replicates.1, filter( CW.clean, EcosystemType2 == "TemperateGrassland"))
summary(fit)
fit = lm(Percent_change ~ Stemp_Avg_C + change_in_T, weights = n_Replicates.1, filter( CW.clean, EcosystemType2 == "TemperateForest"))
summary(fit)
fit = lm(Percent_change ~ Stemp_Avg_C + change_in_T, weights = n_Replicates.1, filter( CW.clean, EcosystemType2 == "TemperateAgriculture"))
summary(fit)
fit = lm(Percent_change ~ Stemp_Avg_C + change_in_T, weights = n_Replicates.1, filter( CW.clean, EcosystemType2 == "Meadow"))
summary(fit)
fit = lm(Percent_change ~ Stemp_Avg_C + change_in_T, weights = n_Replicates.1, filter( CW.clean, EcosystemType2 == "Shrubland"))
summary(fit)
fit = lm(Percent_change ~ Stemp_Avg_C + change_in_T, weights = n_Replicates.1, filter( CW.clean, EcosystemType2 == "Desert"))
summary(fit)
fit = lm(Percent_change ~ Stemp_Avg_C + change_in_T, weights = n_Replicates.1, filter( CW.clean, EcosystemType2 == "TropicalForest"))
summary(fit)

ggplot(CW.clean) +
  aes(x = Stemp_Avg_C, y = Percent_change,
      size = n_Replicates.1, weight = n_Replicates.1) +
  scale_size_continuous(name="Number of replicates",  range = c(0.5,4)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = "gray50", size = 0.5)+
  geom_point(shape = "circle open", alpha = 0.7,color = "grey50") + # grey70color = "grey40"
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = T, aes(weight = n_Replicates.1, color = EcosystemType2,fill = EcosystemType2)) + #color = "black"
  # geom_ribbon(aes(weight = n_Replicates.1, color = EcosystemType),stat = "smooth",
  #             method = "lm",
  #             formula = y ~ x,
  #             se = TRUE,
  #             alpha = 0, # or, use fill = NA
  #             size = 0.25,
  #             linetype = "solid") +
  scale_y_continuous(limits = c(-100, 140), breaks = c(-100, -50, 0, 50, 100, 150)) +
  labs(x = expression(Soil ~ temperature ~ (degree * C)), y = "Percent change in Rs (%)") +
  ggthemes::theme_few() +
  facet_wrap(vars(EcosystemType2), ncol = 4, scales = "free_x") +
  scale_fill_manual(values = safe_colorblind_palette[c(3,4,7,6,8,10,13,15)]) +
  scale_color_manual(values = safe_colorblind_palette[c(3,4,7,6,8,10,13,15)]) +
  theme(strip.text = element_text(size = 14),
        legend.spacing.y = unit(0.5, 'cm'),
        legend.title = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 14),
        axis.text = element_text(color = "black", size = 14),
        axis.title.y = element_text(color = "black", size = 14, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        #legend.position = "none",
        axis.title.x = element_text(color = "black", size = 14, margin = margin(t = 6, r = 0, b = 0, l = 0)),
        panel.border = element_rect(size = 1)) +
  facetted_pos_scales(
    x = list(
      EcosystemType2 == "TemperateForest" ~ scale_x_continuous(limits = c(0,35), breaks = c(0, 10, 20,30)),
      EcosystemType2 == "TropicalForest" ~ scale_x_continuous(limits = c(25,27), breaks = c(25,26,27)),
      EcosystemType2 == "Meadow" ~ scale_x_continuous(limits = c(0,22), breaks = c(0,5,10,15,20))
    )
  )

ggsave("PC_T_AllEco_bw_withtropics.pdf", width = 13, height = 6)
