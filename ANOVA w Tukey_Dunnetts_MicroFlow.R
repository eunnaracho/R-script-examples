library(dplyr)
library(tidyverse)
library(DescTools)
library(writexl)
library(readxl)
library(ggplot2)

data <- read_xlsx("C:/Users/EUCHO/Desktop/4NQO DS paper/4NQO_MN_1.xlsx")

###############################################
#Data summary

#Calculate nuclei-counting bead ratio
data <- data %>%
  mutate(bead_ratio = Nuclei/Beads)

#Calculate average values within each concentration group
Avg_data <- data %>%
  group_by(Dose) %>%
  summarise_all(mean)

#Calculate standard deviation on MN frequency
MN_stdev <- data %>%
  select(Dose, 'MN%', bead_ratio)

colnames(MN_stdev)[colnames(MN_stdev)== "MN%"] <- "MN_ratio"

MN_stdev <- MN_stdev %>%
  group_by(Dose) %>%
  summarise(MN_stdev = sd(MN_ratio, na.rm=TRUE)/sqrt(n()), bead_stdev = sd(bead_ratio, na.rm=TRUE))

#Merge std dev data frame with averaged data
Avg_data <- dplyr::left_join(Avg_data, MN_stdev)

#Calculate relative survival using the n-to-b ratios
bead_ctrl <- Avg_data$Nuclei[1]/Avg_data$Beads[1]
bead_ctrl_stdev <- MN_stdev$bead_stdev[1]

Avg_data <- Avg_data %>%
  mutate(Relative_survival = (Avg_data$Nuclei/Avg_data$Beads)/bead_ctrl*100)

colnames(Avg_data)[colnames(Avg_data)== "MN%"] <- "MN_ratio"
Avg_data$Dose <- as.character(Avg_data$Dose)

relsur_stdev <- as.data.frame(Avg_data$bead_stdev/bead_ctrl_stdev)

######################################################
#Post-hoc Dunnett's Test

MN_dunnetts <- DunnettTest(x=data$`MN%`, g=data$Dose)
 
MN_dunnetts <- data.frame(MN_dunnetts$`0`)
MN_dunnetts$comparisons <- row.names(MN_dunnetts)
MN_dunnetts <- as.data.frame(MN_dunnetts) %>%
  mutate(Fold_Change = 1 + diff)

colnames(MN_dunnetts)[colnames(MN_dunnetts)== "diff"] <- "Fold change difference"
colnames(MN_dunnetts)[colnames(MN_dunnetts)== "lwr.ci"] <- "Lower CI"
colnames(MN_dunnetts)[colnames(MN_dunnetts)== "upr.ci"] <- "Upper  CI"
colnames(MN_dunnetts)[colnames(MN_dunnetts)== "pval"] <- "P-Value"

MN_dunnetts <- MN_dunnetts %>%
  mutate(label = case_when(
    MN_dunnetts$`P-Value` <0.001 ~ "***",
    MN_dunnetts$'P-Value' < 0.01 ~ "**",
    MN_dunnetts$'P-Value' < 0.05 ~ "*",
    TRUE ~ NA_character_
  ))

MN_dunnetts$comparisons[MN_dunnetts$comparisons == "0.0625-0" ] <- "0.0625"
MN_dunnetts$comparisons[MN_dunnetts$comparisons == "0.00781-0" ] <- "0.00781"
MN_dunnetts$comparisons[MN_dunnetts$comparisons == "0.0156-0" ] <- "0.0156"
MN_dunnetts$comparisons[MN_dunnetts$comparisons == "0.0313-0" ] <- "0.0313"

colnames(MN_dunnetts)[colnames(MN_dunnetts) == "comparisons"] <- "Dose"

MN_dunnetts

dun_signif <- select(MN_dunnetts, Dose, label)

write_csv(MN_dunnetts,"C:/Users/EUCHO/Desktop/4NQO DS paper/4NQO_MN_Dunnetts.csv")

######################################################
#ANOVA with post-hoc Tukey test

MN_anova <- aov(formula = data$'MN%' ~ factor(data$Dose), data = data)
summary(MN_anova)

MN_Tukey <- TukeyHSD(MN_anova)

MN_Tukey <- data.frame(MN_Tukey$`factor(data$Dose)`)
MN_Tukey <- as.data.frame(MN_Tukey) %>%
  mutate(Fold_Change = 1 + diff) %>%
  mutate(Comparison = row.names(MN_Tukey))

colnames(MN_Tukey)[colnames(MN_Tukey)== "diff"] <- "Fold change difference"
colnames(MN_Tukey)[colnames(MN_Tukey)== "lwr"] <- "Lower CI"
colnames(MN_Tukey)[colnames(MN_Tukey)== "upr"] <- "Upper  CI"
colnames(MN_Tukey)[colnames(MN_Tukey)== "p.adj"] <- "P-Value"

MN_Tukey <- MN_Tukey %>%
  mutate(label = case_when(
    MN_Tukey$'P-Value' <0.001 ~ "***",
    MN_Tukey$'P-Value' < 0.01 ~ "**",
    MN_Tukey$'P-Value' < 0.05 ~ "*",
    TRUE ~ NA_character_
  ))

MN_Tukey

write_csv(MN_Tukey,"C:/Users/EUCHO/Desktop/4NQO DS paper/4NQO_MN_Tukey.csv")

################################################
#Plot MN and relative survival with statistical significance

Avg_data <- dplyr::left_join(Avg_data, dun_signif) %>%
  mutate(fill = "% RS", colour = "% MN")

plot <- ggplot(data = Avg_data, mapping = aes(x = Dose, y = MN_ratio, group = 1, fill = fill, colour = colour)) + 
  geom_bar(stat = "identity",width = 0.5, position = position_dodge(width = 0.01), fill = "darkorange2") + 
  geom_errorbar(aes(ymin = MN_ratio - MN_stdev, ymax= MN_ratio + MN_stdev), width=.1, colour = "black") +
  labs(x = "Concentration (\U03BCg/ml)", 
       y = "Micronucleus Frequency (%)") +
  geom_path(aes(y = Relative_survival/4), colour = "dodgerblue", size = 1) +
  geom_point(aes(y = Relative_survival/4), colour = "dodgerblue", size = 2) +
  scale_y_continuous(limits = c(0, 25), sec.axis = sec_axis(~.*4, name = "Relative Survival (%)"), expand = c(0,0)) +
  theme_classic() +
  geom_text(aes(label = label), size = 4,  position=position_nudge(y = c(0, 0, 0, 2, 3)), colour = "black") +
  theme(legend.position = "top",
        axis.text.x = element_text(colour = "black", size=10),
        axis.text.y = element_text(colour = "black", size=10),
        axis.title.x = element_text(margin=margin(t=10), size=11),
        axis.title.y = element_text(margin=margin(r=10), size=11)) +
  scale_fill_discrete(name = "") +
  scale_colour_discrete(name = "")

plot

setwd("C:/Users/EUCHO/Desktop/4NQO DS paper")
ggsave(plot, filename = "MN_RS_4NQO_plot2.jpeg",  width = 15, height = 10, units = "cm")


                          
                          