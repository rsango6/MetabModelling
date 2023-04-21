
#### Structure of GEMs, reporterMetabolites, reporter Subnetworks, FBA results, FVA results ... ####

library(tidyverse)
library(readxl)
library(xlsx)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(ggvenn)
library(gridExtra)
library(ggpubr)
library(ggallin)
library(stringr)
library(magrittr)
library(pheatmap)
require(grid)

file.choose()

ControlRXN = read.csv('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/OnlyReactionsFromModels/ControlModel.csv',
                      header = T)
# LPS

LPSRXN = read.csv('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/OnlyReactionsFromModels/LPSModel.csv',
                  header = T)
# IL4

IL4RXN = read.csv('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/OnlyReactionsFromModels/IL4Model.csv',
                  header = T)
reactions = read.csv('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/reactions.csv')


# ResReactionsMatrix = read.csv('ResReactionsMatrix.csv')
# SubsystemCoverage = read.csv('SubsystemCoverage.csv')
# MetabolicTasks = read.csv('res_func.funcComp.matrix.csv')
# 
# ResReactionsMatrixNum = ResReactionsMatrix %>%
#   select(where(is.numeric))

######################### PCoA for Reaction Matrix ######################

# distance.matrix <- dist(scale(t(ResReactionsMatrixNum), center = TRUE, scale = TRUE),
#                         method = "euclidean")
# mds <- cmdscale(distance.matrix, eig = T, x.ret = T)
# mds.var.per <- round(mds$eig/sum(mds$eig)*100, 1)
# mds.values <- mds$points
# mds.data <- data.frame(Sample = rownames(mds.values),
#                        X = mds.values[, 1],
#                        Y = mds.values[, 2])
# ggplot(mds.data, aes(x = X, y = Y, color = Sample)) +
#   geom_point() + theme_bw() +
#   xlab(paste("MDS1 -", mds.var.per[1], "%", sep = "")) +
#   ylab(paste("MDS2 -", mds.var.per[2], "%", sep = "")) +
#   ggtitle("PCoA tINIT GEMs")
# 
# 

########################### Subsystem Coverage ########################

# colors <- colorRampPalette(brewer.pal(8, "YlOrBr"))(7)
# 
# SubsystemCoverageNum = SubsystemCoverage %>%
#   column_to_rownames("SubsystemsIDs") %>%
#   select(where(is.numeric))
# 
# Filtered <- filter_all(SubsystemCoverageNum, any_vars(. > 25)) #for some reason this dplyr function not working
# Filtered = SubsystemCoverageNum[rowSums(abs(SubsystemCoverageNum > 15)) >= 1, ]
# 
# 
# pheatmap(Filtered, 
#          #clustering_distance_rows = sampleDists,
#          #clustering_distance_cols = sampleDists,
#          col = colors,
#          main = "SubsystemCoverage for tINIT models (15% distance from the rowMean)",
#          labels_col = c('LPSGEM', "IL4GEM", "ControlGEM"),
#          angle_col = 45)


#### Performance on Full Metabolic Tasks ( list of 257) #####
# 
# 
# 
# MetabolicTasksNum = MetabolicTasks %>%
#   column_to_rownames("MetabolicTasksIDs") %>%
#   select(where(is.numeric))
#   
# testfilter = MetabolicTasksNum[rowSums(MetabolicTasksNum) == 3 | rowSums(MetabolicTasksNum) == 0, ]
# 
# MetabolicTasksNum = MetabolicTasksNum[!rownames(MetabolicTasksNum) %in% rownames(testfilter),]
# 

#convert wide to long format
#melt <- melt(as.matrix(MetabolicTasksNum))
#plot
# ggplot(melt, aes(Var2, Var1)) +
#   geom_point(aes(colour = c("Fail", "Pass")[value + 1 ])) +
#   scale_fill_identity() + theme_minimal() + theme(axis.text.x = element_text(angle = 45),
#                                                   legend.text = element_text(face = 'bold'),
#                                                   plot.subtitle = element_text(face = "italic", size = 10)) +
#   labs(title = "Metabolic Tasks Performance", subtitle = "Those that differ among conditions", 
#        x = "", y = "", color = "Status")

#### Getting proper media compositions for model constraining ####
setwd("/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/GEM_Comparison_R")
options(scipen=999)
MediaCompDf = read.csv('MediaCompDf.csv', header = T, row.names = 1)

MediaCompDf = MediaCompDf %>%
  mutate(MetPercentage = (Concentration..mg.L. / sum(Concentration..mg.L.)) / 10 ) %>% #changing mg/L to g/L and then taking percentage in one step
  mutate(MetWeight = MetPercentage * 3) %>%
  mutate(FluxConstraint_mmol_per_mouse_day = (MetWeight / Molecular.Weight) * 1000)

write.csv(MediaCompDf, 'MediaCompDf.csv')


#Venn diagram of Constrained GEMs 
setwd("~/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models")


LPS = read.csv('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/LPS_GEM_genes.csv',
               header = T)
IL4 = read.csv('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/IL4_GEM_genes.csv',
               header = T)
Control = read.csv('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/Control_GEM_genes.csv',
                   header = T)

x <- list(
  LPS = as.character(LPS$NAME), 
  IL4 = as.character(IL4$NAME), 
  Control = as.character(Control$NAME)
)

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4) + ggtitle("Genes in GEMs")


#### ReporterMetabolites analysis ####

par(mfrow=c(3,2))

generate_RepMets = function(RepMets, title) {
  
  RepMets = RepMets %>% 
    filter(metPValues < 0.05) %>%
    arrange(desc(metZScore)) %>%
    mutate(Compartment = case_when(
      endsWith(mets, "c") ~ "Cytosol",
      endsWith(mets, "p") ~ "Peroxisome",
      endsWith(mets, "m") ~ "Mitochondria",
      endsWith(mets, "s") ~ "Extracellular",
      endsWith(mets, "l") ~ "Lysosome",
      endsWith(mets, "r") ~ "Endoplasmic_reticulum",
      endsWith(mets, "g") ~ "Golgi_apparatus",
      endsWith(mets, "n") ~ "Nucleus",
      endsWith(mets, "i") ~ "Inner_mitochondria"
    )) %>%
    relocate(Compartment, .before = metZScore) %>%
    subset(!(metNames %in% c('H+','NADP+', 'NADPH', 'NAD+', 'NADH', 'H2O', 'Pi',
                             'O2', 'Na+', 'Fe3+', 'Fe2+', 'Cu2+', 'Ca2+'))) %>%
    select(-mets)
  
  RepMets$Compartment = as.factor(RepMets$Compartment)
  
  myColors <- colorRampPalette(brewer.pal(8, "Dark2"))(9)
  names(myColors) <- levels(RepMets$Compartment)
  custom_colors <- scale_colour_manual(name = "Compartment", values = myColors)
  
  g = ggplot(RepMets[1:20,], aes(x=reorder(metNames, metZScore), y=metZScore)) +
    geom_bar(aes(fill = Compartment), stat="identity") + coord_flip() + theme_light() +
    custom_colors +
    theme(axis.text = element_text(face="bold", size = 12),
          axis.title = element_text(face="bold", size = 14),
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 14),
          plot.title = element_text(size = 20)) +
    labs(x = NULL, title = title)
  return(g)
}

#file.choose()

setwd('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/Plots/RepMets/')
### IL4 ###

# general
IL4 = read.csv("/Users/rokosango/PhD/MetabModelling/RS/RS_IL4_full/RepMetsIL4.csv",
               header = T)

Compartment = str_sub(IL4$mets, -1, -1) # extract compartment letter
Mets = str_sub(IL4$mets, end=-2) # extract Metabolite name

IL4$ModelMetsNames = paste0(IL4$metNames,"[", Compartment, "]")

IL4 = IL4 %>%
  select(ModelMetsNames, metPValues) %>%
  rename(Name = ModelMetsNames) %>%
  rename("p-value" = metPValues)

write.xlsx2(IL4, "GSS_IL4.xlsx")

a = generate_RepMets(IL4, "Top 20 RepMets IL4")

a
# Up 
IL4 = read.csv("/Users/rokosango/PhD/MetabModelling/RS/RS_IL4_full/RepMetsIL4Up.csv",
               header = T)

b = generate_RepMets(IL4, "Top 20 RepMets IL4 Up")



# Down
IL4 = read.csv("/Users/rokosango/PhD/MetabModelling/RS/RS_IL4_full/RepMetsIL4Down.csv",
               header = T)
c = generate_RepMets(IL4, "Top 20 RepMets IL4 Down") 

c

### LPS ###

# general
LPS = read.csv("/Users/rokosango/PhD/MetabModelling/RS/RS_LPS_full/RepMetsLPS.csv",
               header = T)

Compartment = str_sub(LPS$mets, -1, -1) # extract compartment letter 
Mets = str_sub(LPS$mets, end=-2) # extract Metabolite name, I need both for Kiwi package

LPS$ModelMetsNames = paste0(LPS$metNames,"[", Compartment, "]")

LPS = LPS %>%
  select(ModelMetsNames, metPValues) %>%
  rename(Name = ModelMetsNames) %>%
  rename("p-value" = metPValues) # need quotation marks when there is a hyphen in the name

write.xlsx2(LPS, "GSS_LPS.xlsx")


d = generate_RepMets(LPS, "Top 20 RepMets LPS")



# Up
LPSUp = read.csv("/Users/rokosango/PhD/MetabModelling/RS/RS_LPS_full/RepMetsLPSUp.csv",
               header = T)
e = generate_RepMets(LPS, "Top 20 RepMets LPS Up")

e

# Down
LPS = read.csv("/Users/rokosango/PhD/MetabModelling/RS/RS_LPS_full/RepMetsLPSDown.csv",
               header = T)
f = generate_RepMets(LPS, "Top 20 RepMets LPS Down")

f

# combine plots

ggarrange(b, c, e, f,
             common.legend = T, legend = "bottom")

ggarrange(b,c, common.legend = T, legend = 'top')

ggarrange(e, f, common.legend = T, legend = 'top')


#### Exchange reactions ####
options(scipen = 1)

# Control

ControlRXN = read.csv('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/OnlyReactionsFromModels/ControlModel.csv',
                   header = T)
# LPS

LPSRXN = read.csv('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/OnlyReactionsFromSimplifiedModels/Simplified_LPS_rxns_only.csv',
                   header = T)
# IL4

IL4RXN = read.csv('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/OnlyReactionsFromSimplifiedModels/Simplified_IL4_rxns_only.csv',
               header = T)


# generate_barplot = function(x, titlep, titlec) {
#   x$X. = NULL
#   
#   x = x %>%
#     mutate(RXNNAMES =gsub("\\[s]\\ <=> ", "", x$EQUATION)) %>%
#     relocate(RXNNAMES, .after = ID) %>%
#     filter(SUBSYSTEM == "Exchange/demand reactions") %>%
#     filter(!RXNNAMES %in%  c('HCO3-', 'Pi', 'NADP+', 'NAD+', 
#                              'NH4+', 'Fe2+', 'O2-', 'H+')) %>%
#     arrange(desc(Flux))
#   #filter(Flux < -1)
#   #arrange(desc(Flux)) #%>%
#   #mutate(Activity = cut(Flux, breaks = c(0, -100, -500, -1001),
#   #labels = c("High", "Medium", "Low"))) %>%
#   #relocate(Activity, .after = RXNNAMES)
#   
#   p = ggplot(x[x$Flux > 1e-6,], aes(x=reorder(RXNNAMES, -Flux), y=Flux)) +
#     geom_bar(stat="identity") + coord_flip() + theme_light() +
#     labs(x = NULL, title = titlep)
#   
#   c = ggplot(x[x$Flux < -1e-6,], aes(x=reorder(RXNNAMES, Flux), y=Flux)) +
#     geom_bar(stat="identity") + coord_flip() + theme_light() +
#     labs(x = NULL, title = titlec)
#   
#   
#   g = ggarrange(p, c, 
#                 #labels = c("Low", "Medium", "High"), vjust = c(1.5,1.5,0.5),
#                 ncol = 2, nrow = 1)
#   return(g)
# }
# 
# generate_barplot(Control, 'Control GEM pFBA Production', 'Control GEM pFBA Consumption')
# generate_barplot(LPS, 'LPS GEM pFBA Production', 'LPS GEM pFBA Consumption')
# generate_barplot(IL4, 'IL4 GEM pFBA Production', 'IL4 GEM pFBA Consumption')

#### Alternative Bounds ####

file.choose(
)

AltB = read.csv('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/AlternativeBounds.csv',
                header = T)



ggplot(data=AltB, aes(x=reorder(Bounds, -Obj_Value), y=Obj_Value, group=Model)) +
  geom_line(aes(color = Model)) +
  geom_point(aes(color = Model)) + 
  scale_y_continuous(breaks = seq(0, 7, by = 0.5)) + 
  labs(y = "Biomass Flux", x = "Alternative Bounds", title = "Model Robustness to Extreme Values") +
  theme_bw()


#### FVA ####

setwd('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces')
# IL4
# IL4 = read.csv('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/FVAIL4.csv')
# IL4$Equation = IL4RXN$EQUATION[match(IL4$ReactionID, IL4RXN$ID)]
# 
# IL4 = IL4 %>%
#   filter(!MinFlux == 0 & !MaxFlux == 0) %>%
#   mutate(FluxDiff = abs(MinFlux - MaxFlux)) %>%
#   relocate(FluxDiff, .after = MaxFlux) %>%
#   filter(FluxDiff < 1) %>%
#   arrange(FluxDiff)
# 
# i = arrange(as.data.frame(table(IL4$Compartment)), desc(Freq))
# 

# LPS
# LPS = read.csv('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/NewFVAResults/FVALPS.csv')
# LPS$Equation = LPSRXN$EQUATION[match(LPS$ReactionID, LPSRXN$ID)]
# 
# LPS = LPS %>%
#   filter(!MinFlux == 0 & !MaxFlux == 0) %>%
#   mutate(FluxDiff = abs(MinFlux - MaxFlux)) %>%
#   relocate(FluxDiff, .after = MaxFlux) %>%
#   filter(FluxDiff < 1) %>%
#   arrange(FluxDiff)
# 
# l = arrange(as.data.frame(table(LPS$Compartment)), desc(Freq))


# Control

# Control = read.csv('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/FVACONTROL.csv')
# Control$Equation = ControlRXN$EQUATION[match(Control$ReactionID, ControlRXN$ID)] 
# 
# Control = Control %>%
#   filter(!MinFlux == 0 & !MaxFlux == 0) %>%
#   mutate(FluxDiff = abs(MinFlux - MaxFlux)) %>%
#   relocate(FluxDiff, .after = MaxFlux) %>%
#   filter(FluxDiff < 1) %>%
#   arrange(FluxDiff)
# 
# 
# c = arrange(as.data.frame(table(Control$Compartment)), desc(Freq))
# 
# x <- list(
#   LPS = as.character(l$Var1), 
#   IL4 = as.character(i$Var1), 
#   Control = as.character(c$Var1)
# )
# 
# ggvenn(
#   x, 
#   fill_color = c("#E41A1C", "#377EB8", "#4DAF4A"),
#   stroke_size = 0.5, set_name_size = 4) + ggtitle("FVA - biomass - Similarity by Subsystem")


#### Examine predicted fluxes for all reactions when setting Biomass as objF ####

# file.choose()
# 
# CFlux = read.xlsx2('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults/ControlBiomassFlux.xlsx',
#                    sheetIndex = 1)
# LPSFlux = read.xlsx2('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults/LPSBiomassFlux.xlsx',
#                      sheetIndex = 1)
# IL4Flux = read.xlsx2('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults/IL4BiomassFlux.xlsx',
#                      sheetIndex = 1)
# 
# 
# df = merge(LPSFlux, IL4Flux, by = "ReactionID", all = TRUE) %>%
#      merge(., CFlux, by = "ReactionID", all = TRUE) %>%
#     filter(!FluxControl == FluxLPS & !FluxLPS == FluxIL4 & !FluxIL4 == FluxControl)
# 
# df$ReactionName = reactions$rxnRecon3DID[match(df$Reaction, reactions$rxns)]
# df$Subsystem = LPSRXN$SUBSYSTEM[match(df$ReactionID, LPSRXN$ID)]
# df$Subsystem = as.factor(df$Subsystem)
# df$Equation = LPSRXN$EQUATION[match(df$ReactionID, LPSRXN$ID)]
# df = df %>% arrange(Subsystem)
# 
# write.xlsx(df, "OptimFlux.xlsx", col.names = T)
# 
# ### check if the fluxes are significant statistically only in LPS and IL4 comparison ###
# 
# OnlyLpsIL4 = merge(IL4Flux, LPSFlux, by = "ReactionID", all = TRUE)
# 
# OnlyLpsIL4 = na.omit(OnlyLpsIL4)
# 
# OnlyLpsIL4$FluxIL4 = as.numeric(OnlyLpsIL4$FluxIL4)
# OnlyLpsIL4$FluxLPS = as.numeric(OnlyLpsIL4$FluxLPS)
# 
# OnlyLpsIL4$Pvals = t.test(OnlyLpsIL4$FluxIL4, OnlyLpsIL4$FluxLPS)$p.value
# OnlyLpsIL4$AdjPvals = p.adjust(OnlyLpsIL4$Pvals, method = 'fdr')
# 
# table(OnlyLpsIL4$Pvals < 0.05)
# 
# #### end ###
# 
#  df1 = df %>%
#    filter(!FluxControl == FluxLPS & !FluxLPS == FluxIL4 & !FluxIL4 == FluxControl) %>%
#    dplyr::filter( Subsystem == "Exchange/demand reactions")
#  
# 
# write.xlsx(df1, "ModelExchangeFluxes.xlsx", sheetName = "Sheet1", 
#            col.names = TRUE, row.names = TRUE, append = FALSE)
# 
# 
# dfNotExtreme = df1 %>%
#   filter(!rowSums(.[2:4]) >= 30) %>%
#   filter(!abs(FluxControl - FluxLPS) < 0.05 &
#          !abs(FluxControl - FluxIL4) < 0.05 &
#          !abs(FluxIL4 - FluxLPS) < 0.05)
# 
# dfNotExtreme[20,5] = "LacCer pool[s]"
# 
# mat2 = as.matrix(cbind(dfNotExtreme$FluxLPS,dfNotExtreme$FluxIL4, dfNotExtreme$FluxControl))
# rownames(mat2) = dfNotExtreme$ReactionName
# colnames(mat2) = c("FluxLPS", "FluxIL4", "FluxControl")
# annotation_col2 = data.frame(
#   Model = c("FluxLPS", "FluxIL4", "FluxControl")
# )
# rownames(annotation_col2) = annotation_col2$Model
# 
# pheatmap(mat2,
#          color = rev(colorRampPalette(brewer.pal(8, "Spectral"))(50)),
#          annotation_col = annotation_col2,
#          main = "ExRxns: Non-extreme values",
#          cluster_rows = F,
#          show_rownames = T,
#          show_colnames = F,
#          cellwidth = 22,
#          cellheight = 10,
#          display_numbers = T,
#          fontsize_row = 7,
#          number_color = "black",
#          angle_col = 45
# )
# 
# write.xlsx(dfNotExtreme, "ExRxnsNotExtreme.xlsx", sheetName = "Sheet1", 
#            col.names = TRUE, row.names = TRUE, append = FALSE)
# 
# dfExtreme = df1 %>%
#   filter(rowSums(.[2:4]) >= 30)
# 
# mat3 = as.matrix(cbind(dfExtreme$FluxLPS,dfExtreme$FluxIL4, dfExtreme$FluxControl))
# rownames(mat3) = dfExtreme$ReactionName
# colnames(mat3) = c("FluxLPS", "FluxIL4", "FluxControl")
# annotation_col3 = data.frame(
#   Model = c("FluxLPS", "FluxIL4", "FluxControl")
# )
# rownames(annotation_col3) = annotation_col3$Model
# 
# pheatmap(mat3,
#          color = rev(colorRampPalette(brewer.pal(8, "Spectral"))(50)),
#          annotation_col = annotation_col3,
#          main = "ExRxns: Extreme values",
#          cluster_rows = F,
#          show_rownames = T,
#          show_colnames = F,
#          display_numbers = T,
#          fontsize_number = 10,
#          number_color = "black",
#          cellwidth = 35,
#          cellheight = 22,
#          angle_col = 45
# )
# 
# 
# write.xlsx(dfExtreme, "ExRxnsExtreme.xlsx", sheetName = "Sheet1", 
#            col.names = TRUE, row.names = TRUE, append = FALSE)
# 
# 
# 
# 
# ### box plots for each subsystem and then saved as a separate pdf ###
# 
# df$Subsystem = str_replace_all(df$Subsystem, " ", "_")
# df$Subsystem = str_replace_all(df$Subsystem, "/", "_")
# 
# for (i in unique(df$Subsystem)) {
#   
#   #print(i)
#   test = df %>% filter(Subsystem == i) %>%
#     dplyr::select(!c(ReactionID, ReactionName, Equation)) %>%
#     mutate(FluxLPS = as.numeric(FluxLPS)) %>%
#     mutate(FluxIL4 = as.numeric(FluxIL4)) %>%
#     mutate(FluxControl = as.numeric(FluxControl)) %>%
#     melt(.)
#   
# print(ggplot(test, aes(x=variable, y=value)) +
#     geom_boxplot(aes(fill=variable),
#                  notch = F) + 
#     geom_jitter(shape=16, position=position_jitter(0.2)) +
#     scale_fill_brewer(palette="Dark2") + theme_minimal() +
#     xlab("") +
#     ylab("") +
#     ggtitle(test$Subsystem[1]))
#   
#     ggsave(paste0(i,".pdf"))
#   
# }

# df1$rowMeans = apply(df1[c('FluxControl', 'FluxLPS', 'FluxIL4')], 1, mean)
# 
# # subset(df1, complete.cases(df1))
# rownames(df1) = df1$ReactionID
# 
# mat = as.matrix(df1[c(2,3,4)])
# rownames(mat) = df1$ReactionID

# keep those reactions where at least for one model the reaction(i) diverges more than 25% from the mean
# more interesting to explore this than absolute numbers
# testfilter = mat[rowSums(abs(mat / df1$rowMeans > 0.25)) == 3 | rowSums(abs(mat / df1$rowMeans > 0.25)) == 0, ]
# 
# mat.filt = mat[!rownames(mat) %in% rownames(testfilter),]
# 
# df1.filt = df1[rownames(mat.filt),]
# 
# 
# annotation_row = data.frame(
#   Subsystems = as.factor(df1.filt$Subsystem)
# )
# annotation_col = data.frame(
#   Models = as.factor(colnames(mat.filt))
# )
# rownames(annotation_col) = colnames(mat.filt)
# 
# rownames(annotation_row) = df1.filt$ReactionName
# rownames(mat.filt) = df1.filt$ReactionName
# annotation_row$Subsystems = _all(annotation_str_replacerow$Subsystems, "/", "_")
# annotation_row$Subsystems = str_replace_all(annotation_row$Subsystems, " ", "_")
# annotation_row$Subsystems = str_replace_all(annotation_row$Subsystems, ",", "")
# annotation_row$Subsystems = str_replace_all(annotation_row$Subsystems,"3", "three")
# annotation_row$Subsystems = str_replace_all(annotation_row$Subsystems,"9", "nine")
# annotation_row$Subsystems = str_replace_all(annotation_row$Subsystems,
#                             'Beta_oxidation_of_unsaturated_fatty_acids_(n-nine)_(mitochondrial)',
#                             'Beta_oxidation_of_unsaturated_fatty_acids_(n-nine, mitochondrial)')
# annotation_row$Subsystems = str_replace_all(annotation_row$Subsystems,"\\(", "")
# annotation_row$Subsystems = str_replace_all(annotation_row$Subsystems,"\\)", "")
# annotation_row$Subsystems = str_replace_all(annotation_row$Subsystems,"-", "_")
# 
# unique(df1.filt$Subsystem)
# nb.cols <- 27
# mycolors <- colorRampPalette(brewer.pal(8, "Accent"))(nb.cols)
# 
# display.brewer.all()
# 
# ann_colors = list(
#   Models = c(FluxControl = "#7570B3", FluxIL4 = "#E7298A", FluxLPS = "#66A61E"),
#   Subsystems = c(Glycolysis___Gluconeogenesis = "#7FC97F",
#                  Butanoate_metabolism = "#8FC195",
#                  Purine_metabolism = "#A0BAAC",  
#                  Pyrimidine_metabolism = "#B1B3C3",
#                  Nucleotide_metabolism = "#C2AFCE",
#                  Arginine_and_proline_metabolism = "#D3B4B9",
#                  Glycine_serine_and_threonine_metabolism = "#E4B9A3",
#                  Oxidative_phosphorylation = "#F5BD8E",
#                  ROS_detoxification = "#FDC988",
#                  Fatty_acid_activation_cytosolic = "#FDDA8E",
#                  Fatty_acid_activation_endoplasmic_reticular = "#FEEB93",
#                  Omega_three_fatty_acid_metabolism = "#FEFC98",
#                  Carnitine_shuttle_cytosolic = "#D1DD9E",
#                  Carnitine_shuttle_mitochondrial = "#9BB5A4",
#                  Carnitine_shuttle_peroxisomal = "#658DAA",
#                  Carnitine_shuttle_endoplasmic_reticular = "#3F67AE",
#                  Beta_oxidation_of_unsaturated_fatty_acids_n_nine_mitochondrial = "#704BA0",
#                  Cholesterol_metabolism = "#A22E93",
#                  Sphingolipid_metabolism = "#D31286",
#                  Glycerophospholipid_metabolism = "#EA0C72",
#                  Folate_metabolism = "#F4D055",
#                  Retinol_metabolism = "#DD2456",
#                  Transport_reactions = "#CF3C3A",
#                  Fatty_acid_oxidation = "#C2541E", 
#                  Fatty_acid_biosynthesis = "#AD5D26",
#                  Pool_reactions = "#95603B",
#                  Exchange_demand_reactions = "#7D6350"))
# 
# 
# 
# pheatmap(mat.filt,
#          color = colorRampPalette(brewer.pal(8, "BrBG"))(20),
#          annotation_row = annotation_row,
#          annotation_col = annotation_col,
#          main = "Biomass Flux Differences - 25% Diff from Mean",
#          cluster_rows = F,
#          show_rownames = T,
#          fontsize_row = 5,
#          show_colnames = F,
#          cellwidth = 15,
#          angle_col = 45,
#          annotation_colors = ann_colors
#          )
# 
# write.csv(df1.filt, "BiomassFiltered.csv")


FluxBarplots = function(x, title) {
  ggplot(x, aes(x=ReactionName, y=value, fill=variable)) +
    geom_bar(stat="identity", position = position_dodge(), alpha = 0.75) +theme_minimal() +
    theme(axis.text = element_text(face="bold", size = 12),
          axis.title = element_text(face="bold", size = 14)) +
    labs(x = NULL, y = "Fluxes", title = title)
  
}

# Interesting subsystems with Biomass as ObjF
ppp = melt(df[df$Subsystem == "Pentose phosphate pathway",])
nuc = melt(df[df$Subsystem == "Nucleotide metabolism",])
glc = melt(df[df$Subsystem == "Glycolysis / Gluconeogenesis",])
ox = melt(df[df$Subsystem == "Oxidative phosphorylation",])
pur = melt(df[df$Subsystem == "Purine metabolism",])
pyr =  melt(df[df$Subsystem == "Pyrimidine metabolism",])

fa = df %>%
  filter(df$Subsystem == "Fatty acid activation (cytosolic)" |
         df$Subsystem == "Fatty acid activation (endoplasmic reticular)")
fa = fa[rowSums(abs(fa > 2)) >= 4, ]
fa = melt(fa)

tca = melt(df[df$Subsystem == "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism",])

ppp = FluxBarplots(ppp, "Pentose phosphate pathway")
nuc = FluxBarplots(nuc, "Nucleotide metabolism")
glc = FluxBarplots(glc, "Glycolysis / Gluconeogenesis")
ox =  FluxBarplots(ox, "Oxidative Phosphorylation")
pyr = FluxBarplots(pyr, "Pyrimidine metabolism")
pur = FluxBarplots(pur, "Purine metabolism")
fa  =  FluxBarplots(fa, "Fatty Acid Activation")
tca = FluxBarplots(tca, "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism")

g = ggarrange(ppp, nuc, glc, 
              ox, pur, pyr,
              ncol = 3, nrow = 2,
              common.legend = T)
g

annotate_figure(g, left = text_grob("Fluxes", face = "bold", rot = 90))



#### Single Reaction Deletion Analysis ####

# LPS - 1) Nitric Oxide, 2) Succinate


# NODel = read.csv("/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/NODel.csv",
#                  header = T)
#   
# NODel = NODel %>%
#   filter(NOhasEffect == 1)
# 
# NODel$ReconID = reactions$rxnRecon3DID[match(NODel$NOdelRxn, reactions$rxns)]
# NODel$Equation = LPSRXN$EQUATION[match(NODel$NOdelRxn, LPSRXN$ID)]
# NODel$Subsystem = LPSRXN$SUBSYSTEM[match(NODel$NOdelRxn, LPSRXN$ID)]


# LPSDel = read.csv('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults/LPSBiomassDel.csv')
# LPSDel$ReconID = reactions$rxnRecon3DID[match(LPSDel$delRxn, reactions$rxns)]
# LPSDel$Equation = LPSRXN$EQUATION[match(LPSDel$delRxn, LPSRXN$ID)]
# LPSDel$Subsystem = LPSRXN$SUBSYSTEM[match(LPSDel$delRxn, LPSRXN$ID)]
# LPSDel$GENE.ASSOCIATION = LPSRXN$GENE.ASSOCIATION
# 
# IL4Del = read.csv('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults/IL4BiomassDel.csv')
# IL4Del$ReconID = reactions$rxnRecon3DID[match(IL4Del$delRxn, reactions$rxns)]
# IL4Del$Equation = IL4RXN$EQUATION[match(IL4Del$delRxn, IL4RXN$ID)]
# IL4Del$Subsystem = IL4RXN$SUBSYSTEM[match(IL4Del$delRxn, IL4RXN$ID)]
# IL4Del$GENE.ASSOCIATION = IL4RXN$GENE.ASSOCIATION
# 
# 
# CtrlDel = read.csv('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults/ControlBiomassDel.csv')
# CtrlDel$ReconID = reactions$rxnRecon3DID[match(CtrlDel$delRxn, reactions$rxns)]
# CtrlDel$Equation = ControlRXN$EQUATION[match(CtrlDel$delRxn, ControlRXN$ID)]
# CtrlDel$Subsystem = ControlRXN$SUBSYSTEM[match(CtrlDel$delRxn, ControlRXN$ID)]
# CtrlDel$GENE.ASSOCIATION = ControlRXN$GENE.ASSOCIATION
# 
# LPSDel = LPSDel %>%
#   filter(grRatio == 0)
# 
# Del = function(x, title) {
#   
#     x = x %>%
#     filter(grRatio < 0.9) %>%
#     arrange(grRateKO)
#   
#   df = arrange(as.data.frame(table(x$Subsystem)), desc(Freq))
#   
#   ggplot(df, aes(x = reorder(Var1, Freq), y = Freq)) +
#     geom_bar(stat="identity", fill = "steelblue") +
#     theme_minimal() + coord_flip() + labs(x = "", y = title) +
#     theme(axis.text = element_text(face = 'bold',
#                                    size = 10))
#   
#   
#   
# }
#   
# 
# Del(IL4Del, "Essential Reaction in each Subsystem - Biomass - IL4Model")
# Del(CtrlDel, "Essential Reaction in each Subsystem - Biomass - ControlModel")
 #### Biomass Vs NO ObjF Flux Comparison ####

# NO non zero fluxes - 45/7623
# Biomass non zero fluxes - 1212/7623


# Conclusion: can't really use a objective like NO alone since it
# is neglecting important physiological needs of the cell
# but I can still plot them to see what changes:

#### Plot Flux Heatmaps -> Nitric Oxide, Succinate ####

# NO: pre-process
# RAVEN_NO = read.csv('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/LPS_NO_Obj.csv')
# 
# RAVEN_NO$ReconID = reactions$rxnRecon3DID[match(RAVEN_NO$ReactionID, reactions$rxns)]
# RAVEN_NO$Equation = LPSRXN$EQUATION[match(RAVEN_NO$ReactionID, LPSRXN$ID)]
# RAVEN_NO$Subsystem = LPSRXN$SUBSYSTEM[match(RAVEN_NO$ReactionID, LPSRXN$ID)]
# 
# RAVEN_NO = RAVEN_NO %>% 
#   filter(Flux > 0 | Flux < 0) %>%
#   mutate_at(vars(Equation), ~replace(., is.na(.), c("NO[c] <=> NO[s]",
#                                                     "NO[s] <=> "))) %>%
#   mutate_at(vars(Subsystem), ~replace(., is.na(.), c("Exchange/demand reactions")))

# new thing learned: if there are more than one NA in a single variable,
# make a vector of new values and each respective values is then added to variable
# in a step-wise manner where NA were previously. if only one value is being iterated
# over the variable, then no need for creating a vector.

# Succinate 1:  pre-process
# 
# RAVEN_Succ1 = read.csv("/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/FluxResults/Succ1.csv")
# 
# RAVEN_Succ1$ReconID = reactions$rxnRecon3DID[match(RAVEN_Succ1$ReactionID, reactions$rxns)]
# RAVEN_Succ1$Equation = LPSRXN$EQUATION[match(RAVEN_Succ1$ReactionID, LPSRXN$ID)]
# RAVEN_Succ1$Subsystem = LPSRXN$SUBSYSTEM[match(RAVEN_Succ1$ReactionID, LPSRXN$ID)]
# 
# RAVEN_Succ1 = RAVEN_Succ1 %>%
#   filter(Flux > 0 | Flux < 0) %>%
#   mutate_at(vars(Equation), ~replace(., is.na(.), c("6 H+[c] + coproporphyrin I[c] <=> coproporphyrinogen I[c]",
#                                                     "coproporphyrin I[s] <=> ",
#                                                     "ATP[c] + H2O[c] + coproporphyrin I[c] => ADP[c] + H+[c] + Pi[c] + coproporphyrin I[s]"))) %>%
#   mutate_at(vars(Subsystem), ~replace(., is.na(.), c("Transport reactions",
#                                                      "Exchange/demand reactions",
#                                                      "Transport reactions")))
# 
# # Succinate 2: GDP -> GTP
# 
# Raven_Succ2 = read.csv('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/Succ2.csv')
# 
# Raven_Succ2$ReconID = reactions$rxnRecon3DID[match(Raven_Succ2$ReactionID, reactions$rxns)]
# Raven_Succ2$Equation = LPSRXN$EQUATION[match(Raven_Succ2$ReactionID, LPSRXN$ID)]
# Raven_Succ2$Subsystem = LPSRXN$SUBSYSTEM[match(Raven_Succ2$ReactionID, LPSRXN$ID)]
# 
# Raven_Succ2 = Raven_Succ2 %>%
#   filter(Flux > 0 | Flux < 0)
# 
# 
# #  Succinate 3: FADH
# 
# Raven_Succ3 = read.csv('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/Succ3.csv')
# 
# Raven_Succ3$ReconID = reactions$rxnRecon3DID[match(Raven_Succ3$ReactionID, reactions$rxns)]
# Raven_Succ3$Equation = LPSRXN$EQUATION[match(Raven_Succ3$ReactionID, LPSRXN$ID)]
# Raven_Succ3$Subsystem = LPSRXN$SUBSYSTEM[match(Raven_Succ3$ReactionID, LPSRXN$ID)]
# 
# Raven_Succ3 = Raven_Succ3 %>%
#   filter(Flux > 0 | Flux < 0)
# 
# HeatmapPlot(RAVEN_NO, "NO pFBA")
# HeatmapPlot(RAVEN_Succ1, "Succinate (Semialdehyde) Optim")
# HeatmapPlot(Raven_Succ2, "Succinate (GDP -> GTP) Optim")
# HeatmapPlot(Raven_Succ3, "Succinate (FADH) Optim ")
# 
# HeatmapPlot = function(x, title) {
#   
#   x$Subsystem = str_replace_all(x$Subsystem, "/", "_")
#   x$Subsystem = str_replace_all(x$Subsystem, " ", "_")
#   x$Subsystem = str_replace_all(x$Subsystem, ",", "")
#   
#   annotation_row = data.frame(
#   Subsystem = factor(x$Subsystem)
#   )
#   row.names(annotation_row) = x$Equation
#   
#   m = as.matrix(x$Flux)
#   row.names(m) = x$Equation
#   
#   
#   ht_opt(
#     heatmap_column_title_gp = gpar(fontsize = 15),
#     legend_labels_gp = gpar(fontsize = 15,
#                             fontface = 'bold'),
#     legend_border = "black",
#     heatmap_border = TRUE,
#     annotation_border = TRUE
#   )
#   
#   row_ha <- rowAnnotation(df = annotation_row,
#                           col = list(Subsystem = c(Arginine_and_proline_metabolism ='#66C2A5', 
#                                                    Exchange_demand_reactions =  '#FC8D62',
#                                                    Folate_metabolism = '#8DA0CB',
#                                                    Glycolysis___Gluconeogenesis = '#E78AC3',
#                                                    Miscellaneous = '#A6D854', 
#                                                    ROS_detoxification = '#FFD92F',
#                                                    Transport_reactions = '#E5C494',
#                                                    Alanine_aspartate_and_glutamate_metabolism = '#377EB8',
#                                                    Bile_acid_biosynthesis = '#E41A1C',
#                                                    Glutathione_metabolism = '#4DAF4A',
#                                                    Glycine_serine_and_threonine_metabolism = '#984EA3',
#                                                    Oxidative_phosphorylation = '#FF7F00',
#                                                    Peptide_metabolism = '#FFFF33',
#                                                    Porphyrin_metabolism = '#A65628',
#                                                    Purine_metabolism = '#B3CDE3',
#                                                    Pyruvate_metabolism = '#7FC97F',
#                                                    ROS_detoxification = '#FDC086',
#                                                    Tricarboxylic_acid_cycle_and_glyoxylate_dicarboxylate_metabolism = '#F0027F',
#                                                    Tyrosine_metabolism = '#FFFF99'
#                           )),
#                           show_annotation_name = F,
#                           annotation_legend_param = list(legend_gp = gpar(fontsize = 20,
#                                                                           fontface = 'bold'))
#                           
#                           
#   )
#   
#   ht_list = Heatmap(as.matrix(x$Flux), name = "Flux",
#                     column_title = title,
#                     width = unit(0.5, "cm"),
#                     row_labels = x$Equation,
#                     left_annotation = row_ha,
#                     row_names_gp = gpar(fontsize = 14,
#                                         fontface = "bold"),
#                     column_title_gp = gpar(fontsize = 20, fontface = "bold"),
#                     heatmap_legend_param = list(
#                       legend_gp = gpar(fontsize = 15))
#   )
#   
#   draw(ht_list, heatmap_legend_side="left")
# }




#### Single Gene Deletion ####

# file.choose()
# 
# IL4DelGenes = read_xlsx('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults/IL4GeneDel.xlsx')
# 
# 
# IL4DelGenes = IL4DelGenes %>%
#   filter(GrowthRatio == 0) %>%
#   arrange(desc(NumReactions))
# 
# 
# 
# LPSDelGenes = read_xlsx('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults/LPSGeneDel.xlsx')
# 
# 
# LPSDelGenes = LPSDelGenes %>%
#   filter(GrowthRatio == 0) %>%
#   arrange(desc(NumReactions))
# 
# 
# 
# ControlDelGenes = read_xlsx('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults/ControlGeneDel.xlsx')
# 
# 
# ControlDelGenes = ControlDelGenes %>%
#   filter(GrowthRatio == 0) %>%
#   arrange(desc(NumReactions))
# 
# 
# IL4DelGenes$Genes %in% ControlDelGenes$Genes # all IL4 Model essential genes are shared with Control Model
# 
# LPSDelGenes$Genes %in% ControlDelGenes$Genes # LPS-specific: Fosl1, https://www.genecards.org/cgi-bin/carddisp.pl?gene=FOSL1
# 
# ControlDelGenes$Genes %in% IL4DelGenes$Genes # Control-specific: Cyp2r1, https://medlineplus.gov/genetics/gene/cyp2r1/
# 
# LPSDelGenes$Genes %in% IL4DelGenes$Genes # LPS-specific: Fosl1
# 
# IL4DelGenes$Genes[!IL4DelGenes$Genes %in% LPSDelGenes$Genes] # IL4-specific genes  compared to LPS:
# 
#" Elovl1"  "Msmo1"   "Rrm1"    "Rrm2"    "Rrm2b"   "Dhcr24"  "Nsdhl"   "Cyp51"   "Akr1a1"  "Ebp"     "Fdft1"   "Slc29a1" "Coasy"  
# "Hsd17b7" "Cmpk2"   "Hmgcs1"  "Lss"     "Slc37a4" "Sqle"    "Gnpat"   "Ppcs"    "Slc35d1" "Ugdh"  


# Entrez Gene Summary for FOSL1 Gene 
# The Fos gene family consists of 4 members: FOS, FOSB, FOSL1, and FOSL2. These genes encode leucine zipper proteins that can dimerize with proteins of the JUN family, thereby forming the transcription factor complex AP-1. As such, the FOS proteins have been implicated as regulators of cell proliferation, differentiation, and transformation. Several transcript variants encoding different isoforms have been found for this gene. [provided by RefSeq, Jul 2014]

# GeneCards Summary for FOSL1 Gene
# FOSL1 (FOS Like 1, AP-1 Transcription Factor Subunit) is a Protein Coding gene. Diseases associated with FOSL1 include Human T-Cell Leukemia Virus Type 1 and Krabbe Disease. Among its related pathways are Immune response IL-2 activation and signaling pathway and IL-1 Family Signaling Pathways. Gene Ontology (GO) annotations related to this gene include DNA-binding transcription factor activity and RNA polymerase II transcription regulatory region sequence-specific DNA binding. An important paralog of this gene is FOSL2.


## Fosl1 LPS vs Control DESEq2 log2 fold change: -2.54, adj.p-value: 4e-04


# write.csv(IL4DelGenes, 'IL4DelGenes.csv')
# write.csv(LPSDelGenes, 'LPSDelGenes.csv')
# write.csv(ControlDelGenes, 'ControlDelGenes.csv')


#### ATP yield LPS v IL4 v Control (for all objective values is max = 1000) ####
# 
# LPS_ATP = read.csv('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/FluxResults/LPS_ATP.csv')
# IL4_ATP = read.csv('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/FluxResults/IL4_ATP.csv')
# Control_ATP = read.csv('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/FluxResults/Control_ATP.csv')
# 
# 
# LPS_ATP$ReconID = reactions$rxnRecon3DID[match(LPS_ATP$ReactionID, reactions$rxns)]
# LPS_ATP$Equation = LPSRXN$EQUATION[match(LPS_ATP$ReactionID, LPSRXN$ID)]
# LPS_ATP$Subsystem = LPSRXN$SUBSYSTEM[match(LPS_ATP$ReactionID, LPSRXN$ID)]
# LPS_ATP = LPS_ATP %>% filter(LPS_ATP > 1e-6 | LPS_ATP < -1e-6) %>%
#   select(ReactionID, LPS_ATP) %>% arrange(desc(LPS_ATP))
# 
# 
# 
# IL4_ATP$ReconID = reactions$rxnRecon3DID[match(IL4_ATP$ReactionID, reactions$rxns)]
# IL4_ATP$Equation = IL4RXN$EQUATION[match(IL4_ATP$ReactionID, IL4RXN$ID)]
# IL4_ATP$Subsystem = IL4RXN$SUBSYSTEM[match(IL4_ATP$ReactionID, IL4RXN$ID)]
# IL4_ATP = IL4_ATP %>% filter(IL4_ATP > 1e-6 | IL4_ATP < -1e-6) %>%
#   select(ReactionID, IL4_ATP) %>% arrange(desc(IL4_ATP))
# 
# 
# Control_ATP$ReconID = reactions$rxnRecon3DID[match(Control_ATP$ReactionID, reactions$rxns)]
# Control_ATP$Equation = ControlRXN$EQUATION[match(Control_ATP$ReactionID, ControlRXN$ID)]
# Control_ATP$Subsystem = ControlRXN$SUBSYSTEM[match(Control_ATP$ReactionID, ControlRXN$ID)]
# Control_ATP = Control_ATP %>% filter(Control_ATP > 1e-6 | Control_ATP < -1e-6) %>%
#   select(ReactionID, Control_ATP) %>% arrange(desc(Control_ATP))
# 
# 
# ATPdf = merge(IL4_ATP, LPS_ATP, by = "ReactionID", all = TRUE) %>%
#   merge(., Control_ATP, by = "ReactionID", all = TRUE)
# 
# reactions = read.csv('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/reactions.csv')
# 
# ATPdf$ReactionName = reactions$rxnRecon3DID[match(ATPdf$Reaction, reactions$rxns)]
# ATPdf$Subsystem = ControlRXN$SUBSYSTEM[match(ATPdf$ReactionID, ControlRXN$ID)]
# ATPdf$Equation = ControlRXN$EQUATION[match(ATPdf$ReactionID, ControlRXN$ID)]
# 
# ATPdf$Subsystem = as.factor(ATPdf$Subsystem) # 133 subsystems
# 
# ATPdf = ATPdf %>% arrange(Subsystem)
# 
# # df1 = df %>%
# #   filter(!FluxControl == 0 & !FluxLPS == 0 & !FluxIL4 == 0) %>%
# #   filter(!FluxControl == 1000 & !FluxLPS == 1000 & !FluxIL4 == 1000) %>%
# #   filter(!FluxControl == -1000 & !FluxLPS == -1000 & !FluxIL4 == -1000)
# # df1$rowMeans = apply(df1[c('FluxControl', 'FluxLPS', 'FluxIL4')], 1, mean)
# rownames(ATPdf) = ATPdf$ReactionName
# 
# matATP = as.matrix(ATPdf[c(2,3,4)])
# rownames(matATP) = ATPdf$ReactionName
# 
# annotation_rowATP = data.frame(
#   Subsystems = factor(ATPdf$Subsystem, ordered = T)
# )
# annotation_colATP = data.frame(
#   Models = as.factor(colnames(matATP))
# )
# rownames(annotation_colATP) = colnames(matATP)
# 
# rownames(annotation_rowATP) = ATPdf$ReactionName
# 
# annotation_rowATP$Subsystems = str_replace_all(annotation_rowATP$Subsystems, "/", "_")
# annotation_rowATP$Subsystems = str_replace_all(annotation_rowATP$Subsystems, " ", "_")
# annotation_rowATP$Subsystems = str_replace_all(annotation_rowATP$Subsystems, ",", "")
# annotation_rowATP$Subsystems = str_replace_all(annotation_rowATP$Subsystems,"3", "three")
# annotation_rowATP$Subsystems = str_replace_all(annotation_rowATP$Subsystems,"9", "nine")
# annotation_rowATP$Subsystems = str_replace_all(annotation_rowATP$Subsystems,
#                                             'Beta_oxidation_of_unsaturated_fatty_acids_(n-nine)_(mitochondrial)',
#                                             'Beta_oxidation_of_unsaturated_fatty_acids_(n-nine, mitochondrial)')
# annotation_rowATP$Subsystems = str_replace_all(annotation_rowATP$Subsystems,"\\(", "")
# annotation_rowATP$Subsystems = str_replace_all(annotation_rowATP$Subsystems,"\\)", "")
# annotation_rowATP$Subsystems = str_replace_all(annotation_rowATP$Subsystems,"-", "_")
# 
# levels(annotation_rowATP$Subsystems)
# 
# 
# nb.cols <- 15
# mycolors <- colorRampPalette(brewer.pal(8, "Accent"))(nb.cols)
# 
# display.brewer.pal(n = 8, name = 'Set1')
# 
# display.brewer.all()
# 
# ann_colors = list(
#   Models = c(Control_ATP = "#7570B3", IL4_ATP = "#E7298A", LPS_ATP = "#66A61E"),
#   Subsystems = c(Beta_oxidation_of_even_chain_fatty_acids_mitochondrial = "#7FC97F",
#                  Bile_acid_biosynthesis = "#9EBBA9",
#                  Carnitine_shuttle_cytosolic = "#BEAED4",  
#                 Carnitine_shuttle_mitochondrial = "#DDB7AD",
#                 Exchange_demand_reactions = "#FDC086",
#                 Fatty_acid_activation_cytosolic = "#FEDF8F",
#                   Glutathione_metabolism = "#FFFF99",
#                 Glycolysis___Gluconeogenesis = "#9BB5A4",
#                 Nucleotide_metabolism = "#386CB0",
#                   Oxidative_phosphorylation = "#933797",
#                   ROS_detoxification = "#F0027F",
#                   Purine_metabolism = "#D72E4A",
#                   Pyrimidine_metabolism = "#BF5B17",
#                   Transport_reactions = "#92603E",
#                  Tricarboxylic_acid_cycle_and_glyoxylate_dicarboxylate_metabolism = "#666666"))
# 
# 
# 
# pheatmap(matATP,
#          color = colorRampPalette(brewer.pal(8, "BrBG"))(20),
#          annotation_row = annotation_rowATP,
#          annotation_col = annotation_colATP,
#          main = "ATP",
#          cluster_rows = F,
#          show_rownames = T,
#          fontsize_row = 5,
#          show_colnames = F,
#          cellwidth = 15,
#          angle_col = 45,
#          annotation_colors = ann_colors
# )
# 
# 
# write.csv(ATPdf, "ATPdf.csv")

#### NO optimization - normal FBA)
# 
# NO_Normal_FBA = read.csv(
#   '/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/NO_Normal_FBA.csv')
# 
# 
# NO_Normal_FBA$ReconID = reactions$rxnRecon3DID[match(NO_Normal_FBA$ReactionID, reactions$rxns)]
# NO_Normal_FBA$Equation = LPSRXN$EQUATION[match(NO_Normal_FBA$ReactionID, LPSRXN$ID)]
# NO_Normal_FBA$Subsystem = LPSRXN$SUBSYSTEM[match(NO_Normal_FBA$ReactionID, LPSRXN$ID)]


#### Dual NO + Biomass optimization 
# 
# Dual = read.csv(
#   '/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/FluxResults/BMplusNO.csv')
# 
# Dual$ReconID = reactions$rxnRecon3DID[match(Dual$ReactionID, reactions$rxns)]
# Dual$Equation = LPSRXN$EQUATION[match(Dual$ReactionID, LPSRXN$ID)]
# Dual$Subsystem = LPSRXN$SUBSYSTEM[match(Dual$ReactionID, LPSRXN$ID)]
# 
# Dual = Dual %>% filter(Flux > 1e-6 | Flux < -1e-6)

#### Reg Genes - check whether regulatory load genes from Fastcormics publication are present in our datasets already ####
# 
# 
# ## instead of using mapIDs from 'annotate' and 'org.Mm.eg.db' packages, go to: 
# # https://biit.cs.ut.ee/gprofiler/convert to for converting names of the genes
# # since these functions yieleded NA's for Mouse genome
# library(annotate)
# library(org.Mm.eg.db)
# 
# file.choose()
# LPSgenes = read.xlsx2('/Users/rokosango/PhD/RNA-seq/csv_xslx_files/shrink_results_LPS.xlsx',
#                     sheetIndex = 1)
# # LPSgenes$X = substr(LPSgenes$X, 1, 18)
# rownames(LPSgenes) = LPSgenes$X
# 
# Orthologs = read.xlsx2('/Users/rokosango/PhD/RNA-seq/HighRegLoadGenes/RegGenes/gProfiler_55MouseOrthologues.xls',
#                        sheetIndex = 1)
# rownames(Orthologs) = Orthologs$EnsemblMouse
# 
# LPS = LPSgenes[rownames(Orthologs),]
# LPS = na.omit(LPS)
# 
# LPS = LPS %>% arrange(padj)
# 
# pvalue = LPS$padj
# 
# annotation_row = data.frame(
#   Padj = LPS$padj
# )
# row.names(annotation_row) = LPS$symbol
# 
# mLPS = as.matrix(LPS$log2FoldChange)
# row.names(mLPS) = LPS$symbol
# colnames(mLPS) = "L2FC"
# 
# row_ha <- rowAnnotation(df = annotation_row,
#                         #col = list(Pvalue = RColorBrewer::brewer.pal(name = "Oranges", n = 9)),
#                         show_annotation_name = F,
#                         annotation_legend_param = list(legend_gp = gpar(fontsize = 20,
#                                                                         fontface = 'bold'))
#                         
#                         
# )
# 
# ht_list = Heatmap(mLPS, name = "L2FC",
#                   column_title = "LPS high regulatory load genes",
#                   width = unit(0.5, "cm"),
#                   row_labels = LPS$symbol,
#                   left_annotation = row_ha,
#                   row_names_gp = gpar(fontsize = 14,
#                                       fontface = "bold"),
#                   column_title_gp = gpar(fontsize = 20, fontface = "bold"),
#                   heatmap_legend_param = list(
#                     legend_gp = gpar(fontsize = 15))
# )
# 
# draw(ht_list, heatmap_legend_side="left")
# 
# IL4genes = read.csv('/Users/rokosango/PhD/RNA-seq/csv_xslx_files/shrink_results_IL4.csv')
# IL4genes$X = substr(IL4genes$X, 1, 18)
# rownames(IL4genes) = IL4genes$X
# 
# IL4 = IL4genes[rownames(Orthologs),]
# IL4 = na.omit(IL4)
# 
# IL4 = IL4 %>% arrange(padj)
# 
# 
# pvalue = IL4$padj
# 
# annotation_row = data.frame(
#   Padj = IL4$padj
# )
# row.names(annotation_row) = IL4$symbol
# 
# mIL4 = as.matrix(IL4$log2FoldChange)
# row.names(mIL4) = IL4$symbol
# colnames(mIL4) = "L2FC"
# row_ha <- rowAnnotation(df = annotation_row,
#                         #col = list(Pvalue = RColorBrewer::brewer.pal(name = "Oranges", n = 9)),
#                         show_annotation_name = F,
#                         annotation_legend_param = list(legend_gp = gpar(fontsize = 20,
#                                                                         fontface = 'bold'))
#                         
#                         
# )
# 
# ht_list = Heatmap(m, name = "L2FC",
#                   column_title = "IL4 high regulatory load genes",
#                   width = unit(0.5, "cm"),
#                   row_labels = IL4$symbol,
#                   left_annotation = row_ha,
#                   row_names_gp = gpar(fontsize = 14,
#                                       fontface = "bold"),
#                   column_title_gp = gpar(fontsize = 20, fontface = "bold"),
#                   heatmap_legend_param = list(
#                     legend_gp = gpar(fontsize = 15))
# )
# 
# draw(ht_list, heatmap_legend_side="left")


# Cross-check reg genes with the reactions where they occur and whether this reactions are essential for the model
# first, extract GENE ASSOCIATION from the LPSRXN (spreadsheet of all reactions)
# second, use grepl in a double for loop to obtain a logical df of participating genes in each reaction
# third, subset logical df for only TRUE cells
# fourth, find those rxns in LPSDel dataframe
# 
# df = data.frame()
# 
# for (i in 1:length(LPSRXN$GENE.ASSOCIATION)) {
#  df[i,1] = LPSRXN$GENE.ASSOCIATION[i]
#   
# }
# names(df) = 'GENE.ASSOCIATION'
# 
# 
# df1 = data.frame()
# 
# for (i in 1:length(rownames(mLPS))) {
#   for (j in 1:length(df$GENE.ASSOCIATION)) {
#     df1[j,i] = grepl(rownames(mLPS)[i], df$GENE.ASSOCIATION[j])
#   }
# }
# rownames(df1) = LPSRXN$ID
# colnames(df1) = rownames(mLPS)
# 
# df2 = df1 %>% 
#   filter_at(vars(colnames(df1)), any_vars( . == "TRUE")) %>%
#   mutate(ID = rownames(df2)) %>%
#   relocate(ID, .before = Glrx)

# subset LPSDel dataframe by the names of the reactions in df2
# 
# 
# LPSRegGenesRxns = LPSDel[LPSDel$delRxn %in% df2$ID,]
# 
# 
# LPSRegGenesRxns = LPSRegGenesRxns %>%
#   arrange(grRateKO)
# setwd("~/PhD/RNA-seq/HighRegLoadGenes")
# 
# write.csv(LPSRegGenesRxns, "LPSRegGenesRxns.csv")
# write.xlsx2(mLPS, "RegGenesLPS.xls")

## check if high regulated genes are in Control Model

ControlGEMGenes = read.xlsx2('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults/ControlGeneDel.xlsx',
                             sheetIndex = 1)
ControlGEMGenes = ControlGEMGenes %>%
  filter(!GrowthRatio == 1) # growth ratio of 1 means no change, hence we look for those that change the ratio in
# ControlGEM, and see if these are in regulated genes from LPS condition

rownames(mLPS)[!rownames(mLPS) %in% ControlGEMGenes$Genes]

# have a list of resulting genes:

# Glrx"    "Slc5a3"  "Dgkz"    "Itpkb"   "Lta4h"   "Slc38a2" "Hk2"     "Gns"     "Pi4k2a"  "Gpx1"    "Pgls"    "Slc20a1" "Pfkfb4" 
#"Slc16a3" "Mgll"    "Hk3"     "Pde8a"   "Glul"    "Acp5"    "Slc23a2" "Slc6a6"  "Ado"     "Pla2g7"  "Acp2"    "Mthfs"   "Slc1a2" 
# Fth1"    "Plpp3"   "Naglu"   "Dbi"     "Prdx1"   "Abcc3"  
#look if they're in essential reactions:

#"Plpp3" It is not
## now do the same for IL4

IL4df = data.frame()

for (i in 1:length(IL4RXN$GENE.ASSOCIATION)) {
  IL4df[i,1] = IL4RXN$GENE.ASSOCIATION[i]
  
}
names(IL4df) = 'GENE.ASSOCIATION'


IL4df1 = data.frame()

for (i in 1:length(rownames(mIL4))) {
  for (j in 1:length(IL4df$GENE.ASSOCIATION)) {
    IL4df1[j,i] = grepl(rownames(mIL4)[i], IL4df$GENE.ASSOCIATION[j])
  }
}
rownames(IL4df1) = IL4RXN$ID
colnames(IL4df1) = rownames(mIL4)

IL4df2 = IL4df1 %>% 
  filter_at(vars(colnames(IL4df1)), any_vars( . == "TRUE")) %>%
  mutate(ID = rownames(IL4df2)) %>%
  relocate(ID, .before = Glul)

# subset IL4Del dataframe by the names of the reactions in IL4df2


IL4RegGenesRxns = IL4Del[IL4Del$delRxn %in% IL4df2$ID,]


IL4RegGenesRxns = IL4RegGenesRxns %>%
  arrange(grRateKO)
setwd("~/PhD/RNA-seq/HighRegLoadGenes")

write.csv(IL4RegGenesRxns, "IL4RegGenesRxns.csv")
write.csv(mIL4, "RegGenesIL4.csv")





#### FEA Flux Enrichment Analysis ####

file.choose()

LPSFEA = read.csv('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults/FEA_LPS_Subsystems.csv')

LPSFEA = LPSFEA %>% 
  filter(Adjusted.P.value < 0.05) %>%
  mutate(RxnRatio = Enriched.set.size / Total.set.size)

ggplot(LPSFEA, aes(x = RxnRatio, y = reorder(Group, -Adjusted.P.value), size = Enriched.set.size, fill = Adjusted.P.value)) +
  geom_point(shape=21, color="black") +
  scale_size(range = c(2, 10)) + 
  scale_fill_gradient(
    low = "red",
    high = "blue",
    aesthetics = "fill")  + 
  labs(y = "", title = "LPS FEA") + 
  theme(axis.text = element_text(face = 'bold',
                                 size = 10))



  
IL4FEA = read.csv('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults/FEA_IL4_Subsystems.csv')

IL4FEA = IL4FEA %>% 
  filter(Adjusted.P.value < 0.05) %>%
  mutate(RxnRatio = Enriched.set.size / Total.set.size)

ggplot(IL4FEA, aes(x = RxnRatio, y = reorder(Group, -Adjusted.P.value), size = Enriched.set.size, fill = Adjusted.P.value)) +
  geom_point(shape=21, color="black") +
  scale_size(range = c(2, 10)) + 
  scale_fill_gradient(
    low = "red",
    high = "blue",
    aesthetics = "fill")  + 
  labs(y = "", title = "IL4 FEA") + 
  theme(axis.text = element_text(face = 'bold',
                                 size = 10))



ControlFEA = read.csv('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults/FEA_Control_Subsystems.csv')

ControlFEA = ControlFEA %>% 
  filter(Adjusted.P.value < 0.05) %>%
  mutate(RxnRatio = Enriched.set.size / Total.set.size)

ggplot(ControlFEA, aes(x = RxnRatio, y = reorder(Group, -Adjusted.P.value), size = Enriched.set.size, fill = Adjusted.P.value)) +
  geom_point(shape=21, color="black") +
  scale_size(range = c(2, 10)) + 
  scale_fill_gradient(
    low = "red",
    high = "blue",
    aesthetics = "fill")  + 
  labs(y = "", title = "Control FEA") + 
  theme(axis.text = element_text(face = 'bold',
                                 size = 10))

comb = rbind(LPSFEA, IL4FEA, ControlFEA)
comb$Model = rep(c("LPS", "IL4", "Control"), times = c(17,12,13))

comb$Model = as.factor(comb$Model)
# filter1 = IL4FEA$Group[!IL4FEA$Group %in% LPSFEA$Group]
# IL4FEA$Group[!IL4FEA$Group %in% ControlFEA$Group] # no character
# filter2 = LPSFEA$Group[!LPSFEA$Group %in% IL4FEA$Group]
# 
# IL4FEA = IL4FEA %>%
#   dplyr::filter(Group %in% filter1)
# 
# LPSFEA = LPSFEA %>%
#   dplyr::filter(Group %in% filter2)
# 
# 
# Combined = rbind(IL4FEA, LPSFEA)
# Combined$Model = rep(c("IL4", "LPS"), times = c(4,16))


# ggplot(Combined, aes(x = RxnRatio, y = reorder(Group, -Adjusted.P.value), 
#                      size = Enriched.set.size, fill = Adjusted.P.value,
#                      color = Model)) +
#   geom_point(shape=21, stroke = 3) +
#   scale_size(range = c(2, 10)) + 
#   scale_fill_gradient(
#     low = "#F7FCB9",
#     high = "#41AB5D",
#     aesthetics = "fill")  + 
#   labs(y = "", title = "Enriched Subsystems") + 
#   theme(axis.text = element_text(face = 'bold',
#                                  size = 10))


comb = comb %>%
  filter(!Group %in% c('Exchange/demand reactions', 
                           'Miscellaneous',
                           'Aminoacyl-tRNA biosynthesis'))

ggplot(comb, aes(x = RxnRatio, y = reorder(Group, -Adjusted.P.value), 
                 size = Enriched.set.size, 
                 fill = Adjusted.P.value)) +
  #color = Comparison)) +
  geom_point(shape=21, stroke = 1) +
  scale_size(range = c(2, 6)) + 
  scale_fill_gradient(
    low = "#D6604D", 
    high = "#4393C3", 
    aesthetics = "fill")  + 
  facet_wrap(vars(Model)) +
  labs(y = "", x = 'Reaction Ratio') +
  theme(axis.text = element_text(face = 'bold',
                               size = 15),
        axis.text.x = element_text(angle = 30,
                                   size = 15),
        legend.text = element_text(face = 'bold',
                                   size = 15),
        legend.title = element_text(face = 'bold',
                                    size = 15),
        strip.text.x = element_text(
          size = 15, color = "black", face = "bold.italic"))


#### Curated FBA ####

file.choose()

options(scipen = 999)

# curated fba: includes ExMet, Media Comp, manual curation to reduce amount of sink reactions,
# while keeping a functional objective (biomass).
# in percentages for each model, we have: 66% in LPS, 18% in IL4, 20% in Control GEM.
# another proof that IL4 and Control are similar, with LPS needing more sinks and always having higher biomass value,
# probably because they are needed to support aggregation of biomass precursors (proteins, dna, lipids, etc.) much 
# more than in [IL4, Control] scenarios

path = '/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults/'

LPSGapFilled = read.xlsx2(paste0(path, "LPSFluxTable.xlsx"),
                          sheetIndex = 1)
IL4GapFilled = read.xlsx2(paste0(path, "IL4FluxTable.xlsx"),
                          sheetIndex = 1)
CtrlGapFilled = read.xlsx2(paste0(path, "ControlFluxTable.xlsx"),
                          sheetIndex = 1)

prepFluxTables = function(FluxTable, allReactions, modelSpecificRxns) {
  
  FluxTable = FluxTable %>%
    dplyr::select(ReactionID, Flux)
  
  
  FluxTable$ReactionName = allReactions$rxnRecon3DID[match(FluxTable$Reaction, allReactions$rxns)]
  FluxTable$Flux = as.numeric(FluxTable$Flux)
  FluxTable$Subsystem = modelSpecificRxns$SUBSYSTEM[match(FluxTable$ReactionID, modelSpecificRxns$ID)]
  FluxTable$Subsystem = as.factor(FluxTable$Subsystem)
  FluxTable$Equation = modelSpecificRxns$EQUATION[match(FluxTable$ReactionID, modelSpecificRxns$ID)]
  
  return(FluxTable)
  
}

LPSGapFilled = prepFluxTables(LPSGapFilled, reactions, LPSRXN)
IL4GapFilled = prepFluxTables(IL4GapFilled, reactions, IL4RXN)
CtrlGapFilled = prepFluxTables(CtrlGapFilled, reactions, ControlRXN)


LPSGapFilled[6870, c("Subsystem", "Equation")] = c("Exchange/demand reactions", "NO[s] <=>")
LPSGapFilled[6869, c("Subsystem", "Equation")] = c("Transport reactions", "NO[c] <=> NO[s]")
LPSGapFilled[5153, "ReactionName"] = "NOS2"



# biomass comparison

FluxIL4Biomass = IL4GapFilled %>%
  filter(Equation == "biomass[s] <=> ")
FluxLPSBiomass  = LPSGapFilled %>%
  filter(Equation == "biomass[s] <=> ")
FluxCtrlBiomass = CtrlGapFilled %>%
  filter(Equation == "biomass[s] <=> ")

BiomassTot = rbind(FluxLPSBiomass,
                   FluxIL4Biomass,
                   FluxCtrlBiomass)
BiomassTot$Model = c("LPS", "IL4", "Control")

ggplot(BiomassTot, aes(x=Model, y=Flux, fill= Model))+
  geom_bar(stat="identity", color="black") +
  scale_fill_manual(values=c("LPS" = "#56B4E9", # explicitly map color to the factor
                             "IL4" = "#E69F00",
                             "Control" = "#984EA3")) + # explicitly change legend text
  theme_minimal() +
  theme(axis.text = element_text(face = "bold", size = 15),
        axis.title = element_text(face = "bold", size = 15),
        legend.title = element_text(face = "bold", size = 15),
        legend.text = element_text(face = "bold", size = 15)) +
  labs(y = "Biomass Growth", x = NULL)


#### top 10 producing mets in every GEM ####

Top10LPS = LPSGapFilled %>%
  filter(Subsystem == "Exchange/demand reactions") %>%
  arrange(desc(Flux)) %>%
  subset(!(Equation %in% c('biomass[s] <=>', 'H+[s] <=> ','NADP+[s] <=> ', 'NADPH[s] <=> ', 'H2O[s] <=> ', 
                           'Pi[s] <=> ', 'Na+[s] <=> ', 'Fe3+[s] <=> ', 'Fe2+[s] <=> ', 
                           'Cu2+[s] <=> ', 'Ca2+[s] <=> '))) %>%
  top_n(10, Flux)

Top10IL4 = IL4GapFilled %>%
  filter(Subsystem == "Exchange/demand reactions") %>%
  arrange(desc(Flux)) %>%
  subset(!(Equation %in% c('H+[s] <=> ','NADP+[s] <=> ', 'NADPH[s] <=> ', 'H2O[s] <=> ', 
                           'Pi[s] <=> ', 'Na+[s] <=> ', 'Fe3+[s] <=> ', 'Fe2+[s] <=> ', 
                           'Cu2+[s] <=> ', 'Ca2+[s] <=> ')))

Top10IL4 = Top10IL4[1:10,]

Top10Ctrl = CtrlGapFilled %>%
  filter(Subsystem == "Exchange/demand reactions") %>%
  arrange(desc(Flux)) %>%
  subset(!(Equation %in% c('biomass[s] <=>', 'H+[s] <=> ','NADP+[s] <=> ', 'NADPH[s] <=> ', 'H2O[s] <=> ', 
                           'Pi[s] <=> ', 'Na+[s] <=> ', 'Fe3+[s] <=> ', 'Fe2+[s] <=> ', 
                           'Cu2+[s] <=> ', 'Ca2+[s] <=> '))) %>%
  top_n(10, Flux)

Total = rbind(Top10LPS,
              Top10IL4,
              Top10Ctrl)
 
Total$Model = rep(c("LPS", "IL4", "Control"), each = 10
                  )
Total$Model = factor(Total$Model, levels=c("LPS", "IL4", "Control"))
#Total = Total %>%
#  filter(!ReactionID == "MAR07108")

ggplot(Total, aes(x=reorder(Equation, Flux), y=sqrt(Flux), fill=Model))+
  geom_bar(stat="identity", color="black") +
  coord_flip() +
  scale_fill_manual(values=c("#56B4E9", "#E69F00", "#984EA3")) +
  theme_minimal() +
  theme(axis.text = element_text(face = "bold", size = 13),
        axis.title = element_text(face = "bold", size = 15),
        legend.title = element_text(face = "bold", size = 15),
        legend.text = element_text(face = "bold", size = 15)) +
  labs(x = NULL, y = expression(sqrt(Flux)))

LPSBMfba = 0.001325928
IL4BMfba = 0.001565171
CtrlBMfba = 0.001063008


#### RelativeContribution ####

FunRelativeContrib = function(LPSModel, IL4Model, ControlModel, title) {
  
LPSFluxSum = LPSGapFilled  %>%
  dplyr::group_by(Subsystem) %>%
  dplyr::summarize(AbsSum = sum(abs(Flux), na.rm = T))

IL4FluxSum = IL4Model %>%
  dplyr::group_by(Subsystem) %>%
  dplyr::summarize(AbsSum = sum(abs(Flux), na.rm = T))

CtrlFluxSum = ControlModel %>%
  dplyr::group_by(Subsystem) %>%
  dplyr::summarize(AbsSum = sum(abs(Flux), na.rm = T))

CommonSubsystems = intersect(intersect(LPSFluxSum$Subsystem, IL4FluxSum$Subsystem),
                             CtrlFluxSum$Subsystem)

LPSFluxSum = LPSFluxSum  %>%
  dplyr::filter(Subsystem %in% CommonSubsystems) %>%
  dplyr::filter(!Subsystem == "NA")

IL4FluxSum = IL4FluxSum %>%
  dplyr::filter(Subsystem %in% CommonSubsystems) %>%
  dplyr::filter(!Subsystem == "NA")

CtrlFluxSum = CtrlFluxSum %>%
  dplyr::filter(Subsystem %in% CommonSubsystems) %>%
  dplyr::filter(!Subsystem == "NA")

TotalFlux = data.frame(Subsystem = LPSFluxSum$Subsystem,
                       LPSFlux = LPSFluxSum$AbsSum,
                       IL4Flux = IL4FluxSum$AbsSum,
                       CtrlFlux = CtrlFluxSum$AbsSum)

TotalFlux$SumAllModels = apply(TotalFlux[c(2,3,4)], MARGIN = 1, sum)

TotalFlux = TotalFlux[!rowSums(TotalFlux[c(2,3,4)]) <= 1e-2,] #if the sum across different GEMs is < 1e-2
TotalFluxRelContribution = TotalFlux %>%
  mutate(LPSRelContribution = (LPSFlux * 1) / SumAllModels) %>%
  mutate(IL4RelContribution = (IL4Flux * 1) / SumAllModels) %>%
  mutate(CtrlRelContribution = (CtrlFlux * 1) / SumAllModels) %>%
  dplyr::select(Subsystem, LPSRelContribution, IL4RelContribution, CtrlRelContribution)

TotalFluxRelContribution = na.omit(TotalFluxRelContribution)

TotalFluxRelContribution = TotalFluxRelContribution %>%
  mutate(Order = case_when(
    LPSRelContribution > IL4RelContribution & LPSRelContribution > CtrlRelContribution ~ "1",
    IL4RelContribution > LPSRelContribution & IL4RelContribution > CtrlRelContribution ~ "2",
    CtrlRelContribution > LPSRelContribution & CtrlRelContribution > IL4RelContribution ~ "3"
  ))

TotalFluxRelContribution = na.omit(TotalFluxRelContribution)

meltedDf = melt(TotalFluxRelContribution)



meltedDf = meltedDf %>%
  arrange(Order)

meltedDf$Order = as.numeric(meltedDf$Order)


ggplot(meltedDf, aes(x=reorder(Subsystem, Order), y=value, fill=variable))+
  geom_bar(position='stack', stat="identity", color="black", linewidth=0.7) +
  coord_flip() +
  scale_fill_manual(values=c("LPSRelContribution" = "#56B4E9", # explicitly map color to the factor
                             "IL4RelContribution" = "#E69F00",
                             "CtrlRelContribution" = "#984EA3"),
                    labels = c('LPS', 'IL4', 'Control')) + # explicitly change legend text
  theme_minimal() +
  theme(axis.text = element_text(face = "bold", size = 11),
        axis.title = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(face = "bold", size = 10),
        title = element_text(face = "bold")) +
  labs(y = "Relative Contribution", x = NULL, title = title) +
  guides(fill=guide_legend(title="Models:"))

}


FunRelativeContrib(LPSGapFilled, IL4GapFilled, CtrlGapFilled, "Wildtype")

# print model-specific distribution of fluxes per subsystem

bpFluxes = function(FluxDistr, title) {
  
FluxSumarize = FluxDistr %>% 
  group_by(Subsystem) %>%
  summarise(
    FluxSum = sum(abs(Flux))
  ) %>%
  filter(FluxSum > 0.005)


test = FluxDistr %>% 
    select(c(Flux, Subsystem)) %>%
    mutate(Flux = as.numeric(Flux)) %>%
    filter(Subsystem %in% FluxSumarize$Subsystem) %>%
    melt(.)
  
ggplot(test, aes(x=Subsystem, y=value)) +
          geom_boxplot(aes(fill=variable),
                       notch = F) + 
          coord_flip() +
          geom_jitter(shape=16, position=position_jitter(0.2)) +
          scale_fill_brewer(palette="Dark2") + theme_minimal() +
          labs(x = "", y = "Flux", title = title)


}

FluxDistr = rbind(LPSGapFilled, IL4GapFilled, CtrlGapFilled)

FluxDistr$Model = rep(c("LPS", "IL4", "Control"), times = c(7950, 7950, 7668))

FluxSumarize = FluxDistr %>% 
  group_by(Subsystem)
FluxSumarize$Order = TotalFluxRelContribution$Order[match(FluxSumarize$Subsystem,
                                                          TotalFluxRelContribution$Subsystem)]
FluxSumarize = FluxSumarize %>%
  filter(!is.na(Order)) %>% #remove NAs
  filter(!Flux == 0) %>%
  filter(!Subsystem %in% c("Transport reactions", "Exchange/demand reactions"))

FluxSumarize$Order = as.numeric(FluxSumarize$Order)

ggplot(FluxSumarize, aes(x=Flux, y=reorder(Subsystem, Order), fill=Model)) +
  geom_boxplot() +
  #xlim(-2.5,2.5) +
  #geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_manual(values=c("LPS" = "#56B4E9", # explicitly map color to the factor
                             "IL4" = "#E69F00",
                             "Control" = "#984EA3")) +
  theme_minimal() +
  theme(axis.text = element_text(face = "bold", size = 10),
        axis.title = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(face = "bold", size = 10)) +
  labs(x = "Flux", y = "", title = "Sorted Boxplot for all 3 Models")


# bpFluxes(CtrlGapFilled, "Control Flux Activity")
# bpFluxes(IL4GapFilled, "IL4 Flux Activity")
# bpFluxes(LPSGapFilled, "LPS Flux Activity")


#### Creation of Flux Maps #####
# First needed LPS Interaction File extracted from revised ExtractedMetaboliteGSN.m #
# Filter it by flux-containing reactions as calculated via FBA
# 
# LPSInteractionFile = read.table('/Users/rokosango/PhD/Modelling/MATLAB_scripts_workspaces/LPSReactionInteraction.txt',
#                                 sep = ",", col.names = c("Reaction1", "Reaction2"))
# 
# 
# NoZeroFluxesLPS = LPSGapFilled %>%
#   filter(!abs(Flux) < 1e-3) %>%
#   filter(Subsystem %in% c("Glycolysis / Gluconeogenesis",
#                           "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism",
#                           "Pentose phosphate pathway",
#                           "Oxidative phosphorylation",
#                           "Arginine and proline metabolism"
#                           ))
# 
# 
# NoZeroFluxesIL4 = IL4GapFilled %>%
#   filter(!abs(Flux) < 1e-3)
# 
# 
# NoZeroFluxesControl = CtrlGapFilled %>%
#   filter(!abs(Flux) < 1e-3)
# 
# 
# LPS_RSN = LPSInteractionFile$Reaction1 %in% NoZeroFluxesLPS$ReactionID
# 
# #how many reactions from FBA non-zero rxns in interaction file, reaction 1 column?
# table(LPS_RSN) #FALSE    TRUE 
#                #2045220  231341 
# test = LPSInteractionFile[LPS_RSN, ]
# LPS_RSN2 = LPSInteractionFile$Reaction2 %in% NoZeroFluxesLPS$ReactionID 
# 
# 
# #how many reactions from FBA non-zero rxns in interaction file, reaction 2 column?
# table(LPS_RSN2) #FALSE    TRUE 
#                 #2062353  214208
# test2 = LPSInteractionFile[LPS_RSN2, ]
# intersect(test, test2) #22906
# LPS_RSN_Combined = rbind(test, test2)
# write.table(LPS_RSN_Combined, 'LPS_RSN_Combined.txt')
# 


### idea 1: look at the stacked barplot of subsystems, 
### idea 2: look at select interesting pathways in central carbon metabolism,
# then the same pathways in the other GEMs 
### then decide which systems to 
# ### include in the flux map configuration
# 
# 
# NoZeroFluxesLPS = LPSGapFilled %>%
#   filter(!abs(Flux) < 1e-3) %>%
#   filter(Subsystem %in% c("Glycolysis / Gluconeogenesis"
#                           #"Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism",
#                           #"Pentose phosphate pathway"
#                           #"Oxidative phosphorylation",
#                           #"Arginine and proline metabolism",
#                           #"Fatty acid biosynthesis"
#   ))
# 
# 
# LPS_RSN = LPSInteractionFile$Reaction1 %in% NoZeroFluxesLPS$ReactionID
# 
# #how many reactions from FBA non-zero rxns in interaction file, reaction 1 column?
# table(LPS_RSN) #FALSE    TRUE 
# #2045220  231341 
# test = LPSInteractionFile[LPS_RSN, ]
# LPS_RSN2 = LPSInteractionFile$Reaction2 %in% NoZeroFluxesLPS$ReactionID 
# 
# 
# #how many reactions from FBA non-zero rxns in interaction file, reaction 2 column?
# table(LPS_RSN2) #FALSE    TRUE 
# #2062353  214208
# test2 = LPSInteractionFile[LPS_RSN2, ]
# intersect(test, test2) #22906
# LPS_RSN_Combined = rbind(test, test2)
# write.table(LPS_RSN_Combined, 'LPS_RSN_Combined2.txt')
# 
# 
# write.xlsx2(NoZeroFluxesLPS, "NoZeroFluxesLPS.xlsx")
# 
# 

### conclusion: too many reactions, need to condense somehow the map
## maybe groupby by the subsystem and then take mean 

# LPSFluxSumarize = LPSGapFilled %>% 
#   group_by(Subsystem) %>%
#   summarise(
#     FluxSumLPS = sum(abs(Flux))
#   ) %>%
#   filter(FluxSumLPS > 0.005)
# 
# IL4FluxSumarize = IL4GapFilled %>% 
#   group_by(Subsystem) %>%
#   summarise(
#     FluxSumIL4 = sum(abs(Flux))
#   ) %>%
#   filter(FluxSumIL4 > 0.005)
# 
# LPSGEM = LPSFluxSumarize[LPSFluxSumarize$Subsystem %in% intersect(LPSFluxSumarize$Subsystem, IL4FluxSumarize$Subsystem),]
# IL4GEM = IL4FluxSumarize[IL4FluxSumarize$Subsystem %in% intersect(IL4FluxSumarize$Subsystem, LPSFluxSumarize$Subsystem),]
# 
# Combined = data.frame(Subsystem = LPSGEM$Subsystem,
#                      LPSFluxSum = LPSGEM$FluxSumLPS, 
#                      IL4FluxSum = IL4GEM$FluxSumIL4)



FluxDistr = rbind(LPSGapFilled, IL4GapFilled, CtrlGapFilled)

FluxDistr$Model = rep(c("LPS", "IL4", "Control"), times = c(7950, 7950, 7668))

# FluxSumarize = FluxDistr %>% 
#   group_by(Subsystem) %>%
#   summarise(
#     FluxSum = sum(abs(Flux))
#   ) %>%
#   filter(FluxSum > 0.005)

# FluxSumarize = FluxDistr %>% 
#   group_by(Subsystem)


#### ANOVA for flux values ####
# equal variance testing in the case of Glyco
# test = FluxDistr %>%
#   filter(Subsystem == "Glycolysis / Gluconeogenesis")

#Bartlett test of homogeneity of variances

#data:  Flux by Model
#Bartlett's K-squared = 50.26, df = 2, p-value = 1.22e-11


# test = test %>%
#   mutate(FluxSquared = Flux**2)
# 
# bartlett.test(Flux ~ Model, test)
# 
# 
# ac = aov(FluxSquared ~ Model, test)
# summary(ac)
# plot(ac)



# subsystems.list = list()
# SubsystemNames = c()
# 
#   
# for (i in unique(FluxSumarize$Subsystem)){
#   df = FluxSumarize %>%
#     filter(Subsystem == i) %>%
#     mutate(FluxSquared = Flux ** 2) #squared because there are negative data points which then interfers with ANOVA
#   
#   column <- df$Subsystem[1]
#   
#   ac = aov(FluxSquared ~ Model, df)
#   
#   result <- summary(ac) # summarize each ANOVA in a table
#   hsd <- TukeyHSD(ac)
#   
#   print(column)
#   print(result)
#   print(hsd)
#   
#   subsystems.list[[i]] = hsd
# }
# 
# non_na_subs = names(subsystems.list)[41:129]
# subsystems.list = subsystems.list[non_na_subs]


# setwd('/Users/rokosango/PhD/Modelling/Flux_Reaction_Network')
# 
# na_sub_names = read.xlsx2('na_sub_names.xlsx', sheetIndex = 1)
# 
# subsystems.list = subsystems.list[!names(subsystems.list) %in% na_sub_names$Name]
# 
# for (i in names(subsystems.list)) {
#   if (any(subsystems.list[[i]][["Model"]][,4] < 0.1)) {
#     SubsystemNames[i] = subsystems.list[[i]]
#   }
# }
# 
# test = data.frame(values = unlist(SubsystemNames))
# 
# test$Group = rownames(test)
# 
# 
# test1 = test %>%
#   filter(str_detect(Group, c("10","11","12"))) %>%
#   mutate(Comparison = rep(c("IL4-Control", "LPS-Control", "LPS-IL4"), times=24)) %>%
#   mutate(Subsystem = rep(names(SubsystemNames), each = 3)) %>%
#   dplyr::select(-Group) %>%
#   dplyr::filter(values < 0.1) %>%
#   arrange(values)

# write.xlsx2(test1, "PostHocSigSubsystems.xlsx")

# bpdf = FluxDistr %>% 
#   dplyr::select(c(Flux, Subsystem, Model)) %>%
#   dplyr::filter(Subsystem %in% names(SubsystemNames)) %>%
#   mutate(Flux = as.numeric(Flux)) %>%
#   filter(Subsystem %in% FluxSumarize$Subsystem)
#   #melt(.)
# 
# order = data.frame(order = unique(test1$Subsystem),
#                    number = 1:24)
# 
# bpdf$Order = order$number[match(bpdf$Subsystem, order$order)]
# 
# 
# bpdf$Model = as.factor(bpdf$Model)
# 
# bpdf = bpdf %>%
#   filter(Subsystem != "Transport reactions") #remove transport rxns as they are not informative
# 
# 
# ggplot(bpdf, aes(x=Flux, y=reorder(Subsystem, -Order), fill=Model)) +
#   geom_bar(stat = 'identity') +
#   #geom_boxplot() +
#   #geom_jitter(shape=16, position=position_jitter(0.2)) +
#   scale_fill_manual(values=c("LPS" = "#56B4E9", # explicitly map color to the factor
#                              "IL4" = "#E69F00",
#                              "Control" = "#984EA3")) +
#   theme_minimal() +
#   labs(x = "Flux", y = "", title = "Sorted Flux Distribution for all GEMs PostHoc ANOVA")


#### Curated FVA ####
#load all 3 datasets containing fluxes for models (IL4RXN etc., plus "reactions" file)
options(scipen = 999)
path = '/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/NewFVAResults/'

fvaprep = function(df, rxnlist, AllKnownRxnsFile) {
  InterestingSubsystems = c("Glycolysis / Gluconeogenesis",
                            "Pentose phosphate pathway",
                            "Arginine and proline metabolism",
                            "Oxidative phosphorylation",
                            "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism")
  
  df$Equation = rxnlist$EQUATION[match(df$ReactionID, rxnlist$ID)]
  df$ReactionName = AllKnownRxnsFile$rxnRecon3DID[match(df$Reaction, AllKnownRxnsFile$rxns)]
  df = df %>%
    #dplyr::filter_if(., is.numeric, any_vars((.) != 0)) %>%
    dplyr::filter(MinFlux != 0 & MaxFlux !=0) %>%
    mutate(FluxDiff = abs(MinFlux - MaxFlux)) %>%
    relocate(FluxDiff, .after = MaxFlux) %>%
    relocate(ReactionName, .after = ReactionID) %>%
    dplyr::filter(Compartment %in% InterestingSubsystems) %>%
    arrange(FluxDiff)
  
  
}

#IL4
IL4 = read.csv(paste0(path, 'FVAIL4.csv'))
IL4 = fvaprep(IL4, IL4RXN, reactions)

#i = arrange(as.data.frame(table(IL4$Compartment)), desc(Freq))


# LPS
LPS = read.csv(paste0(path, 'FVALPS.csv'))
LPS = fvaprep(LPS, LPSRXN, reactions)

#l = arrange(as.data.frame(table(LPS$Compartment)), desc(Freq))


# Control
Control = read.csv(paste0(path, 'FVACONTROL.csv'))
Control = fvaprep(Control, ControlRXN, reactions)

#c = arrange(as.data.frame(table(Control$Compartment)), desc(Freq))

# LPS Knockout 

KO_LPS = read.csv(paste0(path, 'FVA_KO_LPS.csv'))

#x <- list(
#   LPS = as.character(l$Var1), 
#   IL4 = as.character(i$Var1), 
#   Control = as.character(c$Var1)
# )

# ggvenn(
#   x, 
#   fill_color = c("#E41A1C", "#377EB8", "#4DAF4A"),
#   stroke_size = 0.5, set_name_size = 4) + ggtitle("FVA - biomass - Similarity by Subsystem")
# 


#barplot


FVAdf = rbind(LPS[1:10, ], IL4[1:10, ], Control[1:10, ])
FVAdf$Model = rep(c("LPS", "IL4", "Control"), each = 10)
FVAdf$Model = factor(FVAdf$Model, levels = c("LPS", "IL4", "Control"))
FVAdf$Compartment = str_replace_all(FVAdf$Compartment, 
                                    "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism",
                                    "TCA cycle & glyoxylate/dicarboxylate")


FVAdf[12, 7] = "H2O [c] + arginine [c] ⇒ ornithine [c] + urea [c]"
FVAdf[14, 7] = "1-pyrroline-5-carboxylate [m] + H+ [m] + H2O [m] ⇔ L-glutamate 5-semialdehyde [m]"


ggplot(FVAdf, aes(x = Compartment, fill=Model))+
  geom_bar(stat = "count", position = position_dodge()) +
  coord_flip() +
  scale_y_continuous(breaks = seq(0, 7, by = 1)) +
  scale_fill_manual(values=c("LPS" = "#56B4E9", # explicitly map color to the factor
                             "IL4" = "#E69F00",
                             "Control" = "#984EA3")) +
  theme(axis.text = element_text(face = "bold", size = 8),
        axis.title = element_text(face = "bold", size = 8),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(face = "bold", size = 10)) +
  labs(x = NULL, y = "Reaction Counts", title = "FVA top 10 reactions smallest flux range")



# another idea: one plot per subsystem - 
# in each plot bars are flux values, while error bars are FVA min and max.
# in each plot one bar per model

fvabarplotprep = function(fba, fva) {
  
  df = fba %>%
    dplyr::filter(ReactionID %in% fva$ReactionID)
  
  df = df %>%
    mutate(FVAMin = fva$MinFlux[match(df$ReactionID, fva$ReactionID)]) %>%
    mutate(FVAMax = fva$MaxFlux[match(df$ReactionID, fva$ReactionID)]) %>%
    mutate(FluxDiff = fva$FluxDiff[match(df$ReactionID, fva$ReactionID)]) %>%
    relocate(FVAMin, .after = Flux) %>%
    relocate(FVAMax, .after = FVAMin) %>%
    relocate(FluxDiff, .after = FVAMax) %>%
    dplyr::filter(FluxDiff < 10) # get rid of reactions that have min and max FLux -1000 and 1000, respectively
  
#if(fba == TRUE) {print("this is IL4")}
    #for (i in 23:26) {
     # df[i, 7] = "Arginine and proline metabolism"
     # }
  
  return(df)
  
}


IL4df = fvabarplotprep(IL4GapFilled, IL4)
LPSdf = fvabarplotprep(LPSGapFilled, LPS)
Ctrldf = fvabarplotprep(CtrlGapFilled, Control)


for (i in 23:26) {
 IL4df[i, 7] = "Arginine and proline metabolism"
 }


FVACombined = rbind(IL4df, LPSdf, Ctrldf)

FVACombined$Model = rep(c("IL4", "LPS", "Control"), times = c(dim(IL4df)[1], 
                                                              dim(LPSdf)[1],
                                                              dim(Ctrldf)[1]))


FVACombined$ReactionName = str_replace_all(FVACombined$ReactionName, 
                                        "r0145",
                                        "NOS2")

FVACombined = FVACombined %>%
  mutate(Order = case_when(
    Subsystem == "Glycolysis / Gluconeogenesis" ~ "1",
    Subsystem == "Pentose phosphate pathway" ~ "2",
    Subsystem == "Arginine and proline metabolism" ~ "3",
    Subsystem == "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism" ~ "4",
    Subsystem == "Oxidative phosphorylation" ~ "5")) %>%
  mutate(Order = as.numeric(Order))

FVACombined$Subsystem = str_replace_all(FVACombined$Subsystem, 
                "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism",
                "TCA")
FVACombined$Subsystem = str_replace_all(FVACombined$Subsystem, 
                                        "Glycolysis / Gluconeogenesis",
                                        "Glyco")
FVACombined$Subsystem = str_replace_all(FVACombined$Subsystem, 
                                        "Pentose phosphate pathway",
                                        "PPP")
FVACombined$Subsystem = str_replace_all(FVACombined$Subsystem, 
                                        "Arginine and proline metabolism",
                                        "Arginine/Proline")
FVACombined$Subsystem = str_replace_all(FVACombined$Subsystem, 
                                        "Oxidative phosphorylation",
                                        "OxPhos")



ggplot(FVACombined, aes(x = reorder(ReactionName, Order), y = Flux,
                        fill=Model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values=c("LPS" = "#56B4E9", # explicitly map color to the factor
                             "IL4" = "#E69F00",
                             "Control" = "#984EA3")) +
  geom_errorbar(aes(ymin=FVAMin, ymax=FVAMax), width=.2,
                position=position_dodge(.9)) +
  facet_grid(rows = vars(Subsystem),
             scales = "free_y",
             space = "free_y",
             switch = "y") +
  coord_flip() +
  theme(axis.text = element_text(face = "bold", size = 10),
        axis.title = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "bold", size = 14),
        strip.placement = "outside",                      # Place facet labels outside x axis labels.
        strip.background = element_rect(fill = "white"),
        strip.text.y = element_text(
          size = 9, color = "black", face = "bold.italic")) +  # Make facet label background white.
  labs(x = NULL, y = "Flux")
#### Targeted subsystems ####



LPSGlyco = LPSGapFilled %>%
  filter(Subsystem == "Glycolysis / Gluconeogenesis")
#write.xlsx2(LPSGlyco, "LPSGlyco.xlsx")


IL4Glyco = IL4GapFilled %>%
  filter(Subsystem == "Glycolysis / Gluconeogenesis")
#write.xlsx2(IL4Glyco, "IL4Glyco.xlsx")


ControlGlyco = CtrlGapFilled %>%
  filter(Subsystem == "Glycolysis / Gluconeogenesis")
#write.xlsx2(ControlGlyco, "ControlGlyco.xlsx")


GlycoNodes = read.csv("GlycolysisPathwayNodes.csv")

GlycoNodes = GlycoNodes %>%
  relocate(Control_Flux, .after = Equation) %>%
  mutate(WhichModel = case_when(
    abs(LPS_Flux) > abs(IL4_Flux) & abs(LPS_Flux) > abs(Control_Flux) ~ "LPS",
    abs(IL4_Flux) > abs(LPS_Flux) & abs(IL4_Flux) > abs(Control_Flux) ~ "IL4",
    abs(Control_Flux) > abs(IL4_Flux) & abs(Control_Flux) > abs(LPS_Flux) ~ "Control"
  ))

write.csv(GlycoNodes, "GlycoNodesNewInfo.csv")




#### CheckTasksGenes ####

EssGenesPrep = function(GeneVector, GeneMatrix) {
  
    Tasks = '/Users/rokosango/PhD/MetabModelling/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.xlsx'
    Tasks = read.xlsx2(Tasks, sheetIndex = 1)
    TasksFiltered = Tasks$DESCRIPTION[Tasks$DESCRIPTION != ""]
  
    Path = '/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/EssentialGenes/'
    GeneVector = read.csv(paste0(Path, GeneVector), col.names = "Gene")
    EssentialGenesDf = read.csv(paste0(Path, GeneMatrix))
    EssentialGenesDf$X0 = NULL
    colnames(EssentialGenesDf) = TasksFiltered
    rownames(EssentialGenesDf) = GeneVector$Gene
    
    return(EssentialGenesDf)
}


# LPS #
LPS = EssGenesPrep('LPSGapFilledModelGenes.csv', 'LPSessentialGenes.csv')
IL4 = EssGenesPrep('IL4GapFilledModelGenes.csv', 'IL4essentialGenes.csv')
Control = EssGenesPrep('ControlGapFilledModelGenes.csv', 'CtrlessentialGenes.csv')

colSums(LPS) 

GenesToCheck = c("Cyp27a1", "Phgdh", "Vmat2", "ATGL", "Arg1", "Alox15",
                 "IRG1", "LDLR", "Setdb2", "Sphk")

matches <- grep(paste(GenesToCheck,collapse="|"), 
                rownames(Control), value=TRUE)
matches

LPS = LPS %>%
  filter(rownames(LPS) %in% matches)
IL4 = IL4 %>%
  filter(rownames(IL4) %in% matches)
Control = Control %>%
  filter(rownames(Control) %in% matches)

### checkTasksGenes not working, try singleGeneDeletion.m

LPSBM = 0.001325928
IL4BM = 0.001565171
CtrlBM = 0.001063008

Path = '/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/EssentialGenes/'

LPSGapFilledEssGenes = read.csv(paste0(Path, 'LPSGapFilledEssGenes.csv'))

IL4GapFilledEssGenes = read.csv(paste0(Path, 'IL4GapFilledEssGenes.csv'))

CtrlGapFilledEssGenes = read.csv(paste0(Path, 'CtrlGapFilledEssGenes.csv'))

GenesToCheck = c("Cyp27a1", "Phgdh", "Vmat2", "ATGL", "Arg1", "Alox15",
                 "IRG1", "LDLR", "Setdb2", "Sphk")

matches <- grep(paste(GenesToCheck,collapse="|"), 
                       LPSGapFilledEssGenes$Var4, value=TRUE)
matches

LPSGapFilledEssGenes = LPSGapFilledEssGenes %>%
  filter(Var4 %in% matches)
IL4GapFilledEssGenes = IL4GapFilledEssGenes %>%
  filter(Var4 %in% matches)
CtrlGapFilledEssGenes = CtrlGapFilledEssGenes %>%
  filter(Var4 %in% matches)



# new idea: get fba results from knockout models for LPS, IL4 and Control
#### knock-out Fluxes ####
setwd('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults')
#LPS knockout
LPSGeneDelFluxSolution = read.xlsx('LPSGeneDelFluxSolution.xlsx',
                                   sheetIndex = 1)
IL4GeneDelFluxSolution = read.xlsx('IL4GeneDelFluxSolution.xlsx',
                                   sheetIndex = 1)
CtrlGeneDelFluxSolution = read.xlsx('CtrlGeneDelFluxSolution.xlsx',
                                   sheetIndex = 1)

#for each knockout flux, make a flux table with reaction information etc.

GeneDelFlux = function(df, gene) {
  rownames(df) = df$Row
  
  df$Row = NULL
  
  df = df %>%
    dplyr::select(all_of(gene))
  
  names(df)[1] = "Flux"
  df = df %>%
    mutate(ReactionID = rownames(.))
  
  df$ReactionName = reactions$rxnRecon3DID[match(df$Reaction, reactions$rxns)]
  df$Flux = as.numeric(df$Flux)
  df$Subsystem = LPSRXN$SUBSYSTEM[match(df$ReactionID, LPSRXN$ID)]
  df$Subsystem = as.factor(df$Subsystem)
  df$Equation = LPSRXN$EQUATION[match(df$ReactionID, LPSRXN$ID)]
  
  return(df)
  
}

#make a list for each model and store each respect Knockout and Flux Distribution
LPS_Knockout_List <- list()

for (i in names(LPSGeneDelFluxSolution)[2:8]) {
  
  name =  paste("LPS_KO", i, sep = "_")
  assign(name, 
         LPSGeneDelFluxSolution[i])
  
  LPS_Knockout_List[[i]] = GeneDelFlux(LPSGeneDelFluxSolution, i)
}

IL4_Knockout_List <- list()

for (i in names(IL4GeneDelFluxSolution)[2:8]) {
  
  name =  paste("IL4_KO", i, sep = "_")
  assign(name, 
         IL4GeneDelFluxSolution[i])
  
  IL4_Knockout_List[[i]] = GeneDelFlux(IL4GeneDelFluxSolution, i)
}

Ctrl_Knockout_List <- list()

for (i in names(CtrlGeneDelFluxSolution)[2:8]) {
  
  name =  paste("Ctrl_KO", i, sep = "_")
  assign(name, 
         CtrlGeneDelFluxSolution[i])
  
  Ctrl_Knockout_List[[i]] = GeneDelFlux(CtrlGeneDelFluxSolution, i)
}

for (i in names(LPS_Knockout_List)) {
  FunRelativeContrib(LPS_Knockout_List[[i]], IL4_Knockout_List[[i]], 
                     Ctrl_Knockout_List[[i]], paste0(i, "knockout"))
  
}

FunRelativeContrib(LPS_Knockout_List[["Cyp27a1"]], IL4_Knockout_List[["Cyp27a1"]], 
                   Ctrl_Knockout_List[["Cyp27a1"]], "Cyp27a1 knockout")

FunRelativeContrib(LPS_Knockout_List[["Phgdh"]], IL4_Knockout_List[["Phgdh"]], 
                   Ctrl_Knockout_List[["Phgdh"]], "Phgdh knockout")

FunRelativeContrib(LPS_Knockout_List[["Arg1"]], IL4_Knockout_List[["Arg1"]], 
                   Ctrl_Knockout_List[["Arg1"]], "Arg1 knockout")

FunRelativeContrib(LPS_Knockout_List[["Alox15"]], IL4_Knockout_List[["Alox15"]], 
                   Ctrl_Knockout_List[["Alox15"]], "Alox15 knockout")

FunRelativeContrib(LPS_Knockout_List[["Setdb2"]], IL4_Knockout_List[["Setdb2"]], 
                   Ctrl_Knockout_List[["Setdb2"]], "Setdb2 knockout")

FunRelativeContrib(LPS_Knockout_List[["Sphk1"]], IL4_Knockout_List[["Sphk1"]], 
                   Ctrl_Knockout_List[["Sphk1"]], "Sphk1 knockout")

FunRelativeContrib(LPS_Knockout_List[["Sphk2"]], IL4_Knockout_List[["Sphk2"]], 
                   Ctrl_Knockout_List[["Sphk2"]], "Sphk2 knockout")


#check contributions within each model, i.e. LPS WT against all 8 knockouts


FunRelativeContribWithAllKnockouts = function(Model1, Model2, Model3, Model4,
                                              Model5, Model6, Model7, Model8,
                                                                      title) {
  
  models = list(Model1, Model2, Model3, Model4, Model5, Model6, Model7, Model8)                                            
  names(models) = c("WT", "Cyp27a1", "Phgdh", "Arg1", "Alox15", "Setdb2", "Sphk1", "Sphk2")
  
  CommonSubsystems = Reduce(intersect, list(models[["WT"]]$Subsystem,
                                            models[["Cyp27a1"]]$Subsystem,
                                            models[["Phgdh"]]$Subsystem,
                                            models[["Arg1"]]$Subsystem,
                                            models[["Alox15"]]$Subsystem,
                                            models[["Setdb2"]]$Subsystem,
                                            models[["Sphk1"]]$Subsystem,
                                            models[["Sphk2"]]$Subsystem))
  

for (i in names(models)) {
        
    models[[i]] = models[[i]] %>%
          dplyr::group_by(Subsystem) %>%
          dplyr::summarize(AbsSum = sum(abs(Flux), na.rm = T))
        
        models[[i]] = models[[i]] %>%
        dplyr::filter(!Subsystem == "NA") %>%
        dplyr::filter(Subsystem %in% CommonSubsystems)
}

  options(scipen = 999)
  
  TotalFlux = data.frame(Subsystem = models[["WT"]]$Subsystem,
                         WT = models[["WT"]][["AbsSum"]],
                         Cyp27a1 = models[["Cyp27a1"]][["AbsSum"]],
                         Phgdh = models[["Phgdh"]][["AbsSum"]],
                         Arg1 = models[["Arg1"]][["AbsSum"]],
                         Alox15 = models[["Alox15"]][["AbsSum"]],
                         Setdb2 =models[["Setdb2"]][["AbsSum"]],
                         Sphk1 = models[["Sphk1"]][["AbsSum"]],
                         Sphk2 = models[["Sphk2"]][["AbsSum"]])
  

  TotalFlux$SumAllModels = rowSums(TotalFlux[2:9])
  
  TotalFlux = TotalFlux[!rowSums(TotalFlux[c(2:8)]) <= 1e-2,] #if the sum across different GEMs is < 1e-2
  TotalFluxRelContribution = TotalFlux %>%
    mutate(WT_RelContribution = (WT * 1) / SumAllModels) %>%
    mutate(Cyp27a1_RelContribution = (Cyp27a1 * 1) / SumAllModels) %>%
    mutate(Phgdh_RelContribution = (Phgdh * 1) / SumAllModels) %>%
    mutate(Arg1_RelContribution = (Arg1 * 1) / SumAllModels) %>%
    mutate(Alox15_RelContribution = (Alox15 * 1) / SumAllModels) %>%
    mutate(Setdb2_RelContribution = (Setdb2 * 1) / SumAllModels) %>%
    mutate(Sphk1_RelContribution = (Sphk1 * 1) / SumAllModels) %>%
    mutate(Sphk2_RelContribution = (Sphk2 * 1) / SumAllModels) %>%
    dplyr::select(Subsystem, WT_RelContribution, Cyp27a1_RelContribution, Phgdh_RelContribution,
           Arg1_RelContribution, Alox15_RelContribution, Setdb2_RelContribution,
           Sphk1_RelContribution, Sphk2_RelContribution)
  
  TotalFluxRelContribution = na.omit(TotalFluxRelContribution)
  
  meltedDf = melt(TotalFluxRelContribution)
    
  meltedDf$Subsystem = as.factor(meltedDf$Subsystem)
  
  ggplot(meltedDf, aes(x=Subsystem, y=value, fill=variable))+
    geom_bar(position='stack', stat="identity", color="black") +
    coord_flip() +
    #scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1.1)) +
    scale_fill_manual(values=c("WT_RelContribution" = "#1B9E77", # explicitly map color to the factor
                               "Cyp27a1_RelContribution" = "#D95F02",
                               "Phgdh_RelContribution" = "#7570B3",
                               "Arg1_RelContribution" = "#E7298A",
                               "Alox15_RelContribution" = "#66A61E",
                               "Setdb2_RelContribution" = "#E6AB02",
                               "Sphk1_RelContribution" = "#A6761D",
                               "Sphk2_RelContribution" = "#666666"),
                labels = c('WT', 'Cyp27a1_KO', 'Phgdh_KO', 'Arg1_KO',
                'Alox15_KO', 'Setdb2_KO', 'Sphk1_KO', 'Sphk2_KO')) + # explicitly change legend text
    theme_minimal() +
    theme(axis.text = element_text(face = "bold", size = 11),
          axis.title = element_text(face = "bold", size = 10),
          legend.title = element_text(face = "bold", size = 10),
          legend.text = element_text(face = "bold", size = 10),
          title = element_text(face = "bold")) +
    labs(y = "Relative Contribution", x = NULL, title = title) +
    guides(fill=guide_legend(title="Models:"))
  
}



FunRelativeContribWithAllKnockouts(LPSGapFilled, LPS_Knockout_List[["Cyp27a1"]],
                                   LPS_Knockout_List[["Phgdh"]], LPS_Knockout_List[["Arg1"]],
                                   LPS_Knockout_List[["Alox15"]], LPS_Knockout_List[["Setdb2"]],
                                   LPS_Knockout_List[["Sphk1"]], LPS_Knockout_List[["Sphk2"]],
                                   "LPS Flux Predictions")

FunRelativeContribWithAllKnockouts(IL4GapFilled, IL4_Knockout_List[["Cyp27a1"]],
                                   IL4_Knockout_List[["Phgdh"]], IL4_Knockout_List[["Arg1"]],
                                   IL4_Knockout_List[["Alox15"]], IL4_Knockout_List[["Setdb2"]],
                                   IL4_Knockout_List[["Sphk1"]], IL4_Knockout_List[["Sphk2"]],
                                   "IL4 Flux Predictions")

FunRelativeContribWithAllKnockouts(CtrlGapFilled, Ctrl_Knockout_List[["Cyp27a1"]],
                                   Ctrl_Knockout_List[["Phgdh"]], Ctrl_Knockout_List[["Arg1"]],
                                   Ctrl_Knockout_List[["Alox15"]], Ctrl_Knockout_List[["Setdb2"]],
                                   Ctrl_Knockout_List[["Sphk1"]], Ctrl_Knockout_List[["Sphk2"]],
                                   "Control Flux Predictions")




######  LPS Knockout Fluxes: choose: #####
# Bile acid biosynthesis, 
# Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism
# Pyrimidine metabolism
# Heme synthesis

InterestingSubs = c("Bile acid biosynthesis", 
                    "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism",
                    "Pyrimidine metabolism", "Heme synthesis")

ProcessPlotFun = function(df1, df2, InterestingSubs, gene) {
  
   df1 = LPSGapFilled %>%
    dplyr::filter(Subsystem %in% InterestingSubs) %>%
    dplyr::group_by(Subsystem) %>%
    dplyr::summarize(AbsMean = mean(abs(Flux), na.rm = T),
              SDMean = sd(Flux, na.rm = T))
   
  df2 = LPS_Knockout_List[[gene]] %>%
    dplyr::filter(Subsystem %in% InterestingSubs) %>%
    dplyr::group_by(Subsystem) %>%
    dplyr::summarize(AbsMean = mean(abs(Flux), na.rm = T),
              SDMean = sd(Flux, na.rm = T))
  
  combined = rbind(df1, df2)
  combined$Model = rep(c("WT", "Cyp27a1"), each = 4)
  combined$Model = as.factor(combined$Model)
  combined$Subsystem = str_wrap(combined$Subsystem, width = 10)
  
  ggplot(combined, aes(x = Subsystem, y = AbsMean,
                 fill=Model)) +
    geom_bar(stat = "identity", color = "black", 
             position = position_dodge()) +
    scale_fill_manual(values=c("WT" = "#1B9E77", #D95F02
                               'Cyp27a1' = "#D95F02"), # 1B9E77
                      labels = c('Cyp27a1_KO', 'WT')) +
    geom_errorbar(aes(ymin=(AbsMean - SDMean), 
                      ymax=(AbsMean + SDMean)), 
                      width=.2,
                      position=position_dodge(.9)) +
    #coord_flip() +
    labs(x = NULL, y = "Flux") +
    theme_bw() +
    theme(axis.text = element_text(face = "bold", size = 11),
          axis.title = element_text(face = "bold", size = 10),
          legend.title = element_text(face = "bold", size = 10),
          legend.text = element_text(face = "bold", size = 10))
    
  
}

ProcessPlotFun(LPSGapFilled, LPS_Knockout_List[["Cyp27a1"]],
               InterestingSubs, 'Cyp27a1')

### getting p vals for subsystems ###
gr = unique(df1$Subsystem)
res = matrix(NA, nrow = length(gr), ncol = 1)

for (i in unique(df1$Subsystem)) {
  
  x = wilcox.test(df1$Flux[df1$Subsystem == i],
                  df2$Flux[df2$Subsystem == i])
  
  res[i] = x$p.value
  
}

#Pyrimidine metabolism 
#0.6821914 
#Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism 
#0.6351590 
#Bile acid biosynthesis 
#0.6125045 
#Heme synthesis 
#0.6170751 


####
#don't think I can use statistical testing because it gives weird results
####

#### focus on specific subsystems and reactions within:
#### Bile acid biosynthesis ####


SpecFluxesFun = function(df, df2, WhichSubsystem, title) {
      #needs FVA, Knockout and FBA dataframes.
  
      df = LPSGapFilled %>%
        dplyr::filter(Subsystem %in% WhichSubsystem) %>%
        dplyr::filter(!Flux == 0)

      df = df %>%
        dplyr::mutate(FVAMin = LPS$MinFlux[match(df$ReactionID, LPS$ReactionID)]) %>%
        dplyr::mutate(FVAMax = LPS$MaxFlux[match(df$ReactionID, LPS$ReactionID)]) %>%
        dplyr::mutate(FVAHalfPoint = (FVAMin + FVAMax) / 2) %>%
        relocate(FVAMin, .after = Flux) %>%
        relocate(FVAMax, .after = FVAMin)


      df2 = LPS_Knockout_List[["Cyp27a1"]] %>%
        dplyr::filter(Subsystem %in% WhichSubsystem) %>%
        dplyr::filter(ReactionID %in% df$ReactionID)
      
      df2 = df2 %>%
        dplyr::mutate(FVAMin = KO_LPS$MinFlux[match(df2$ReactionID, KO_LPS$ReactionID)]) %>%
        dplyr::mutate(FVAMax = KO_LPS$MaxFlux[match(df2$ReactionID, KO_LPS$ReactionID)]) %>%
        dplyr::mutate(FVAHalfPoint = (FVAMin + FVAMax) / 2) %>%
        relocate(FVAMin, .after = Flux) %>%
        relocate(FVAMax, .after = FVAMin)
      
      combined = rbind(df, df2)
      combined$Model = rep(c("WT", "Cyp27a1"), each = dim(df)[1])
      #combined$ReactionID = factor(combined$ReactionID, levels = unique(rev(sort(combined$ReactionID))))
      #uncomment above line for displaying plots in alphanumerical order.
      #however, uncommenting it  makes the density plot `arrangedplots` perform weirdly!
      
      fba = ggplot(combined, aes(x = ReactionID, y = Flux,
                           fill=Model)) +
        geom_bar(stat = "identity", color = "black", 
                 position = position_dodge()) +
        scale_fill_manual(values=c("WT" = "#1B9E77",
                                   'Cyp27a1' = "#D95F02"),
                          labels = c('Cyp27a1_KO', 'WT')) +
        coord_flip() +
        labs(x = NULL, y = "Flux", title = 'pFBA') +
        theme_bw() +
        theme(axis.text = element_text(face = "bold", size = 11),
              axis.title = element_text(face = "bold", size = 10),
              legend.title = element_text(face = "bold", size = 10),
              legend.text = element_text(face = "bold", size = 10),
              plot.title = element_text(face = "bold"))
      
      fva = ggplot(combined, aes(x = ReactionID, y = FVAHalfPoint, 
                                 color=Model, group=Model)) +
        geom_errorbar(aes(ymin=FVAMin, ymax=FVAMax), width=1,
                      position=position_dodge(0.5)) +
        scale_color_manual(values=c("WT" = "#1B9E77",
                                    'Cyp27a1' = "#D95F02"),
                           labels = c('Cyp27a1_KO', 'WT')) +
        labs(x = NULL, y = "Flux", title = 'FVA') +
        coord_flip() +
        theme_bw() +
        theme(axis.text = element_text(face = "bold", size = 11),
              axis.title = element_text(face = "bold", size = 10),
              #axis.text.y = element_blank(),
              legend.title = element_text(face = "bold", size = 10),
              legend.text = element_text(face = "bold", size = 10),
              plot.title = element_text(face = "bold"))
      
      fig = ggarrange(fba, fva,
                      ncol = 2, nrow = 1,
                      common.legend = T,
                      legend = "bottom")
      
      fig = annotate_figure(fig, top = text_grob(title, face = "bold"))
      
      return_list = list(fig, combined)
      names(return_list) = c("CombinedFigure", "CombinedDataset")
      
      return(return_list)
      
}

BileAcid = SpecFluxesFun(LPSGapFilled, LPS_Knockout_List[["Cyp27a1"]],
                     "Bile acid biosynthesis", 
                     "Bile acid biosynthesis")

PyrMet = SpecFluxesFun(LPSGapFilled, LPS_Knockout_List[["Cyp27a1"]],
                       "Pyrimidine metabolism", 
                       "Pyrimidine metabolism")

HemeSynth = SpecFluxesFun(LPSGapFilled, LPS_Knockout_List[["Cyp27a1"]],
                          "Heme synthesis", 
                          "Heme synthesis")
TCAMet = SpecFluxesFun(LPSGapFilled, LPS_Knockout_List[["Cyp27a1"]],
                       "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism", 
                       "Tricarboxylic acid cycle")


#### Random sampling by ACHR ####

LPSSampledModelRxnID = read.table('~/PhD/MetabModelling/MATLAB_scripts_workspaces/SamplingResults/LPSModelSamplingRxns.txt', 
                                  col.names = "Reaction")

IL4SampledModelRxnID = read.table('~/PhD/MetabModelling/MATLAB_scripts_workspaces/SamplingResults/IL4ModelSamplingRxns.txt', 
                                  col.names = "Reaction")

KO_LPSSampledModelRxnID = read.table('~/PhD/MetabModelling/MATLAB_scripts_workspaces/SamplingResults/KO_LPSModelSamplingRxns.txt',
                                  col.names = "Reaction")

#CommonRxns = intersect(IL4SampledModelRxnID$Reaction,
#                       LPSSampledModelRxnID$Reaction)

CommonRxns = intersect(KO_LPSSampledModelRxnID$Reaction,
                       LPSSampledModelRxnID$Reaction)

filter = function(df) {
  df = df %>%
    dplyr::filter(rownames(df) %in% CommonRxns) %>%
    arrange(rownames(.))
  return(df)
}
                       


#LPS

path = '~/PhD/MetabModelling/MATLAB_scripts_workspaces/SamplingResults/LPS/'


# load up the data: reaction IDs & flux matrices
#find reactions that are shared across LPS & IL4
#arrange rows in both conditions (10 frames per model)
# for each output file (10), find mean of each sampled flux
# explore densities of those reactions (side by side densities per model)
#focus on reactions in glycolysis, tca, oxphos, ppp and arginine pathways
# use Kruskal–Wallis test to find significant reactions
LPSdata_files <- list.files(paste0(path))  # Identify file names
LPSdata_files


for(i in 1:length(LPSdata_files)) {                              # Head of for-loop
  assign(paste0("LPSSamplingFile_", i),                                   # Read and store data frames
         read.csv(paste0(path,
                          LPSdata_files[i]),
                   header = F, 
                   col.names = paste0(rep("sample", times = 1000), rep(1:1000)),
                   row.names =  LPSSampledModelRxnID$Reaction))
}


LPS_list = list(LPSSamplingFile_1, LPSSamplingFile_2, LPSSamplingFile_3, LPSSamplingFile_4,
                LPSSamplingFile_5, LPSSamplingFile_6, LPSSamplingFile_7, LPSSamplingFile_8,
                LPSSamplingFile_9, LPSSamplingFile_10)

SampleFileNames = paste0("LPSSamplingFile_10", 1:10)

for (i in 1:length(LPS_list)) {
  
  LPS_list[[i]] = filter(LPS_list[[i]])
  names(LPS_list)[i] = SampleFileNames[i]
  
}

LPSSamplingMeanDf = Reduce(`+`, LPS_list)/length(LPS_list)

#LPSSamplingMeanDf = Reduce(`+`, mget(paste0("LPSSamplingFile_", 1:10)))/10

#IL4


path = '~/PhD/MetabModelling/MATLAB_scripts_workspaces/SamplingResults/IL4/'

IL4data_files <- list.files(paste0(path))  # Identify file names
IL4data_files

for(i in 1:length(IL4data_files)) {                              # Head of for-loop
  assign(paste0("IL4SamplingFile_", i),                                   # Read and store data frames
         read.csv(paste0(path,
                         IL4data_files[i]),
                  header = F, 
                  col.names = paste0(rep("sample", times = 1000), rep(1:1000)),
                  row.names = IL4SampledModelRxnID$Reaction))
}


IL4_list = list(IL4SamplingFile_1, IL4SamplingFile_2, IL4SamplingFile_3, IL4SamplingFile_4,
                IL4SamplingFile_5, IL4SamplingFile_6, IL4SamplingFile_7, IL4SamplingFile_8,
                IL4SamplingFile_9, IL4SamplingFile_10)

SampleFileNames = paste0("IL4SamplingFile_", 1:10)

for (i in 1:length(IL4_list)) {
  
  IL4_list[[i]] = filter(IL4_list[[i]])
  names(IL4_list)[i] = SampleFileNames[i]
  
}

IL4SamplingMeanDf = Reduce(`+`, IL4_list)/length(IL4_list)

#IL4SamplingMeanDf = Reduce(`+`, mget(paste0("IL4SamplingFile_", 1:10)))/10

# LPS Knockout

path = '~/PhD/MetabModelling/MATLAB_scripts_workspaces/SamplingResults/KO_LPS/'

KO_LPSdata_files <- list.files(paste0(path))  # Identify file names
KO_LPSdata_files

for(i in 1:length(KO_LPSdata_files)) {                              # Head of for-loop
  assign(paste0("KO_LPSSamplingFile_", i),                                   # Read and store data frames
         read.csv(paste0(path,
                         KO_LPSdata_files[i]),
                  header = F, 
                  col.names = paste0(rep("sample", times = 1000), rep(1:1000)),
                  row.names = KO_LPSSampledModelRxnID$Reaction))
}

KO_LPS_list = list(KO_LPSSamplingFile_1, KO_LPSSamplingFile_2, KO_LPSSamplingFile_3, KO_LPSSamplingFile_4,
                   KO_LPSSamplingFile_5, KO_LPSSamplingFile_6, KO_LPSSamplingFile_7, KO_LPSSamplingFile_8,
                   KO_LPSSamplingFile_9, KO_LPSSamplingFile_10)

SampleFileNames = paste0("KO_LPSSamplingFile_", 1:10)

for (i in 1:length(KO_LPS_list)) {
  
  KO_LPS_list[[i]] = filter(KO_LPS_list[[i]])
  names(KO_LPS_list)[i] = SampleFileNames[i]
  
}

#do element-wise addition and divide by the number of the list elements. i.e., get the mean of
#each sampled flux for each cell across 10 list elements
KO_LPSSamplingMeanDf = Reduce(`+`, KO_LPS_list)/length(KO_LPS_list)

#concat different sampling models to compare significance
concat_sampling_dfs = function(sampling_mean_df1, 
                               sampling_mean_df2,
                               sampling_df1_names,
                               sampling_df2_names,
                               model,
                               InterestingSubsystems) {
  
  options(scipen = 4)
  
  df = cbind(sampling_mean_df1,
             sampling_mean_df2)
  
  pvals <- apply(df, 1, function(x) {
    wilcox.test(x[1:1000], x[1001:2000])$p.value
  })
  
  adjpvals = p.adjust(pvals, method = "bonferroni")
  adjpvals = adjpvals[adjpvals < 1e-5]
  
  df = df[names(adjpvals),]
  names(df) = c(paste0(sampling_df1_names, 1:1000), paste0(sampling_df2_names, 1:1000))
  df$Subsystem = model$Subsystem[match(rownames(df),
                                              model$ReactionID)]
  df$ReactionName= model$ReactionName[match(rownames(df),
                                                   model$ReactionID)]
  
  df = df %>% 
    relocate(Subsystem, .before = LPS_Sample_1) %>%
    relocate(ReactionName, .after = Subsystem) %>%
    arrange(Subsystem) %>%
    mutate(AdjPvals = adjpvals[match(names(adjpvals),
                                     rownames(df))]) %>%
    relocate(AdjPvals, .after = ReactionName) %>%
    mutate(ReactionID = rownames(.)) %>%
    relocate(ReactionID, .after = AdjPvals) %>%
    dplyr::filter(Subsystem %in% InterestingSubsystems)
  
}

InterestingSubsystems = c("Bile acid biosynthesis",
                          "Heme synthesis",
                          "Pyrimidine metabolism",
                          "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism")


combined_mean_dfs = concat_sampling_dfs(LPSSamplingMeanDf, KO_LPSSamplingMeanDf,
                         "LPS_Sample_", "KO_LPS_Sample_",
                         LPSGapFilled, InterestingSubsystems)

#function for generating density plots from significant
#reactions sampled from 5 different subsystems

densityplotsfun = function(sampling_mean_df1, 
                           sampling_mean_df2,
                           reactionID, Model1Factor,
                           Model2Factor,
                           FBAModel1, FBAModel2,
                           title) {
  
  combined = melt(rbind(sampling_mean_df1[reactionID, ],
                        sampling_mean_df2[reactionID, ]))
  
  combined$Model = rep(c(Model1Factor, Model2Factor))
  
  m = combined$value[seq(1, length(combined$value), 2)] #every step-wise first sampling model value
  n = combined$value[seq(2, length(combined$value), 2)] #every step-wise second sampling model value
  
  Model1Bar = FBAModel1 %>%
    dplyr::filter(ReactionID == reactionID) %>%
    select(Flux)
  
  Model2Bar = FBAModel2 %>%
    dplyr::filter(ReactionID == reactionID) %>%
    select(Flux)
  
  
  ggplot(combined, aes(x=value, fill=Model)) + 
    geom_density(color="black", alpha=0.9) + 
    scale_x_continuous(trans = pseudolog10_trans) +
    #geom_vline(aes(xintercept=mean(m)),
    #           color="#1B9E77", linetype="dashed", linewidth=1) +
    #geom_vline(aes(xintercept=Model1Bar[1, 1]),
    #         color="#1B9E77", linetype="solid", linewidth=1) +
    #geom_vline(aes(xintercept=mean(n)),
    #           color="#D95F02", linetype="dashed", linewidth=1) +
    #geom_vline(aes(xintercept=Model2Bar[1, 1]),
    #           color="#D95F02", linetype="solid", linewidth=1) +
    scale_fill_manual(values=c("WT" = "#1B9E77", # explicitly map color to the factor
                               "Cyp27a1_KO" = "#D95F02")) +
    theme_bw() +
    theme(axis.text = element_text(face = "bold", size = 8),
          axis.title = element_text(face = "bold", size = 8),
          legend.title = element_text(face = "bold", size = 10),
          legend.text = element_text(face = "bold", size = 10)) +
    labs(x = NULL, y = NULL,  title = title)
  
}

densityplotsfun(LPSSamplingMeanDf, KO_LPSSamplingMeanDf,
                "MAR04145", "WT", "Cyp27a1_KO", LPSGapFilled,
                LPS_Knockout_List[["Cyp27a1"]], "MAR04141")

#use this densityplotsfun to get sampling results from rxns that are changing
#between conditions in SpecFluxesFun

unique_rxns = factor(unique_rxns, levels = unique(rev(sort(unique_rxns))))

unique_rxns = TCAMet[[2]]$ReactionID


for (i in 1:length(unique_rxns)){
          
           plotlist[[i]] = print(densityplotsfun(LPSSamplingMeanDf,
                                        KO_LPSSamplingMeanDf,
                                        unique_rxns[i], 
                                        'WT', 'Cyp27a1_KO', 
                                        LPSGapFilled, LPS_Knockout_List[["Cyp27a1"]],
                                        unique_rxns[i]))
}


arrangedplots = function(sampling_mean_df1, 
          sampling_mean_df2,
          subsystem_name,
          subsystem_list, #subsystem = "HemeSynth", "TCAMet", "PyrMet" from SpecFluxesFun. No "BileAcid"
          model1Factor,
          model2Factor,
          FBAModel1, 
          FBAModel2,
          nrow) { 
  
          plotlist = list()
  
          unique_rxns = unique(subsystem_list[[2]]$ReactionID)
          
          
    for (i in 1:length(unique_rxns)){
    
    plotlist[[i]] = print(densityplotsfun(sampling_mean_df1,
                                          sampling_mean_df2,
                                          unique_rxns[i], 
                                          model1Factor, model2Factor, 
                                          FBAModel1, FBAModel2,
                                          unique_rxns[i]))
  }
  
  fig = ggarrange(plotlist=plotlist,
                  labels = letters[1:length(plotlist)],
                  ncol = 2, nrow = nrow,
                  common.legend = T)
  
  
  fig
  #require(grid)
  anotfig = annotate_figure(fig, 
                            left = textGrob("Density", rot = 90, vjust = 1, gp = gpar(cex = 0.8)),
                            bottom = textGrob("Sampled Flux", gp = gpar(cex = 0.8)),
                            top = textGrob(subsystem_name))
  
  return(anotfig)
  
}


arrangedplots(LPSSamplingMeanDf,KO_LPSSamplingMeanDf,"TCA Cycle",
              TCAMet, "WT", "Cyp27a1_KO", LPSGapFilled,
              LPS_Knockout_List[["Cyp27a1"]], 4)

arrangedplots(LPSSamplingMeanDf,KO_LPSSamplingMeanDf,"Pyrimidine Metabolism",
              PyrMet, "LPS", "Cyp27a1_KO", LPSGapFilled,
              LPS_Knockout_List[["Cyp27a1"]], 7)

arrangedplots(LPSSamplingMeanDf,KO_LPSSamplingMeanDf, "Heme synthesis",
              HemeSynth, "LPS", "Cyp27a1_KO", LPSGapFilled,
              LPS_Knockout_List[["Cyp27a1"]], 1)



#from concatenated LPS + IL4 sampling dataset find reactions 
#that are also found in FVA dataset

FVAfiltSamplingDf = df %>%
  dplyr::filter(rownames(.) %in% FVACombined$ReactionID)

#PYK - MAR04358 - pyruvate kinase reaction:
#it is not in FVA dataset as it is only captured in LPS (not in Control or IL4) but with wide flux range {1, 1000}
# hence it is not in the FBA-FVA figure. but still it would be nice to show last step in
# glycolysis where ATP is being made:

FVAfiltSamplingDf['MAR04358', ] = df['MAR04358', ] #add reaction PYK
FVAfiltSamplingDf['MAR04379', ] = df['MAR04379', ] #add reaction PFK
FVAfiltSamplingDf['MAR04368', ] = df['MAR04368', ] #add reaction PGK
FVAfiltSamplingDf['MAR04394', ] = df['MAR04394', ] #add reaction HEX1

FVAfiltSamplingDf = FVAfiltSamplingDf[!FVAfiltSamplingDf$ReactionName %in% "r0173",]  #remove reaction r0173 (peroxisomal lactate dehydrogenase)


df[!(row.names(df) %in% c("1","2")),]

arrangedplots = function(x, subsystem, nrow) {
 
   plotlist = list()
   
   filtered = x %>%
     dplyr::filter(Subsystem == subsystem)
  
  for (i in 1:dim(filtered)[1]) {
    
    plotlist[[i]] = print(densityplotsfun(rownames(filtered)[i], filtered[i, "ReactionName"]))
    
  }
   
   fig = ggarrange(plotlist=plotlist,
                   labels = letters[1:length(plotlist)],
                   ncol = 2, nrow = nrow, common.legend = T)
   
   fig
   #require(grid)
   anotfig = annotate_figure(fig, 
                left = textGrob("Density", rot = 90, vjust = 1, gp = gpar(cex = 0.8)),
                bottom = textGrob("Sampled Flux", gp = gpar(cex = 0.8)),
                 top = textGrob(subsystem))
   
   return(anotfig)
}

FVAfiltSamplingDf %>% 
  group_by(Subsystem) %>%
  summarise(no_rows = length(Subsystem))


#Subsystem                                                        no_rows
#<fct>                                                                  <int>
#1 Arginine and proline metabolism                                        6
#2 Glycolysis / Gluconeogenesis                                           7
#3 Oxidative phosphorylation                                              2
#4 Pentose phosphate pathway                                              8
#5 Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism       2

arrangedplots(FVAfiltSamplingDf, "Arginine and proline metabolism", 3)
arrangedplots(FVAfiltSamplingDf, "Glycolysis / Gluconeogenesis", 5)
arrangedplots(FVAfiltSamplingDf, "Pentose phosphate pathway", 4)
arrangedplots(FVAfiltSamplingDf, 
              "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism", 2)
arrangedplots(FVAfiltSamplingDf, "Oxidative phosphorylation", 2)



#this is just to compare means and variances LPS vs IL4
LPSSamplingMeanDf$MeanLPS = apply(LPSSamplingMeanDf, 1, mean)
IL4SamplingMeanDf$MeanIL4 = apply(IL4SamplingMeanDf, 1, mean)
LPSSamplingMeanDf$VarLPS = apply(LPSSamplingMeanDf, 1, var)
IL4SamplingMeanDf$VarIL4 = apply(LPSSamplingMeanDf, 1, var)

SamplingDfCombined = data.frame(cbind(
                           rownames(LPSSamplingMeanDf),
                           LPSSamplingMeanDf$MeanLPS,
                           IL4SamplingMeanDf$MeanIL4,
                           LPSSamplingMeanDf$VarLPS,
                           IL4SamplingMeanDf$VarIL4))

names(SamplingDfCombined) = c(
  "ReactionID", "MeanLPS", "MeanIL4", "VarLPS", "VarIL4"
)


rownames(SamplingDfCombined) = rownames(LPSSamplingMeanDf)




