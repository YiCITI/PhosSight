library(tidyverse)
library(data.table)



if (FALSE){
  phosphors <- fread('E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/PhosphoRS/PhosphoRS.txt')
  auto_rt <- fread('E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/autoRT_Results/tf_prediction/phospho.prediction.tsv')
  pDeep <- fread('E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/pDeep3_Results/pDeep3PredictionResults.Phospho.txt')
  features <- fread('E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/Features/features.PhosphoRS.txt')
  output <- 'E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/Features/Features.Localization.entropy.txt'
  VariableMods = '1,Oxidation,M,15.994919,1;2,Phospho,S,79.966331,2;3,Phospho,T,79.966331,2;4,Phospho,Y,79.966331,2' # Identification modifications
  FixedMods = 'null' # Identification modifications
  phosSight <- fread('E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/phosSight_Results/pth_prediction/PhosSight.Predict.Phospho.txt')
} else {
  args <- commandArgs(T)
  phosphors <- fread(args[1])
  auto_rt <- fread(args[2])
  pDeep <- fread(args[3])
  features <- fread(args[4])
  output <- args[5]
  VariableMods <- args[6]
  FixedMods <- args[7]
  phosSight <- fread(args[8])
}

VariableModsInfo = unlist(strsplit(VariableMods,';',fixed = TRUE))
if (FixedMods!='null'){
  FixedModsInfo = unlist(strsplit(FixedMods,';',fixed = TRUE))
}

colnames(phosphors)[2] <- 'Title'
features$Link <- NULL
phosphors$Link <- NULL

# ------------ step 1: phosphors_autort ------------

colnames(auto_rt)[1] <- 'AutoRTSequence'

auto_rt$IsoformSequence <- auto_rt$AutoRTSequence
if (FixedMods!='null'){
  for (i in 1:length(FixedModsInfo)){
    tmp = FixedModsInfo[i]
    tmp2 = unlist(strsplit(tmp,',',fixed = TRUE))
    number = tmp2[1]
    name = tmp2[2]
    aa = tmp2[3]
    mass = tmp2[4]
    sym = tmp2[5]
    
    if (aa == 'AnyN-term'){
      next
    }
    
    auto_rt$IsoformSequence <- gsub(aa,number,auto_rt$IsoformSequence)
  }
  #auto_rt$IsoformSequence <- str_replace_all(auto_rt$AutoRTSequence,c("C" = "5"))
}else{
  auto_rt$IsoformSequence <- auto_rt$AutoRTSequence
}

colnames(auto_rt)[3] <- 'Title'
auto_rt$rtError <- abs(auto_rt$y_pred - auto_rt$y)
auto_rt$rtRatio <- ifelse(auto_rt$y_pred>=auto_rt$y, auto_rt$y/auto_rt$y_pred, auto_rt$y_pred/auto_rt$y)
auto_rt$y <- NULL
auto_rt$y_pred <- NULL
auto_rt$rtError <- NULL

auto_rt2 <- auto_rt %>%
  group_by(Title) %>%
  slice(which.max(rtRatio))

colnames(auto_rt2)[4] <- 'MaxRatio'
auto_rt2$AutoRTSequence <- NULL
auto_rt2$rtRatio <- NULL
auto_rt2$IsoformSequence <- NULL
auto_rt2$rtError <- NULL

auto_rt3 <- merge(auto_rt2, auto_rt, all.x=TRUE, by='Title')

phosphors_autort <- left_join(phosphors, auto_rt3, by = c("Title" = "Title", "IsoformSequence" = "IsoformSequence")) %>% select(Title, IsoformSequence, rtRatio, MaxRatio, Isoform.Score)
colnames(phosphors_autort)[5] <- 'PhosphoRSScore'
phosphors_autort$CombinedScoreAutoRT <- phosphors_autort$PhosphoRSScore * (phosphors_autort$rtRatio/phosphors_autort$MaxRatio)

phosphors_autort$AutoRTProb1 <- 10^(phosphors_autort$CombinedScoreAutoRT/10)

tmp <- phosphors_autort[ ,list(AutoRTProb2=sum(AutoRTProb1)), by=Title]
phosphors_autort <- left_join(phosphors_autort,tmp, by = "Title")
phosphors_autort$AutoRTProb <- phosphors_autort$AutoRTProb1/phosphors_autort$AutoRTProb2
colnames(phosphors_autort)[2] <- 'IsoformSequence_AutoRT'

# ------------ step 2: phosphors_pDeep3 ------------

colnames(pDeep)[1] <- 'Title'
colnames(pDeep)[3] <- 'IsoformModification'
pDeep$Charge <- NULL

pDeep$dot_product <- NULL
pDeep$sr_dot_product <- NULL
pDeep$spectral_contrast_angle <- NULL
pDeep$pearson_correlation <- NULL
pDeep$unweighted_entropy <- NULL

pDeep[is.na(pDeep$entropy),]$entropy <- 0

pDeep2 <- pDeep %>%
  group_by(Title) %>%
  slice(which.max(entropy))

colnames(pDeep2)[4] <- 'Max_entropy'
pDeep2$Peptide <- NULL
pDeep2$IsoformModification <- NULL


pDeep3 <- merge(pDeep2, pDeep, all.x=TRUE, by='Title')

phosphors_pDeep <- left_join(phosphors, pDeep3, by = c("Title" = "Title", "IsoformModification" = "IsoformModification")) %>% select(Title, entropy, Max_entropy, Isoform.Score,IsoformSequence)
colnames(phosphors_pDeep)[4] <- 'PhosphoRSScore'
phosphors_pDeep$CombinedScorepDeep <- phosphors_pDeep$PhosphoRSScore * (phosphors_pDeep$entropy/phosphors_pDeep$Max_entropy)
phosphors_pDeep$pDeepProb1 <- 10^(phosphors_pDeep$CombinedScorepDeep/10)

tmp <- phosphors_pDeep[ ,list(pDeepProb2=sum(pDeepProb1)), by=Title]
phosphors_pDeep <- left_join(phosphors_pDeep,tmp, by = "Title")

phosphors_pDeep$pDeepProb <- phosphors_pDeep$pDeepProb1/phosphors_pDeep$pDeepProb2
phosphors_pDeep[is.na(phosphors_pDeep$pDeepProb),]$pDeepProb <- 0
colnames(phosphors_pDeep)[5] <- 'IsoformSequence_pDeep'

# ------------ step 3: phosphors_autort_pDeep ------------

phosphors_autort_pDeep <- left_join(phosphors_autort, phosphors_pDeep, by = c("Title" = "Title", "IsoformSequence_AutoRT" = "IsoformSequence_pDeep", "PhosphoRSScore" = "PhosphoRSScore")) %>% select(Title, rtRatio, MaxRatio, entropy, Max_entropy, PhosphoRSScore, IsoformSequence_AutoRT)
phosphors_autort_pDeep$CombinedScore_autort_pDeep <- phosphors_autort_pDeep$PhosphoRSScore * (phosphors_autort_pDeep$rtRatio/phosphors_autort_pDeep$MaxRatio) * (phosphors_autort_pDeep$entropy/phosphors_autort_pDeep$Max_entropy)

phosphors_autort_pDeep$autort_pDeep_Prob1 <- 10^(phosphors_autort_pDeep$CombinedScore_autort_pDeep/10)

tmp <- phosphors_autort_pDeep[ ,list(autort_pDeep_Prob2=sum(autort_pDeep_Prob1)), by=Title]
phosphors_autort_pDeep <- left_join(phosphors_autort_pDeep,tmp, by = "Title")

phosphors_autort_pDeep$autort_pDeep_Prob <- phosphors_autort_pDeep$autort_pDeep_Prob1/phosphors_autort_pDeep$autort_pDeep_Prob2
phosphors_autort_pDeep[is.na(phosphors_autort_pDeep$autort_pDeep_Prob),]$autort_pDeep_Prob <- 0
colnames(phosphors_autort_pDeep)[7] <- 'IsoformSequence_autort_pDeep'

# ------------ step 4: phosphors_phosSight ------------

colnames(phosSight)[1] <- 'Title'
phosSight$peptide <- NULL
colnames(phosSight)[3] <- 'Detectability'

phosSight2 <- phosSight %>%
  group_by(Title) %>%
  slice(which.max(Detectability))

colnames(phosSight2)[3] <- 'MaxDetectability'
phosSight2$IsoformSequence <- NULL

phosSight3 <- merge(phosSight2, phosSight, all.x=TRUE, by='Title')

phosphors_phosSight <- left_join(phosphors, phosSight3, by = c("Title" = 
"Title", "IsoformSequence" = "IsoformSequence")) %>% select(Title, IsoformSequence, Detectability, MaxDetectability, Isoform.Score)

colnames(phosphors_phosSight)[5] <- 'PhosphoRSScore'

phosphors_phosSight$CombinedScorePhosSight <- phosphors_phosSight$PhosphoRSScore * (phosphors_phosSight$Detectability/phosphors_phosSight$MaxDetectability)

phosphors_phosSight$PhosSightProb1 <- 10^(phosphors_phosSight$CombinedScorePhosSight/10)

tmp <- phosphors_phosSight[ ,list(PhosSightProb2=sum(PhosSightProb1)), by=Title]
phosphors_phosSight <- left_join(phosphors_phosSight,tmp, by = "Title")
phosphors_phosSight$PhosSightProb <- phosphors_phosSight$PhosSightProb1/phosphors_phosSight$PhosSightProb2
colnames(phosphors_phosSight)[2] <- 'IsoformSequence_PhosSight'

# ------------ step 5: phosphors_autort_pDeep3_phosSight ------------

phosphors_autort_pDeep_phosSight <- left_join(phosphors_autort_pDeep, phosphors_phosSight, by = c("Title" = "Title", "IsoformSequence_autort_pDeep" = "IsoformSequence_PhosSight", "PhosphoRSScore" = "PhosphoRSScore")) %>% select(Title, rtRatio, MaxRatio, entropy, Max_entropy, Detectability, MaxDetectability, PhosphoRSScore, IsoformSequence_autort_pDeep)

phosphors_autort_pDeep_phosSight$CombinedScore_autort_pDeep_phosSight <- phosphors_autort_pDeep_phosSight$PhosphoRSScore * (phosphors_autort_pDeep_phosSight$rtRatio/phosphors_autort_pDeep_phosSight$MaxRatio) * (phosphors_autort_pDeep_phosSight$entropy/phosphors_autort_pDeep_phosSight$Max_entropy) * (phosphors_autort_pDeep_phosSight$Detectability/phosphors_autort_pDeep_phosSight$MaxDetectability)

phosphors_autort_pDeep_phosSight$autort_pDeep_phosSight_Prob1 <- 10^(phosphors_autort_pDeep_phosSight$CombinedScore_autort_pDeep_phosSight/10)


tmp <- phosphors_autort_pDeep_phosSight[ ,list(autort_pDeep_phosSight_Prob2=sum(autort_pDeep_phosSight_Prob1)), by=Title]
phosphors_autort_pDeep_phosSight <- left_join(phosphors_autort_pDeep_phosSight,tmp, by = "Title")

phosphors_autort_pDeep_phosSight$autort_pDeep_phosSight_Prob <- phosphors_autort_pDeep_phosSight$autort_pDeep_phosSight_Prob1/phosphors_autort_pDeep_phosSight$autort_pDeep_phosSight_Prob2
phosphors_autort_pDeep_phosSight[is.na(phosphors_autort_pDeep_phosSight$autort_pDeep_phosSight_Prob),]$autort_pDeep_phosSight_Prob <- 0
colnames(phosphors_autort_pDeep_phosSight)[9] <- 'IsoformSequence_autort_pDeep_phosSight'

# ------------------------------------------------------------

# Update features
features$PhosphoLabel <- ifelse(grepl('Phospho',features$Modification), 1, 0)
feature_phos = features[features$PhosphoLabel==1,]
feature_nonphos = features[features$PhosphoLabel==0,]

phosphors2 <- phosphors %>%
  group_by(Title) %>%
  slice(which.max(Isoform.Probability))

phosphors_autort2 <- phosphors_autort %>%
  group_by(Title) %>%
  slice(which.max(AutoRTProb))

phosphors_pDeep2 <- phosphors_pDeep %>%
  group_by(Title) %>%
  slice(which.max(pDeepProb))

phosphors_phosSight2 <- phosphors_phosSight %>%
  group_by(Title) %>%
  slice(which.max(PhosSightProb))

phosphors_autort_pDeep2 <- phosphors_autort_pDeep %>%
  group_by(Title) %>%
  slice(which.max(autort_pDeep_Prob))

phosphors_autort_pDeep_phosSight2 <- phosphors_autort_pDeep_phosSight %>%
  group_by(Title) %>%
  slice(which.max(autort_pDeep_phosSight_Prob))


phosphors2$Spectrum.ID <- NULL
phosphors2$Spectrum.PrecursorCharge <- NULL
phosphors2$Spectrum.ActivationType <- NULL
phosphors2$Peptide.ID <- NULL
phosphors2$Peptide.Sequence <- NULL
phosphors2$Peptide.SitePrediction <- NULL
phosphors2$Isoform.ID <- NULL
phosphors2$Isoform.Sites <- NULL
phosphors2$Isoform.Score <- NULL
phosphors2$Isoform.Probability <- NULL
phosphors2$IsoformModification <- NULL
colnames(phosphors2)[2] <- 'IsoformSequence_PhosphoRS'

phosphors_autort2$rtRatio <- NULL
phosphors_autort2$MaxRatio <- NULL
phosphors_autort2$PhosphoRSScore <- NULL
phosphors_autort2$CombinedScoreAutoRT <- NULL
phosphors_autort2$AutoRTProb1 <- NULL
phosphors_autort2$AutoRTProb2 <- NULL

phosphors_pDeep2$entropy <- NULL
phosphors_pDeep2$Max_entropy <- NULL
phosphors_pDeep2$PhosphoRSScore <- NULL
phosphors_pDeep2$CombinedScorepDeep <- NULL
phosphors_pDeep2$pDeepProb1 <- NULL
phosphors_pDeep2$pDeepProb2 <- NULL

phosphors_autort_pDeep2$rtRatio <- NULL
phosphors_autort_pDeep2$MaxRatio <- NULL
phosphors_autort_pDeep2$entropy <- NULL
phosphors_autort_pDeep2$Max_entropy <- NULL
phosphors_autort_pDeep2$PhosphoRSScore <- NULL
phosphors_autort_pDeep2$CombinedScore_autort_pDeep <- NULL
phosphors_autort_pDeep2$autort_pDeep_Prob1 <- NULL
phosphors_autort_pDeep2$autort_pDeep_Prob2 <- NULL

phosphors_phosSight2$Detectability <- NULL
phosphors_phosSight2$MaxDetectability <- NULL
phosphors_phosSight2$PhosphoRSScore <- NULL
phosphors_phosSight2$CombinedScorePhosSight <- NULL
phosphors_phosSight2$PhosSightProb1 <- NULL
phosphors_phosSight2$PhosSightProb2 <- NULL

phosphors_autort_pDeep_phosSight2$rtRatio <- NULL
phosphors_autort_pDeep_phosSight2$MaxRatio <- NULL
phosphors_autort_pDeep_phosSight2$entropy <- NULL
phosphors_autort_pDeep_phosSight2$Max_entropy <- NULL
phosphors_autort_pDeep_phosSight2$PhosphoRSScore <- NULL
phosphors_autort_pDeep_phosSight2$Detectability <- NULL
phosphors_autort_pDeep_phosSight2$MaxDetectability <- NULL
phosphors_autort_pDeep_phosSight2$CombinedScore_autort_pDeep_phosSight <- NULL
phosphors_autort_pDeep_phosSight2$autort_pDeep_phosSight_Prob1 <- NULL
phosphors_autort_pDeep_phosSight2$autort_pDeep_phosSight_Prob2 <- NULL

feature_phos2 <- left_join(phosphors2, feature_phos, by="Title")
feature_phos_autort <- left_join(phosphors_autort2, feature_phos2, by="Title")
feature_phos_autort_pdeep <- left_join(phosphors_pDeep2, feature_phos_autort, by="Title")
feature_phos_autort_pdeep2 <- left_join(phosphors_autort_pDeep2, feature_phos_autort_pdeep, by="Title")
feature_phos_autort_pdeep_phosSight <- left_join(phosphors_phosSight2, feature_phos_autort_pdeep2, by="Title")

# cat("=== DEBUG: feature_phos_autort_pdeep_phosSight columns ===\n")
# cat("feature_phos_autort_pdeep_phosSight columns:", paste(colnames(feature_phos_autort_pdeep_phosSight), collapse=", "), "\n")
# cat("feature_phos_autort_pdeep_phosSight dimensions:", dim(feature_phos_autort_pdeep_phosSight), "\n")

feature_phos_autort_pdeep_phosSight_combined <- left_join(phosphors_autort_pDeep_phosSight2, feature_phos_autort_pdeep_phosSight, by="Title")
feature_autort_pdeep_phosSight_combined <- merge(feature_nonphos, feature_phos_autort_pdeep_phosSight_combined, all = TRUE)

# cat("=== DEBUG: feature_phos_autort_pdeep_phosSight_combined columns ===\n")
# cat("feature_phos_autort_pdeep_phosSight_combined columns:", paste(colnames(feature_phos_autort_pdeep_phosSight_combined), collapse=", "), "\n")
# cat("feature_phos_autort_pdeep_phosSight_combined dimensions:", dim(feature_phos_autort_pdeep_phosSight_combined), "\n")

feature_autort_pdeep_phosSight_combined$AutoRTProb <- ifelse(is.na(feature_autort_pdeep_phosSight_combined$AutoRTProb), 0, feature_autort_pdeep_phosSight_combined$AutoRTProb)
feature_autort_pdeep_phosSight_combined$pDeepProb <- ifelse(is.na(feature_autort_pdeep_phosSight_combined$pDeepProb), 0, feature_autort_pdeep_phosSight_combined$pDeepProb)
feature_autort_pdeep_phosSight_combined$PhosSightProb <- ifelse(is.na(feature_autort_pdeep_phosSight_combined$PhosSightProb), 0, feature_autort_pdeep_phosSight_combined$PhosSightProb)
feature_autort_pdeep_phosSight_combined$autort_pDeep_Prob <- ifelse(is.na(feature_autort_pdeep_phosSight_combined$autort_pDeep_Prob), 0, feature_autort_pdeep_phosSight_combined$autort_pDeep_Prob)
feature_autort_pdeep_phosSight_combined$autort_pDeep_phosSight_Prob <- ifelse(is.na(feature_autort_pdeep_phosSight_combined$autort_pDeep_phosSight_Prob), 0, feature_autort_pdeep_phosSight_combined$autort_pDeep_phosSight_Prob)

feature_autort_pdeep_phosSight_combined$Mod_Sequence_for_phosphoRS2 <- feature_autort_pdeep_phosSight_combined$Mod_Sequence_for_phosphoRS
for (i in 1:length(VariableModsInfo)){
  tmp = VariableModsInfo[i]
  tmp2 = unlist(strsplit(tmp,',',fixed = TRUE))
  number = tmp2[1]
  name = tmp2[2]
  aa = tmp2[3]
  mass = tmp2[4]
  sym = tmp2[5]
  if (name!='Phospho'){
        feature_autort_pdeep_phosSight_combined$Mod_Sequence_for_phosphoRS2 <- gsub(paste0(aa,sym),number,feature_autort_pdeep_phosSight_combined$Mod_Sequence_for_phosphoRS2)
  }
}
if (FixedMods!='null'){
  for (i in 1:length(FixedModsInfo)){
    tmp = FixedModsInfo[i]
    tmp2 = unlist(strsplit(tmp,',',fixed = TRUE))
    number = tmp2[1]
    name = tmp2[2]
    aa = tmp2[3]
    mass = tmp2[4]
    sym = tmp2[5]
    
    if (aa == 'AnyN-term'){
      next
    }
    
    if (name!='Phospho'){
      feature_autort_pdeep_phosSight_combined$Mod_Sequence_for_phosphoRS2 <- gsub(paste0(aa,sym),number,feature_autort_pdeep_phosSight_combined$Mod_Sequence_for_phosphoRS2)
    }
  }
}
#feature_autort_pdeep_phosSight_combined$Mod_Sequence_for_phosphoRS2 <- str_replace_all(feature_autort_pdeep_phosSight_combined$Mod_Sequence_for_phosphoRS,c("M1" = "1", "C3" = "5"))
  
feature_autort_pdeep_phosSight_combined$IsoformSequence_PhosphoRS <- ifelse(is.na(feature_autort_pdeep_phosSight_combined$IsoformSequence_PhosphoRS), feature_autort_pdeep_phosSight_combined$Mod_Sequence_for_phosphoRS2, feature_autort_pdeep_phosSight_combined$IsoformSequence_PhosphoRS)
feature_autort_pdeep_phosSight_combined$IsoformSequence_AutoRT <- ifelse(is.na(feature_autort_pdeep_phosSight_combined$IsoformSequence_AutoRT), feature_autort_pdeep_phosSight_combined$Mod_Sequence_for_phosphoRS2, feature_autort_pdeep_phosSight_combined$IsoformSequence_AutoRT)
feature_autort_pdeep_phosSight_combined$IsoformSequence_pDeep <- ifelse(is.na(feature_autort_pdeep_phosSight_combined$IsoformSequence_pDeep), feature_autort_pdeep_phosSight_combined$Mod_Sequence_for_phosphoRS2, feature_autort_pdeep_phosSight_combined$IsoformSequence_pDeep)
feature_autort_pdeep_phosSight_combined$IsoformSequence_PhosSight <- ifelse(is.na(feature_autort_pdeep_phosSight_combined$IsoformSequence_PhosSight), feature_autort_pdeep_phosSight_combined$Mod_Sequence_for_phosphoRS2, feature_autort_pdeep_phosSight_combined$IsoformSequence_PhosSight)
feature_autort_pdeep_phosSight_combined$IsoformSequence_autort_pDeep <- ifelse(is.na(feature_autort_pdeep_phosSight_combined$IsoformSequence_autort_pDeep), feature_autort_pdeep_phosSight_combined$Mod_Sequence_for_phosphoRS2, feature_autort_pdeep_phosSight_combined$IsoformSequence_autort_pDeep)
feature_autort_pdeep_phosSight_combined$IsoformSequence_autort_pDeep_phosSight <- ifelse(is.na(feature_autort_pdeep_phosSight_combined$IsoformSequence_autort_pDeep_phosSight), feature_autort_pdeep_phosSight_combined$Mod_Sequence_for_phosphoRS2, feature_autort_pdeep_phosSight_combined$IsoformSequence_autort_pDeep_phosSight)
feature_autort_pdeep_phosSight_combined$Mod_Sequence_for_phosphoRS2 <- NULL
write.table(feature_autort_pdeep_phosSight_combined, output, row.names=F, quote=F, sep="\t")