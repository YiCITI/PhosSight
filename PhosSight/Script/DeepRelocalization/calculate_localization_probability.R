library(tidyverse)
library(data.table)

args <- commandArgs(T)
phosphors <- fread(args[1])
auto_rt <- fread(args[2])
pDeep <- fread(args[3])
features <- fread(args[4])
PhosSight <- fread(args[5])
output <- args[6]
VariableMods <- args[7]
FixedMods <- args[8]

if (FALSE){
phosphors <- fread('PhosphoRS.txt')
auto_rt <- fread('phospho.prediction.tsv')
pDeep <- fread('pDeep3PredictionResults.Phospho.txt')
features <- fread('features.PhosphoRS.txt')
PhosSight <- fread('predictions.csv')
output <- 'Features.Localization.entropy.PhosSight.txt'
VariableMods = '1,Oxidation,M,15.994919,1;2,Phospho,S,79.966331,2;3,Phospho,T,79.966331,2;4,Phospho,Y,79.966331,2' # Identification modifications
FixedMods = 'null' # Identification modifications
}

VariableModsInfo = unlist(strsplit(VariableMods,';',fixed = TRUE))
if (FixedMods!='null'){
  FixedModsInfo = unlist(strsplit(FixedMods,';',fixed = TRUE))
}

colnames(phosphors)[2] <- 'Title'
features$Link <- NULL
phosphors$Link <- NULL

# AutoRT
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

# pDeep3
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


# AutoRT_pDeep3
phosphors_autort_pDeep <- left_join(phosphors_autort, phosphors_pDeep, by = c("Title" = "Title", "IsoformSequence_AutoRT" = "IsoformSequence_pDeep", "PhosphoRSScore" = "PhosphoRSScore")) %>% select(Title, rtRatio, MaxRatio, entropy, Max_entropy, PhosphoRSScore, IsoformSequence_AutoRT)
phosphors_autort_pDeep$CombinedScore_autort_pDeep <- phosphors_autort_pDeep$PhosphoRSScore * (phosphors_autort_pDeep$rtRatio/phosphors_autort_pDeep$MaxRatio) * (phosphors_autort_pDeep$entropy/phosphors_autort_pDeep$Max_entropy)
phosphors_autort_pDeep$autort_pDeep_Prob1 <- 10^(phosphors_autort_pDeep$CombinedScore_autort_pDeep/10)

tmp <- phosphors_autort_pDeep[ ,list(autort_pDeep_Prob2=sum(autort_pDeep_Prob1)), by=Title]
phosphors_autort_pDeep <- left_join(phosphors_autort_pDeep,tmp, by = "Title")

phosphors_autort_pDeep$autort_pDeep_Prob <- phosphors_autort_pDeep$autort_pDeep_Prob1/phosphors_autort_pDeep$autort_pDeep_Prob2
phosphors_autort_pDeep[is.na(phosphors_autort_pDeep$autort_pDeep_Prob),]$autort_pDeep_Prob <- 0
colnames(phosphors_autort_pDeep)[7] <- 'IsoformSequence_autort_pDeep'


# PhosSight

colnames(PhosSight)[1] <- 'IsoformSequence'
colnames(PhosSight)[2] <- 'detectability'
phosphors_PhosSight <- left_join(phosphors_autort_pDeep, PhosSight, by = c('IsoformSequence_autort_pDeep' = 'IsoformSequence')) %>% select(Title, rtRatio, MaxRatio, entropy, Max_entropy, detectability, PhosphoRSScore, IsoformSequence_autort_pDeep)

colnames(phosphors_PhosSight)[8] <- 'IsoformSequence'


phosphors_PhosSight2 <- phosphors_PhosSight %>%
  group_by(Title) %>%
  slice(which.max(detectability))  

colnames(phosphors_PhosSight2)[6] <- 'Maxdetectability'
# phosphors_PhosSight2$detectability <- NULL
phosphors_PhosSight2$rtRatio <- NULL
phosphors_PhosSight2$MaxRatio <- NULL
phosphors_PhosSight2$entropy <- NULL
phosphors_PhosSight2$Max_entropy <- NULL
phosphors_PhosSight2$PhosphoRSScore <- NULL
phosphors_PhosSight2$IsoformSequence <- NULL

phosphors_PhosSight3 = merge(phosphors_PhosSight, phosphors_PhosSight2, all.x=TRUE, by='Title')

# cat('ok1\n')
# write.table(phosphors_PhosSight3, 'test.txt', row.names=F, quote=F, sep="\t")
# cat('ok2\n')


phosphors_PhosSight3$CombinedScorePhosSight <- phosphors_PhosSight3$PhosphoRSScore * (phosphors_PhosSight3$detectability/phosphors_PhosSight3$Maxdetectability)
phosphors_PhosSight3$CombinedScoreAll <- phosphors_PhosSight3$PhosphoRSScore * (phosphors_PhosSight3$rtRatio/phosphors_PhosSight3$MaxRatio) * (phosphors_PhosSight3$entropy/phosphors_PhosSight3$Max_entropy) * (phosphors_PhosSight3$detectability/phosphors_PhosSight3$Maxdetectability)

phosphors_PhosSight3$PhosSightProb1 <- 10^(phosphors_PhosSight3$CombinedScorePhosSight/10)
phosphors_PhosSight3$AllProb1 <- 10^(phosphors_PhosSight3$CombinedScoreAll/10)


setDT(phosphors_PhosSight3)



tmp1 <- phosphors_PhosSight3[ ,list(PhosSightProb2=sum(PhosSightProb1)), by=Title]
tmp2 <- phosphors_PhosSight3[ ,list(AllProb2=sum(AllProb1)), by=Title]



phosphors_PhosSight3 <- left_join(phosphors_PhosSight3,tmp1, by = "Title")
phosphors_PhosSight3 <- left_join(phosphors_PhosSight3,tmp2, by = "Title")

phosphors_PhosSight3$PhosSightProb <- phosphors_PhosSight3$PhosSightProb1/phosphors_PhosSight3$PhosSightProb2
phosphors_PhosSight3$AllProb <- phosphors_PhosSight3$AllProb1/phosphors_PhosSight3$AllProb2

phosphors_PhosSight3[is.na(phosphors_PhosSight3$PhosSightProb),]$PhosSightProb <- 0
phosphors_PhosSight3[is.na(phosphors_PhosSight3$AllProb),]$AllProb <- 0



phosphors_PhosSight = phosphors_PhosSight3 

all = phosphors_PhosSight3

# write.table(phosphors_autort_pDeep_PhosSight, 'test.txt', row.names=F, quote=F, sep="\t")







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

phosphors_autort_pDeep2 <- phosphors_autort_pDeep %>%
  group_by(Title) %>%
  slice(which.max(autort_pDeep_Prob))

phosphors_PhosSight2 <- phosphors_PhosSight %>%
  group_by(Title) %>%
  slice(which.max(PhosSightProb))

all2 <- all %>%
  group_by(Title) %>%
  slice(which.max(AllProb))

# write.table(phosphors_autort_pDeep2, 'test.txt', row.names=F, quote=F, sep="\t")

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

colnames(phosphors_PhosSight2)[8] <- 'IsoformSequence_PhosSight'
phosphors_PhosSight2$rtRatio <- NULL
phosphors_PhosSight2$MaxRatio <- NULL
phosphors_PhosSight2$entropy <- NULL
phosphors_PhosSight2$Max_entropy <- NULL
phosphors_PhosSight2$PhosphoRSScore <- NULL
phosphors_PhosSight2$CombinedScorePhosSight <- NULL
phosphors_PhosSight2$CombinedScoreAll <- NULL
phosphors_PhosSight2$PhosSightProb1 <- NULL
phosphors_PhosSight2$PhosSightProb2 <- NULL
phosphors_PhosSight2$AllProb1 <- NULL
phosphors_PhosSight2$AllProb2 <- NULL
phosphors_PhosSight2$detectability <- NULL
phosphors_PhosSight2$Maxdetectability <- NULL
phosphors_PhosSight2$AllProb <- NULL

# write.table(phosphors_PhosSight2, 'phosphors_PhosSight2.txt', row.names=F, quote=F, sep="\t")

colnames(all2)[8] <- 'IsoformSequence_All'
all2$rtRatio <- NULL
all2$MaxRatio <- NULL
all2$entropy <- NULL
all2$Max_entropy <- NULL
all2$PhosphoRSScore <- NULL
all2$CombinedScorePhosSight <- NULL
all2$CombinedScoreAll <- NULL
all2$PhosSightProb1 <- NULL
all2$PhosSightProb2 <- NULL
all2$AllProb1 <- NULL
all2$AllProb2 <- NULL
all2$detectability <- NULL
all2$Maxdetectability <- NULL
all2$PhosSightProb <- NULL

# write.table(all2, 'all2.txt', row.names=F, quote=F, sep="\t")

cat("OK!\n")

######################################################################################################################################################################################

feature_phos2 <- left_join(phosphors2, feature_phos, by="Title")
feature_phos_autort <- left_join(phosphors_autort2, feature_phos2, by="Title")
feature_phos_autort_pdeep <- left_join(phosphors_pDeep2, feature_phos_autort, by="Title")
feature_phos_autort_pdeep_combined <- left_join(phosphors_autort_pDeep2, feature_phos_autort_pdeep, by="Title")


# feature_autort_pdeep_combined <- merge(feature_nonphos, feature_phos_autort_pdeep_combined, all = TRUE)

# feature_autort_pdeep_combined$AutoRTProb <- ifelse(is.na(feature_autort_pdeep_combined$AutoRTProb), 0, feature_autort_pdeep_combined$AutoRTProb)
# feature_autort_pdeep_combined$pDeepProb <- ifelse(is.na(feature_autort_pdeep_combined$pDeepProb), 0, feature_autort_pdeep_combined$pDeepProb)
# feature_autort_pdeep_combined$autort_pDeep_Prob <- ifelse(is.na(feature_autort_pdeep_combined$autort_pDeep_Prob), 0, feature_autort_pdeep_combined$autort_pDeep_Prob)

feature_phosphors_pdeep_combined_PhosSight <- left_join(phosphors_PhosSight2, feature_phos_autort_pdeep_combined, by="Title")



feature_all <- left_join(all2, feature_phosphors_pdeep_combined_PhosSight, by="Title")



####################################################################################################################
feature_autort_pdeep_combined <- merge(feature_nonphos, feature_phos_autort_pdeep_combined, all = TRUE)
############################################################################################################
feature_all <- merge(feature_nonphos, feature_all, all = TRUE)
############################################################################################################

feature_all$AutoRTProb <- ifelse(is.na(feature_all$AutoRTProb), 0, feature_all$AutoRTProb)
feature_all$pDeepProb <- ifelse(is.na(feature_all$pDeepProb), 0, feature_all$pDeepProb)
feature_all$autort_pDeep_Prob <- ifelse(is.na(feature_all$autort_pDeep_Prob), 0, feature_all$autort_pDeep_Prob)
feature_all$PhosSightProb <- ifelse(is.na(feature_all$PhosSightProb), 0, feature_all$PhosSightProb)
feature_all$AllProb <- ifelse(is.na(feature_all$AllProb), 0, feature_all$AllProb)




#######################################################################################################################################################################################



feature_autort_pdeep_combined$Mod_Sequence_for_phosphoRS2 <- feature_autort_pdeep_combined$Mod_Sequence_for_phosphoRS
for (i in 1:length(VariableModsInfo)){
  tmp = VariableModsInfo[i]
  tmp2 = unlist(strsplit(tmp,',',fixed = TRUE))
  number = tmp2[1]
  name = tmp2[2]
  aa = tmp2[3]
  mass = tmp2[4]
  sym = tmp2[5]
  if (name!='Phospho'){
    feature_autort_pdeep_combined$Mod_Sequence_for_phosphoRS2 <- gsub(paste0(aa,sym),number,feature_autort_pdeep_combined$Mod_Sequence_for_phosphoRS2)
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
      feature_autort_pdeep_combined$Mod_Sequence_for_phosphoRS2 <- gsub(paste0(aa,sym),number,feature_autort_pdeep_combined$Mod_Sequence_for_phosphoRS2)
    }
  }
}


feature_all$Mod_Sequence_for_phosphoRS2 <- feature_all$Mod_Sequence_for_phosphoRS
for (i in 1:length(VariableModsInfo)) {
    tmp  <- VariableModsInfo[i]
    tmp2 <- unlist(strsplit(tmp, ',', fixed = TRUE))
    number <- tmp2[1]; name <- tmp2[2]; aa <- tmp2[3]; sym <- tmp2[5]

    if (name != 'Phospho') {
        feature_all$Mod_Sequence_for_phosphoRS2 <- 
            gsub(paste0(aa, sym), number, feature_all$Mod_Sequence_for_phosphoRS2)
    }
}
if (FixedMods != 'null') {
    for (i in 1:length(FixedModsInfo)) {
        tmp  <- FixedModsInfo[i]
        tmp2 <- unlist(strsplit(tmp, ',', fixed = TRUE))
        number <- tmp2[1]; name <- tmp2[2]; aa <- tmp2[3]; sym <- tmp2[5]

        if (aa != 'AnyN-term' && name != 'Phospho') {
            feature_all$Mod_Sequence_for_phosphoRS2 <- 
                gsub(paste0(aa, sym), number, feature_all$Mod_Sequence_for_phosphoRS2)
        }
    }
}



#feature_autort_pdeep_combined$Mod_Sequence_for_phosphoRS2 <- str_replace_all(feature_autort_pdeep_combined$Mod_Sequence_for_phosphoRS,c("M1" = "1", "C3" = "5"))
  
feature_autort_pdeep_combined$IsoformSequence_PhosphoRS <- ifelse(is.na(feature_autort_pdeep_combined$IsoformSequence_PhosphoRS), feature_autort_pdeep_combined$Mod_Sequence_for_phosphoRS2, feature_autort_pdeep_combined$IsoformSequence_PhosphoRS)
feature_autort_pdeep_combined$IsoformSequence_AutoRT <- ifelse(is.na(feature_autort_pdeep_combined$IsoformSequence_AutoRT), feature_autort_pdeep_combined$Mod_Sequence_for_phosphoRS2, feature_autort_pdeep_combined$IsoformSequence_AutoRT)
feature_autort_pdeep_combined$IsoformSequence_pDeep <- ifelse(is.na(feature_autort_pdeep_combined$IsoformSequence_pDeep), feature_autort_pdeep_combined$Mod_Sequence_for_phosphoRS2, feature_autort_pdeep_combined$IsoformSequence_pDeep)
feature_autort_pdeep_combined$IsoformSequence_autort_pDeep <- ifelse(is.na(feature_autort_pdeep_combined$IsoformSequence_autort_pDeep), feature_autort_pdeep_combined$Mod_Sequence_for_phosphoRS2, feature_autort_pdeep_combined$IsoformSequence_autort_pDeep)
feature_autort_pdeep_combined$Mod_Sequence_for_phosphoRS2 <- NULL
##################################################################################
feature_all$IsoformSequence_PhosphoRS <- ifelse(is.na(feature_all$IsoformSequence_PhosphoRS), feature_all$Mod_Sequence_for_phosphoRS2, feature_all$IsoformSequence_PhosphoRS)
feature_all$IsoformSequence_AutoRT <- ifelse(is.na(feature_all$IsoformSequence_AutoRT), feature_all$Mod_Sequence_for_phosphoRS2, feature_all$IsoformSequence_AutoRT)
feature_all$IsoformSequence_pDeep <- ifelse(is.na(feature_all$IsoformSequence_pDeep), feature_all$Mod_Sequence_for_phosphoRS2, feature_all$IsoformSequence_pDeep)
feature_all$IsoformSequence_autort_pDeep <- ifelse(is.na(feature_all$IsoformSequence_autort_pDeep), feature_all$Mod_Sequence_for_phosphoRS2, feature_all$IsoformSequence_autort_pDeep)
feature_all$IsoformSequence_PhosSight <- ifelse(is.na(feature_all$IsoformSequence_PhosSight), feature_all$Mod_Sequence_for_phosphoRS2, feature_all$IsoformSequence_PhosSight)
feature_all$IsoformSequence_All <- ifelse(is.na(feature_all$IsoformSequence_All), feature_all$Mod_Sequence_for_phosphoRS2, feature_all$IsoformSequence_All)
feature_all$Mod_Sequence_for_phosphoRS2 <- NULL

############################################################################################
write.table(feature_all, output, row.names=F, quote=F, sep="\t")




