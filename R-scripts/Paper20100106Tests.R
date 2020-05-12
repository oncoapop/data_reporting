### Table 4

tdf <- read.csv("Table4TUNELPosNucleiData.csv")

tdf$MouseIDf <- as.factor(tdf$MouseID)

mf <- glm(cbind(NucPos, NucNeg) ~ Genotype + MouseIDf, data = tdf, family = binomial())
## Obs 23 outlier

mf2 <- glm(cbind(NucPos, NucNeg) ~ Genotype + MouseIDf, data = tdf, subset = NucPos != 84, family = binomial())
## QQ plot reasonable

mr2 <- glm(cbind(NucPos, NucNeg) ~ MouseIDf, data = tdf, subset = NucPos != 84, family = binomial())
mrr2 <- glm(cbind(NucPos, NucNeg) ~ 1, data = tdf, subset = (NucPos != 84), family = binomial())

anova(mr2, mf2)

lmef1 <- glmer(cbind(NucPos, NucNeg) ~ Genotype + (1|MouseIDf), data = tdf, family = binomial, REML = FALSE)
lmef2 <- glmer(cbind(NucPos, NucNeg) ~ Genotype + (1|MouseIDf), data = tdf, subset = (NucPos != 84), family = binomial, REML = FALSE)
lmer1 <- glmer(cbind(NucPos, NucNeg) ~ 1 + (1|MouseIDf), data = tdf, family = binomial, REML = FALSE)
lmer2 <- glmer(cbind(NucPos, NucNeg) ~ 1 + (1|MouseIDf), data = tdf, subset = (NucPos != 84), family = binomial, REML = FALSE)

anova(lmer1, lmef1)
anova(lmer2, lmef2)

tdf$PropPos <- tdf$NucPos / tdf$NucleiTotal
tdf$asp <- asin(sqrt(tdf$PropPos))
boxplot(PropPos ~ Genotype, data = tdf)
boxplot(asp ~ Genotype, data = tdf)
t.test(asp ~ Genotype, data = tdf)
t.test(PropPos ~ Genotype, data = tdf)
t.test(asp ~ Genotype, data = tdf, subset = NucPos != 84)
t.test(PropPos ~ Genotype, data = tdf, subset = NucPos != 84)

maovf1 <- aov(asp ~ Genotype + Error(MouseIDf), data = tdf)
maovf1
summary(maovf1)


twdf <- reshape(tdf, v.names = c("NucleiTotal", "NucPos", "NucNeg"), timevar = "Rep",
                idvar =  "MouseID", drop = c("PropPos", "asp"), direction = "wide")
twdf$NucNeg <- twdf$NucNeg.1 + twdf$NucNeg.2 + twdf$NucNeg.3 + twdf$NucNeg.4
twdf$NucPos <- twdf$NucPos.1 + twdf$NucPos.2 + twdf$NucPos.3 + twdf$NucPos.4
twdf$NucleiTotal <- twdf$NucleiTotal.1 + twdf$NucleiTotal.2 + twdf$NucleiTotal.3 + twdf$NucleiTotal.4
twdf
twdf$PropPos <- twdf$NucPos/twdf$NucleiTotal
twdf$asp <- asin(sqrt(twdf$PropPos))
t.test(twdf$asp[twdf$Genotype == "WT"], twdf$asp[twdf$Genotype == "KO"])

###------------------------------------------------------------------------------------
### qPCR


## Need a copy of SM utility functions
require("MOutils")


if ( !file.exists("./Plots") ) dir.create("./Plots")
require(lattice)

### qPCR test plate data
qtdf <-
  read.table("./qPCR_Data/DYTestis10Jul08.sdm-Result_Data.csv",
###  read.table("./qPCR_Data/DYTestis16Jul2008.sdm-Result_Data.csv",
             header = TRUE, row.names = NULL,
             na.strings = c("Undetermined", ""),
             stringsAsFactors = FALSE, sep = ",")


genodf <- data.frame(Sample = c("T1", "T2", "T3", "T4", "T5", "T6"),
                     SampleID = c("T1", "T2", "T3", "T4", "T5", "T6"),
                     Genotype = c("WT", "WT", "WT", "KO", "KO", "KO"),
             stringsAsFactors = FALSE)

### Double check that sample names in qtdf match sample names in genodf.
if ( !all.equal(sort(unique(qtdf$Sample)),  sort(unique(genodf$Sample))) ) {
  msg <- paste("Sample variable values in qtdf do not match Sample variable values in genodf\n")
  stop(msg)
}
qtdf <- merge(qtdf, genodf, all = TRUE)

### qtdf$Well <- as.numeric(qtdf$WellNo)
qtdf$PlateID <- "Plate1"
qtdf$GeneID <- qtdf$Detector
qtdf$SubjID <- qtdf$Sample
qtdf$RepID <- with(qtdf, paste(PlateID, SubjID, GeneID, sep = "|"))

qtdf$sGenotype <- qtdf$Genotype

qtdf$repn <- rep(NA_integer_, nrow(qtdf))
for ( repi in unique(qtdf$RepID) ) {
  posidxp <- qtdf$RepID %in% repi
  qtdf$repn[posidxp] <- seq(sum(posidxp))
}  
qtdf$repc <- paste("Rep", as.character(qtdf$repn))
## need unique id for each obs.
qtdf$id <- with(qtdf, paste(PlateID, Sample, repc, sep = "|"))


### ### Well position info - 384 well plate
### qtdf$Row384n <- (((qtdf$Well - 1) %/% 24) + 1)
### qtdf$Row384c <- LETTERS[qtdf$Row384n]
### qtdf$Col384n <- (((qtdf$Well - 1) %% 24) + 1)
### qtdf$Col384c <- leftPad0(qtdf$Col384n, 2)
### qtdf$Well384c <- paste(qtdf$Row384c, qtdf$Col384c, sep = "")

### More RQ Mangler funk:
### 2010 Jan 07:  NOTE:  Well384c does NOT match "Pos"
### from RQ manager.  RQ Manager "Pos" is just a sequential
### count of wells, no matter what they are.  The above scheme
### only works if all wells on the plate are used.
### Thus the only way to get proper well position is to make
### sure RQ Mangler puts in the well, via "Plate-Centric" export.
###

qtdf$Row384c <- substring(qtdf$Pos, 1, 1)
qtdf$Col384c <- leftPad0(substring(qtdf$Pos, 2), 2)
qtdf$Well384c <- paste(qtdf$Row384c, qtdf$Col384c, sep = "")
qtdf$Row384n <- match(toupper(qtdf$Row384c), LETTERS)
qtdf$Col384n <- as.numeric(qtdf$Col384c)

qtalldf <- qtdf

weakDetectionidxs <- which(qtalldf$Ct >= 36.0)
qtdf <- qtalldf
is.na(qtdf[weakDetectionidxs, "Ct"]) <- TRUE
is.na(qtdf[weakDetectionidxs, "Ctmed"]) <- TRUE


### Calculate median Ct value in each set of technical replicates
CtMedVals <- with(qtdf, tapply(Ct, RepID, function(x) { median(x) } , simplify = FALSE) )
midxs <- match(qtdf$RepID, names(CtMedVals))
### all.equal(names(CtMedVals)[midxs], qtdf$RepID) ## TRUE
qtdf$Ctmed <- unlist(CtMedVals[midxs])

qpcrplot <- function(qdf, gtarget, gcontrol, panel.data = NULL)
{
  require("lattice")
  dat <- merge(qdf[toupper(qdf$GeneID) == toupper(gtarget),],
               qdf[toupper(qdf$GeneID) == toupper(gcontrol),
                   c("Ct", "Ctmed", "id")], by = 'id')

  dat$SubjIDf <- as.factor(dat$SubjID)
  nSubj <- length(unique(dat$SubjID))
  if ( nSubj <= 8 ) { keyncol <- 2 } else { keyncol <- 4 }
  try(xyplot((Ct.x - Ctmed.y) ~ Ct.y | Genotype, groups = SubjIDf, data = dat,
###              xlab = gcontrol,
###              ylab = paste(as.character(dat$GeneID[1]), " - median( ", gcontrol, " )", sep = ""),
             xlab = strsplit(gcontrol, split = "_UPL")[[1]][1],
             ylab = paste(strsplit(as.character(dat$GeneID[1]), split = "_")[[1]][1],
               " - median( ", strsplit(gcontrol, split = "_")[[1]][1], " )", sep = ""),
             main = "MLL5",
             layout = c(2, 1), ## No page count if using grid package layouts.
             pch = c(1:nSubj), col = c(1:nSubj), pwd = 2,
             panel = function(..., panel.dat = panel.data) {
               panel.xyplot(...)
               if (!is.null(panel.dat)) {
                 panel.abline(h = panel.dat[packet.number(), 1], lwd = 1, col = "blue")
                 panel.abline(h = panel.dat[packet.number(), 2], lwd = 2, col = "blue")
                 panel.abline(h = panel.dat[packet.number(), 3], lwd = 1, col = "blue")
               }
             },               
             key = list(text = list(levels(dat$SubjIDf)),
               points = list(pch = 1:nSubj, col = c(1:nSubj)), columns = keyncol)
             ))
}


GeneNames <- sort(unique(as.character(qtdf$GeneID)))
GeneNamesControl <- c("Eef1a1_UPL31")
GeneNamesExp <- setdiff(GeneNames, GeneNamesControl)



### Plot WT, KiSS KO, GPR54 KO    Testosterone vs None
### Same axis

pdf('./Plots/qPCRplotsMLL510Jul2008v02.pdf', width=10, height=8)
require("lattice")
require("grid")

for ( i in seq(((length(GeneNamesExp) - 1) %/% 4) + 1) ) {

  grid.newpage()

  pushViewport(viewport(layout = grid.layout(2, 2)))

  gni <- GeneNamesExp[(4*i) - 3]
  if ( !is.na(gni) ) {
    ## UL
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
    qp1 <- qpcrplot(qdf = qtdf, gtarget = gni, gcontrol = GeneNamesControl[1])
    if ( !(class(qp1) == "try-error") ) {
      print(qp1, newpage = FALSE)
    } else {
      msg <- paste("Plot not available:", gni, GeneNamesControl[1])
      grid.text(msg)
    }
    popViewport()
  }
  
  gni <- GeneNamesExp[(4*i) - 2]
  if ( !is.na(gni) ) {
    ## UR
    pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 1))
    qp1 <- qpcrplot(qdf = qtdf, gtarget = gni, gcontrol = GeneNamesControl[1])
    if ( !(class(qp1) == "try-error") ) {
      print(qp1, newpage = FALSE)
    } else {
      msg <- paste("Plot not available:", gni, GeneNamesControl[1])
      grid.text(msg)
    }
    popViewport()
  }

  gni <- GeneNamesExp[(4*i) - 1]
  if ( !is.na(gni) ) {
    ## LL
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2))
    qp1 <- qpcrplot(qdf = qtdf, gtarget = gni, gcontrol = GeneNamesControl[1])
    if ( !(class(qp1) == "try-error") ) {
      print(qp1, newpage = FALSE)
    } else {
      msg <- paste("Plot not available:", gni, GeneNamesControl[1])
      grid.text(msg)
    }
    popViewport()
  }

  gni <- GeneNamesExp[(4*i) - 0]
  if ( !is.na(gni) ) {
    ## LR
    pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 2))
    qp1 <- qpcrplot(qdf = qtdf, gtarget = gni, gcontrol = GeneNamesControl[1])
    if ( !(class(qp1) == "try-error") ) {
      print(qp1, newpage = FALSE)
    } else {
      msg <- paste("Plot not available:", gni, GeneNamesControl[1])
      grid.text(msg)
    }
    popViewport()
  }

}
dev.off()



### The problem is that the mean doesn't act like this: technically, you need
### to multiply it by the Jacobian of the transformation. Fortunately, this
### is a standard result, so the mean of e^x is mu + 0.5*s^2, where
### mu and s^2 are the mean and variance on the un-transformed scale. 
qPCRModel <- function(qdf, gtarget, gcontrol, genotype = "KO", cref = "WT", alphalevel = 0.05)
{

  require(lme4)
  
  plotPanelgroupNames <- c("WT",  "KO")
  plotPanelestNames <- c("FCest", "95% LCL", "95% UCL")
  groupNames <- c("KO/WT")
  estNames <- c("FCest", "FC 95% LCL", "FC 95% UCL",
                 "Abs FCest", "Abs FC Direction", "Abs FC 95% LCL", "Abs FC 95% UCL")
  qpcrout.names <- c("GeneTarget", "Loading Control", "FC Description",
                     "Technical Replicates Fligner-Killeen:med chi-squared", "TR F-K:df", "TR F-K:p-value",
                     "Biological Replicates Fligner-Killeen:med chi-squared", "BR F-K:df", "BR F-K:p-value",
                     "qPCR GeneTarget Chisq", "Chi Df", "Pr(>Chisq)", "BenHoch Adj Pvalue",
                     paste(rep(groupNames, each = length(estNames)), estNames),
                     paste(rep(plotPanelgroupNames, each = length(plotPanelestNames)), plotPanelestNames)
                     )
  qpcrout <- rep(NA_character_, length(qpcrout.names))
  qpcrout[1:3] <- c(gtarget, gcontrol, paste(genotype, "relative to", cref))
  names(qpcrout) <- qpcrout.names

  ## Model gtarget adjusting for gcontrol
  ## Get subset of qPCR data
  qsdf <- qdf[((toupper(qdf$Detector) == toupper(gtarget) |
                toupper(qdf$Detector) == toupper(gcontrol)) &
               ((qdf$Genotype == genotype) |
                (qdf$Genotype == cref))), ]

  qsdf$Genotype <- relevel(as.factor(qsdf$Genotype), ref = cref)
  qsdf$Detector <- relevel(as.factor(toupper(qsdf$Detector)), ref = toupper(gcontrol))

  msg <- paste("\n\n###--- ", genotype, ":  Modeling ", gtarget, " for loading control ", gcontrol, "   ---###\n\n", sep = "") #
  cat(msg)
  
  ## Want loading controls only from plates containing gene target
  gPlates <- unique(qsdf$PlateID[toupper(qsdf$Detector) == toupper(gtarget)])
  nPlates <- length(gPlates)
  qsdf <- qsdf[qsdf$PlateID %in% gPlates, ]
  if (length(table(as.character(qsdf$Detector))) <= 1) {
    ## Missing loading control or other data
    msg <- paste("\nNOTE: Insufficient data to model ",
                 gtarget, " for ", gcontrol, "\n\n",
                 sep = "")
    warning(msg)
    print(table(qsdf$Detector))
    return(qpcrout)
  }
  qsdf$RepIDf <- as.factor(paste(qsdf$RepID, qsdf$Genotype, sep = "|"))
  qsdf$SubjID <- as.factor(qsdf$Sample)

  ## full model with interaction
  fim.reml <- try(lmer(Ct ~ Detector * Genotype + (1|SubjID), data = qsdf, REML = TRUE))

  if ( class(fim.reml) == "try-error" ) {
    msg <- paste("\nNOTE: fim.reml fit failed: ",
                 gtarget, " for ", gcontrol, "\n\n",
                 sep = "")
    warning(msg)
    return(qpcrout)
  }
  
  if (!all(sapply(coef(fim.reml)[[1]], function(x) length(unique(x)) == 1)[-1])) {
    warning("NOTE:  Not constant estimates for all cases")
    print(coef(fim.reml)[[1]])
  }
  ## Estimates
  ## WT T- :
  ones <- function(which, length) {
    one <- rep(0, length)
    one[which] <- 1
    one
  }
  Estsmat <- cbind(ones(c(2), 4), ## WT dCt
                   ones(c(2, 4), 4), ## KO dCt
                   ones(c(4), 4) ## (KO FC)/(WT FC) ddCt
                   )
  FCestidxs <- c(3)
  EstslFC <- as.vector(t(Estsmat) %*% unlist(coef(fim.reml)[[1]][1, ]))
  StdDevslFC <- sqrt(diag(t(Estsmat) %*% vcov(fim.reml) %*% Estsmat))
  CImult <- qnorm(1.0 - alphalevel/2.0)
  LCLlFC <- as.vector(EstslFC - CImult * StdDevslFC)
  UCLlFC <- as.vector(EstslFC + CImult * StdDevslFC)
  EstsFC <- 2^(-1.0 * EstslFC[FCestidxs])
  LCLFC <- 2^(-1.0 * UCLlFC[FCestidxs])
  UCLFC <- 2^(-1.0 * LCLlFC[FCestidxs])
  absEstsFC <- ifelse(EstsFC < 1.0, 1.0/EstsFC, EstsFC)
  absDirectionFC <- ifelse(EstsFC < 1.0, "Down", "Up")
  absLCLFC <- ifelse(EstsFC < 1.0, 1.0/UCLFC, LCLFC)
  absUCLFC <- ifelse(EstsFC < 1.0, 1.0/LCLFC, UCLFC)

  ## full model with interaction
  fim.ml <- try(lmer(Ct ~ Detector * Genotype + (1|SubjID), data = qsdf, REML = FALSE)) ## method = "ML"))

  if ( class(fim.ml) == "try-error" ) {
    msg <- paste("\nNOTE: fim.ml fit failed: ",
                 gtarget, " for ", gcontrol, "\n\n",
                 sep = "")
    warning(msg)
    return(qpcrout)
  } 
  ## full additive model, no interaction
  fam.ml <- lmer(Ct ~ Detector + Genotype + (1|SubjID), data = qsdf, REML = FALSE) ## method = "ML")    

  if ( class(fam.ml) == "try-error" ) {
    msg <- paste("\nNOTE: fam.ml fit failed: ",
                 gtarget, " for ", gcontrol, "\n\n",
                 sep = "")
    warning(msg)
    return(qpcrout)
  } 
  ## Reduced model, no difference between genotypes.
  rrm.ml <- lmer(Ct ~ Detector + (1|SubjID), data = qsdf, REML = FALSE) ## method = "ML")

  if ( class(rrm.ml) == "try-error" ) {
    msg <- paste("\nNOTE: rrm.ml fit failed: ",
                 gtarget, " for ", gcontrol, "\n\n",
                 sep = "")
    warning(msg)
    return(qpcrout)
  } 

  ## Interaction anova:  Any difference between genotypes in the gene target
  gtAnova <- anova(fim.ml, fam.ml)
  print(gtAnova)
  ## full model summary
  fim.ml.su <- summary(fim.ml)
  print(fim.ml.su)
###   fim.ml.su.coefs <- fim.ml.su@coefs[nrow(fim.ml.su@coefs), ]
  fim.reml.su <- summary(fim.reml)
  print(fim.reml.su)
###   fim.reml.su.coefs <- fim.reml.su@coefs[nrow(fim.reml.su@coefs), c("Estimate", "Std. Error")]
  ## Row 2 has main effect estimate - save to use in plots
  fim.reml.maineff.su.coefs <- fim.reml.su@coefs[2, c("Estimate", "Std. Error")]
  ## Need to test for constant variance.
  ## fligner.test
  ## Need to test technical replicates for non-constant variance.
  ## If technical replicates okay, test means (i.e. test biological
  ##   rep variance) across genotype groups.
  ## technical rep constant variance:
  trfktestout <- try(fligner.test(Ct ~ RepIDf, data = qsdf, subset = !(toupper(Detector) == toupper(gcontrol))))
  if ( class(trfktestout) == "try-error" )  {
    trfktestout <- list(statistic = "failed", parameter = NA, p.value = NA)
  }
  qsadf <- aggregate(qsdf[toupper(qsdf$Detector) != toupper(gcontrol), ]$Ct,
                     by = list(qsdf[toupper(qsdf$Detector) != toupper(gcontrol), ]$RepIDf), FUN = mean)
  names(qsadf) <- c("RepID", "meanCt")
  qsadf$genotype <- NA_character_
### Need to fix this genotype code as RepID differs here
  qsadf$genotype[grep("KO", qsadf$RepID)] <- "KO"
  qsadf$genotype[grep("WT", qsadf$RepID)] <- "WT"
  qsadf$genotypef <- as.factor(qsadf$genotype)
  ## Biological rep var test
  brfktestout <- try(fligner.test(meanCt ~ genotypef, data = qsadf))
  if ( class(brfktestout) == "try-error" )  {
    brfktestout <- list(statistic = "failed", parameter = NA, p.value = NA)
  }
###   ## exp scale mean ratio estimate and 95% Wald CI
###   names(exp.scale) <- c("FCest", "95% LCL", "95% UCL")
  ## Need also GKO and KKO estimates
  ## Genotype relative to WT
###   genotype.est <- as.vector(rbind(EstsFC, LCLFC, UCLFC, absEstsFC, absDirectionFC, absLCLFC, absUCLFC))
  genotype.T.Xn.ests <- as.vector(rbind(EstsFC, LCLFC, UCLFC, absEstsFC, absDirectionFC, absLCLFC, absUCLFC))
  names(genotype.T.Xn.ests) <- paste(rep(groupNames, each = length(estNames)), estNames)
  plotMeansCIs <- as.vector(rbind(EstslFC[-FCestidxs], LCLlFC[-FCestidxs], UCLlFC[-FCestidxs]))

  qpcrout <- c(gtarget, gcontrol, paste(genotype, "relative to", cref),
               trfktestout$statistic, trfktestout$parameter, trfktestout$p.value,
               brfktestout$statistic, brfktestout$parameter, brfktestout$p.value,
               gtAnova[2, "Chisq"], gtAnova[2, "Chi Df" ], gtAnova[2, "Pr(>Chisq)"], "NAforBHhere",
               genotype.T.Xn.ests, plotMeansCIs
               )
  names(qpcrout) <- qpcrout.names
  print(qpcrout)
  qpcrout
}


### Modeling
### i <- 1; j <- 1; k <- 1


modelResults <- NULL

Genotypes <- c("KO")
date()
for ( k in seq(along = Genotypes) ) {
  for ( j in seq(along = GeneNamesControl) ) {
    for ( i in seq(along = GeneNamesExp) ) {
      fit <- qPCRModel(qdf = qtdf, gtarget = GeneNamesExp[i], gcontrol = GeneNamesControl[j], genotype = Genotypes[k])
      modelResults <- rbind(modelResults, fit)
    }
  }
}
date()
dimnames(modelResults) <- list(NULL, dimnames(modelResults)[[2]])
### Adjust p-values for multiple comparisons
### modelResults <-  cbind(modelResults[, 1:12], rep("", nrow(modelResults)), modelResults[, 13:15])
### dimnames(modelResults)[[2]][13] <- "BenHochAdjPvalue"
### GKO
idxp <- modelResults[, "FC Description"] == "KO relative to WT"
modelResults[idxp, "BenHoch Adj Pvalue"] <-
  as.character(p.adjust(as.numeric(modelResults[idxp, "Pr(>Chisq)"]), method = "BH"))

write.csv(modelResults, file = "qPCR.MLL5.16Jul2008.modelResults.csv")
