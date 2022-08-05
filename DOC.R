# filter and format SEER DOC 
load('./IR.SEER.RData')
SEER <- IR.SEER$IR.SEER
SEER1 <- SEER[which(SEER$cancer == "Breast" & SEER$Age.Interval != "0 <= Age < 1"), 
              c("Age.Interval","gender","Prob..Of.Dying.From.Other.Causes.In.Interval")]
colnames(SEER1) <- c("age","sex","deathPen")
SEER1$age <- rep(1:95, 2)
SEER1$sex <- rep(c("F","M"), each = 95)
saveRDS(SEER1, file = "../DOC.rds")
