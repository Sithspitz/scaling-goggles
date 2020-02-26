### Loop Troubleshooting Two, Fixing a Hole ###
# Where the rain gets in: aka making sure that the log2FC values work #

mypackages <- c("GSEABase", "GSVA", "Biobase","RColorBrewer", "scales", "dplyr", "gtools", "stats")

lapply(mypackages, library, character.only = T)
source("C:/Users/rbuch/Local Documents/PhD/Local_R_Repos/functions/functions.R")

setwd("C:/Users/rbuch/DataShare/Volcano Plot HALLMARK/Troubleshooting Metabolic GSEA/Loop Troubleshooting/Fixing Loop/")

# Pre-processing

df <- read.csv("fixing_csv_test_one.csv", row.names = 1, header = T)
df2 <- as.matrix(df)
data_names <- df

mt_position_notes <- c("MT are from col 1:10, WT are from col 11:20")

# Detecting if both median group comparison values are negative numbers
## Change the comparison columns

neg_value_detector <- vector()
for (i in 1:nrow(df2)) {
  print(rownames(data_names[i, ]))
  neg_value_detector[length(neg_value_detector)+1] <- rownames(data_names[i, ])
  if (median(df2[i, 1:10]) < 0 && median(df2[i, 11:20]) < 0) {
    neg_value_detector[length(neg_value_detector)+1] <- print("FLAG")
  } else {
    neg_value_detector[length(neg_value_detector)+1] <- print("ok")
  }
}

neg_value_detector_df <- as.data.frame(neg_value_detector)
setwd("C:/Users/rbuch/DataShare/Volcano Plot HALLMARK/Troubleshooting Metabolic GSEA/Loop Troubleshooting/Fixing Loop/")
write.csv(neg_value_detector_df, "neg_value_detector_test.csv", row.names = F)

# Willcox tests, p value, q value and fold change calculation
## Change the comparison columns and FDR number to the number of GS being investigated

test_output <- vector()
for (i in 1:nrow(df2)) {
  print(rownames(data_names[i, ]))
  wil_tested <- wilcox.test(df2[i, 1:10], df2[i, 11:20])
  q_val <- as.numeric(wil_tested$p.value)
  q_val <- p.adjust(q_val, method = "fdr", n = 4)
  test_output[length(test_output)+1] <- rownames(data_names[i, ])
  test_output[length(test_output)+1] <- wil_tested$p.value
  test_output[length(test_output)+1] <- q_val
  if (median(df2[i, 1:10]) > median(df2[i, 11:20])) {
    test_output[length(test_output)+1] <- print ("Up")
  } else {
    test_output[length(test_output)+1] <- print ("Down")
  }
  f_change <- foldchange(median(df2[i, 1:10]), median(df2[i, 11:20]))
  log_fold <- foldchange2logratio(f_change, base = 2)  
  if (median(df2[i, 1:10]) < 0 && median(df2[i, 11:20]) < 0 && log_fold > 0)  {
  log_fold <- log_fold*-1
  test_output[length(test_output)+1] <- log_fold
  log_fold <- 1
  }
  if (median(df2[i, 1:10]) < 0 && median(df2[i, 11:20]) < 0 && log_fold < 0) {
    log_fold <- abs(log_fold)
    test_output[length(test_output)+1] <- log_fold
  }
  if (median(df2[i, 1:10]) > 0 | median(df2[i, 11:20]) > 0) {
    test_output[length(test_output)+1] <- log_fold
  }
}

test_output_df <- as.data.frame(test_output)
setwd("C:/Users/rbuch/DataShare/Volcano Plot HALLMARK/Troubleshooting Metabolic GSEA/Loop Troubleshooting/Fixing Loop/")
write.csv(test_output_df, "test_output.csv", row.names = F)

