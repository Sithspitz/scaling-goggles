# scaling-goggles
GSEA Master Scripts with integrated FDR and log2FC

# Background
When trying to make the metabolic analysis volcano plots, I found an issue with the direction of pathway regulation and log2FC. 

I.e. With the HALLMARK gene signatures, GS47, and GS49 were called as downregulated by the loop script, whilst the log2FC showed an upregulation. And GS50 was visa versa.

Discovered that the reason this was is because the SSGSEA values are negatives which throws up problems calculating log2 fold change.

# This Repo
After troubleshooting, the comprehensive yet generic scripts for calculating statistical difference between MT and WT SSGSEA scores (by Wilcoxon test) is housed here.

The R file 'fixing_loop.R' found in '~\GSEA Analysis Willcoxon Test Master Script\' HAS TO BE MODIFIED EACH TIME, but can calculate wilcoxon tested difference, p value, FDR adjusted q value and log2FC.

The biggest change from the previous version is that log2FC can be calculated from two groups of samples both with median values which are < 0. FDR adjusted q value is also integrated into this loop to negate further R scripts being needed.
The generic loop here has the output as 'test_output'

Also the 'neg_value_detector' output provides a CSV file with FLAGS if the median values of that GS comparison are both negative so you can go and check the data if appropriate.

