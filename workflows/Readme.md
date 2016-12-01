# Cancer Genomics Workflows
These are some example workflows that will get you started using the Cancer Genomics Toolkit
For details, the full toolkit is described [here](http://biorxiv.org/content/early/2016/11/26/089631)

###Calling somatic SNVs using an ensemble approach
![ScreenShot](ensemble_caller_workflow.png)
This workflow demonstrates our implementation of a simple voting-based ensemble variant calling approach using Galaxy. Here, four somatic SNV callers are run on the same tumour/normal pair (RADIA, SomaticSniper, MutationSeq and Strelka). The vcf-formatted output of each tool is passed to ensemble_vcf, which in turn outputs variants detected by a user-specified proportion of these callers (e.g. >50%). The last step in this workflow is to annotate the remaining variants using vcf2maf. 

###Taking a batch of screenshots in IGV for variants affecting specific loci
![ScreenShot](igv_screenshot.png)

###Using OncodriveFM to detect significantly mutated genes and visualizations
![ScreenShot](oncodrivefm_gene_discovery.png)

###Running GISTIC on CNV results and generating Oncocircos visualization
![ScreenShot](gistic_snv_workflow.png)

###Running OncodriveClust and Visualizing Mutations in Lollipop Plots
![ScreenShot](oncodriveclust_workflow.png)
