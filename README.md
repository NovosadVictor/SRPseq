# SRPseq

Splicing Regulation Prediction from RNA-seq data. Please cite this paper if you are using SRPseq in your work:

Nersisyan S, Novosad V, Tonevitsky A. SRPseq: Splicing Regulation Prediction from RNA-seq data. bioRxiv. 2021 Aug 04. doi: [10.1101/2021.08.03.454798](https://doi.org/10.1101/2021.08.03.454798).


## Table of Contents  
[Introduction](#introduction)  
[Installation](#installation)  
[Tutorial](#tutorial)  
[Running SRPseq](#running-srpseq)  

# Introduction
<img align="right" width="400px" src="https://github.com/NovosadVictor/SRPseq/blob/dev/static/flowchart.png?raw=true">
<div>
<p>In general words, the core idea of our approach consists in the sequentiality of the splicing process. Specifically, by focusing on a single pre-RNA of a certain gene, our aim is to model the splicing process as a set of specific steps where each step is dependent on the previous steps and the “environment” of the current step. In particular, we assign each step to the spliceosome decision on whether to include or to exclude a considered exon. Importantly, each step has its probability which depends on the previously included/excluded exons and the expression levels of splicing factors related to the considered exon. The latter means that it is believed that if the splicing factor is an enhancer for a specific exon in a specific tissue, then the higher its expression the higher the probability of the exon inclusion, and otherwise for the silencing splicing factors. Briefly, a pipeline is implemented as follows:</p>
<ol>
  <li><i>Construction of isoform-exon tree:</i> Divide groups of isoforms by sequence of common exons.</li>
  <li><i>List of RBPs with matching motifs:</i> For each node of the tree set a list of RBPs which motifs are present on the node's exon.</li>
  <li><i>Model training:</i>Train linear model based on the specified RNA-seq data for each node of the tree.</li>
  <li><i>Evaluation:</i> evaluate each model and make summary of all steps.</li>
</ol>

Input data can consist from different tissue-samples and datasets, and each dataset should be labeled by one of the following types:
<ol>
<li><i>Training set:</i> samples from training datasets will be used for tuning the regression models. At least one such dataset is required.</li>
<li><i>Validation (test) set:</i> performance of models evaluated on the validation sets. At least one such dataset is required.</li>
</ol>

If not specified, the whole data randomly split into training and validation datasets with 3 to 1 ratio.

</div>

# Installation

### Prerequisites:
Make sure you have installed all of the following prerequisites on your development machine:
  - python3.6+  
  - pip3


### SRPseq installation:  
`pip3 install srpseq`

# Tutorial

In this section we illustrate the main functionality of SRPseq, and
together with that show how to reproduce the results present in 
the manuscript.

<details>
  <summary>The CD44 example</summary>
  
  We illustrate SRPseq basics by using the CD44 gene isoforms, default TCGA and the default Attract + SpliceAid-F datasets. All necessary data for this example can be found in [`tutorial/CD44`](https://github.com/NovosadVictor/SRPseq/blob/main/tutorial/CD44) directory.  
 
  Configuration and output files for this example are
  located in [`tutorial/CD44`](https://github.com/NovosadVictor/SRPseq/blob/main/tutorial/CD44) folder, input data can be downloaded [here](https://eduhseru-my.sharepoint.com/:f:/g/personal/snersisyan_hse_ru/EihEOok4stJFnCjGxqL14qgBSqxLzUy7hBThWp4_jE3HKw?e=bges1q). The microarray data are split into independent training, filtration and validation sets.
 
  The only `"gene": "CD44"` and `"tissue_specific": true` option was used in ([`config.json`](https://github.com/NovosadVictor/SRPseq/blob/main/tutorial/CD44/config.json)):
  
  `srpseq build -c {path_to_config.json}`

  Here we review the results:
  
  - [`CD44_isoforms.png`](https://github.com/NovosadVictor/SRPseq/blob/main/tutorial/CD44/results/CD44_isoforms.png)
  
  7 of 8 CD44 isoforms were used for the analysis.
  The structure of these isoforms, constant and variable exons are shown in the image.
  
  - [`isoforms_tree.png`](https://github.com/NovosadVictor/SRPseq/blob/main/tutorial/CD44/results/isoforms_tree.png)

  Based on the structure of CD44 isoforms, constructed isoform-exon tree is fully imbalanced, meaning that each variable exon separate only one isoform from the remaining group.
  In total, 5 of the variable exons were divider exons (exons 3, 6, 7, 12, 14) and, hence, 5 linear models were trained.
 	
  As in the previous toy example, let us take a closer look to the single gene signature
  (see [`config_for_summary_classifiers.json`](https://github.com/NovosadVictor/SRPseq/blob/main/tutorial/breast_cancer/config_for_summary_classifiers.json)). The following output files were not
  covered in the toy example:
  - [`feature_importances.pdf`](https://github.com/NovosadVictor/SRPseq/blob/main/tutorial/breast_cancer/results_summary_classifiers/feature_importances.pdf) (contains coefficients of the linear SVM classifier)
  - [`model.pkl`](https://github.com/NovosadVictor/SRPseq/blob/main/tutorial/breast_cancer/results_summary_classifiers/model.pkl) (Python-pickled classifier and pre-processor objects)

  ExhauFS also allows one to evaluate constucted classifiers on time-to-event data.
  For example, let us evaluate the same ten-gene signature on 
  additional RNA-seq TCGA-BRCA dataset. To do that, we should include to desired feature 
  subset and pickled model path to the configuration file ([`config_for_km_plot.json`](https://github.com/NovosadVictor/SRPseq/blob/main/tutorial/breast_cancer/config_for_km_plot.json)).
  The analysis could be done by running

  `exhaufs km_plot -c config_for_km_plot.json`
  
  This will generate the Kaplan-Meier plot ([`results_km_plot/KM_TCGA-BRCA_Validation.pdf`](https://github.com/NovosadVictor/SRPseq/blob/main/tutorial/breast_cancer/results_km_plot/KM_TCGA-BRCA_Validation.pdf)).
</details>

# Running SRPseq

## Step 1: data preparation

Before running the tool, you should prepare optional two tsv tables containing expression data of RNA-binding proteins (or all genes expression levels), expression levels of considered gene isoforms and a table with RNA-binding proteins motifs. RNA-seq tables should contain log2 transformed FPKM values associated with samples (rows) and genes / isoforms (columns) with optional columns correspoding to sample tissue and dataset type:

<details>
  <summary>Example</summary>
  
  |            | ESRP1     | QKI       | [Tissue]  | [Dataset.Type] |
  | ---------- | --------- | --------- | --------- | --------       |
  | Sample 1   | 17.17     | 365.1     | TCGA-COAD | Training       |
  | Sample 2   | 56.99     | 123.9     | TCGA-COAD | Validation     |
  | ...        |           |           |           |                |
  | Sample 98  | 22.22     | 123.4     | TCGA-BRCA | Training       |
  | Sample 99  | 23.23     | 567.8     | TCGA-BRCA | Training       |
  | ...        |           |           | | |
  | Sample 511 | 10.82     | 665.8     | TCGA-READ | Validation     |
  | Sample 512 | 11.11     | 200.2     | TCGA-READ | Validation     |
</details>


Table with RNA-binding protein motifs should consist of one column with the gene name of RBP and a second column with RBP motif nucleotide sequences.
<details>
  <summary>Example</summary>
  
  |      RBP      | Motif |
  | ---------- | ----- |
  | ESRP1   | AGGGAU     |
  | ESRP1   | UGGGAAU     |
  | ...        |       |
  | QKI | ACACACUAACCU     |
  | QKI | ACUUAU     |
</details>

## Step 2: creating configuration file

Configuration file is a json file containing all customizable parameters for the pipeline.  

<details>
  <summary>Available parameters</summary> 

  🔴!NOTE! - All paths to files / directories can be either relative to the configuration file directory or absolute paths 
  * `rbp_data_path`
      Optional path to a tsv table containing expression levels of RBPs (by default, combined [TCGA](https://doi.org/10.1038%2Fng.2764) dataset is used).

  * `isoforms_data_path`
      Optional path to a tsv table containing expression levels of selected gene isoforms (by default, combined [TCGA](https://doi.org/10.1038%2Fng.2764) dataset is used).

  * `rbps_path`
      Optional path to a tsv table containing list of RBPs and their motifs (by default, [Attract](https://doi.org/10.1093/database/baw035) and [SpliceAid-F](https://doi.org/10.1093/nar/gks997) datasets are used).

  * `output_dir`
      Path to directory for output files. If it doesn't exist, it will be created.

  * `gene`  
      Gene name for splicing analysis

  * `rbps_tresh_mean`  
      Optional threshold value for expression median of RBPs for them to be considered in the analysis (RBPs with the median expression value lowe than the specified threshold are excluded).

  * `rbps_tresh_var`  
      Optional threshold value for expression variance of RBPs for them to be considered in the analysis (RBPs with the expression variance lowe than the specified threshold are excluded).

  * `isoforms_tresh_mean`  
      Optional threshold value for expression median of isoforms for them to be considered in the analysis (isoforms with the median expression value lowe than the specified threshold are excluded).

  * `isoforms_tresh_var`
      Optional threshold value for expression variance of isoforms for them to be considered in the analysis (isoforms with the expression variance lowe than the specified threshold are excluded).

  * `n_processes`
      Number of processes / threads to run on.
  
  * `random_state`
      Random seed (set to an arbitrary integer for reproducibility).

</details>

## Step 3: running the pipeline

When input data and configuration file are ready,  
the  pipeline can be executed as follows -  
```bash
srpseq build -c <config_file>
```


This will generate multiple files and directory in the specified output folder:
* `isoforms_tree.png`: image of the constructed isoform-exon tree for the specified gene isoforms.
* `{gene}_isoforms.png`: image of the exon-intron structure of all gene isoforms.
* `scores/{score}.png`: directory with plots for accuracy scores on the validation set for all isoforms.
* `tree`: tree-like directory containing results for each tree node.
  * `transcirpts.json`: file containing list of node transcripts, list of parent trancripts and the number of the divider exon.
  * `scores.json`: accuracy scores for the node model.
  * `coefs.csv`: table with the node model coefficients. 
  * `df.csv`: table containing data used for the model fitting (expression levels of used RBPs, isoforms ratio) and the model predicted values. 
  * `[left]`: directory with the results for the left node. 
  * `[right]`: directory with the results for the right node.
  * `[by_tissue]`: directory with coefficients tables for tissues.