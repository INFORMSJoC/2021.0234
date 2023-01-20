[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Data for the Empirical Analysis of the Paper["A Bayesian Semi-supervised Approach to Keyword Extraction with Only Positive and Unlabeled Data"](https://doi.org/) 

This archive is distributed in association with the [INFORMS Journal on Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

This repository includes the data used in the empirical analysis of the paper
[A Bayesian Semi-supervised Approach to Keyword Extraction with Only Positive and Unlabeled Data](https://doi.org/) by G. Wang and Y. Cheng and Y. Xia and Q. Ling and X. Wang.

## Cite

To cite this software, please cite the [paper](https://doi.org/10.1287/ijoc.2021.0234) using its DOI and the software itself, using the following DOI.

[![DOI](https://zenodo.org/badge/568171078.svg)](https://zenodo.org/badge/latestdoi/568171078)

Below is the BibTex for citing this version of the repository.

```
@article{BSS,
  author =        {G. Wang and Y. Cheng and Y. Xia and Q. Ling and X. Wang},
  publisher =     {INFORMS Journal on Computing},
  title =         {{BSS} Version v2021.0234},
  year =          {2022},
  doi =           {10.5281/zenodo.7348935},
  note =          {available for download at https://github.com/INFORMSJoC/2021.0234},
}  
```

## Datasets
--------------------------------------------------------------------------


[text files: /data_JOC]

- short article data from Hulth 2003 :  500 articles, 300 (No.201-No.500) of them are used: ./short_articles/pre_process
all short articles data files are tokenized and pos-taggered.

- short article data ground truths from Hulth 2003 : ./short_article/keywords

Hulth 2003 dataset reference: 

Improved automatic keyword extraction given more linguistic knowledge Anette Hulth. In Proceedings of EMNLP 2003. p. 216-223.

- long article example with ground truth: ./C-42.txt
the article is tokenized, pos-taggered and also stemmed. 

- long article example without ground truth: ./EM_Preprocessed.txt
the article is tokenized, pos-taggered and also stemmed. 

- Amazon review example: ./amazon_pre.txt
the review is tokenized, pos-taggered and also stemmed.

- long article data from Semeval 2010: ./pre_process
- observed keywords from Semeval 2010: ./pre_process_author_truth
- ground truth from Semeval 2010: ./pre_process_reader_truth
all articles are tokenized, pos-taggered and also stemmed.

Semeval 2010 dataset reference:

Kim SN, Medelyan O, Kan MY, Baldwin T (2010) Semeval-2010 task 5 : Automatic keyphrase extraction from scientific articles. Proceedings of the 5th International Workshop on Semantic Evaluation, 21-26.

The original dataset is available: https://semeval2.fbk.eu/semeval2.php?location=download&task_id=6&datatype=test

--------------------------------------------------------------------------
## Code for analysis

[text files: /code_JOC]

Code implementation for BSS keyword extraction method, and its comparison with textrank(TR) and Semi-supervised(SS).

code for 300 short articles:

- ./main_graph.R: main function to run short article examples for Hulth data.

Output - result: overall result for short article examples. To view results, please read comments of this code.

- ./Keyphrase_functions.R: all user-defined functions and packages required, this file is called by running main_graph.R

- ./Component_MCMC.cpp: this is Rcpp file for Component MH method, this file is loaded by Keyphrase_functions.R

code for long article with ground truth:

- ./Long_article_with_truth.R: methods comparison for long articles with the ground truth. The output would be Table 5, but without keywords assembled, which was done manually.

- ./keyphrase_examples.R: comparison for data examples without ground truth. There are two examples: 1. a famous statistics paper. 2. A sample of Amazon review.
The output would be Table 8 and Table 10,respectively. However, adjacent keywords were grouped manually once again. 

- ./semeval_author_obs.R: Experiments on long articles data from Semeval 2010.



