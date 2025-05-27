# Aniti-Malarial-Activity-Prediction-Project
This Project aims to predict and classify the anti-malarial activity of compounds using machine learning techniques.

## Table of Contents  
- [Project Overview](#project-overview)
- [Live Project](#live-project)  
- [Dataset Information](#dataset-information)
- [Project Structure](#project-structure)  
- [Setup Instructions](#setup-instructions)  
  - [Prerequisites](#prerequisites)  
  - [Download and Installation](#download-and-installation)  
- [Featurisation](#featurisation)  
- [Model Building](#model-building)  
- [Model Evaluation](#model-evaluation)  
- [Results and Analysis](#results-and-analysis)
  - [Model Performance Summary](#model-performance-summary)  
  - [Areas of Improvement](#areas-of-improvement)  
- [References](#references)

---

## Project Overview

The antimalaria prediction project is design to implement machine learning algorithms in the early phases of hit identification for potent antimalarial compounds focusing mainly on Plasmodium falciparum. Considering the vast chemical library out there on demand, assessing drug likeliness becomes challenging. This machine learning algorithm will learn from structural features of tested antimalarial compounds to predict novel compounds with potential activity. This is a classification model that will only take the SMILES notation of compounds as input and then predict whether the compound will be active against plasmodium falciparum or not. We will employ the model to screen chemical databases unrelated to antimalarial discovery in a bid for novel hits. These hits will be further subjected to docking and molecular dynamics protocols. 

## Live Project

This project is live as a web app and can be accessed here.

Click on the link above and enter the SMILES notation of any chemical compound to get a prediction! You can as well upload strings of SMILES notation in a CSV or TSV file for batch prediction.

## Dataset Information

This project uses the activity dataset for compounds with reported assays against Plasmodium falciparum from the [Chembl](https://www.ebi.ac.uk/chembl/) and [Pubchem](https://pubchem.ncbi.nlm.nih.gov/) databases.

These data were retrieved via an API executed in python, cleaned and classified as either Active or Inactive based on similar criteria used in Lin, M., Cai, J., Wei, Y., Peng, X., Luo, Q., Li, B., Chen, Y., & Wang, L. (2024). MalariaFlow: A comprehensive deep learning platform for multistage phenotypic antimalarial drug discovery. European journal of medicinal chemistry, 277, 116776 (https://doi.org/10.1016/j.ejmech.2024.116776)

## Project Structure
```bash
├── data/
│   ├──
│   ├──
│
│
│
│
│    └──
├──  models/
│
│
│
│
├──  notebooks/
│
│
│
│
├──  scripts/
│
│
├──  environment.yml
└──  read.md
