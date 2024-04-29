# GAN-ML
Generative based data augmentation for ACPs


#GAN-ML: Generative based data augmentation for anticancer peptide prediction using inter-
pretable machine learning

## Overview

This project aims to [develop a generative based data augmentation method as a nascent approach to Sampling based and feature based augmentation].

## Dataset

- **Description**: [Main Dataset] Dataset The dataset used in this project is sourced from [AntiCP2.0] and consists of [861] POSITVE samples and [861] NEGATIVE samples  in fasta format. Additionally, this project consist of alterante dataset sourced from AntiCP 2.0 which consists of 970 positive sample and 970 negative samples in fasta format.
- **Features**: The dataset includes the following features:
  - Main Dataset: Positive samples : ACP (Anticancer peptides without antimicrobial activity) and Negative Samples : AMP(Antimicrobial peptides without Anticancer Activity (ACPs).
  - Alternate Dataset: Positive samples : ACP (Anticancer peptides (with or without antimicrobial activity) and Negative Samples :Non-Secretory Protein samples from [UniProt](https://www.uniprot.org/).
  - ...
- **Data Preprocessing**: Prior to modeling, the dataset underwent preprocessing steps, including cleaning, normalization/scaling, and feature engineering specially MRMD based method, SFS Method, Boruta Method, and Shapely Method for Feature Selection .

## Feature Encoding

- **Feature Encoding Techniques**: Peptide were encoded using AAC, AAE, DPC, AAI, GTPC. CKSAAGP, PAAC, CTD, PSAAC, ASDC, QSO, Comprehensive feature .
- **Feature Scaling**: Numerical features were scaled using Min-Max Scaling and Normalization method to ensure uniformity across features.

## GAN Augmentation

- **GAN Architecture**: The Generative Adversarial Network (GAN) architecture used for data augmentation consisted of DC-GAN architecture where peptides were encoded using a PC6 (Physiochemcial Encoding method).
- **Data Augmentation Process**: The GAN was trained on the dataset to generate synthetic data samples specially Anticancer peptides, which were then combined with the original dataset for training.
- **Benefits of GAN Augmentation**: GAN augmentation helped in increasing the diversity of the dataset and improving model generalization.

## Comparison

- **Evaluation Metrics**: Model performance was evaluated using balanced and imbalanced case metrics, including accuracy, precision, sensitivity, specificity, MCC, G-mean, AUC-ROC and F1-score.
- **Experimental Setup**: Models were trained using DC-GAN Architecture. The dataset was split into training/validation/test sets in the ratio 80:20.
- 
- **Discussion**: Analysis of results revealed Generative based data augmentation surpasses Oversampling based approach and information based data augmentation, highlighting the strengths GAN-ML over other state of art method in ACPs prediction.

## Usage

- **Installation**: To install dependencies, run `pip install -r requirements.txt`.
- **Usage Instructions**: Encode peptides, run `fea_extract.py`. generates a peptide encoding in Data Frame format.
- **Example Code**: Below is an example code snippet for generatingh encoding:

```python
# from fea_extract import read_fasta,insert_AAC,insert_DPC,insert_CKSAAGP,insert_CTD,insert_PAAC,insert_AAI,insert_GTPC,insert_QSO,insert_AAE,insert_PSAAC,insert_word2int,insert_ASDC
