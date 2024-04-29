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
- **Data Preprocessing**: Prior to modeling, the dataset underwent preprocessing steps, including cleaning, normalization, and feature engineering.

## Feature Encoding

- **Feature Encoding Techniques**: Categorical features were encoded using [technique], while numerical features were [scaled/normalized/etc.].
- **Feature Scaling**: Numerical features were scaled using [scaling method] to ensure uniformity across features.

## GAN Augmentation

- **GAN Architecture**: The Generative Adversarial Network (GAN) architecture used for data augmentation consisted of [describe architecture].
- **Data Augmentation Process**: The GAN was trained on the dataset to generate synthetic data samples, which were then combined with the original dataset for training.
- **Benefits of GAN Augmentation**: GAN augmentation helped in increasing the diversity of the dataset and improving model generalization.

## Comparison

- **Evaluation Metrics**: Model performance was evaluated using [metrics], including accuracy, precision, recall, and F1-score.
- **Experimental Setup**: Models were trained using [architectures], with hyperparameters tuned using [technique]. The dataset was split into training/validation/test sets in the ratio [ratio].
- **Results**: [Table or description of results, including performance metrics and visualizations].
- **Discussion**: Analysis of results revealed [insights], highlighting the strengths and weaknesses of each approach.

## Usage

- **Installation**: To install dependencies, run `pip install -r requirements.txt`.
- **Usage Instructions**: To train the model, run `python train.py`. For inference, use `python predict.py`.
- **Example Code**: Below is an example code snippet for training the model:

```python
# Insert example code here
