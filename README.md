[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17162664.svg)](https://doi.org/10.5281/zenodo.17162664)


<h1>DNALongBench: A Benchmark Suite for Long-Range DNA Prediction Tasks</h1>

<h2>Introduction</h2>

DNALongBench is a benchmark of realistic and biologically meaningful genomic DNA prediction tasks that require long-range sequence input and involve long-range dependencies. There are five tasks in our DNALongBench.

![image](./Figure1.v3.png)

<h2>Data Download</h2>

| LR Tasks                       | LR Type                          | Input Length | Output Shape                                    | # Samples               | Metric    |
|--------------------------------|----------------------------------|--------------|-------------------------------------------------|-------------------------|-----------|
| Enhancer-target Gene           | Binary Classification            | 450,000      | 1                                               | 2,602                   | AUROC     |
| eQTL                           | Binary Classification            | 450,000      | 1                                               | 31,282                  | AUROC     |
| Contact Map                    | Binned (2048bp) 2D Regression    | 1,048,576    | 99,681                                          | 7,840                   | SCC & PCC |
| Regulatory Sequence Activity   | Binned (128bp) 1D Regression     | 196,608      | Human: (896, 5,313)<br>Mouse: (896, 1,643)      | Human: 38,171<br>Mouse: 33,521 | PCC       |
| Transcription Initiation Signal | Nucleotide-wise 1D Regression    | 100,000      | (100,000, 10)                                   | 100,000\*               | PCC       |

The data for each task could be downloaded via the following link. Alternatively, you could download the data from [Box](https://cmu.box.com/s/cyn3tqfej3v4tg4xwv1god3jemq7916y).

<h3>Regulatory Sequence Activity Prediction</h3>

<h4>Data Link</h4>

The data can be downloaded at <a href="https://doi.org/10.7910/DVN/MNUEZR">Regulatory Sequence Activity Prediction</a>. 

<h4>Data Details</h4>

We provide the sequences.bed, statistics.json, hg38.ml.fa.fai and hg38.ml.fa.gz files, and reformulate these data into train*/valid*/test*.tfr. The only files you needed are the corresponding tfr files.

<h3>Transcription Initiation Signal Prediction</h3>

<h4>Data Link</h4>

The data can be downloaded at <a href="https://doi.org/10.7910/DVN/VXQKWO">Transcription Initiation Signal Prediction</a>. 

<h4>Data Details</h4>

We provide all the correspoding bed files.

<h3>Enhancer-Target Gene Prediction</h3>

<h4>Data Link</h4>

The data can be downloaded at <a href="https://doi.org/10.7910/DVN/CTEQXX">Enhancer-Target Gene Prediction</a>. 

<h4>Data Details</h4>

The sequences, fa and metrics data are provided.

<h3>Contact Map Prediction</h3>

<h4>Data Link</h4>

The data can be downloaded at <a href="https://doi.org/10.7910/DVN/AZM25S">Contact Map Prediction</a>. 

<h4>Data Details</h4>

We provide the well-split train/valid/test files.

<h3>eQTLP</h3>

<h4>Data Link</h4>

The data can be downloaded at <a href="https://doi.org/10.7910/DVN/YUP2G5">eQTL</a>. 

<h4>Data Details</h4>

The corresponding bed files are provided. 


<!--
### 1. [Regulatory Sequence Activity Prediction](https://dataverse.harvard.edu/privateurl.xhtml?token=4c6b250c-26fc-412a-b3e1-bc15f1332f0c)

### 2. [Transcription Initiation Signal Prediction](https://dataverse.harvard.edu/privateurl.xhtml?token=9810103a-b8b8-4a4d-95c4-b26b6e153446)

### 3. [Enhancer-Target Gene Prediction](https://dataverse.harvard.edu/privateurl.xhtml?token=c238c0dd-528f-4d04-a3c8-0ff1eee1d651)

### 4. [Contact Map Data](https://dataverse.harvard.edu/privateurl.xhtml?token=a990b515-d76e-4b63-ba74-5c78c469ae53)

### 5. [eQTL Data](https://dataverse.harvard.edu/privateurl.xhtml?token=93d446a5-9c75-44bf-be1c-7622563c48d0)
-->


<h2>Experiments</h2>

We've provided the performance of three types of models, which are Expert Model, a lightweight CNN baseline, and a finetuned DNA foundation model (HyenaDNA, Caduceus-Ph and Caduceus-PS). We'll introduce below how to run these models by taking the task of Enhancer-Target Gene Prediction (ETGP) as an example.

| Model |   Expert Model   |  CNN   |  HyenaDNA  |  Caduceus-Ph  |   Caduceus-PS    |   
|:---------------|:---------:|:---------:|:---------:|:---------:|:----------:|
| ETGP        |   **0.926**   |  0.797   |   0.828    |   0.826    |   0.821    |   

<Download Code>

Following the commands below to download our code:

```ruby
conda create -n dnalongbench python=3.9 -y 
conda activate dnalongbench

git clone https://github.com/ma-compbio/DNALONGBENCH.git
pip install .
```

<Load Data>
Use the following Python code to load data for a specific task:

```python
import dnalongbench
from dnalongbench.utils import load_data
train_loader, valid_loader, test_loader = load_data(root=root, task_name = 'contact_map_prediction', subset='HFF', batch_size=16)
```

We also provide data loaders for each task in scripts/data_loaders.ipynb.

<h2>CNN</h2>

Please refer to experiments/CNN/README.md.

<h2>HyenaDNA</h2>

<h3>Environment Setup</h3>

We used the official code of HyenaDNA. The environment setup can be found at <a href="https://github.com/HazyResearch/hyena-dna?tab=readme-ov-file#dependencies">HyenaDNA Enviroment Setup</a>.

Be careful if you would like to use flash attention. Sometimes there are some issues when installing flash attention. We recommend first setup the environment, then activate the enviroment, and finally install flash attention inside the environment. 

<h3>Training & Inference</h3>

For RSAP and TISP tasks, please refer to the readme under experiments/HyenaDNA_RSAP_TISP/README.md 

```ruby
experiments/HyenaDNA/HyenaDNA_RSAP_TISP/README.md
```

For CMP, eQTLP and ETGP tasks, please refer to experiments/HyenaDNA_ETGP_CMP_eQTLP/README.md

```ruby
experiments/HyenaDNA/HyenaDNA_ETGP_CMP_eQTLP/README.md
```

<h2>Caduceus</h2>

<h3>Environment Setup</h3>

We used the official environment provided by Caduceus.

To get started, create a conda environment containing the required dependencies.

```bash
cd experiments/Caduceus/Caduceus_CMP_eQTLP_ETGP

conda env create -f caduceus_env.yml
```

Activate the environment.

```bash
conda activate caduceus_env
```

<h3>Training & Inference</h3>

For RSAP and TISP tasks, please refer to experiments/Caduceus_RSAP_TISP/README.md 

```ruby
experiments/Caduceus_RSAP_TISP/README.md
```

For CMP, eQTLP and ETGP tasks, please refer to experiments/Caduceus_CMP_eQTLP_ETGP/README.md

```ruby
experiments/Caduceus/Caduceus_CMP_eQTLP_ETGP/README.md
```

<h2>GENERator</h2>

<h3>Environment Setup</h3>

We used the official environment setup of GENERator: https://github.com/GenerTeam/GENERator.git.


```bash
git clone https://github.com/GenerTeam/GENERator.git
cd GENERator
pip install -r requirements.txt
```

<h3>Training & Inference</h3>

Please refer to experiments/GENERator/README.md.

```ruby
experiments/GENERator/README.md
```

<h2>Evo2</h2>

<h3>Environment Setup</h3>

We used the official environment setup of Evo2.

Evo 2 is based on [StripedHyena 2](https://github.com/Zymrael/vortex) which requires python>=3.11. Evo 2 uses [Transformer Engine](https://github.com/NVIDIA/TransformerEngine) FP8 for some layers which requires an H100 (or other GPU with compute capability ≥8.9). We are actively investigating ways to avoid this requirement.

To install Evo 2 for inference or generation, please clone and install from GitHub. We recommend using a new conda environment with python>=3.11.

```bash
git clone --recurse-submodules git@github.com:ArcInstitute/evo2.git
cd evo2
pip install .
```

<h3>Training & Inference</h3>

Please refer to experiments/Evo2/README.md.

```ruby
experiments/Evo2/README.md
```

<h2>Citation</h2>
If you find our work helpful, please consider citing our paper.

```
@inproceedings{chengdna,
  title={DNALongBench: A Benchmark Suite for Long-Range DNA Prediction Tasks},
  author={Cheng, Wenduo and Song, Zhenqiao and Zhang, Yang and Wang, Shike and Wang, Danqing and Yang, Muyu and Li, Lei and Ma, Jian}
}
```

