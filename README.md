# multiScaleFootprinting

![image](https://user-images.githubusercontent.com/44768711/193921131-c7a9f8ab-d123-4689-b62b-d71fdf2abd43.png)

### 1. Overview

In this Github repo we present the multi-scale footprinting framework in our Hu et al. paper. In general, the algorithm takes ATAC-Seq (bulk or single cells) data as input, and try to detect DNA-protein interaction across spatial scales. More specifically, the model first internally corrects the sequence preference of Tn5, and then use a statistical model to calculate footprint score for each position within enhancers and promoters. The process is performed with a range of footprint kernel sizes, capturing DNA-binding proteins of different sizes and shapes. Conceptually, this procedure is similar to wavelet analysis where we decompose the input signal at each location across scales. 

<img src="https://user-images.githubusercontent.com/44768711/193936026-b49715d8-7ec9-4e23-8aa9-330c1f93f2e7.png" width="350" align="left">

The multi-scale footprint pattern at any genomic location delineates local chromatin structure and can be used to infer TF and nucleosome binding. We have shown in our paper that multi-scale footprints can be used as input to neural network models to predict TF binding, even for TFs that do not leave visible footprints on their own.

Additionally, we have implemented the infrastructure for generating pseudo-bulks using single cell data, as well as running multi-scale footprinting using the pseudo-bulked data. This provides us with the unique opportunity to track chromatin structure dynamics across pseudo-time.

<br>

### 2. Key Components

* Correction of Tn5 insertion bias and obtain single-base pair resolution chromatin accessibility

* Calling footprint signal across spatial scales to resolve local chromatin structure

* Tracking nucleosome binding/eviction/sliding across pseudo-time

* Infer TF binding within cis-regulatory elements (CREs)

* Segmentation of CREs into substructures

### 3. Vignettes and Tutorials

Tutorials for running multi-scale footprinting on example data can be found [here][tutorial]

[tutorial]:https://github.com/HYsxe/PRINT/blob/main/analyses/BMMCTutorial/BMMCVignette.pdf

Before running the tutorial, please download the pre-computed bias files from https://zenodo.org/record/7121027#.ZCbw4uzMI8N and put it in the data/shared/precomputedTn5Bias folder.

### 4. References

Hu et al., Multi-scale chromatin footprinting reveals wide-spread encoding of CRE substructures

### 5. Installation

Currently the framework can be installed by cloning the github repo. 

### 6. Support

If you have any questions, please feel free to open an issue. You are also welcome to email me at yanhu@g.harvard.edu. We appreciate everyone's contribution!


### 7. Comming soon:

Currently the tool is implemented in R, which doesn't handle certein computation in the most efficient way. We are working on a ultra-fast python package which will be released soon. Stay tuned!

