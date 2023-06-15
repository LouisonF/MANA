# MANA
MANA: mMoA identification Assisted by modelling and Network Analysis
This repository contains code and a test case associated with the article named : 

    New genome scale network modelling and mining workflow to detect metabolic modulations induced by the exposure to xenobiotics

The workflow presented in this article aims at improving our understanding of the metabolic Mechanism of Action and can be divided in three steps:
1. Condition-specific metabolic network modelling with partial enumeration from gene expression data
2. Identify Differentially Activated Reactions (DAR) from the modelised sets of condition-specific metabolic networks
3. Network anaylsis to extract minimal subnetworks représentative of the chemical's mMoA

Each step of the workflow is performed by a jupyter notebook:
* **partial_enumeration.ipynb**
* **dars_calculation.ipynb**
* **analysis.ipynb**

Properties and parameters for the workflow are stored in a unique file to update in order to change parameters such as compound, dose, time, etc:
* **props.properties**

The package source code is contained in the mana folder and can be installed as a python module with the following command (require pip) launched from the git repository:

<code>pip install .</code>
## Requirements:

* Python3.9X
* Java 11
* Met4j 1.2.2 jar (stored in this repository)
* CPLEX 12.10 or newer (not required for the test case)

## Launching the test case:

As a test case, we compared Primary Human Hepatocytes (PHH) exposed to Amiodarone during 24 hours to its controls.
Launching **master_notebook.ipnyb** will perform the complete workflow on this test case. The notebook will pause near the end of the workflow waiting for you to provide the desired number of clusters at the hierarchical clustering step.


