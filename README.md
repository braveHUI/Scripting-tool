# Scripting-tool
用于做数据分析的项目
Step 1, download stable version from [Scripting-tool](https://github.com/braveHUI/Scripting-tool.git)

>Step 2, conda create Scripting-tool

```Bash
conda create -n Scripting-tool python=3
source activate Scripting-tool
pip install -r requirement.txt
```
>Step 3, bwa project is to perform a pan-genome comparison of the packages that need to be installed

```Bash
conda install -c conda-forge -c bioconda -c defaults prokka
conda install install roary
```
>Step 4, install taxonkit: 

```Bash
conda activate taxonkit
```
