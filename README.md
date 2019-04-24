# Installation recommendations

At the time of submission, it is recommended to install songbird within a [qiime2-2019.1 environment](https://docs.qiime2.org/2019.1/install/native/#install-qiime-2-within-a-conda-environment).
The installation of songbird would look like
```
source activate qiime2-2019.1
pip install songbird
```

The scripts used to run songbird can be found in the scripts directory

# Data files

All data files are stored under the data directory.

- oral_deblurred.biom : the deblurred microbe abundance table.
- oral_metadata.txt : sample metadata information for the saliva time series dataset
- oral-collapsed-table.qza : the qiime2 artifact of collapsed microbial abundances in the saliva time series dataset.  This used taxonomy collapsing and closed ref picking.
- oral-collapsed-table.biom : the biom file of collapsed microbial abundances in the saliva time series dataset. This was extracted from oral-collapsed-table.qza
- byrd_skin_table.biom : The microbe abundances from the shotgun dataset in the atopic dermitatis dataset in Byrd et al 2017.
- byrd_final_microbes.txt: Microbes that were highlighted to be interesting Byrd et al
- byrd_metadata.txt : Sample metadata in Byrd et al
- read_counts_BacteriaMalassezia_62subjectsUPto3mismatchs_ADonly.txt : M. globosa and P. acnes counts in Leung et al
- soils/ : directory containing biom tables and sample metadata to reproduce the Central Parks analysis.

Note that the qiime2 artifacts contain all of the processing information in the providence, so it is possible to retrace all of the qiime2 commands that were performed in order to reproduce those artifacts.

# Analysis notebooks

The analysis notebooks can be found under the ipynb folder
