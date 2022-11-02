
---

<div align="center">    
 
# KA-Search: Rapid and exhaustive sequence identity search of known antibodies

---
    
by
Tobias H. Olsen $^{1,\dagger}$, Brennan A. Kenyon $^{1,\dagger}$, Iain H. Moal $^{2}$ and Charlotte M. Deane $^{1,3}$
    
$^{1}$ Oxford Protein Informatics Group, Department of Statistics, University of Oxford, Oxford, United Kingdom  
$^{2}$ GSK Medicines Research Centre, GlaxoSmithKline plc, Stevenage, United Kingdom  
$^{3}$ Exscientia plc, Oxford, United Kingdom  
$^{\dagger}$ These authors contributed equally to this work and share first authorship  
    
</div>


<!---<div style="text-align:center"><img src="data/tool_comparison.png" width="800"/></div> 

*Speed and sensitivity comparison between KA-Search and three commonly used protein sequence identity search tools. Speed was measured by the time to search 10 million sequence with a single query and sensitivity by how often the tools returned the closest or the closest within the top-100 match for 100 heavy chains against the same 10 million sequences.*
--->

## Abstract
Antibodies with similar amino acid sequences, especially across their complementary-determining regions, often share properties. Finding that an antibody of interest has a similar sequence to naturally expressed antibodies in healthy or diseased repertoires is a powerful approach for the prediction of antibody properties, such as immunogenicity or antigen specificity. However, as the number of available antibody sequences is now in the billions and continuing to grow, repertoire mining for similar sequences has become increasingly computationally expensive. Existing approaches are limited by either being low-throughput, non-exhaustive, not antibody-specific, or only searching against entire chain sequences. Therefore, there is a need for a specialized tool, optimized for a rapid and exhaustive search of any antibody region against all known antibodies, to better utilize the full breadth of available repertoire sequences.

We introduce Known Antibody Search (KA-Search), a tool that allows for rapid search of billions of antibody sequences by sequence identity across either the whole chain, the CDRs, or a user defined antibody region. We show KA-Search in operation on the ~2.4 billion antibody sequences available in the OAS database. KA-Search can be used to find the most similar sequences from OAS within 30 minutes using 5 CPUs. We give examples of how KA-Search can be used to obtain new insights about an antibody of interest. KA-Search is freely available at https://github.com/oxpig/kasearch.


-----------

# Software implementation

KA-Search is freely available and can be installed with pip.

~~~.sh
    pip install kasearch
~~~

or directly from github.

~~~.sh
    pip install -U git+https://github.com/oxpig/kasearch
~~~


**NB:** You need to manually install a version of [ANARCI](https://github.com/oxpig/ANARCI) in the same environment.

----------

# Download pre-aligned data to search against

This list contains the download links for the paper version of the pre-aligned OAS and any future releases, ready for KA-Search. 

**NB**: The following datasets are large, you should therefore ensure you have enough space before trying to download them.

- [OAS-aligned](http://opig.stats.ox.ac.uk/webapps/ngsdb/kasearch_aligned_oas/paper_aligned_oas_sep2022.tar) (Paper version), a pre-aligned version of OAS, from September 2022, with 2.4 billion sequences taking up ~63GB. 
- [OAS-aligned-small](https://zenodo.org/record/7079547/files/oasdb_small.tar) (Paper version), a pre-aligned version of OAS, from September 2022, with 144 million sequences taking up ~4.4GB. 


After downloading, extract the pre-aligned dataset with "tar -xf downloaded_file.tar". Give the extacted dataset path when initiating KA-Search to search against it. See how to do this by following the KA-Search notebook guide below.


---------

# KA-Search guide

**A Jupyter notebook** showcasing KA-Search can be found [here](https://github.com/oxpig/kasearch/blob/main/notebooks/examples.ipynb). 

KA-Search can also be run using the following [Colab](https://colab.research.google.com/github/TobiasHeOl/kasearch/blob/main/notebooks/KAsearch_colab.ipynb).

---------



### Citation

```
@article{Olsen2022,
  title={KA-Search: Rapid and exhaustive sequence identity search of known antibodies},
  author={Tobias H. Olsen, Brennan A. Kenyon, Iain H. Moal and Charlotte M. Deane},
  journal={bioRxiv},
  doi={},
  year={2022}
}
```  
