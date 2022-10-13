
---

<div align="center">    
 
# KA-Search: Rapid and exhaustive sequence identity search of known antibodies

---
    
by
Tobias H. Olsen\,$^{1,\dagger}$, Brennan A. Kenyon\,$^{1,\dagger}$, Iain H. Moal\,$^{2}$ and Charlotte M. Deane$^{1,3}$
    
$^{1}$ Oxford Protein Informatics Group, Department of Statistics, University of Oxford, Oxford, United Kingdom  
$^{2}$ GSK Medicines Research Centre, GlaxoSmithKline plc, Stevenage, United Kingdom  
$^{3}$ XXXX, Exscientia plc, Oxford, United Kingdom  
$^{\dagger}$ These authors contributed equally to this work and share first authorship  
    
</div>


<div style="text-align:center"><img src="data/tool_comparison.png" width="800"/></div>

*Speed and sensitivity comparison between KA-Search and three commonly used protein sequence identity search tools. Speed was measured by the time to search 10 million sequence with a single query and sensitivity by how often the tools returned the closest or the closest within the top-100 match for 100 heavy chains against the same 10 million sequences.*




## Abstract
Antibodies with similar amino acid sequences, especially in the complementary-determining regions, often share certain properties. Finding that an antibody of interest has a similar variable domain sequence to naturally expressed antibodies in healthy or diseased repertoires can therefore be a powerful approach for the prediction of antibody properties, such as immunogenicity or antigen specificity. However, as the number of available antibody repertoire sequences has exceeded billions and is continuing to grow, repertoire mining for highly similar sequences has become increasingly computationally expensive. Existing approaches are limited by either being low-throughput, non-exhaustive, not antibody-specific, or only searching against entire chain sequences. Therefore, there is a need for a specialized tool, optimized for a rapid and exhaustive search of any antibody region against all known antibodies, to better utilize the full number of available repertoire sequences.

We introduce Known Antibody Search (KA-Search), a tool that allows for rapid search of the 2.4 billion antibodies in the Observed Antibody Space (OAS) database by sequence identity across either the whole chain, the CDRs, or a user defined antibody region. KA-Search can be used to find the most similar sequences from OAS within 30 minutes using 5 CPUs. We demonstrate how KA-Search can be used to obtain new insights about an antibody of interest. KA-Search is freely available at https://github.com/oxpig/kasearch.

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

# KA-Search guide

**A Jupyter notebook** showcasing KA-Search can be found [here](https://github.com/oxpig/kasearch/blob/main/notebooks/examples.ipynb). 


---------



### Citation   
Work in preparation
