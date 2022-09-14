
---

<div align="center">    
 
# KA-Search: Rapid and exhaustive sequence identity search of known antibodies


</div>


Antibodies with similar amino acid sequences, especially in the complementary-determining regions (CDRs), often share certain properties. It is often powerful to compare the sequence of an antibody of interest against natural antibody repertoires, as finding similar antibodies in nature can indicate likely specificity or immunogenicity. However, as the number of available antibody repertoire sequences has exceeded a billion and is continuing to grow, repertoire mining for highly similar sequences has become increasingly computationally expensive. Existing approaches are limited by either being low-throughput, non-exhaustive, not antibody-specific, or only searching against entire chain sequences. Therefore, there is a need for a specialized tool, optimized for a rapid and exhaustive search of any antibody region against all known antibodies, to better utilize the full number of available repertoire sequences.

Here, we introduce Known Antibody Search (KA-Search), a tool that allows for rapid search of the 2.4 billion antibodies in the Observed Antibody Space (OAS) database by sequence identity across either the whole chain, the CDRs, or a user defined antibody region. KA-Search can be used to find the most similar sequences from OAS within 20 minutes using 2 CPUs. We demonstrate how KA-Search can be used to obtain new insights about an antibody of interest. KA-Search is freely available at https://github.com/oxpig/kasearch.

-----------

# Install KA-Search

KA-Search is freely available and can be installed with pip.

~~~.sh
    pip install kasearch
~~~

or directly from github.

~~~.sh
    pip install -U git+https://github.com/oxpig/kasearch
~~~


Additionally, you need to install a version of [ANARCI](https://github.com/oxpig/ANARCI) in the same environment.

----------

# KA-Search guide

**A Jupyter notebook** showcasing KA-Search can be found [here](https://github.com/oxpig/kasearch/blob/main/notebooks/examples.ipynb). 


---------



### Citation   
Work in preparation
