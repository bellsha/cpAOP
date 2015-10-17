# cpAOP
Code used to generate cpAOP network
This repository containes code used in the paper "Integrating publically-available data to generate computationally-predicted adverse outcome pathways for fatty liver" by Shannon Bell, Michelle Angrish, Charles E. Wood, Stephen W. Edwards.

A break down of the workflow is as follows:

Data sources

TG GATES data was downloaded from ftp://ftp.biosciencedbc.jp/archive/open-tggates/LATEST/

EPA ToxCast data is available at http://epa.gov/comptox/toxcast/data.html (accessed 11/2014)

Reactome pathway information (V53) was obtained from http://www.reactome.org/pages/download-data/ (accessed 7/2015)

Code
1)gene2biospace.R <https://github.com/bellsha/cpAOP/blob/master/gene2biospace.R> is used to go from the dataframe of differentailly expressed probes (output of CombineArrayDE,  <https://github.com/bellsha/TGGATESProc/blob/master/CombineArrayDE.R>), and using probe to entreze gene id mappings (ProbeAnnotation.R, <https://github.com/bellsha/TGGATESProc/blob/master/ProbeAnnotation.R>) and rat reactome pathway to uniprot protien mapping (ReactomeCalssv2.R, <https://github.com/bellsha/Reactome2Network/blob/master/ReactomeClassv2.R>) generate a table of enriched Reactome pathways for each chemical treatment, using the hypergeometric distribution <https://github.com/bellsha/cpAOP/blob/master/HyperEnrich.R>.

2)TC2TGGateschemmap.R <https://github.com/bellsha/cpAOP/blob/master/TC2TGGateschemmap.R> maps the matching chemicals from toxcast dataset to the tggates dataset based on the chemcial name to produce an edgelist

3)rules2net.R <https://github.com/bellsha/cpAOP/blob/master/rules2net.R> takes the output of the gene2biospace.R mapping (dataframe of chemical treatments and the differentially expressed reactome pathways) along with the descretized phenotype data (from the pathology and lab workflows https://github.com/bellsha/TGGATESProc/blob/master/TGGATESpathologyprep.R and https://github.com/bellsha/TGGATESProc/blob/master/TGGATESlabprep.R) and generates association rules that are convereted to edge lists. Edge lists are generate fro Phenotype-Phenotype connections, Reactome-Reactome connections, Reactome-Chemical, and Reactome-Phenotype

4)NodeEdgeLabel.R <https://github.com/bellsha/cpAOP/blob/master/NodeEdgeLabel.R> takes all the edge lists generated in srules2net.R along with the other processing steps to make a master none/edge dataframe for cytoscape along with a label file for annotating the nodes.

5)cpAOPextraction.R <https://github.com/bellsha/cpAOP/blob/master/cpAOPextraction.R> is code to extract the cpAOP subnetworks
