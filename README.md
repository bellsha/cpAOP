# cpAOP
Code used to generate cpAOP network
This repository containes code used in the paper "Integrating publically-available data to generate computationally-predicted adverse outcome pathways for fatty liver" by Shannon Bell, Michelle Angrish, Charles E. Wood, Stephen W. Edwards.

A break down of the workflow is as follows:

Data sources

TG GATES data was downloaded from ftp://ftp.biosciencedbc.jp/archive/open-tggates/LATEST/

EPA ToxCast data is available at http://epa.gov/comptox/toxcast/data.html (accessed 11/2014)

Reactome pathway information (V53) was obtained from http://www.reactome.org/pages/download-data/ (accessed 7/2015)

Code
The following scripts are run in order. Paths internal to the code may need to be changed to reflect location on the local drive.

1)gene2biospace.R <https://github.com/bellsha/cpAOP/blob/master/gene2biospace.R> is used to go from the dataframe of differentailly expressed probes (output of CombineArrayDE,  <https://github.com/bellsha/TGGATESProc/blob/master/CombineArrayDE.R>), and using probe to entreze gene id mappings (ProbeAnnotation.R, <https://github.com/bellsha/TGGATESProc/blob/master/ProbeAnnotation.R>) and rat reactome pathway to uniprot protien mapping (ReactomeCalssv2.R, <https://github.com/bellsha/Reactome2Network/blob/master/ReactomeClassv2.R>) generate a table of enriched Reactome pathways for each chemical treatment, using the hypergeometric distribution <https://github.com/bellsha/cpAOP/blob/master/HyperEnrich.R>. Output file used in the networks: https://github.com/bellsha/cpAOP_Supplementary/blob/master/PathEnrichReactomeDataFrame.txt

2)TC2TGGateschemmap.R <https://github.com/bellsha/cpAOP/blob/master/TC2TGGateschemmap.R> maps the matching chemicals from toxcast dataset to the tggates dataset based on the chemcial name to produce an edgelist (TC2TGgatesChemMap20150204.txt)

3)rules2net.R <https://github.com/bellsha/cpAOP/blob/master/rules2net.R> takes the output of the gene2biospace.R mapping (dataframe of chemical treatments and the differentially expressed reactome pathways, https://github.com/bellsha/cpAOP_Supplementary/blob/master/PathEnrichReactomeDataFrame.txt) along with the descretized phenotype data from the pathology (https://github.com/bellsha/TGGATESProc/blob/master/TGGATESpathologyprep.R, https://github.com/bellsha/cpAOP_Supplementary/blob/master/open_tggates_pathologyDESC.txt) and lab workflows ( https://github.com/bellsha/TGGATESProc/blob/master/TGGATESlabprep.R, https://github.com/bellsha/cpAOP_Supplementary/blob/master/LiverLabDescData201411.txt) and generates association rules that are convereted to edge lists along with edge lists from direct experimental observations. Edge lists are generate from Phenotype-Phenotype connections, Reactome-Reactome connections, Reactome-Chemical, and Reactome-Phenotype. Edgelists are combined in (https://github.com/bellsha/cpAOP_Supplementary/blob/master/NetworkTable.txt)

4)NodeEdgeLabel.R <https://github.com/bellsha/cpAOP/blob/master/NodeEdgeLabel.R> takes all the edge lists generated in rules2net.R along with the other processing steps to make a master none/edge dataframe for cytoscape along with a label file for annotating the nodes. This is an intermediate step, final 

5)cpAOPextraction.R <https://github.com/bellsha/cpAOP/blob/master/cpAOPextraction.R> is code to extract the cpAOP subnetworks. major products are: https://github.com/bellsha/cpAOP_Supplementary/blob/master/NetworkTable.txt; and all 3 networks have been uploaded into https://github.com/bellsha/cpAOP_Supplementary/blob/master/cpAOPNetworkFile.cys
