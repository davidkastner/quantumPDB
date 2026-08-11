The find_backbone_atoms: true line in the config.yaml script allows you to output a list of backbone atoms for future QM calculations.

To run this script you will need to use:

qp run -c config.yaml


This will save a *_backbone.txt file containing a list of indices of the Ca, C, N and O atoms in your cluster. An example has been provided in this directory. 
