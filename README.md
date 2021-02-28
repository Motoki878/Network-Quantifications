# Description
   These codes reproduces results/figures in the report, which 
   expeted to be published as Kajiwara et al. (2021).
  
--------------
# Requirments
    
  Matlab, Statistical and Signal Processing Toolbox and Statistics and Machine Learning Toolbox are necessary to run these code.
  
----------------
# To generate figure 5b
  
  addpath(genpath('./'));
  fig5b
  
----------------
# Description of fig5b.m
 
  This program products a figure showing degree histgrams, histgrams of number of connections.

-----------------
# To generate figures 5c and 5e
  
   addpath(genpath('./'));
   fig5ce
  
----------------
# Description of fig5ce.m

  This program products two figures showing the histgram of firing rate for excitatory and inhibitory neurons, and the histgram of connectivity strenghs for excitatory and inhibitory neurons.

-----------------
# To generate figure 5d
  
   addpath(genpath('./'));
   fig5d
  
----------------
# Description of fig5d.m

  This program products tje figure showing the total number of identified excitatory and inhibitory neurons within cortex, and their ralative ratios.

-----------------
# To generate figure 6a and 6b
  
   addpath(genpath('./'));
   fig6ab
  
----------------
# Description of fig6ab.m

  This program products two figures showing the difference of averaged k-core values for all cortical excitatory and inhibitory neurons, and showing  the difference of averaged k-core values for excitatory and inhibitory neurons within individual cortical layers.

-----------------
# To generate figure 7a and 7b
  
   addpath(genpath('./'));
   fig7ab
  
----------------
# Description of fig7ab.m

   This program products two figures showing the difference of ratios of FVS  (Feedback Vertex Set) for all cortical excitatory or inhibitory neurons, and showing the difference of ratios of FVS for excitatory or inhibitory neurons within individual cortical layers.

-----------------
# To generate basic data for drawing figure 7c
  
   addpath(genpath('./'));
   fig7c
  
----------------
# Description of fig7c.m

   This program products data set for comparing between high K-Core nodes and FVS (Feedback Vertex Set) of nodes. The data was utilized to drowing the Venn diagram.
   
------------------
   # Reference
   If you use this code, cite this following article: 
   Kajiwara et al. (2021) under review.

