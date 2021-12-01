# Genetic architecture of dispersal and local adaptation drives accelerating range expansions.
Code for:

Jhelam N. Deshpande and Emanuel A. Fronhofer: Genetic architecture of dispersal and local adaptation drives accelerating range expansions.

Abstract: 

Contemporary evolution has the potential to significantly alter biotic responses to global change, including range expansion dynamics and biological invasions. However, predictive models often make highly simplifying assumptions about the genetic architecture underlying relevant traits. This can be problematic since genetic architecture defines evolvability, that is, evolutionary rates, and higher order evolutionary processes, which determine whether evolution will be able to keep up with environmental change or not. Therefore, we here study the impact of the genetic architecture of dispersal and local adaptation, two central traits of high relevance for range expansion dynamics, on the speed and variability of range expansions into an environmental gradient, such as temperature. In our theoretical model we assume that dispersal and local adaptation traits result from the products of two non-interacting gene-regulatory networks (GRNs). We compare our model to simpler quantitative genetics models and show that in the GRN model, range expansions are accelerated, faster and more variable. Increased variability implies that these evolutionary changes reduce predictability. We further find that acceleration in the GRN model is primarily driven by an increase in the rate of local adaptation to novel habitats which results from greater sensitivity to mutation (decreased robustness) and increased gene expression. Our results highlight how processes at microscopic scales, here, within genomes, can impact the predictions of large scale, macroscopic phenomena such as range expansions by modulating the rate of evolution.

Description:

In each folder there is the associated cpp code range_expansion.cpp, the parameter input file and the Makefile corresponding to the models described below.

The code range_expansion.cpp in GRN_model assumes that two non-interacting GRNs encode dispersal and local adaptation

The code range_expansion.cpp in single _locus_model assumes that one locus each encodes dispersal and local adaptation

The code range_expansion.cpp in GRN_local_adaptation_single_locus_dispersal assumes that a GRN encodes local adaptation and dispersal is encoded by a single quantitative locus

The code range_expansion.cpp in single_locus_local_adaptation_GRN_dispersal assumes that a single locus encodes local adaptation and dispersal is encoded by a GRN
