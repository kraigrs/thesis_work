#!/usr/bin/sh

##### Kellis network #####

# regulators

perl chromosomal_regulators.pl ../../S1_integrative_regulatory_networks/flynet_supervised_0.6.txt ../../FlyBase_Fields_download.txt 2L > ../../Kellis_chromosomal_regulators_2L.txt

perl chromosomal_regulators.pl ../../S1_integrative_regulatory_networks/flynet_supervised_0.6.txt ../../FlyBase_Fields_download.txt 2R > ../../Kellis_chromosomal_regulators_2R.txt

perl chromosomal_regulators.pl ../../S1_integrative_regulatory_networks/flynet_supervised_0.6.txt ../../FlyBase_Fields_download.txt 3L > ../../Kellis_chromosomal_regulators_3L.txt

perl chromosomal_regulators.pl ../../S1_integrative_regulatory_networks/flynet_supervised_0.6.txt ../../FlyBase_Fields_download.txt 3R > ../../Kellis_chromosomal_regulators_3R.txt

perl chromosomal_regulators.pl ../../S1_integrative_regulatory_networks/flynet_supervised_0.6.txt ../../FlyBase_Fields_download.txt X > ../../Kellis_chromosomal_regulators_X.txt

# targets

perl chromosomal_targets.pl ../../S1_integrative_regulatory_networks/flynet_supervised_0.6.txt ../../FlyBase_Fields_download.txt 2L > ../../Kellis_chromosomal_targets_2L.txt

perl chromosomal_targets.pl ../../S1_integrative_regulatory_networks/flynet_supervised_0.6.txt ../../FlyBase_Fields_download.txt 2R > ../../Kellis_chromosomal_targets_2R.txt

perl chromosomal_targets.pl ../../S1_integrative_regulatory_networks/flynet_supervised_0.6.txt ../../FlyBase_Fields_download.txt 3L > ../../Kellis_chromosomal_targets_3L.txt

perl chromosomal_targets.pl ../../S1_integrative_regulatory_networks/flynet_supervised_0.6.txt ../../FlyBase_Fields_download.txt 3R > ../../Kellis_chromosomal_targets_3R.txt

perl chromosomal_targets.pl ../../S1_integrative_regulatory_networks/flynet_supervised_0.6.txt ../../FlyBase_Fields_download.txt X > ../../Kellis_chromosomal_targets_X.txt

##### expression data #####

# regulators

perl chromosomal_regulators_subset.pl ../../S1_integrative_regulatory_networks/flynet_supervised_0.6.txt ../../FlyBase_Fields_download.txt ../../ASE_genes.txt 2L > ../../Dmel_chromosomal_regulators_2L.txt

perl chromosomal_regulators_subset.pl ../../S1_integrative_regulatory_networks/flynet_supervised_0.6.txt ../../FlyBase_Fields_download.txt ../../ASE_genes.txt 2R > ../../Dmel_chromosomal_regulators_2R.txt

perl chromosomal_regulators_subset.pl ../../S1_integrative_regulatory_networks/flynet_supervised_0.6.txt ../../FlyBase_Fields_download.txt ../../ASE_genes.txt 3L > ../../Dmel_chromosomal_regulators_3L.txt

perl chromosomal_regulators_subset.pl ../../S1_integrative_regulatory_networks/flynet_supervised_0.6.txt ../../FlyBase_Fields_download.txt ../../ASE_genes.txt 3R > ../../Dmel_chromosomal_regulators_3R.txt

perl chromosomal_regulators_subset.pl ../../S1_integrative_regulatory_networks/flynet_supervised_0.6.txt ../../FlyBase_Fields_download.txt ../../ASE_genes.txt X > ../../Dmel_chromosomal_regulators_X.txt

# targets

perl chromosomal_targets_subset.pl ../../S1_integrative_regulatory_networks/flynet_supervised_0.6.txt ../../FlyBase_Fields_download.txt ../../ASE_genes.txt 2L > ../../Dmel_chromosomal_targets_2L.txt

perl chromosomal_targets_subset.pl ../../S1_integrative_regulatory_networks/flynet_supervised_0.6.txt ../../FlyBase_Fields_download.txt ../../ASE_genes.txt 2R > ../../Dmel_chromosomal_targets_2R.txt

perl chromosomal_targets_subset.pl ../../S1_integrative_regulatory_networks/flynet_supervised_0.6.txt ../../FlyBase_Fields_download.txt ../../ASE_genes.txt 3L > ../../Dmel_chromosomal_targets_3L.txt

perl chromosomal_targets_subset.pl ../../S1_integrative_regulatory_networks/flynet_supervised_0.6.txt ../../FlyBase_Fields_download.txt ../../ASE_genes.txt 3R > ../../Dmel_chromosomal_targets_3R.txt

perl chromosomal_targets_subset.pl ../../S1_integrative_regulatory_networks/flynet_supervised_0.6.txt ../../FlyBase_Fields_download.txt ../../ASE_genes.txt X > ../../Dmel_chromosomal_targets_X.txt