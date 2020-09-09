# pygmes is a wraper for GeneMark-ES

This project aims at making GeneMark-ES more usable. By default 
GeneMark-ES has some problems with fragmented and incomplete genomes.

For metagenomic analysis this is an issue, thus I developed pygmes.

## Status
Currently this project is very much WIP, so do expect problems and
incomplete documentation.

It is just on pip already as it is a suggested package for EukCC (https://github.com/Finn-Lab/EukCC/)

## Initial Models used
Currently in the two step algorythm pygmes uses these 6 models 
as a initial set. This might change in the future

GCA_000409445.2.mod: Aureobasidium pullulans (Fungi)
GCA_001243155.1.mod: Sporisorium scitamineum (smut fungi)
GCF_000004695.1.mod: Dictyostelium discoideum AX4 (Amoebozoa)
GCF_000151545.1.mod: Saprolegnia parasitica (Stramenopiles)
GCF_000523435.1.mod: Bipolaris zeicola (Fungi)
GCF_000956335.1.mod: Plasmodium fragile (Alveolata)

I am contemplating reducing the fungi to one model, as they usually are not
the issue, but rather the SAR clade.

