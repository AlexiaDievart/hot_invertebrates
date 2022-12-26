# Repository description
Infrared pictures, data and R scripts - Dievart et al. (2023)

## General information

**Citation**: 

**Principal investigator**: [Alexia M. A. DIEVART](https://scholar.google.com/citations?user=1CQgX5kAAAAJ&hl=fr&oi=ao), Coastal Research Group, Rhodes University (South Africa)

**Co-investigators**:
* [Christopher D. McQUAID](https://scholar.google.com/citations?user=uNl9g6wAAAAJ&hl=fr&oi=ao)
* [Gerardo I. ZARDI](https://scholar.google.com/citations?user=s8019k0AAAAJ&hl=fr&oi=ao)
* [Katy R. NICASTRO](https://scholar.google.com/citations?user=UUOXLPcAAAAJ&hl=fr&oi=ao)
* [Pierre W. FRONEMAN](https://scholar.google.com/citations?user=G5tEQu4AAAAJ&hl=fr&oi=ao)

**Corresponding investigator**: Alexia M. A. DIEVART, alexia.dievart@hotmail.fr

**Keywords**: infrared thermography; ecosystem engineers; *Perna perna*; ecosystem functioning; desiccation stress; heat stress; parasitism; mutualism; invertebrate communities

**Fundings**: This study was funded by the National Research Foundation (NRF) of South Africa, grant number 64801 and further supported by ANR SAN22202 to K.R.N. This project has received funding from the European Union's Horizon 2020 Research and Innovation programme, under the Marie Skłodowska Curie, grant agreement No 101034329. G.I.Z. is the recipient of the WINNINGNor-mandy Program supported by the Normandy Region. 

### Abstract
Mussel beds form important intertidal matrices that provide thermal buffering to associated in-faunal communities, especially under stressful environmental conditions. Mussel shells are often parasitized by photoautotrophic euendoliths, which have indirect beneficial thermoregulatory effects on both solitary and aggregated mussels by increasing the albedo of the shell. We investigated whether euendolithic infestation of artificial mussel beds (*Perna perna*) influences the body temperatures of 4 associated macroinvertebrate species during simulated periods of emersion, using non-invasive infrared thermography and the captured shell temperature as a proxy. Shell temperatures of the limpet [*Scutellastra granularis*](https://fockfish.wordpress.com/2022/10/07/granular-limpet-scutellastra-granularis/) and the chiton *Acanthochitona garnoti* were higher in non-infested than in infested mussel beds, indicating that endolith-induced improvements of humidity and temperature in mussel beds could benefit associated invertebrate organisms. However, this was not the case for the limpet *Helcion pectunculus* or the topshell *Oxystele antoni*. This beneficial thermal buffering offered by euendolithic infestation of the mussel beds was relevant only if the organism is under heat stress. With global climate change, the indirect beneficial effect of euendolithic infestation for invertebrate communities associated with mussel beds may mitigate intertidal local extinction events triggered by marine heatwaves.

## Data collection and curation

### Geographic location
The series of experiments were conducted on a flat roof of the Zoology and Ento-mology Department of Rhodes University, in Grahamstown, during austral summer (April and May 2021 and February 2022), on sunny days, with little to no overcast and low wind speeds and high sun elevation (10 am – 3 pm) to maximize the like-lihood of extreme temperature and desiccation stress.

### Experimental dates
#### Body temperature assessment
* _Acanthochitona garnoti_ (ACAGAR): 9 February 2022, 10 February 2022
* _Helcion pectunculus_ (HELPEC): 13 April 2021, 15 April 2021
* _Oxystele antoni_ (OXYANT): 13 April 2021, 15 April 2021
* _Scutellastra granularis_ (SCUGRA): 14 April 2021, 18 April 2021
#### Desiccation assessment
* _Acanthochitona garnoti_: 1 May 2021
* _Helcion pectunculus_: 17 May 2021
* _Oxystele antoni_: 18 April 2021
* _Scutellastra granularis_: 2 May 2021

#### Softwares
* [**IRSoft 4.8**](https://www.testo.com/fr-FR/produit/thermography-irsoft) was used to calculate the body temperatures of biomimetic mussels and invertebrates from infrared pictures.
* [**R and R Studio**](https://www.R-project.org/) was used to run the statistical analyses.
* **R packages**: _gratia_, _mgcv_.  

## Data and file overview

### Data files

* **infrared_data.csv** includes body temperatures (in °C) of each invertebrate investigated on each experimental date.
* **desic_data.csv** includes wet and dry weights (in mg) of each invertebrate investigated on each experimental date.

### R scripts

* **infared_ms_mussels.R** includes the statistical analyses for all biomimetic mussels on each experimental date.
* **infrared_ms_acagar.R** includes the statistical analyses for the chiton _Acanthochitona garnoti_ on each experimental date.
* **infrared_ms_oxyant+helpec.R** includes the statistical analyses for the topshell _Oxystele antoni_ and the limpet _Helcion pectunculus_ on each experimental date, these invertebrates being investigated simultaneously.
* **infrared_ms_scugra.R** includes the statistical analyses for the limpet _Scutellastra granularis_ on each experimental date.
* **desic_ms.R** includes the statistical analyses for all invertebrates on each experimental date.

### Infrared pictures

All infrared pictures are accessible in a file for each experimental date (YYMMDD) and each species (species code), and further organizes in a 'RAW PICTURES' file and a 'ANALYZED PICTURES' file.

e.g. '220209_ACAGAR' for infrared pictures of _Acanthochitona garnoti_ on 9 February 2022

Raw and analyzed infrared pictures are labelled for each time point (T+xx), each treatment (infested or clean) and each replicate (1, 2 or 3).

e.g. 'T+15_Clean beds_2' for the infrared picture taken 15 minutes after the beginning of the experiment, in non-infested mussel beds, in replicate n°2.



