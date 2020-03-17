# Requirement Inputs and Default Parameters from `findSpanin.pl` script that has now been obsoleted
## INPUT : 
* Genbank File
## PARAMETERS :
* strand : +, -, both
* start codons: ATG, GTG and TTG
* isp_min: minimal length of the ORF, measured in AAs --> default 60
* isp_nterm_mindist: minimal distance to first AA of TMD, measured in AA --> default 10
* isp_nterm_maxdist: maximum distance to first AA of TMD, measured in AA --> default 30
* osp_min: minimal length of the ORF, measured in AAs --> default 30
* osp_signal_mindist: minimal distance to first AA of Lipobox, measured in AA --> default 10
* osp_signal_maxdist: maximum distance to first AA of Lipobox, measured in AA --> default 30
* Use LipoRy (?)
* max_isp_osp_distance: maximum distance between the END of the isp, and the beginning of the osp, measured in AA --> default 10