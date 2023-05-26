Creator "JGraphT GML Exporter"
Version 1
graph
[
	label ""
	directed 1
	node
	[
		id 1
		label "R_GLNS"
		Reversible "false"
		Name "glutamine synthetase"
	]
	node
	[
		id 2
		label "R_GLUPRT"
		Reversible "false"
		Name "glutamine phosphoribosyldiphosphate amidotransferase"
	]
	node
	[
		id 3
		label "R_GF6PTA"
		Reversible "false"
		Name "glutamine-fructose-6-phosphate transaminase"
	]
	node
	[
		id 4
		label "R_G6PDA"
		Reversible "false"
		Name "glucosamine-6-phosphate deaminase"
	]
	node
	[
		id 5
		label "R_PRAGSr"
		Reversible "true"
		Name "phosphoribosylglycinamide synthase"
	]
	node
	[
		id 6
		label "R_AIRCr_PRASCS"
		Reversible "true"
		Name "phosphoribosylaminoimidazole carboxylase / phosphoribosylaminoimidazolesuccinocarboxamide synthase"
	]
	node
	[
		id 7
		label "R_ADSL2"
		Reversible "false"
		Name "adenylosuccinate lyase"
	]
	node
	[
		id 8
		label "R_ASPCTr"
		Reversible "true"
		Name "aspartate carbamoyltransferase (reversible)"
	]
	node
	[
		id 9
		label "R_PHACCOAGLNAC"
		Reversible "false"
		Name "Phenylacetyl-CoA:L-glutamine alpha-N-phenylacetyltransferase"
	]
	node
	[
		id 10
		label "R_PACCOAL"
		Reversible "false"
		Name "phenylacetate-CoA ligase"
	]
	node
	[
		id 11
		label "R_r0545"
		Reversible "false"
		Name "Phenylacetaldehyde:NAD+ oxidoreductase Phenylalanine metabolism / Styrene degradation EC:1.2.1.5 EC:1.2.1.39"
	]
	node
	[
		id 12
		label "R_GARFT"
		Reversible "true"
		Name "phosphoribosylglycinamide formyltransferase"
	]
	node
	[
		id 13
		label "R_PRFGS"
		Reversible "false"
		Name "phosphoribosylformylglycinamidine synthase"
	]
	node
	[
		id 14
		label "R_RE2514E"
		Reversible "true"
		Name "RE2514"
	]
	node
	[
		id 15
		label "R_RE2513E"
		Reversible "false"
		Name "RE2513"
	]
	node
	[
		id 16
		label "R_GLNB0AT3tc"
		Reversible "false"
		Name "glutamine transport by B0AT3"
	]
	node
	[
		id 17
		label "R_r0666"
		Reversible "true"
		Name "2-(Formamido)-N1-(5-phosphoribosyl)acetamidine cyclo-ligase (ADP-forming) Purine metabolism EC:6.3.3.1"
	]
	node
	[
		id 18
		label "R_ACGAMtly"
		Reversible "false"
		Name "N-acetyl-glucosamine lysosomal efflux"
	]
	node
	[
		id 19
		label "R_ACGAMK"
		Reversible "false"
		Name "N-acetylglucosamine kinase"
	]
	node
	[
		id 20
		label "R_AGDC"
		Reversible "false"
		Name "N-acetylglucosamine-6-phosphate deacetylase"
	]
	node
	[
		id 21
		label "R_CBPS"
		Reversible "false"
		Name "carbamoyl-phosphate synthase (glutamine-hydrolysing)"
	]
	edge
	[
		source 1
		target 2
		label "M_gln_L_c"
		Formula "C5H10N2O3"
		Name "L-glutamine"
	]
	edge
	[
		source 3
		target 1
		label "M_glu_L_c"
		Formula "C5H8NO4"
		Name "L-glutamate(1-)"
	]
	edge
	[
		source 4
		target 3
		label "M_f6p_c"
		Formula "C6H11O9P"
		Name "D-fructose 6-phosphate(2-)"
	]
	edge
	[
		source 2
		target 5
		label "M_pram_c"
		Formula "C5H11NO7P"
		Name "5-phospho-beta-D-ribosylaminium(1-)"
	]
	edge
	[
		source 6
		target 7
		label "M_25aics_c"
		Formula "C13H15N4O12P"
		Name "SAICAR(4-)"
	]
	edge
	[
		source 6
		target 8
		label "M_asp_L_c"
		Formula "C4H6NO4"
		Name "L-aspartate(1-)"
	]
	edge
	[
		source 1
		target 9
		label "M_gln_L_c"
		Formula "C5H10N2O3"
		Name "L-glutamine"
	]
	edge
	[
		source 10
		target 9
		label "M_phaccoa_c"
		Formula "C29H38N7O17P3S"
		Name "phenylacetyl-CoA(4-)"
	]
	edge
	[
		source 11
		target 10
		label "M_pac_c"
		Formula "C8H7O2"
		Name "phenylacetate"
	]
	edge
	[
		source 12
		target 13
		label "M_fgam_c"
		Formula "C8H13N2O9P"
		Name "N(2)-formyl-N(1)-(5-phospho-D-ribosyl)glycinamide(2-)"
	]
	edge
	[
		source 14
		target 15
		label "M_CE4633_e"
		Formula "ClHO"
		Name "hypochlorous acid"
	]
	edge
	[
		source 15
		target 16
		label "M_cl_e"
		Formula "Cl"
		Name "chloride"
	]
	edge
	[
		source 16
		target 2
		label "M_gln_L_c"
		Formula "C5H10N2O3"
		Name "L-glutamine"
	]
	edge
	[
		source 5
		target 12
		label "M_gar_c"
		Formula "C7H14N2O8P"
		Name "N(1)-(5-phospho-D-ribosyl)glycinamide(1-)"
	]
	edge
	[
		source 13
		target 17
		label "M_fpram_c"
		Formula "C8H15N3O8P"
		Name "2-formamido-N(1)-(5-O-phosphonato-D-ribosyl)acetamidine"
	]
	edge
	[
		source 18
		target 19
		label "M_acgam_c"
		Formula "C8H15NO6"
		Name "aldehydo-N-acetyl-D-glucosamine"
	]
	edge
	[
		source 19
		target 20
		label "M_acgam6p_c"
		Formula "C8H14NO9P"
		Name "N-acetyl-D-glucosamine 6-phosphate(2-)"
	]
	edge
	[
		source 20
		target 4
		label "M_gam6p_c"
		Formula "C6H13NO8P"
		Name "D-glucosamine 6-phosphate"
	]
	edge
	[
		source 21
		target 8
		label "M_cbp_c"
		Formula "CH2NO5P"
		Name "carbamoyl phosphate(2-)"
	]
	edge
	[
		source 6
		target 17
		label "M_air_c"
		Formula "C8H12N3O7P"
		Name "5-amino-1-(5-phosphonato-D-ribosyl)imidazol-3-ium"
	]
]
