#V3.30
#C data file created using the SS_writedat function in the R package r4ss
#C should work with SS version: 
#C file write time: 2020-02-28 13:30:55
#
1 #_styr
30 #_endyr
1 #_nseas
12 #_months_per_seas
2 #_Nsubseasons
1.00001 #_spawn_month
-1 #_Nsexes
12 #_Nages
1 #_Nareas
2 #_Nfleets
#_fleetinfo
#_type	surveytiming	area	units	need_catch_mult	fleetname
1	-1	1	1	0	FISHERY1	#_1
3	 1	1	2	0	SURVEY1 	#_2
#_Catch data
#_year	season	fleet	catch	catch_se
    1	1	1	 159.997	0.00999975	#_1         
    2	1	1	 467.376	0.00999975	#_2         
    3	1	1	 753.296	0.00999975	#_3         
    4	1	1	 996.263	0.00999975	#_4         
    5	1	1	 767.041	0.00999975	#_5         
    6	1	1	1309.001	0.00999975	#_6         
    7	1	1	1285.032	0.00999975	#_7         
    8	1	1	2493.337	0.00999975	#_8         
    9	1	1	1330.870	0.00999975	#_9         
   10	1	1	1521.980	0.00999975	#_10        
   11	1	1	1616.846	0.00999975	#_11        
   12	1	1	1618.984	0.00999975	#_12        
   13	1	1	1110.641	0.00999975	#_13        
   14	1	1	1525.800	0.00999975	#_14        
   15	1	1	1499.795	0.00999975	#_15        
   16	1	1	1286.021	0.00999975	#_16        
   17	1	1	2265.296	0.00999975	#_17        
   18	1	1	1600.446	0.00999975	#_18        
   19	1	1	1455.780	0.00999975	#_19        
   20	1	1	1320.323	0.00999975	#_20        
   21	1	1	1620.750	0.00999975	#_21        
   22	1	1	1078.585	0.00999975	#_22        
   23	1	1	1636.124	0.00999975	#_23        
   24	1	1	1203.726	0.00999975	#_24        
   25	1	1	1114.864	0.00999975	#_25        
   26	1	1	 973.763	0.00999975	#_26        
   27	1	1	 927.564	0.00999975	#_27        
   28	1	1	1212.746	0.00999975	#_28        
   29	1	1	 872.856	0.00999975	#_29        
   30	1	1	1278.179	0.00999975	#_30        
-9999	0	0	   0.000	0.00000000	#_terminator
#_CPUE_and_surveyabundance_observations
#_Units:  0=numbers; 1=biomass; 2=F; >=30 for special types
#_Errtype:  -1=normal; 0=lognormal; >0=T
#_SD_Report: 0=no sdreport; 1=enable sdreport
#_Fleet	Units	Errtype	SD_Report
1	1	0	1	#_FISHERY1
2	0	0	1	#_SURVEY1 
#
#_CPUE_data
#_year	seas	index	obs	se_log
    1	1	2	1.338711	0.198042	#_1         
    2	1	2	1.822476	0.198042	#_2         
    3	1	2	1.399601	0.198042	#_3         
    4	1	2	1.415053	0.198042	#_4         
    5	1	2	1.231354	0.198042	#_5         
    6	1	2	1.571295	0.198042	#_6         
    7	1	2	1.614848	0.198042	#_7         
    8	1	2	1.294297	0.198042	#_8         
    9	1	2	1.160441	0.198042	#_9         
   10	1	2	1.067379	0.198042	#_10        
   11	1	2	1.084138	0.198042	#_11        
   12	1	2	0.968867	0.198042	#_12        
   13	1	2	0.853934	0.198042	#_13        
   14	1	2	1.079735	0.198042	#_14        
   15	1	2	1.266543	0.198042	#_15        
   16	1	2	0.918003	0.198042	#_16        
   17	1	2	0.928080	0.198042	#_17        
   18	1	2	0.840902	0.198042	#_18        
   19	1	2	0.736248	0.198042	#_19        
   20	1	2	0.707781	0.198042	#_20        
   21	1	2	0.828417	0.198042	#_21        
   22	1	2	0.580373	0.198042	#_22        
   23	1	2	0.750720	0.198042	#_23        
   24	1	2	0.616462	0.198042	#_24        
   25	1	2	0.578659	0.198042	#_25        
   26	1	2	0.544929	0.198042	#_26        
   27	1	2	0.472636	0.198042	#_27        
   28	1	2	0.578869	0.198042	#_28        
   29	1	2	0.503233	0.198042	#_29        
   30	1	2	0.576452	0.198042	#_30        
-9999	0	0	0.000000	0.000000	#_terminator
0 #_N_discard_fleets
#_discard_units (1=same_as_catchunits(bio/num); 2=fraction; 3=numbers)
#_discard_errtype:  >0 for DF of T-dist(read CV below); 0 for normal with CV; -1 for normal with se; -2 for lognormal
#
#_discard_fleet_info
#
#_discard_data
#
#_meanbodywt
0 #_use_meanbodywt
 #_DF_for_meanbodywt_T-distribution_like
#
#_population_length_bins
2 # length bin method: 1=use databins; 2=generate from binwidth,min,max below; 3=read vector
1 # binwidth for population size comp
10 # minimum size in the population (lower edge of first bin and size at age 0.00)
89 # maximum size in the population (lower edge of last bin)
0 #_use_lencomp
12 #_N_agebins
#
#_agebin_vector
1 2 3 4 5 6 7 8 9 10 11 12 #_agebin_vector
#
#_ageing_error
1 #_N_ageerror_definitions
#_age0	age1	age2	age3	age4	age5	age6	age7	age8	age9	age10	age11	age12
-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	#_1
 0	 0	 0	 0	 0	 0	 0	 0	 0	 0	 0	 0	 0	#_2
#
#_age_info
#_mintailcomp	addtocomp	combine_M_F	CompressBins	CompError	ParmSelect	minsamplesize
0	1e-04	1	0	0	0	1	#_FISHERY1
0	1e-04	1	0	0	0	1	#_SURVEY1 
1 #_Lbin_method: 1=poplenbins; 2=datalenbins; 3=lengths
 #_combine males into females at or below this bin number
#_Yr	Seas	FltSvy	Gender	Part	Ageerr	Lbin_lo	Lbin_hi	Nsamp	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
    1	1	1	0	0	1	-1	-1	200	0.045	0.065	0.155	0.100	0.095	0.145	0.100	0.045	0.050	0.040	0.010	0.150	#_1         
    2	1	1	0	0	1	-1	-1	200	0.090	0.130	0.105	0.080	0.085	0.095	0.080	0.080	0.040	0.040	0.030	0.145	#_2         
    3	1	1	0	0	1	-1	-1	200	0.060	0.130	0.135	0.135	0.115	0.075	0.050	0.055	0.035	0.055	0.020	0.135	#_3         
    4	1	1	0	0	1	-1	-1	200	0.085	0.085	0.140	0.110	0.090	0.060	0.065	0.070	0.040	0.045	0.055	0.155	#_4         
    5	1	1	0	0	1	-1	-1	200	0.095	0.075	0.145	0.175	0.050	0.075	0.090	0.080	0.045	0.045	0.025	0.100	#_5         
    6	1	1	0	0	1	-1	-1	200	0.095	0.150	0.095	0.100	0.155	0.075	0.055	0.045	0.050	0.030	0.040	0.110	#_6         
    7	1	1	0	0	1	-1	-1	200	0.090	0.125	0.150	0.125	0.100	0.130	0.040	0.050	0.035	0.040	0.030	0.085	#_7         
    8	1	1	0	0	1	-1	-1	200	0.125	0.110	0.185	0.140	0.080	0.105	0.070	0.035	0.010	0.030	0.030	0.080	#_8         
    9	1	1	0	0	1	-1	-1	200	0.050	0.140	0.120	0.145	0.135	0.065	0.055	0.070	0.030	0.020	0.020	0.150	#_9         
   10	1	1	0	0	1	-1	-1	200	0.095	0.105	0.190	0.120	0.140	0.085	0.065	0.050	0.040	0.025	0.005	0.080	#_10        
   11	1	1	0	0	1	-1	-1	200	0.055	0.170	0.145	0.170	0.125	0.100	0.060	0.030	0.030	0.015	0.020	0.080	#_11        
   12	1	1	0	0	1	-1	-1	200	0.080	0.150	0.115	0.115	0.125	0.105	0.080	0.080	0.010	0.040	0.030	0.070	#_12        
   13	1	1	0	0	1	-1	-1	200	0.120	0.100	0.160	0.150	0.075	0.125	0.055	0.045	0.050	0.015	0.020	0.085	#_13        
   14	1	1	0	0	1	-1	-1	200	0.120	0.160	0.095	0.165	0.120	0.060	0.090	0.035	0.030	0.035	0.015	0.075	#_14        
   15	1	1	0	0	1	-1	-1	200	0.120	0.135	0.125	0.110	0.165	0.085	0.050	0.055	0.020	0.020	0.030	0.085	#_15        
   16	1	1	0	0	1	-1	-1	200	0.090	0.130	0.165	0.115	0.115	0.125	0.080	0.050	0.040	0.025	0.020	0.045	#_16        
   17	1	1	0	0	1	-1	-1	200	0.135	0.130	0.170	0.195	0.065	0.070	0.075	0.065	0.040	0.015	0.000	0.040	#_17        
   18	1	1	0	0	1	-1	-1	200	0.060	0.195	0.210	0.110	0.135	0.070	0.025	0.045	0.035	0.040	0.035	0.040	#_18        
   19	1	1	0	0	1	-1	-1	200	0.105	0.125	0.230	0.170	0.100	0.060	0.060	0.040	0.020	0.025	0.020	0.045	#_19        
   20	1	1	0	0	1	-1	-1	200	0.090	0.185	0.180	0.145	0.140	0.075	0.055	0.010	0.030	0.020	0.005	0.065	#_20        
   21	1	1	0	0	1	-1	-1	200	0.210	0.150	0.145	0.100	0.130	0.110	0.060	0.035	0.020	0.010	0.005	0.025	#_21        
   22	1	1	0	0	1	-1	-1	200	0.095	0.250	0.190	0.150	0.075	0.095	0.030	0.040	0.025	0.010	0.000	0.040	#_22        
   23	1	1	0	0	1	-1	-1	200	0.115	0.215	0.255	0.145	0.080	0.045	0.060	0.035	0.010	0.020	0.000	0.020	#_23        
   24	1	1	0	0	1	-1	-1	200	0.125	0.185	0.285	0.150	0.080	0.050	0.025	0.025	0.040	0.010	0.005	0.020	#_24        
   25	1	1	0	0	1	-1	-1	200	0.245	0.150	0.160	0.175	0.135	0.065	0.025	0.020	0.020	0.000	0.000	0.005	#_25        
   26	1	1	0	0	1	-1	-1	200	0.110	0.240	0.195	0.160	0.110	0.080	0.030	0.025	0.020	0.015	0.005	0.010	#_26        
   27	1	1	0	0	1	-1	-1	200	0.170	0.180	0.235	0.180	0.080	0.040	0.055	0.015	0.030	0.005	0.010	0.000	#_27        
   28	1	1	0	0	1	-1	-1	200	0.175	0.190	0.155	0.175	0.095	0.085	0.030	0.025	0.015	0.015	0.005	0.035	#_28        
   29	1	1	0	0	1	-1	-1	200	0.215	0.245	0.185	0.100	0.100	0.060	0.030	0.030	0.020	0.005	0.005	0.005	#_29        
   30	1	1	0	0	1	-1	-1	200	0.170	0.320	0.150	0.120	0.090	0.075	0.030	0.020	0.005	0.015	0.000	0.005	#_30        
    1	1	2	0	0	1	-1	-1	200	0.055	0.105	0.185	0.090	0.075	0.090	0.065	0.050	0.065	0.035	0.065	0.120	#_31        
    2	1	2	0	0	1	-1	-1	200	0.095	0.150	0.110	0.145	0.115	0.040	0.080	0.075	0.035	0.035	0.040	0.080	#_32        
    3	1	2	0	0	1	-1	-1	200	0.055	0.200	0.110	0.115	0.115	0.060	0.080	0.040	0.040	0.045	0.015	0.125	#_33        
    4	1	2	0	0	1	-1	-1	200	0.065	0.120	0.145	0.105	0.110	0.120	0.045	0.065	0.040	0.040	0.045	0.100	#_34        
    5	1	2	0	0	1	-1	-1	200	0.100	0.095	0.135	0.155	0.085	0.095	0.070	0.055	0.035	0.025	0.035	0.115	#_35        
    6	1	2	0	0	1	-1	-1	200	0.045	0.155	0.120	0.115	0.185	0.050	0.085	0.025	0.060	0.055	0.020	0.085	#_36        
    7	1	2	0	0	1	-1	-1	200	0.065	0.215	0.155	0.120	0.085	0.060	0.070	0.050	0.035	0.030	0.035	0.080	#_37        
    8	1	2	0	0	1	-1	-1	200	0.060	0.160	0.190	0.100	0.085	0.090	0.055	0.065	0.040	0.020	0.040	0.095	#_38        
    9	1	2	0	0	1	-1	-1	200	0.055	0.225	0.145	0.155	0.140	0.055	0.040	0.045	0.040	0.015	0.015	0.070	#_39        
   10	1	2	0	0	1	-1	-1	200	0.065	0.170	0.245	0.115	0.120	0.065	0.050	0.040	0.030	0.025	0.015	0.060	#_40        
   11	1	2	0	0	1	-1	-1	200	0.080	0.165	0.200	0.135	0.090	0.095	0.075	0.020	0.045	0.035	0.000	0.060	#_41        
   12	1	2	0	0	1	-1	-1	200	0.055	0.225	0.175	0.140	0.110	0.040	0.075	0.060	0.025	0.035	0.010	0.050	#_42        
   13	1	2	0	0	1	-1	-1	200	0.070	0.160	0.205	0.110	0.100	0.095	0.055	0.050	0.040	0.015	0.005	0.095	#_43        
   14	1	2	0	0	1	-1	-1	200	0.100	0.215	0.125	0.115	0.145	0.080	0.065	0.045	0.030	0.000	0.015	0.065	#_44        
   15	1	2	0	0	1	-1	-1	200	0.075	0.220	0.140	0.085	0.125	0.070	0.075	0.045	0.010	0.035	0.045	0.075	#_45        
   16	1	2	0	0	1	-1	-1	200	0.070	0.180	0.220	0.125	0.110	0.095	0.075	0.045	0.015	0.010	0.005	0.050	#_46        
   17	1	2	0	0	1	-1	-1	200	0.090	0.260	0.160	0.160	0.090	0.085	0.040	0.035	0.010	0.035	0.000	0.035	#_47        
   18	1	2	0	0	1	-1	-1	200	0.070	0.245	0.260	0.125	0.075	0.050	0.045	0.010	0.040	0.025	0.020	0.035	#_48        
   19	1	2	0	0	1	-1	-1	200	0.080	0.165	0.250	0.160	0.080	0.080	0.035	0.045	0.030	0.015	0.025	0.035	#_49        
   20	1	2	0	0	1	-1	-1	200	0.060	0.240	0.195	0.160	0.100	0.060	0.070	0.035	0.015	0.015	0.015	0.035	#_50        
   21	1	2	0	0	1	-1	-1	200	0.095	0.145	0.205	0.170	0.160	0.080	0.040	0.025	0.020	0.015	0.010	0.035	#_51        
   22	1	2	0	0	1	-1	-1	200	0.145	0.275	0.215	0.125	0.065	0.060	0.035	0.020	0.030	0.000	0.010	0.020	#_52        
   23	1	2	0	0	1	-1	-1	200	0.110	0.285	0.265	0.100	0.080	0.055	0.030	0.010	0.015	0.005	0.015	0.030	#_53        
   24	1	2	0	0	1	-1	-1	200	0.115	0.205	0.270	0.210	0.045	0.050	0.020	0.040	0.015	0.010	0.005	0.015	#_54        
   25	1	2	0	0	1	-1	-1	200	0.165	0.190	0.170	0.170	0.135	0.030	0.080	0.020	0.010	0.020	0.005	0.005	#_55        
   26	1	2	0	0	1	-1	-1	200	0.060	0.270	0.210	0.145	0.120	0.100	0.030	0.020	0.005	0.020	0.000	0.020	#_56        
   27	1	2	0	0	1	-1	-1	200	0.120	0.140	0.335	0.145	0.095	0.090	0.030	0.010	0.020	0.000	0.010	0.005	#_57        
   28	1	2	0	0	1	-1	-1	200	0.160	0.300	0.205	0.145	0.070	0.040	0.025	0.020	0.005	0.010	0.005	0.015	#_58        
   29	1	2	0	0	1	-1	-1	200	0.175	0.290	0.210	0.115	0.085	0.050	0.035	0.020	0.020	0.000	0.000	0.000	#_59        
   30	1	2	0	0	1	-1	-1	200	0.070	0.345	0.195	0.115	0.080	0.075	0.025	0.045	0.025	0.010	0.005	0.010	#_60        
-9999	0	0	0	0	0	 0	 0	  0	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	#_terminator
#
#_MeanSize_at_Age_obs
0 #_use_MeanSize_at_Age_obs
0 #_N_environ_variables
0 #_N_sizefreq_methods
0 #_do_tags
0 #_morphcomp_data
0 #_use_selectivity_priors
#
999