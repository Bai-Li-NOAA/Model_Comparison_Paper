#V3.30
#C data file created using the SS_writedat function in the R package r4ss
#C should work with SS version: 
#C file write time: 2020-02-28 13:29:43
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
    1	1	1	 159.762	0.00999975	#_1         
    2	1	1	 468.033	0.00999975	#_2         
    3	1	1	 758.946	0.00999975	#_3         
    4	1	1	 990.130	0.00999975	#_4         
    5	1	1	 767.028	0.00999975	#_5         
    6	1	1	1326.167	0.00999975	#_6         
    7	1	1	1281.600	0.00999975	#_7         
    8	1	1	2476.124	0.00999975	#_8         
    9	1	1	1320.008	0.00999975	#_9         
   10	1	1	1547.230	0.00999975	#_10        
   11	1	1	1623.328	0.00999975	#_11        
   12	1	1	1618.088	0.00999975	#_12        
   13	1	1	1102.306	0.00999975	#_13        
   14	1	1	1519.661	0.00999975	#_14        
   15	1	1	1504.993	0.00999975	#_15        
   16	1	1	1283.154	0.00999975	#_16        
   17	1	1	2266.571	0.00999975	#_17        
   18	1	1	1603.724	0.00999975	#_18        
   19	1	1	1447.516	0.00999975	#_19        
   20	1	1	1321.588	0.00999975	#_20        
   21	1	1	1637.645	0.00999975	#_21        
   22	1	1	1086.240	0.00999975	#_22        
   23	1	1	1637.109	0.00999975	#_23        
   24	1	1	1212.816	0.00999975	#_24        
   25	1	1	1116.902	0.00999975	#_25        
   26	1	1	 968.195	0.00999975	#_26        
   27	1	1	 925.269	0.00999975	#_27        
   28	1	1	1206.002	0.00999975	#_28        
   29	1	1	 884.653	0.00999975	#_29        
   30	1	1	1265.093	0.00999975	#_30        
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
    1	1	2	1.513064	0.198042	#_1         
    2	1	2	1.479974	0.198042	#_2         
    3	1	2	1.563473	0.198042	#_3         
    4	1	2	1.489011	0.198042	#_4         
    5	1	2	1.375234	0.198042	#_5         
    6	1	2	1.502784	0.198042	#_6         
    7	1	2	1.508886	0.198042	#_7         
    8	1	2	1.329874	0.198042	#_8         
    9	1	2	1.349042	0.198042	#_9         
   10	1	2	1.355432	0.198042	#_10        
   11	1	2	1.151666	0.198042	#_11        
   12	1	2	1.051211	0.198042	#_12        
   13	1	2	0.985124	0.198042	#_13        
   14	1	2	1.005661	0.198042	#_14        
   15	1	2	1.066340	0.198042	#_15        
   16	1	2	0.797810	0.198042	#_16        
   17	1	2	0.921759	0.198042	#_17        
   18	1	2	0.867186	0.198042	#_18        
   19	1	2	0.905707	0.198042	#_19        
   20	1	2	0.758909	0.198042	#_20        
   21	1	2	0.724204	0.198042	#_21        
   22	1	2	0.713515	0.198042	#_22        
   23	1	2	0.615650	0.198042	#_23        
   24	1	2	0.610563	0.198042	#_24        
   25	1	2	0.572367	0.198042	#_25        
   26	1	2	0.500754	0.198042	#_26        
   27	1	2	0.565227	0.198042	#_27        
   28	1	2	0.550593	0.198042	#_28        
   29	1	2	0.539817	0.198042	#_29        
   30	1	2	0.574637	0.198042	#_30        
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
    1	1	1	0	0	1	-1	-1	200	0.050	0.085	0.160	0.135	0.130	0.090	0.065	0.040	0.055	0.015	0.015	0.160	#_1         
    2	1	1	0	0	1	-1	-1	200	0.085	0.100	0.090	0.145	0.115	0.095	0.065	0.050	0.035	0.050	0.030	0.140	#_2         
    3	1	1	0	0	1	-1	-1	200	0.030	0.140	0.130	0.110	0.090	0.105	0.060	0.060	0.060	0.030	0.040	0.145	#_3         
    4	1	1	0	0	1	-1	-1	200	0.035	0.095	0.105	0.130	0.085	0.095	0.065	0.075	0.050	0.060	0.030	0.175	#_4         
    5	1	1	0	0	1	-1	-1	200	0.130	0.115	0.135	0.155	0.100	0.055	0.055	0.050	0.030	0.050	0.030	0.095	#_5         
    6	1	1	0	0	1	-1	-1	200	0.100	0.190	0.120	0.110	0.095	0.085	0.045	0.040	0.060	0.020	0.020	0.115	#_6         
    7	1	1	0	0	1	-1	-1	200	0.075	0.100	0.175	0.110	0.090	0.085	0.075	0.040	0.035	0.060	0.030	0.125	#_7         
    8	1	1	0	0	1	-1	-1	200	0.080	0.115	0.150	0.185	0.080	0.050	0.095	0.080	0.030	0.015	0.010	0.110	#_8         
    9	1	1	0	0	1	-1	-1	200	0.110	0.185	0.110	0.165	0.135	0.070	0.035	0.040	0.030	0.015	0.010	0.095	#_9         
   10	1	1	0	0	1	-1	-1	200	0.085	0.125	0.175	0.105	0.120	0.090	0.055	0.015	0.065	0.020	0.010	0.135	#_10        
   11	1	1	0	0	1	-1	-1	200	0.075	0.115	0.140	0.205	0.055	0.065	0.130	0.045	0.025	0.040	0.020	0.085	#_11        
   12	1	1	0	0	1	-1	-1	200	0.115	0.180	0.160	0.135	0.060	0.070	0.085	0.060	0.035	0.020	0.005	0.075	#_12        
   13	1	1	0	0	1	-1	-1	200	0.075	0.155	0.135	0.205	0.115	0.080	0.045	0.050	0.040	0.035	0.005	0.060	#_13        
   14	1	1	0	0	1	-1	-1	200	0.105	0.130	0.110	0.190	0.140	0.095	0.090	0.030	0.020	0.030	0.005	0.055	#_14        
   15	1	1	0	0	1	-1	-1	200	0.115	0.155	0.145	0.155	0.100	0.055	0.035	0.065	0.050	0.035	0.055	0.035	#_15        
   16	1	1	0	0	1	-1	-1	200	0.120	0.110	0.195	0.165	0.105	0.070	0.100	0.040	0.025	0.010	0.015	0.045	#_16        
   17	1	1	0	0	1	-1	-1	200	0.105	0.235	0.155	0.170	0.090	0.060	0.035	0.045	0.015	0.040	0.015	0.035	#_17        
   18	1	1	0	0	1	-1	-1	200	0.095	0.190	0.205	0.115	0.140	0.065	0.065	0.025	0.025	0.015	0.010	0.050	#_18        
   19	1	1	0	0	1	-1	-1	200	0.115	0.105	0.205	0.205	0.120	0.075	0.025	0.035	0.035	0.015	0.015	0.050	#_19        
   20	1	1	0	0	1	-1	-1	200	0.060	0.195	0.200	0.175	0.145	0.065	0.055	0.020	0.035	0.020	0.010	0.020	#_20        
   21	1	1	0	0	1	-1	-1	200	0.190	0.115	0.175	0.135	0.130	0.075	0.040	0.045	0.025	0.035	0.010	0.025	#_21        
   22	1	1	0	0	1	-1	-1	200	0.130	0.220	0.155	0.165	0.105	0.100	0.040	0.020	0.010	0.020	0.015	0.020	#_22        
   23	1	1	0	0	1	-1	-1	200	0.125	0.250	0.250	0.080	0.135	0.065	0.050	0.020	0.005	0.005	0.000	0.015	#_23        
   24	1	1	0	0	1	-1	-1	200	0.175	0.140	0.295	0.140	0.075	0.060	0.045	0.025	0.020	0.005	0.010	0.010	#_24        
   25	1	1	0	0	1	-1	-1	200	0.230	0.145	0.170	0.175	0.105	0.055	0.035	0.015	0.015	0.015	0.005	0.035	#_25        
   26	1	1	0	0	1	-1	-1	200	0.125	0.230	0.175	0.145	0.130	0.085	0.045	0.015	0.005	0.010	0.015	0.020	#_26        
   27	1	1	0	0	1	-1	-1	200	0.140	0.170	0.255	0.130	0.125	0.080	0.050	0.030	0.010	0.000	0.000	0.010	#_27        
   28	1	1	0	0	1	-1	-1	200	0.205	0.205	0.160	0.215	0.085	0.030	0.030	0.025	0.020	0.010	0.000	0.015	#_28        
   29	1	1	0	0	1	-1	-1	200	0.185	0.235	0.155	0.095	0.180	0.050	0.035	0.020	0.010	0.005	0.010	0.020	#_29        
   30	1	1	0	0	1	-1	-1	200	0.135	0.250	0.255	0.105	0.060	0.090	0.055	0.010	0.015	0.015	0.000	0.010	#_30        
    1	1	2	0	0	1	-1	-1	200	0.045	0.140	0.175	0.125	0.100	0.065	0.070	0.045	0.030	0.020	0.015	0.170	#_31        
    2	1	2	0	0	1	-1	-1	200	0.090	0.115	0.130	0.120	0.130	0.055	0.075	0.040	0.050	0.025	0.020	0.150	#_32        
    3	1	2	0	0	1	-1	-1	200	0.070	0.165	0.160	0.135	0.120	0.060	0.065	0.035	0.020	0.025	0.030	0.115	#_33        
    4	1	2	0	0	1	-1	-1	200	0.070	0.110	0.195	0.095	0.085	0.070	0.095	0.085	0.030	0.045	0.010	0.110	#_34        
    5	1	2	0	0	1	-1	-1	200	0.075	0.140	0.120	0.180	0.065	0.115	0.090	0.045	0.040	0.015	0.010	0.105	#_35        
    6	1	2	0	0	1	-1	-1	200	0.110	0.210	0.145	0.055	0.125	0.060	0.055	0.045	0.045	0.035	0.035	0.080	#_36        
    7	1	2	0	0	1	-1	-1	200	0.060	0.145	0.205	0.105	0.100	0.115	0.055	0.025	0.030	0.020	0.015	0.125	#_37        
    8	1	2	0	0	1	-1	-1	200	0.120	0.150	0.160	0.105	0.055	0.065	0.115	0.015	0.030	0.030	0.015	0.140	#_38        
    9	1	2	0	0	1	-1	-1	200	0.095	0.170	0.180	0.130	0.115	0.055	0.050	0.035	0.045	0.020	0.025	0.080	#_39        
   10	1	2	0	0	1	-1	-1	200	0.140	0.155	0.160	0.130	0.120	0.080	0.030	0.025	0.035	0.015	0.030	0.080	#_40        
   11	1	2	0	0	1	-1	-1	200	0.085	0.170	0.140	0.160	0.075	0.120	0.060	0.045	0.035	0.035	0.020	0.055	#_41        
   12	1	2	0	0	1	-1	-1	200	0.055	0.180	0.200	0.125	0.130	0.065	0.070	0.035	0.025	0.030	0.025	0.060	#_42        
   13	1	2	0	0	1	-1	-1	200	0.080	0.185	0.210	0.135	0.070	0.080	0.080	0.040	0.040	0.005	0.015	0.060	#_43        
   14	1	2	0	0	1	-1	-1	200	0.075	0.145	0.145	0.145	0.145	0.105	0.045	0.050	0.030	0.045	0.015	0.055	#_44        
   15	1	2	0	0	1	-1	-1	200	0.045	0.175	0.205	0.125	0.145	0.080	0.030	0.045	0.025	0.035	0.020	0.070	#_45        
   16	1	2	0	0	1	-1	-1	200	0.120	0.150	0.150	0.115	0.095	0.135	0.075	0.045	0.025	0.030	0.010	0.050	#_46        
   17	1	2	0	0	1	-1	-1	200	0.105	0.200	0.150	0.135	0.100	0.055	0.080	0.040	0.015	0.025	0.020	0.075	#_47        
   18	1	2	0	0	1	-1	-1	200	0.120	0.220	0.240	0.080	0.065	0.095	0.035	0.020	0.030	0.010	0.025	0.060	#_48        
   19	1	2	0	0	1	-1	-1	200	0.105	0.170	0.275	0.155	0.105	0.050	0.040	0.015	0.025	0.025	0.005	0.030	#_49        
   20	1	2	0	0	1	-1	-1	200	0.080	0.270	0.145	0.180	0.105	0.045	0.025	0.050	0.035	0.015	0.010	0.040	#_50        
   21	1	2	0	0	1	-1	-1	200	0.165	0.160	0.185	0.105	0.095	0.115	0.055	0.050	0.025	0.015	0.005	0.025	#_51        
   22	1	2	0	0	1	-1	-1	200	0.155	0.285	0.145	0.150	0.080	0.040	0.065	0.025	0.025	0.005	0.000	0.025	#_52        
   23	1	2	0	0	1	-1	-1	200	0.130	0.260	0.265	0.085	0.115	0.045	0.040	0.025	0.010	0.010	0.015	0.000	#_53        
   24	1	2	0	0	1	-1	-1	200	0.100	0.180	0.245	0.175	0.115	0.060	0.040	0.030	0.030	0.000	0.005	0.020	#_54        
   25	1	2	0	0	1	-1	-1	200	0.125	0.210	0.230	0.180	0.095	0.050	0.045	0.030	0.010	0.015	0.000	0.010	#_55        
   26	1	2	0	0	1	-1	-1	200	0.095	0.310	0.185	0.190	0.085	0.065	0.035	0.005	0.005	0.005	0.015	0.005	#_56        
   27	1	2	0	0	1	-1	-1	200	0.075	0.220	0.270	0.175	0.100	0.040	0.040	0.020	0.015	0.010	0.010	0.025	#_57        
   28	1	2	0	0	1	-1	-1	200	0.120	0.225	0.235	0.245	0.060	0.035	0.035	0.005	0.025	0.005	0.000	0.010	#_58        
   29	1	2	0	0	1	-1	-1	200	0.150	0.260	0.240	0.130	0.080	0.040	0.025	0.035	0.010	0.015	0.005	0.010	#_59        
   30	1	2	0	0	1	-1	-1	200	0.105	0.365	0.210	0.110	0.065	0.030	0.050	0.020	0.015	0.010	0.015	0.005	#_60        
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