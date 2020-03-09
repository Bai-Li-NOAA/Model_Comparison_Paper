#V3.30
#C data file created using the SS_writedat function in the R package r4ss
#C should work with SS version: 
#C file write time: 2020-02-28 13:29:35
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
    1	1	1	 161.575	0.00999975	#_1         
    2	1	1	 467.736	0.00999975	#_2         
    3	1	1	 759.117	0.00999975	#_3         
    4	1	1	 996.596	0.00999975	#_4         
    5	1	1	 769.721	0.00999975	#_5         
    6	1	1	1318.206	0.00999975	#_6         
    7	1	1	1282.219	0.00999975	#_7         
    8	1	1	2465.091	0.00999975	#_8         
    9	1	1	1330.456	0.00999975	#_9         
   10	1	1	1521.342	0.00999975	#_10        
   11	1	1	1633.265	0.00999975	#_11        
   12	1	1	1606.588	0.00999975	#_12        
   13	1	1	1105.012	0.00999975	#_13        
   14	1	1	1541.745	0.00999975	#_14        
   15	1	1	1529.437	0.00999975	#_15        
   16	1	1	1278.367	0.00999975	#_16        
   17	1	1	2232.736	0.00999975	#_17        
   18	1	1	1609.506	0.00999975	#_18        
   19	1	1	1465.826	0.00999975	#_19        
   20	1	1	1311.036	0.00999975	#_20        
   21	1	1	1628.566	0.00999975	#_21        
   22	1	1	1083.549	0.00999975	#_22        
   23	1	1	1617.873	0.00999975	#_23        
   24	1	1	1185.079	0.00999975	#_24        
   25	1	1	1120.194	0.00999975	#_25        
   26	1	1	 968.565	0.00999975	#_26        
   27	1	1	 929.371	0.00999975	#_27        
   28	1	1	1200.290	0.00999975	#_28        
   29	1	1	 875.566	0.00999975	#_29        
   30	1	1	1275.095	0.00999975	#_30        
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
    1	1	2	1.482760	0.198042	#_1         
    2	1	2	1.317051	0.198042	#_2         
    3	1	2	1.461116	0.198042	#_3         
    4	1	2	1.320284	0.198042	#_4         
    5	1	2	1.211159	0.198042	#_5         
    6	1	2	1.533895	0.198042	#_6         
    7	1	2	1.221455	0.198042	#_7         
    8	1	2	1.306322	0.198042	#_8         
    9	1	2	1.118493	0.198042	#_9         
   10	1	2	1.197817	0.198042	#_10        
   11	1	2	1.190712	0.198042	#_11        
   12	1	2	1.159013	0.198042	#_12        
   13	1	2	1.095386	0.198042	#_13        
   14	1	2	1.067502	0.198042	#_14        
   15	1	2	0.929270	0.198042	#_15        
   16	1	2	0.949470	0.198042	#_16        
   17	1	2	0.968862	0.198042	#_17        
   18	1	2	1.145569	0.198042	#_18        
   19	1	2	0.784747	0.198042	#_19        
   20	1	2	0.757849	0.198042	#_20        
   21	1	2	0.658952	0.198042	#_21        
   22	1	2	0.703679	0.198042	#_22        
   23	1	2	0.780789	0.198042	#_23        
   24	1	2	0.612883	0.198042	#_24        
   25	1	2	0.532742	0.198042	#_25        
   26	1	2	0.552710	0.198042	#_26        
   27	1	2	0.488881	0.198042	#_27        
   28	1	2	0.469033	0.198042	#_28        
   29	1	2	0.522832	0.198042	#_29        
   30	1	2	0.560446	0.198042	#_30        
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
    1	1	1	0	0	1	-1	-1	200	0.055	0.130	0.120	0.100	0.130	0.060	0.070	0.105	0.045	0.025	0.025	0.135	#_1         
    2	1	1	0	0	1	-1	-1	200	0.090	0.100	0.135	0.120	0.090	0.115	0.070	0.070	0.035	0.020	0.015	0.140	#_2         
    3	1	1	0	0	1	-1	-1	200	0.085	0.165	0.115	0.105	0.120	0.090	0.060	0.050	0.025	0.035	0.015	0.135	#_3         
    4	1	1	0	0	1	-1	-1	200	0.055	0.120	0.155	0.095	0.105	0.075	0.075	0.055	0.035	0.045	0.025	0.160	#_4         
    5	1	1	0	0	1	-1	-1	200	0.045	0.115	0.130	0.135	0.100	0.095	0.040	0.095	0.035	0.030	0.030	0.150	#_5         
    6	1	1	0	0	1	-1	-1	200	0.105	0.180	0.145	0.085	0.115	0.080	0.070	0.065	0.055	0.010	0.010	0.080	#_6         
    7	1	1	0	0	1	-1	-1	200	0.105	0.095	0.140	0.100	0.115	0.125	0.070	0.065	0.015	0.010	0.040	0.120	#_7         
    8	1	1	0	0	1	-1	-1	200	0.090	0.110	0.155	0.150	0.075	0.080	0.110	0.035	0.050	0.030	0.015	0.100	#_8         
    9	1	1	0	0	1	-1	-1	200	0.085	0.100	0.140	0.170	0.165	0.040	0.035	0.075	0.030	0.030	0.010	0.120	#_9         
   10	1	1	0	0	1	-1	-1	200	0.095	0.100	0.165	0.145	0.095	0.130	0.040	0.035	0.045	0.030	0.020	0.100	#_10        
   11	1	1	0	0	1	-1	-1	200	0.115	0.125	0.105	0.145	0.105	0.085	0.090	0.040	0.025	0.020	0.025	0.120	#_11        
   12	1	1	0	0	1	-1	-1	200	0.095	0.195	0.140	0.140	0.095	0.080	0.085	0.050	0.010	0.030	0.015	0.065	#_12        
   13	1	1	0	0	1	-1	-1	200	0.085	0.125	0.150	0.175	0.100	0.100	0.090	0.040	0.055	0.025	0.010	0.045	#_13        
   14	1	1	0	0	1	-1	-1	200	0.090	0.150	0.175	0.185	0.110	0.050	0.030	0.040	0.020	0.030	0.030	0.090	#_14        
   15	1	1	0	0	1	-1	-1	200	0.110	0.125	0.130	0.130	0.160	0.120	0.055	0.050	0.015	0.035	0.015	0.055	#_15        
   16	1	1	0	0	1	-1	-1	200	0.110	0.190	0.185	0.120	0.105	0.130	0.030	0.035	0.030	0.015	0.010	0.040	#_16        
   17	1	1	0	0	1	-1	-1	200	0.120	0.185	0.140	0.125	0.125	0.045	0.050	0.070	0.030	0.025	0.025	0.060	#_17        
   18	1	1	0	0	1	-1	-1	200	0.070	0.165	0.220	0.145	0.125	0.075	0.045	0.040	0.010	0.020	0.030	0.055	#_18        
   19	1	1	0	0	1	-1	-1	200	0.145	0.090	0.240	0.165	0.085	0.065	0.070	0.035	0.020	0.015	0.015	0.055	#_19        
   20	1	1	0	0	1	-1	-1	200	0.110	0.140	0.155	0.195	0.130	0.055	0.040	0.045	0.030	0.020	0.015	0.065	#_20        
   21	1	1	0	0	1	-1	-1	200	0.120	0.180	0.250	0.125	0.115	0.045	0.050	0.015	0.025	0.015	0.005	0.055	#_21        
   22	1	1	0	0	1	-1	-1	200	0.190	0.300	0.155	0.125	0.070	0.060	0.035	0.025	0.020	0.005	0.005	0.010	#_22        
   23	1	1	0	0	1	-1	-1	200	0.125	0.225	0.220	0.095	0.085	0.090	0.060	0.025	0.030	0.015	0.020	0.010	#_23        
   24	1	1	0	0	1	-1	-1	200	0.180	0.220	0.205	0.220	0.075	0.030	0.010	0.030	0.010	0.010	0.000	0.010	#_24        
   25	1	1	0	0	1	-1	-1	200	0.185	0.245	0.165	0.145	0.135	0.035	0.025	0.025	0.005	0.005	0.015	0.015	#_25        
   26	1	1	0	0	1	-1	-1	200	0.110	0.260	0.140	0.160	0.120	0.090	0.040	0.025	0.010	0.020	0.000	0.025	#_26        
   27	1	1	0	0	1	-1	-1	200	0.140	0.185	0.205	0.140	0.090	0.100	0.080	0.010	0.025	0.010	0.000	0.015	#_27        
   28	1	1	0	0	1	-1	-1	200	0.180	0.180	0.170	0.210	0.100	0.065	0.045	0.020	0.010	0.010	0.000	0.010	#_28        
   29	1	1	0	0	1	-1	-1	200	0.235	0.200	0.175	0.125	0.135	0.040	0.030	0.030	0.010	0.015	0.005	0.000	#_29        
   30	1	1	0	0	1	-1	-1	200	0.185	0.295	0.210	0.165	0.060	0.015	0.025	0.020	0.010	0.005	0.005	0.005	#_30        
    1	1	2	0	0	1	-1	-1	200	0.055	0.130	0.115	0.155	0.145	0.085	0.070	0.040	0.055	0.030	0.020	0.100	#_31        
    2	1	2	0	0	1	-1	-1	200	0.075	0.100	0.115	0.130	0.125	0.070	0.035	0.070	0.060	0.055	0.045	0.120	#_32        
    3	1	2	0	0	1	-1	-1	200	0.055	0.190	0.145	0.125	0.085	0.040	0.045	0.060	0.050	0.035	0.015	0.155	#_33        
    4	1	2	0	0	1	-1	-1	200	0.060	0.135	0.220	0.100	0.080	0.060	0.055	0.060	0.050	0.045	0.045	0.090	#_34        
    5	1	2	0	0	1	-1	-1	200	0.070	0.165	0.115	0.145	0.095	0.050	0.060	0.045	0.035	0.035	0.035	0.150	#_35        
    6	1	2	0	0	1	-1	-1	200	0.100	0.190	0.135	0.095	0.125	0.085	0.055	0.065	0.040	0.005	0.020	0.085	#_36        
    7	1	2	0	0	1	-1	-1	200	0.075	0.200	0.170	0.115	0.090	0.080	0.055	0.040	0.020	0.040	0.025	0.090	#_37        
    8	1	2	0	0	1	-1	-1	200	0.080	0.130	0.200	0.135	0.080	0.055	0.065	0.045	0.035	0.030	0.020	0.125	#_38        
    9	1	2	0	0	1	-1	-1	200	0.045	0.190	0.120	0.170	0.145	0.060	0.080	0.045	0.035	0.015	0.015	0.080	#_39        
   10	1	2	0	0	1	-1	-1	200	0.070	0.195	0.200	0.120	0.125	0.095	0.010	0.035	0.050	0.010	0.015	0.075	#_40        
   11	1	2	0	0	1	-1	-1	200	0.065	0.225	0.160	0.155	0.080	0.090	0.050	0.030	0.030	0.020	0.015	0.080	#_41        
   12	1	2	0	0	1	-1	-1	200	0.050	0.145	0.210	0.120	0.145	0.065	0.080	0.045	0.035	0.030	0.020	0.055	#_42        
   13	1	2	0	0	1	-1	-1	200	0.080	0.190	0.230	0.115	0.070	0.085	0.040	0.030	0.065	0.015	0.015	0.065	#_43        
   14	1	2	0	0	1	-1	-1	200	0.115	0.135	0.235	0.125	0.090	0.050	0.070	0.040	0.050	0.035	0.005	0.050	#_44        
   15	1	2	0	0	1	-1	-1	200	0.095	0.190	0.170	0.105	0.095	0.095	0.060	0.055	0.030	0.020	0.025	0.060	#_45        
   16	1	2	0	0	1	-1	-1	200	0.115	0.150	0.215	0.140	0.090	0.065	0.065	0.020	0.040	0.015	0.010	0.075	#_46        
   17	1	2	0	0	1	-1	-1	200	0.100	0.210	0.145	0.195	0.100	0.070	0.060	0.055	0.035	0.005	0.010	0.015	#_47        
   18	1	2	0	0	1	-1	-1	200	0.060	0.170	0.280	0.145	0.085	0.060	0.045	0.045	0.030	0.020	0.010	0.050	#_48        
   19	1	2	0	0	1	-1	-1	200	0.125	0.155	0.250	0.155	0.080	0.075	0.050	0.020	0.010	0.015	0.015	0.050	#_49        
   20	1	2	0	0	1	-1	-1	200	0.100	0.215	0.160	0.175	0.135	0.065	0.050	0.035	0.030	0.015	0.000	0.020	#_50        
   21	1	2	0	0	1	-1	-1	200	0.120	0.185	0.200	0.135	0.090	0.060	0.075	0.050	0.030	0.005	0.010	0.040	#_51        
   22	1	2	0	0	1	-1	-1	200	0.165	0.245	0.150	0.135	0.095	0.085	0.045	0.025	0.015	0.005	0.005	0.030	#_52        
   23	1	2	0	0	1	-1	-1	200	0.115	0.275	0.240	0.115	0.080	0.055	0.060	0.015	0.015	0.020	0.000	0.010	#_53        
   24	1	2	0	0	1	-1	-1	200	0.115	0.175	0.255	0.215	0.060	0.080	0.035	0.030	0.015	0.000	0.010	0.010	#_54        
   25	1	2	0	0	1	-1	-1	200	0.095	0.195	0.210	0.240	0.100	0.065	0.025	0.020	0.020	0.015	0.000	0.015	#_55        
   26	1	2	0	0	1	-1	-1	200	0.100	0.270	0.210	0.160	0.105	0.075	0.020	0.030	0.015	0.000	0.000	0.015	#_56        
   27	1	2	0	0	1	-1	-1	200	0.165	0.205	0.275	0.125	0.095	0.050	0.045	0.005	0.005	0.005	0.010	0.015	#_57        
   28	1	2	0	0	1	-1	-1	200	0.145	0.305	0.185	0.190	0.040	0.075	0.030	0.015	0.005	0.005	0.005	0.000	#_58        
   29	1	2	0	0	1	-1	-1	200	0.155	0.250	0.225	0.130	0.105	0.045	0.025	0.025	0.040	0.000	0.000	0.000	#_59        
   30	1	2	0	0	1	-1	-1	200	0.115	0.310	0.280	0.135	0.060	0.055	0.015	0.010	0.005	0.015	0.000	0.000	#_60        
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