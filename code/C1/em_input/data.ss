##############################################################################
#
#  Model Comparison Project, 2019
#
##############################################################################
#
1    #_styr
30   #_endyr
1    #_nseas
12   #_months/season
2    #_N_subseasons(even     number, minimum is      2)
1.00001    #_spawn_month
-1   #_Ngenders (use -1 for 1 sex setup with SSB multiplied by female_frac parameter)
12   #_Nages=accumulator     age
1    #_N_areas
2    #_Nfleets       (including      surveys)
#_fleet_type:   1=catch fleet;  2=bycatch       only    fleet;  3=survey;       4=ignore
#_survey_timing:        #NAME?  use     of      catch-at-age    to      override        the     month   value   associated      with    a       datum
#_fleet_area:   area    the     fleet/survey    operates        in
#_units of      catch:  1=bio;  2=num   (ignored        for     surveys;        their   units   read    later)
#_catch_mult:   0=no;   1=yes
#_rows  are     fleets
#_fleet_type,   timing, area,   units, need_catch_mult, fleetname
1	-1	1	1	0	FISHERY1
3	1	1	2	0	SURVEY1
#_Catch data:   yr,     seas,   fleet,  catch,  catch_se
#_catch_se:     standard        error   of      log(catch);     can     be      overridden      in      control file    with    detailed        F       input
#
#_FISHERY1
#
-999    1       1       15.33487587     0.049968792
1	1	1	153.3487587	0.049968792
2	1	1	448.2024748	0.049968792
3	1	1	756.0112637	0.049968792
4	1	1	1063.475811	0.049968792
5	1	1	788.3290716	0.049968792
6	1	1	1428.098352	0.049968792
7	1	1	1252.16254	0.049968792
8	1	1	2521.336042	0.049968792
9	1	1	1242.39058	0.049968792
10	1	1	1627.329383	0.049968792
11	1	1	1521.674455	0.049968792
12	1	1	1486.081103	0.049968792
13	1	1	1211.559293	0.049968792
14	1	1	1643.913512	0.049968792
15	1	1	1481.116603	0.049968792
16	1	1	1339.412999	0.049968792
17	1	1	2099.745949	0.049968792
18	1	1	1553.059078	0.049968792
19	1	1	1491.883249	0.049968792
20	1	1	1349.861095	0.049968792
21	1	1	1607.131512	0.049968792
22	1	1	1127.26385	0.049968792
23	1	1	1697.90532	0.049968792
24	1	1	1304.216722	0.049968792
25	1	1	1222.08571	0.049968792
26	1	1	922.4094778	0.049968792
27	1	1	1020.263545	0.049968792
28	1	1	1258.063348	0.049968792
29	1	1	879.982686	0.049968792
30	1	1	1283.866817	0.049968792
-9999	0	0	0	0
#
#_CPUE_and_surveyabundance_observations
#_Units:        0=numbers;      1=biomass;      2=F;    >=30    for     special types
#_Errtype:      -1=normal;      0=lognormal;    >0=T
#_Fleet Units   Errtype SD_Report
#
1       1      0       1       #_FISHERY1
2       0      0       1       #_SURVEY1
#
#_yr    month   fleet   obs     stderr
#
1	1	2	1.636678159	0.1980422
2	1	2	1.361204297	0.1980422
3	1	2	1.282149588	0.1980422
4	1	2	1.254426217	0.1980422
5	1	2	1.963527256	0.1980422
6	1	2	1.83913627	0.1980422
7	1	2	1.585324301	0.1980422
8	1	2	1.503790697	0.1980422
9	1	2	1.251297025	0.1980422
10	1	2	1.506597965	0.1980422
11	1	2	1.171635611	0.1980422
12	1	2	1.049418534	0.1980422
13	1	2	0.957554164	0.1980422
14	1	2	1.125046298	0.1980422
15	1	2	0.988373727	0.1980422
16	1	2	1.170114004	0.1980422
17	1	2	0.7143967	0.1980422
18	1	2	0.872899669	0.1980422
19	1	2	0.993783705	0.1980422
20	1	2	0.796413856	0.1980422
21	1	2	0.602314888	0.1980422
22	1	2	0.575182244	0.1980422
23	1	2	0.85276526	0.1980422
24	1	2	0.566478462	0.1980422
25	1	2	0.748408788	0.1980422
26	1	2	0.542710774	0.1980422
27	1	2	0.523787292	0.1980422
28	1	2	0.545070813	0.1980422
29	1	2	0.40603855	0.1980422
30	1	2	0.68177637	0.1980422
-9999	0	0	0	0
#       <TOADS> Discards        here.
0 #_N_fleets_with_discard
#_discard_units (1=same_as_catchunits(bio/num); 2=fraction;     3=numbers)
#_discard_errtype:      >0      for     DF      of      T-dist(read     CV      below); 0       for     normal  with    CV;     -1      for     normal  with    se;     -2      for     lognormal
#       note,   only    have    units   and     errtype for     fleets  with    discard
#_Fleet units   errtype
#
# Year Month  Fleet  DiscardRate.Sp.Wt.Wgting     CV
#
0 #_use meanbodysize_data (0/1)
#_COND_0 #_DF_for_meanbodysize_T-distribution_like
# note:  use positive partition value for mean body wt, negative partition for mean body length
#_yr month fleet part obs stderr
#  -9999 0 0 0 0 0 # terminator for mean body size data
#
# set up population length bin structure (note - irrelevant if not using size data and using empirical wtatage
2 # length bin method: 1=use databins; 2=generate from binwidth,min,max below; 3=read vector
1 # binwidth for population size comp
10 # minimum size in the population (lower edge of first bin and size at age 0.00)
89 # maximum size in the population (lower edge of last bin)
0 # use length composition data (0/1)
#_mintailcomp: upper and lower distribution for females and males separately are accumulated until exceeding this level.
#_addtocomp:  after accumulation of tails; this value added to all bins
#_males and females treated as combined gender below this bin number
#_compressbins: accumulate upper tail by this number of bins; acts simultaneous with mintailcomp; set=0 for no forced accumulation
#_Comp_Error:  0=multinomial, 1=dirichlet
#_Comp_Error2:  parm number  for dirichlet
#_mintailcomp_addtocomp_combM+F_CompressBins_CompError_ParmSelect_minsamplesize
# 
#
12      #_N_age_bins
#
1       2       3       4       5       6       7       8       9       10      11      12      
#
1       #_N_ageerror_definitions Put in the SDs.  These are the Newport agerrs  using these for the Commercial Trawl.
#
-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1 #Mean_Age
0	0	0	0	0	0	0	0	0	0	0	0	0 #SD
#
#
#_mintailcomp: upper and lower distribution for females and males separately are accumulated until exceeding this level.
#_addtocomp:    after   accumulation    of      tails;  this    value   added   to      all     bins
#_males and     females treated as      combined        gender  below   this    bin     number
#_compressbins: accumulate upper tail by this number of bins; acts simultaneous with mintailcomp; set=0 for no forced accumulation
#_Comp_Error:   0=multinomial   1=dirichlet
#_Comp_Error2:  parm    number  for     dirichlet
#_mintailcomp_addtocomp_combM+F_CompressBins_CompError_ParmSelect_minsamplesize
0	0.0001	1	0	0	0	1  #FISHERY1
0	0.0001	1	0	0	0	1  #SURVEY1
#
1       #_Lbin_method_for_Age_Data:     1=poplenbins;   2=datalenbins;  3=lengths
# sex codes: 0=combined; 1=use female only; 2=use male only; 3=use both as joint sexxlength distribution
# partition codes: (0=combined; 1=discard; 2=retained
#_Fleet 1: Fishing Fleet
#_yr    month   fleet   sex     part    ageerr  Lbin_lo Lbin_hi Nsamp   datavector(female-male) 
1	1	1	0	0	1	-1	-1	200	0.065	0.11	0.12	0.115	0.1	0.065	0.075	0.065	0.07	0.04	0.035	0.14
2	1	1	0	0	1	-1	-1	200	0.11	0.095	0.12	0.145	0.08	0.085	0.055	0.08	0.03	0.045	0.03	0.125
3	1	1	0	0	1	-1	-1	200	0.055	0.14	0.1	0.115	0.14	0.08	0.08	0.07	0.04	0.02	0.02	0.14
4	1	1	0	0	1	-1	-1	200	0.045	0.085	0.145	0.1	0.1	0.08	0.065	0.065	0.035	0.04	0.04	0.2
5	1	1	0	0	1	-1	-1	200	0.095	0.1	0.13	0.18	0.065	0.1	0.065	0.055	0.035	0.03	0.03	0.115
6	1	1	0	0	1	-1	-1	200	0.095	0.16	0.115	0.125	0.11	0.09	0.04	0.045	0.035	0.03	0.03	0.125
7	1	1	0	0	1	-1	-1	200	0.075	0.16	0.165	0.08	0.125	0.105	0.055	0.055	0.055	0.025	0.025	0.075
8	1	1	0	0	1	-1	-1	200	0.09	0.145	0.165	0.145	0.075	0.085	0.05	0.035	0.03	0.025	0.045	0.11
9	1	1	0	0	1	-1	-1	200	0.08	0.145	0.105	0.165	0.13	0.045	0.065	0.065	0.035	0.045	0.035	0.085
10	1	1	0	0	1	-1	-1	200	0.095	0.135	0.195	0.11	0.125	0.085	0.045	0.05	0.03	0.035	0.02	0.075
11	1	1	0	0	1	-1	-1	200	0.13	0.12	0.165	0.14	0.095	0.09	0.06	0.02	0.055	0.035	0.02	0.07
12	1	1	0	0	1	-1	-1	200	0.08	0.15	0.16	0.13	0.105	0.065	0.1	0.055	0.02	0.005	0.03	0.1
13	1	1	0	0	1	-1	-1	200	0.11	0.12	0.205	0.15	0.07	0.1	0.035	0.04	0.05	0.015	0.04	0.065
14	1	1	0	0	1	-1	-1	200	0.07	0.145	0.15	0.19	0.125	0.065	0.075	0.04	0.055	0.015	0.015	0.055
15	1	1	0	0	1	-1	-1	200	0.105	0.165	0.13	0.145	0.1	0.115	0.04	0.05	0.015	0.05	0.03	0.055
16	1	1	0	0	1	-1	-1	200	0.145	0.12	0.21	0.115	0.09	0.08	0.065	0.06	0.025	0.005	0.02	0.065
17	1	1	0	0	1	-1	-1	200	0.09	0.19	0.185	0.1	0.115	0.055	0.095	0.04	0.05	0.02	0.015	0.045
18	1	1	0	0	1	-1	-1	200	0.115	0.16	0.165	0.13	0.115	0.07	0.06	0.035	0.045	0.03	0.01	0.065
19	1	1	0	0	1	-1	-1	200	0.13	0.18	0.19	0.195	0.105	0.05	0.035	0.03	0.015	0.005	0.01	0.055
20	1	1	0	0	1	-1	-1	200	0.09	0.165	0.105	0.225	0.165	0.075	0.06	0.015	0.015	0.025	0.02	0.04
21	1	1	0	0	1	-1	-1	200	0.18	0.115	0.175	0.12	0.125	0.07	0.04	0.07	0.025	0.02	0.01	0.05
22	1	1	0	0	1	-1	-1	200	0.145	0.235	0.17	0.12	0.1	0.07	0.07	0.03	0.01	0.02	0	0.03
23	1	1	0	0	1	-1	-1	200	0.105	0.215	0.255	0.12	0.095	0.07	0.04	0.035	0.035	0.005	0.01	0.015
24	1	1	0	0	1	-1	-1	200	0.095	0.2	0.255	0.17	0.14	0.03	0.035	0.01	0.025	0.005	0.005	0.03
25	1	1	0	0	1	-1	-1	200	0.155	0.21	0.17	0.19	0.11	0.045	0.03	0.025	0.035	0.01	0.005	0.015
26	1	1	0	0	1	-1	-1	200	0.135	0.225	0.155	0.205	0.125	0.05	0.035	0.025	0.02	0.01	0.005	0.01
27	1	1	0	0	1	-1	-1	200	0.17	0.18	0.25	0.145	0.115	0.055	0.04	0.02	0.02	0	0.005	0
28	1	1	0	0	1	-1	-1	200	0.145	0.19	0.165	0.175	0.145	0.045	0.05	0.03	0.015	0.015	0.005	0.02
29	1	1	0	0	1	-1	-1	200	0.205	0.23	0.24	0.115	0.075	0.045	0.04	0.03	0.005	0	0	0.015
30	1	1	0	0	1	-1	-1	200	0.16	0.26	0.205	0.12	0.11	0.08	0.03	0.01	0.005	0.015	0	0.005

#_Fleet 2:Survey
1	1	2	0	0	1	-1	-1	200	0.085	0.115	0.175	0.13	0.085	0.055	0.075	0.065	0.05	0.025	0.04	0.1
2	1	2	0	0	1	-1	-1	200	0.065	0.145	0.13	0.135	0.1	0.08	0.08	0.06	0.03	0.045	0.015	0.115
3	1	2	0	0	1	-1	-1	200	0.09	0.16	0.155	0.095	0.095	0.07	0.045	0.055	0.04	0.06	0.02	0.115
4	1	2	0	0	1	-1	-1	200	0.065	0.095	0.145	0.12	0.13	0.07	0.07	0.075	0.03	0.035	0.02	0.145
5	1	2	0	0	1	-1	-1	200	0.06	0.135	0.145	0.16	0.125	0.045	0.06	0.04	0.04	0.05	0	0.14
6	1	2	0	0	1	-1	-1	200	0.055	0.2	0.135	0.105	0.155	0.075	0.035	0.02	0.055	0.015	0.01	0.14
7	1	2	0	0	1	-1	-1	200	0.065	0.14	0.21	0.08	0.08	0.095	0.07	0.07	0.04	0.05	0.01	0.09
8	1	2	0	0	1	-1	-1	200	0.07	0.145	0.175	0.135	0.09	0.08	0.09	0.055	0.035	0.02	0.025	0.08
9	1	2	0	0	1	-1	-1	200	0.065	0.195	0.12	0.14	0.13	0.06	0.065	0.045	0.035	0.025	0.03	0.09
10	1	2	0	0	1	-1	-1	200	0.06	0.19	0.15	0.19	0.105	0.115	0.025	0.045	0.03	0.015	0.015	0.06
11	1	2	0	0	1	-1	-1	200	0.07	0.17	0.115	0.24	0.09	0.085	0.06	0.015	0.035	0.045	0.02	0.055
12	1	2	0	0	1	-1	-1	200	0.095	0.215	0.155	0.135	0.115	0.035	0.045	0.06	0.015	0.025	0.03	0.075
13	1	2	0	0	1	-1	-1	200	0.07	0.16	0.16	0.175	0.055	0.105	0.035	0.075	0.07	0.01	0.025	0.06
14	1	2	0	0	1	-1	-1	200	0.055	0.17	0.2	0.195	0.125	0.04	0.055	0.035	0.045	0.035	0.01	0.035
15	1	2	0	0	1	-1	-1	200	0.09	0.21	0.14	0.13	0.105	0.09	0.055	0.06	0.025	0.02	0.025	0.05
16	1	2	0	0	1	-1	-1	200	0.085	0.235	0.24	0.13	0.075	0.075	0.045	0.025	0.015	0.02	0.045	0.01
17	1	2	0	0	1	-1	-1	200	0.08	0.15	0.205	0.15	0.155	0.1	0.065	0.03	0.03	0.01	0.005	0.02
18	1	2	0	0	1	-1	-1	200	0.09	0.205	0.195	0.165	0.115	0.075	0.07	0.04	0.005	0.005	0.01	0.025
19	1	2	0	0	1	-1	-1	200	0.09	0.18	0.27	0.16	0.075	0.075	0.02	0.035	0.01	0.01	0.04	0.035
20	1	2	0	0	1	-1	-1	200	0.075	0.175	0.175	0.185	0.1	0.125	0.05	0.04	0.02	0.005	0.02	0.03
21	1	2	0	0	1	-1	-1	200	0.115	0.15	0.24	0.125	0.135	0.065	0.05	0.04	0.005	0.015	0.025	0.035
22	1	2	0	0	1	-1	-1	200	0.11	0.26	0.195	0.14	0.09	0.06	0.045	0.035	0.005	0.01	0.02	0.03
23	1	2	0	0	1	-1	-1	200	0.09	0.28	0.235	0.12	0.07	0.045	0.07	0.03	0.02	0.015	0.01	0.015
24	1	2	0	0	1	-1	-1	200	0.11	0.17	0.295	0.215	0.06	0.055	0.02	0.035	0.015	0.01	0.01	0.005
25	1	2	0	0	1	-1	-1	200	0.09	0.21	0.225	0.2	0.11	0.065	0.03	0.03	0.015	0.01	0.005	0.01
26	1	2	0	0	1	-1	-1	200	0.125	0.205	0.17	0.175	0.115	0.1	0.03	0.02	0.015	0.02	0	0.025
27	1	2	0	0	1	-1	-1	200	0.125	0.215	0.265	0.13	0.095	0.07	0.03	0.02	0.015	0.02	0.005	0.01
28	1	2	0	0	1	-1	-1	200	0.13	0.235	0.225	0.16	0.07	0.05	0.02	0.045	0.02	0.01	0.01	0.025
29	1	2	0	0	1	-1	-1	200	0.205	0.275	0.205	0.09	0.1	0.05	0.03	0.025	0.015	0.005	0	0
30	1	2	0	0	1	-1	-1	200	0.085	0.32	0.255	0.175	0.05	0.045	0.015	0.02	0.01	0.015	0	0.01
-9999	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
#
#
#
#
0 #_Use_MeanSize-at-Age_obs (0/1)
#
#
0 #_N_environ_variables
#Year Variable Value
0 # N sizefreq methods to read
#
0 # do tags (0/1)
#
0 #    morphcomp data(0/1)
#  Nobs  Nmorphs  mincomp
#  yr  seas  type  partition  Nsamp  datavector_by_Nmorphs
#
0  #  Do dataread for selectivity priors(0/1)
# Yr  Seas  Fleet   Age/Size   Bin   selex_prior   prior_sd
# feature not yet implemented
#
999

ENDDATA
