#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
  #include "admodel.h"          // Include AD class definitions
  #include "admb2r.cpp"         // Include R-compatible output functions (needs preceding)
  #include <time.h>
	time_t start,finish;
	long hour,minute,second;	
	double elapsed_time;
	
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <BAM-Sim.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
cout << "Starting Beaufort Assessment Model" << endl;
cout << endl;
cout << "                BAM!" << endl;
cout << endl;
  styr.allocate("styr");
  endyr.allocate("endyr");
  styr_rec_dev.allocate("styr_rec_dev");
  endyr_rec_dev.allocate("endyr_rec_dev");
  endyr_rec_phase1.allocate("endyr_rec_phase1");
  endyr_rec_phase2.allocate("endyr_rec_phase2");
  endyr_selex_phase1.allocate("endyr_selex_phase1");
    nyrs=endyr-styr+1.;
    nyrs_rec=endyr_rec_dev-styr_rec_dev+1.;
  nages.allocate("nages");
  agebins.allocate(1,nages,"agebins");
  nages_agec.allocate("nages_agec");
  agebins_agec.allocate(1,nages_agec,"agebins_agec");
  max_F_spr_msy.allocate("max_F_spr_msy");
  n_iter_spr.allocate("n_iter_spr");
		n_iter_msy=n_iter_spr; 
  selpar_n_yrs_wgted.allocate("selpar_n_yrs_wgted");
  set_BiasCor.allocate("set_BiasCor");
  nyr_survey1_cpue.allocate("nyr_survey1_cpue");
  yrs_survey1_cpue.allocate(1,nyr_survey1_cpue,"yrs_survey1_cpue");
  obs_survey1_cpue.allocate(1,nyr_survey1_cpue,"obs_survey1_cpue");
  survey1_cpue_cv.allocate(1,nyr_survey1_cpue,"survey1_cpue_cv");
  nyr_survey1_agec.allocate("nyr_survey1_agec");
  yrs_survey1_agec.allocate(1,nyr_survey1_agec,"yrs_survey1_agec");
  nsamp_survey1_agec.allocate(1,nyr_survey1_agec,"nsamp_survey1_agec");
  nfish_survey1_agec.allocate(1,nyr_survey1_agec,"nfish_survey1_agec");
  obs_survey1_agec.allocate(1,nyr_survey1_agec,1,nages_agec,"obs_survey1_agec");
  nyr_survey2_cpue.allocate("nyr_survey2_cpue");
  yrs_survey2_cpue.allocate(1,nyr_survey2_cpue,"yrs_survey2_cpue");
  obs_survey2_cpue.allocate(1,nyr_survey2_cpue,"obs_survey2_cpue");
  survey2_cpue_cv.allocate(1,nyr_survey2_cpue,"survey2_cpue_cv");
  nyr_survey2_agec.allocate("nyr_survey2_agec");
  yrs_survey2_agec.allocate(1,nyr_survey2_agec,"yrs_survey2_agec");
  nsamp_survey2_agec.allocate(1,nyr_survey2_agec,"nsamp_survey2_agec");
  nfish_survey2_agec.allocate(1,nyr_survey2_agec,"nfish_survey2_agec");
  obs_survey2_agec.allocate(1,nyr_survey2_agec,1,nages_agec,"obs_survey2_agec");
  styr_fleet1_L.allocate("styr_fleet1_L");
  endyr_fleet1_L.allocate("endyr_fleet1_L");
  obs_fleet1_L.allocate(styr_fleet1_L,endyr_fleet1_L,"obs_fleet1_L");
  fleet1_L_cv.allocate(styr_fleet1_L,endyr_fleet1_L,"fleet1_L_cv");
  nyr_fleet1_agec.allocate("nyr_fleet1_agec");
  yrs_fleet1_agec.allocate(1,nyr_fleet1_agec,"yrs_fleet1_agec");
  nsamp_fleet1_agec.allocate(1,nyr_fleet1_agec,"nsamp_fleet1_agec");
  nfish_fleet1_agec.allocate(1,nyr_fleet1_agec,"nfish_fleet1_agec");
  obs_fleet1_agec.allocate(1,nyr_fleet1_agec,1,nages_agec,"obs_fleet1_agec");
  set_Linf.allocate(1,7,"set_Linf");
  set_K.allocate(1,7,"set_K");
  set_t0.allocate(1,7,"set_t0");
  set_len_cv.allocate(1,7,"set_len_cv");
  set_M_constant.allocate(1,7,"set_M_constant");
  set_steep.allocate(1,7,"set_steep");
  set_log_R0.allocate(1,7,"set_log_R0");
  set_R_autocorr.allocate(1,7,"set_R_autocorr");
  set_rec_sigma.allocate(1,7,"set_rec_sigma");
  set_log_dm_fleet1_ac.allocate(1,7,"set_log_dm_fleet1_ac");
  set_log_dm_survey1_ac.allocate(1,7,"set_log_dm_survey1_ac");
  set_log_dm_survey2_ac.allocate(1,7,"set_log_dm_survey2_ac");
  set_selpar_A50_fleet1_B1.allocate(1,7,"set_selpar_A50_fleet1_B1");
  set_selpar_slope_fleet1_B1.allocate(1,7,"set_selpar_slope_fleet1_B1");
  set_selpar_A50_fleet1_B2.allocate(1,7,"set_selpar_A50_fleet1_B2");
  set_selpar_slope_fleet1_B2.allocate(1,7,"set_selpar_slope_fleet1_B2");
  set_selpar_A50_survey1.allocate(1,7,"set_selpar_A50_survey1");
  set_selpar_slope_survey1.allocate(1,7,"set_selpar_slope_survey1");
  set_selpar_A50_survey2.allocate(1,7,"set_selpar_A50_survey2");
  set_selpar_slope_survey2.allocate(1,7,"set_selpar_slope_survey2");
  set_log_q_survey1.allocate(1,7,"set_log_q_survey1");
  set_log_q_survey2.allocate(1,7,"set_log_q_survey2");
  set_F_init.allocate(1,7,"set_F_init");
  set_log_avg_F_fleet1.allocate(1,7,"set_log_avg_F_fleet1");
  set_log_F_dev_fleet1.allocate(1,3,"set_log_F_dev_fleet1");
  set_log_rec_dev.allocate(1,3,"set_log_rec_dev");
  set_log_Nage_dev.allocate(1,3,"set_log_Nage_dev");
  set_log_F_dev_fleet1_vals.allocate(styr_fleet1_L,endyr_fleet1_L,"set_log_F_dev_fleet1_vals");
  set_log_rec_dev_vals.allocate(styr_rec_dev,endyr_rec_dev,"set_log_rec_dev_vals");
  set_log_Nage_dev_vals.allocate(2,nages,"set_log_Nage_dev_vals");
  set_w_L.allocate("set_w_L");
  set_w_I_survey1.allocate("set_w_I_survey1");
  set_w_I_survey2.allocate("set_w_I_survey2");
  set_w_ac_fleet1.allocate("set_w_ac_fleet1");
  set_w_ac_survey1.allocate("set_w_ac_survey1");
  set_w_ac_survey2.allocate("set_w_ac_survey2");
  set_w_Nage_init.allocate("set_w_Nage_init");
  set_w_rec.allocate("set_w_rec");
  set_w_rec_early.allocate("set_w_rec_early");
  set_w_rec_end.allocate("set_w_rec_end");
  set_w_Ftune.allocate("set_w_Ftune");
  wgtpar_a.allocate("wgtpar_a");
  wgtpar_b.allocate("wgtpar_b");
  prop_f_obs.allocate(1,nages,"prop_f_obs");
  maturity_f_obs.allocate(1,nages,"maturity_f_obs");
  set_M.allocate(1,nages,"set_M");
  spawn_time_frac.allocate("spawn_time_frac");
  SR_switch.allocate("SR_switch");
  set_Ftune_yr.allocate("set_Ftune_yr");
  set_Ftune.allocate("set_Ftune");
  minSS_agec.allocate("minSS_agec");
  age_error.allocate(1,nages,1,nages,"age_error");
  endyr_proj.allocate("endyr_proj");
  styr_regs.allocate("styr_regs");
  Fproj_switch.allocate("Fproj_switch");
  Fproj_mult.allocate("Fproj_mult");
   styr_proj=endyr+1;
  end_of_data_file.allocate("end_of_data_file");
   if(end_of_data_file!=999)
   {
       cout << "*** WARNING: Data File NOT READ CORRECTLY ****" << endl;
       exit(0);  
   }
   else
   {cout << "Data File read correctly" << endl;} 
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  const double Linf_LO=set_Linf(2); const double Linf_HI=set_Linf(3); const double Linf_PH=set_Linf(4);
  const double K_LO=set_K(2); const double K_HI=set_K(3); const double K_PH=set_K(4);
  const double t0_LO=set_t0(2); const double t0_HI=set_t0(3); const double t0_PH=set_t0(4);  
  const double len_cv_LO=set_len_cv(2); const double len_cv_HI=set_len_cv(3); const double len_cv_PH=set_len_cv(4); 
  const double steep_LO=set_steep(2); const double steep_HI=set_steep(3); const double steep_PH=set_steep(4);
  const double log_R0_LO=set_log_R0(2); const double log_R0_HI=set_log_R0(3); const double log_R0_PH=set_log_R0(4);
  const double R_autocorr_LO=set_R_autocorr(2); const double R_autocorr_HI=set_R_autocorr(3); const double R_autocorr_PH=set_R_autocorr(4);
  const double rec_sigma_LO=set_rec_sigma(2); const double rec_sigma_HI=set_rec_sigma(3); const double rec_sigma_PH=set_rec_sigma(4);
  
  const double selpar_A50_survey1_LO=set_selpar_A50_survey1(2); const double selpar_A50_survey1_HI=set_selpar_A50_survey1(3); const double selpar_A50_survey1_PH=set_selpar_A50_survey1(4);
  const double selpar_slope_survey1_LO=set_selpar_slope_survey1(2); const double selpar_slope_survey1_HI=set_selpar_slope_survey1(3); const double selpar_slope_survey1_PH=set_selpar_slope_survey1(4);
  const double selpar_A50_survey2_LO=set_selpar_A50_survey2(2); const double selpar_A50_survey2_HI=set_selpar_A50_survey2(3); const double selpar_A50_survey2_PH=set_selpar_A50_survey2(4);
  const double selpar_slope_survey2_LO=set_selpar_slope_survey2(2); const double selpar_slope_survey2_HI=set_selpar_slope_survey2(3); const double selpar_slope_survey2_PH=set_selpar_slope_survey2(4);
  const double selpar_A50_fleet1_B1_LO=set_selpar_A50_fleet1_B1(2); const double selpar_A50_fleet1_B1_HI=set_selpar_A50_fleet1_B1(3); const double selpar_A50_fleet1_B1_PH=set_selpar_A50_fleet1_B1(4);
  const double selpar_slope_fleet1_B1_LO=set_selpar_slope_fleet1_B1(2); const double selpar_slope_fleet1_B1_HI=set_selpar_slope_fleet1_B1(3); const double selpar_slope_fleet1_B1_PH=set_selpar_slope_fleet1_B1(4);
  const double selpar_A50_fleet1_B2_LO=set_selpar_A50_fleet1_B2(2); const double selpar_A50_fleet1_B2_HI=set_selpar_A50_fleet1_B2(3); const double selpar_A50_fleet1_B2_PH=set_selpar_A50_fleet1_B2(4);
  const double selpar_slope_fleet1_B2_LO=set_selpar_slope_fleet1_B2(2); const double selpar_slope_fleet1_B2_HI=set_selpar_slope_fleet1_B2(3); const double selpar_slope_fleet1_B2_PH=set_selpar_slope_fleet1_B2(4);
  const double log_q_survey1_LO=set_log_q_survey1(2); const double log_q_survey1_HI=set_log_q_survey1(3); const double log_q_survey1_PH=set_log_q_survey1(4);
  const double log_q_survey2_LO=set_log_q_survey2(2); const double log_q_survey2_HI=set_log_q_survey2(3); const double log_q_survey2_PH=set_log_q_survey2(4);
  const double F_init_LO=set_F_init(2); const double F_init_HI=set_F_init(3); const double F_init_PH=set_F_init(4);
  const double log_avg_F_fleet1_LO=set_log_avg_F_fleet1(2); const double log_avg_F_fleet1_HI=set_log_avg_F_fleet1(3); const double log_avg_F_fleet1_PH=set_log_avg_F_fleet1(4);
  const double log_dm_fleet1_ac_LO=set_log_dm_fleet1_ac(2); const double log_dm_fleet1_ac_HI=set_log_dm_fleet1_ac(3); const double log_dm_fleet1_ac_PH=set_log_dm_fleet1_ac(4);
  const double log_dm_survey1_ac_LO=set_log_dm_survey1_ac(2); const double log_dm_survey1_ac_HI=set_log_dm_survey1_ac(3); const double log_dm_survey1_ac_PH=set_log_dm_survey1_ac(4);
  const double log_dm_survey2_ac_LO=set_log_dm_survey2_ac(2); const double log_dm_survey2_ac_HI=set_log_dm_survey2_ac(3); const double log_dm_survey2_ac_PH=set_log_dm_survey2_ac(4);
  //-dev vectors-----------------------------------------------------------------------------------------------------------  
  const double log_F_dev_fleet1_LO=set_log_F_dev_fleet1(1); const double log_F_dev_fleet1_HI=set_log_F_dev_fleet1(2); const double log_F_dev_fleet1_PH=set_log_F_dev_fleet1(3);   
  const double log_rec_dev_LO=set_log_rec_dev(1); const double log_rec_dev_HI=set_log_rec_dev(2); const double log_rec_dev_PH=set_log_rec_dev(3);          
  const double log_Nage_dev_LO=set_log_Nage_dev(1); const double log_Nage_dev_HI=set_log_Nage_dev(2); const double log_Nage_dev_PH=set_log_Nage_dev(3);          
  Linf.allocate(Linf_LO,Linf_HI,Linf_PH,"Linf");
  K.allocate(K_LO,K_HI,K_PH,"K");
  t0.allocate(t0_LO,t0_HI,t0_PH,"t0");
  len_cv_val.allocate(len_cv_LO,len_cv_HI,len_cv_PH,"len_cv_val");
  Linf_out.allocate(1,8,"Linf_out");
  #ifndef NO_AD_INITIALIZE
    Linf_out.initialize();
  #endif
  K_out.allocate(1,8,"K_out");
  #ifndef NO_AD_INITIALIZE
    K_out.initialize();
  #endif
  t0_out.allocate(1,8,"t0_out");
  #ifndef NO_AD_INITIALIZE
    t0_out.initialize();
  #endif
  len_cv_val_out.allocate(1,8,"len_cv_val_out");
  #ifndef NO_AD_INITIALIZE
    len_cv_val_out.initialize();
  #endif
  meanlen_TL.allocate(1,nages,"meanlen_TL");
  #ifndef NO_AD_INITIALIZE
    meanlen_TL.initialize();
  #endif
  wgt_g.allocate(1,nages,"wgt_g");
  #ifndef NO_AD_INITIALIZE
    wgt_g.initialize();
  #endif
  wgt_kg.allocate(1,nages,"wgt_kg");
  #ifndef NO_AD_INITIALIZE
    wgt_kg.initialize();
  #endif
  wgt_mt.allocate(1,nages,"wgt_mt");
  #ifndef NO_AD_INITIALIZE
    wgt_mt.initialize();
  #endif
  wgt_klb.allocate(1,nages,"wgt_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_klb.initialize();
  #endif
  wgt_lb.allocate(1,nages,"wgt_lb");
  #ifndef NO_AD_INITIALIZE
    wgt_lb.initialize();
  #endif
  pred_survey1_agec.allocate(1,nyr_survey1_agec,1,nages_agec,"pred_survey1_agec");
  #ifndef NO_AD_INITIALIZE
    pred_survey1_agec.initialize();
  #endif
  pred_survey1_agec_allages.allocate(1,nyr_survey1_agec,1,nages,"pred_survey1_agec_allages");
  #ifndef NO_AD_INITIALIZE
    pred_survey1_agec_allages.initialize();
  #endif
  ErrorFree_survey1_agec.allocate(1,nyr_survey1_agec,1,nages,"ErrorFree_survey1_agec");
  #ifndef NO_AD_INITIALIZE
    ErrorFree_survey1_agec.initialize();
  #endif
  pred_survey2_agec.allocate(1,nyr_survey2_agec,1,nages_agec,"pred_survey2_agec");
  #ifndef NO_AD_INITIALIZE
    pred_survey2_agec.initialize();
  #endif
  pred_survey2_agec_allages.allocate(1,nyr_survey2_agec,1,nages,"pred_survey2_agec_allages");
  #ifndef NO_AD_INITIALIZE
    pred_survey2_agec_allages.initialize();
  #endif
  ErrorFree_survey2_agec.allocate(1,nyr_survey2_agec,1,nages,"ErrorFree_survey2_agec");
  #ifndef NO_AD_INITIALIZE
    ErrorFree_survey2_agec.initialize();
  #endif
  pred_fleet1_agec.allocate(1,nyr_fleet1_agec,1,nages_agec,"pred_fleet1_agec");
  #ifndef NO_AD_INITIALIZE
    pred_fleet1_agec.initialize();
  #endif
  pred_fleet1_agec_allages.allocate(1,nyr_fleet1_agec,1,nages,"pred_fleet1_agec_allages");
  #ifndef NO_AD_INITIALIZE
    pred_fleet1_agec_allages.initialize();
  #endif
  ErrorFree_fleet1_agec.allocate(1,nyr_fleet1_agec,1,nages,"ErrorFree_fleet1_agec");
  #ifndef NO_AD_INITIALIZE
    ErrorFree_fleet1_agec.initialize();
  #endif
  obs_survey1_cpue_allyr.allocate(styr,endyr,"obs_survey1_cpue_allyr");
  #ifndef NO_AD_INITIALIZE
    obs_survey1_cpue_allyr.initialize();
  #endif
  pred_survey1_cpue_allyr.allocate(styr,endyr,"pred_survey1_cpue_allyr");
  #ifndef NO_AD_INITIALIZE
    pred_survey1_cpue_allyr.initialize();
  #endif
  survey1_cpue_cv_allyr.allocate(styr,endyr,"survey1_cpue_cv_allyr");
  #ifndef NO_AD_INITIALIZE
    survey1_cpue_cv_allyr.initialize();
  #endif
  obs_survey2_cpue_allyr.allocate(styr,endyr,"obs_survey2_cpue_allyr");
  #ifndef NO_AD_INITIALIZE
    obs_survey2_cpue_allyr.initialize();
  #endif
  pred_survey2_cpue_allyr.allocate(styr,endyr,"pred_survey2_cpue_allyr");
  #ifndef NO_AD_INITIALIZE
    pred_survey2_cpue_allyr.initialize();
  #endif
  survey2_cpue_cv_allyr.allocate(styr,endyr,"survey2_cpue_cv_allyr");
  #ifndef NO_AD_INITIALIZE
    survey2_cpue_cv_allyr.initialize();
  #endif
  nsamp_survey1_agec_allyr.allocate(styr,endyr,"nsamp_survey1_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_survey1_agec_allyr.initialize();
  #endif
  nsamp_survey2_agec_allyr.allocate(styr,endyr,"nsamp_survey2_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_survey2_agec_allyr.initialize();
  #endif
  nsamp_fleet1_agec_allyr.allocate(styr,endyr,"nsamp_fleet1_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_fleet1_agec_allyr.initialize();
  #endif
  nfish_survey1_agec_allyr.allocate(styr,endyr,"nfish_survey1_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_survey1_agec_allyr.initialize();
  #endif
  nfish_survey2_agec_allyr.allocate(styr,endyr,"nfish_survey2_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_survey2_agec_allyr.initialize();
  #endif
  nfish_fleet1_agec_allyr.allocate(styr,endyr,"nfish_fleet1_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_fleet1_agec_allyr.initialize();
  #endif
  neff_survey1_agec_allyr_out.allocate(styr,endyr,"neff_survey1_agec_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_survey1_agec_allyr_out.initialize();
  #endif
  neff_survey2_agec_allyr_out.allocate(styr,endyr,"neff_survey2_agec_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_survey2_agec_allyr_out.initialize();
  #endif
  neff_fleet1_agec_allyr_out.allocate(styr,endyr,"neff_fleet1_agec_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_fleet1_agec_allyr_out.initialize();
  #endif
  N.allocate(styr,endyr+1,1,nages,"N");
  #ifndef NO_AD_INITIALIZE
    N.initialize();
  #endif
  N_mdyr.allocate(styr,endyr,1,nages,"N_mdyr");
  #ifndef NO_AD_INITIALIZE
    N_mdyr.initialize();
  #endif
  N_spawn.allocate(styr,endyr,1,nages,"N_spawn");
  #ifndef NO_AD_INITIALIZE
    N_spawn.initialize();
  #endif
  log_Nage_dev.allocate(2,nages,log_Nage_dev_LO,log_Nage_dev_HI,log_Nage_dev_PH,"log_Nage_dev");
  log_Nage_dev_output.allocate(1,nages,"log_Nage_dev_output");
  #ifndef NO_AD_INITIALIZE
    log_Nage_dev_output.initialize();
  #endif
  B.allocate(styr,endyr+1,1,nages,"B");
  #ifndef NO_AD_INITIALIZE
    B.initialize();
  #endif
  totB.allocate(styr,endyr+1,"totB");
  #ifndef NO_AD_INITIALIZE
    totB.initialize();
  #endif
  totN.allocate(styr,endyr+1,"totN");
  #ifndef NO_AD_INITIALIZE
    totN.initialize();
  #endif
  SSB.allocate(styr,endyr,"SSB");
  #ifndef NO_AD_INITIALIZE
    SSB.initialize();
  #endif
  MatFemB.allocate(styr,endyr,"MatFemB");
  #ifndef NO_AD_INITIALIZE
    MatFemB.initialize();
  #endif
  rec.allocate(styr,endyr+1,"rec");
  #ifndef NO_AD_INITIALIZE
    rec.initialize();
  #endif
  prop_f.allocate(1,nages,"prop_f");
  #ifndef NO_AD_INITIALIZE
    prop_f.initialize();
  #endif
  maturity_f.allocate(1,nages,"maturity_f");
  #ifndef NO_AD_INITIALIZE
    maturity_f.initialize();
  #endif
  reprod.allocate(1,nages,"reprod");
  #ifndef NO_AD_INITIALIZE
    reprod.initialize();
  #endif
  log_R0.allocate(log_R0_LO,log_R0_HI,log_R0_PH,"log_R0");
  log_R0_out.allocate(1,8,"log_R0_out");
  #ifndef NO_AD_INITIALIZE
    log_R0_out.initialize();
  #endif
  R0.allocate("R0");
  #ifndef NO_AD_INITIALIZE
  R0.initialize();
  #endif
  steep.allocate(steep_LO,steep_HI,steep_PH,"steep");
  steep_out.allocate(1,8,"steep_out");
  #ifndef NO_AD_INITIALIZE
    steep_out.initialize();
  #endif
  rec_sigma.allocate(rec_sigma_LO,rec_sigma_HI,rec_sigma_PH,"rec_sigma");
  rec_sigma_out.allocate(1,8,"rec_sigma_out");
  #ifndef NO_AD_INITIALIZE
    rec_sigma_out.initialize();
  #endif
  R_autocorr.allocate(R_autocorr_LO,R_autocorr_HI,R_autocorr_PH,"R_autocorr");
  R_autocorr_out.allocate(1,8,"R_autocorr_out");
  #ifndef NO_AD_INITIALIZE
    R_autocorr_out.initialize();
  #endif
  rec_sigma_sq.allocate("rec_sigma_sq");
  #ifndef NO_AD_INITIALIZE
  rec_sigma_sq.initialize();
  #endif
  rec_logL_add.allocate("rec_logL_add");
  #ifndef NO_AD_INITIALIZE
  rec_logL_add.initialize();
  #endif
  log_rec_dev.allocate(styr_rec_dev,endyr_rec_dev,log_rec_dev_LO,log_rec_dev_HI,log_rec_dev_PH,"log_rec_dev");
  log_rec_dev_output.allocate(styr,endyr+1,"log_rec_dev_output");
  #ifndef NO_AD_INITIALIZE
    log_rec_dev_output.initialize();
  #endif
  log_rec_dev_out.allocate(styr_rec_dev,endyr_rec_dev,"log_rec_dev_out");
  #ifndef NO_AD_INITIALIZE
    log_rec_dev_out.initialize();
  #endif
  var_rec_dev.allocate("var_rec_dev");
  #ifndef NO_AD_INITIALIZE
  var_rec_dev.initialize();
  #endif
  sigma_rec_dev.allocate("sigma_rec_dev");
  #ifndef NO_AD_INITIALIZE
  sigma_rec_dev.initialize();
  #endif
  BiasCor.allocate("BiasCor");
  #ifndef NO_AD_INITIALIZE
  BiasCor.initialize();
  #endif
  S0.allocate("S0");
  #ifndef NO_AD_INITIALIZE
  S0.initialize();
  #endif
  B0.allocate("B0");
  #ifndef NO_AD_INITIALIZE
  B0.initialize();
  #endif
  R1.allocate("R1");
  #ifndef NO_AD_INITIALIZE
  R1.initialize();
  #endif
  R_virgin.allocate("R_virgin");
  #ifndef NO_AD_INITIALIZE
  R_virgin.initialize();
  #endif
  SdS0.allocate(styr,endyr,"SdS0");
  #ifndef NO_AD_INITIALIZE
    SdS0.initialize();
  #endif
  log_dm_fleet1_ac.allocate(log_dm_fleet1_ac_LO,log_dm_fleet1_ac_HI,log_dm_fleet1_ac_PH,"log_dm_fleet1_ac");
  log_dm_survey1_ac.allocate(log_dm_survey1_ac_LO,log_dm_survey1_ac_HI,log_dm_survey1_ac_PH,"log_dm_survey1_ac");
  log_dm_survey2_ac.allocate(log_dm_survey2_ac_LO,log_dm_survey2_ac_HI,log_dm_survey2_ac_PH,"log_dm_survey2_ac");
  log_dm_fleet1_ac_out.allocate(1,8,"log_dm_fleet1_ac_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_fleet1_ac_out.initialize();
  #endif
  log_dm_survey1_ac_out.allocate(1,8,"log_dm_survey1_ac_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_survey1_ac_out.initialize();
  #endif
  log_dm_survey2_ac_out.allocate(1,8,"log_dm_survey2_ac_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_survey2_ac_out.initialize();
  #endif
  sel_survey1.allocate(styr,endyr,1,nages,"sel_survey1");
  #ifndef NO_AD_INITIALIZE
    sel_survey1.initialize();
  #endif
  selvec_survey1.allocate(1,nages,"selvec_survey1");
  #ifndef NO_AD_INITIALIZE
    selvec_survey1.initialize();
  #endif
  selpar_A50_survey1.allocate(selpar_A50_survey1_LO,selpar_A50_survey1_HI,selpar_A50_survey1_PH,"selpar_A50_survey1");
  selpar_slope_survey1.allocate(selpar_slope_survey1_LO,selpar_slope_survey1_HI,selpar_slope_survey1_PH,"selpar_slope_survey1");
  selpar_A50_survey1_out.allocate(1,8,"selpar_A50_survey1_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_survey1_out.initialize();
  #endif
  selpar_slope_survey1_out.allocate(1,8,"selpar_slope_survey1_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_survey1_out.initialize();
  #endif
  sel_survey2.allocate(styr,endyr,1,nages,"sel_survey2");
  #ifndef NO_AD_INITIALIZE
    sel_survey2.initialize();
  #endif
  selvec_survey2.allocate(1,nages,"selvec_survey2");
  #ifndef NO_AD_INITIALIZE
    selvec_survey2.initialize();
  #endif
  selpar_A50_survey2.allocate(selpar_A50_survey2_LO,selpar_A50_survey2_HI,selpar_A50_survey2_PH,"selpar_A50_survey2");
  selpar_slope_survey2.allocate(selpar_slope_survey2_LO,selpar_slope_survey2_HI,selpar_slope_survey2_PH,"selpar_slope_survey2");
  selpar_A50_survey2_out.allocate(1,8,"selpar_A50_survey2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_survey2_out.initialize();
  #endif
  selpar_slope_survey2_out.allocate(1,8,"selpar_slope_survey2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_survey2_out.initialize();
  #endif
  sel_fleet1.allocate(styr,endyr,1,nages,"sel_fleet1");
  #ifndef NO_AD_INITIALIZE
    sel_fleet1.initialize();
  #endif
  selvec_fleet1_B1.allocate(1,nages,"selvec_fleet1_B1");
  #ifndef NO_AD_INITIALIZE
    selvec_fleet1_B1.initialize();
  #endif
  selvec_fleet1_B2.allocate(1,nages,"selvec_fleet1_B2");
  #ifndef NO_AD_INITIALIZE
    selvec_fleet1_B2.initialize();
  #endif
  selpar_A50_fleet1_B1.allocate(selpar_A50_fleet1_B1_LO,selpar_A50_fleet1_B1_HI,selpar_A50_fleet1_B1_PH,"selpar_A50_fleet1_B1");
  selpar_slope_fleet1_B1.allocate(selpar_slope_fleet1_B1_LO,selpar_slope_fleet1_B1_HI,selpar_slope_fleet1_B1_PH,"selpar_slope_fleet1_B1");
  selpar_A50_fleet1_B2.allocate(selpar_A50_fleet1_B2_LO,selpar_A50_fleet1_B2_HI,selpar_A50_fleet1_B2_PH,"selpar_A50_fleet1_B2");
  selpar_slope_fleet1_B2.allocate(selpar_slope_fleet1_B2_LO,selpar_slope_fleet1_B2_HI,selpar_slope_fleet1_B2_PH,"selpar_slope_fleet1_B2");
  selpar_A50_fleet1_B1_out.allocate(1,8,"selpar_A50_fleet1_B1_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_fleet1_B1_out.initialize();
  #endif
  selpar_slope_fleet1_B1_out.allocate(1,8,"selpar_slope_fleet1_B1_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_fleet1_B1_out.initialize();
  #endif
  selpar_A50_fleet1_B2_out.allocate(1,8,"selpar_A50_fleet1_B2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_fleet1_B2_out.initialize();
  #endif
  selpar_slope_fleet1_B2_out.allocate(1,8,"selpar_slope_fleet1_B2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_fleet1_B2_out.initialize();
  #endif
  sel_wgted_L.allocate(1,nages,"sel_wgted_L");
  #ifndef NO_AD_INITIALIZE
    sel_wgted_L.initialize();
  #endif
  sel_wgted_tot.allocate(1,nages,"sel_wgted_tot");
  #ifndef NO_AD_INITIALIZE
    sel_wgted_tot.initialize();
  #endif
  pred_survey1_cpue.allocate(1,nyr_survey1_cpue,"pred_survey1_cpue");
  #ifndef NO_AD_INITIALIZE
    pred_survey1_cpue.initialize();
  #endif
  N_survey1.allocate(1,nyr_survey1_cpue,1,nages,"N_survey1");
  #ifndef NO_AD_INITIALIZE
    N_survey1.initialize();
  #endif
  pred_survey2_cpue.allocate(1,nyr_survey2_cpue,"pred_survey2_cpue");
  #ifndef NO_AD_INITIALIZE
    pred_survey2_cpue.initialize();
  #endif
  N_survey2.allocate(1,nyr_survey2_cpue,1,nages,"N_survey2");
  #ifndef NO_AD_INITIALIZE
    N_survey2.initialize();
  #endif
  log_q_survey1.allocate(log_q_survey1_LO,log_q_survey1_HI,log_q_survey1_PH,"log_q_survey1");
  log_q_survey1_out.allocate(1,8,"log_q_survey1_out");
  #ifndef NO_AD_INITIALIZE
    log_q_survey1_out.initialize();
  #endif
  log_q_survey2.allocate(log_q_survey2_LO,log_q_survey2_HI,log_q_survey2_PH,"log_q_survey2");
  log_q_survey2_out.allocate(1,8,"log_q_survey2_out");
  #ifndef NO_AD_INITIALIZE
    log_q_survey2_out.initialize();
  #endif
  L_fleet1_num.allocate(styr,endyr,1,nages,"L_fleet1_num");
  #ifndef NO_AD_INITIALIZE
    L_fleet1_num.initialize();
  #endif
  L_fleet1_mt.allocate(styr,endyr,1,nages,"L_fleet1_mt");
  #ifndef NO_AD_INITIALIZE
    L_fleet1_mt.initialize();
  #endif
  pred_fleet1_L_knum.allocate(styr,endyr,"pred_fleet1_L_knum");
  #ifndef NO_AD_INITIALIZE
    pred_fleet1_L_knum.initialize();
  #endif
  pred_fleet1_L_mt.allocate(styr,endyr,"pred_fleet1_L_mt");
  #ifndef NO_AD_INITIALIZE
    pred_fleet1_L_mt.initialize();
  #endif
  L_total_num.allocate(styr,endyr,1,nages,"L_total_num");
  #ifndef NO_AD_INITIALIZE
    L_total_num.initialize();
  #endif
  L_total_mt.allocate(styr,endyr,1,nages,"L_total_mt");
  #ifndef NO_AD_INITIALIZE
    L_total_mt.initialize();
  #endif
  L_total_knum_yr.allocate(styr,endyr,"L_total_knum_yr");
  #ifndef NO_AD_INITIALIZE
    L_total_knum_yr.initialize();
  #endif
  L_total_mt_yr.allocate(styr,endyr,"L_total_mt_yr");
  #ifndef NO_AD_INITIALIZE
    L_total_mt_yr.initialize();
  #endif
  F_fleet1_prop.allocate("F_fleet1_prop");
  #ifndef NO_AD_INITIALIZE
  F_fleet1_prop.initialize();
  #endif
  F_init_fleet1_prop.allocate("F_init_fleet1_prop");
  #ifndef NO_AD_INITIALIZE
  F_init_fleet1_prop.initialize();
  #endif
  F_temp_sum.allocate("F_temp_sum");
  #ifndef NO_AD_INITIALIZE
  F_temp_sum.initialize();
  #endif
  F_end.allocate(1,nages,"F_end");
  #ifndef NO_AD_INITIALIZE
    F_end.initialize();
  #endif
  F_end_L.allocate(1,nages,"F_end_L");
  #ifndef NO_AD_INITIALIZE
    F_end_L.initialize();
  #endif
  F_end_apex.allocate("F_end_apex");
  #ifndef NO_AD_INITIALIZE
  F_end_apex.initialize();
  #endif
  SSB_msy_out.allocate("SSB_msy_out");
  #ifndef NO_AD_INITIALIZE
  SSB_msy_out.initialize();
  #endif
  F_msy_out.allocate("F_msy_out");
  #ifndef NO_AD_INITIALIZE
  F_msy_out.initialize();
  #endif
  msy_mt_out.allocate("msy_mt_out");
  #ifndef NO_AD_INITIALIZE
  msy_mt_out.initialize();
  #endif
  msy_knum_out.allocate("msy_knum_out");
  #ifndef NO_AD_INITIALIZE
  msy_knum_out.initialize();
  #endif
  B_msy_out.allocate("B_msy_out");
  #ifndef NO_AD_INITIALIZE
  B_msy_out.initialize();
  #endif
  R_msy_out.allocate("R_msy_out");
  #ifndef NO_AD_INITIALIZE
  R_msy_out.initialize();
  #endif
  spr_msy_out.allocate("spr_msy_out");
  #ifndef NO_AD_INITIALIZE
  spr_msy_out.initialize();
  #endif
  F30_dum.allocate("F30_dum");
  #ifndef NO_AD_INITIALIZE
  F30_dum.initialize();
  #endif
  F35_dum.allocate("F35_dum");
  #ifndef NO_AD_INITIALIZE
  F35_dum.initialize();
  #endif
  F40_dum.allocate("F40_dum");
  #ifndef NO_AD_INITIALIZE
  F40_dum.initialize();
  #endif
  F30_out.allocate("F30_out");
  #ifndef NO_AD_INITIALIZE
  F30_out.initialize();
  #endif
  F35_out.allocate("F35_out");
  #ifndef NO_AD_INITIALIZE
  F35_out.initialize();
  #endif
  F40_out.allocate("F40_out");
  #ifndef NO_AD_INITIALIZE
  F40_out.initialize();
  #endif
  SSB_F35_out.allocate("SSB_F35_out");
  #ifndef NO_AD_INITIALIZE
  SSB_F35_out.initialize();
  #endif
  B_F35_out.allocate("B_F35_out");
  #ifndef NO_AD_INITIALIZE
  B_F35_out.initialize();
  #endif
  R_F35_out.allocate("R_F35_out");
  #ifndef NO_AD_INITIALIZE
  R_F35_out.initialize();
  #endif
  L_F35_knum_out.allocate("L_F35_knum_out");
  #ifndef NO_AD_INITIALIZE
  L_F35_knum_out.initialize();
  #endif
  L_F35_mt_out.allocate("L_F35_mt_out");
  #ifndef NO_AD_INITIALIZE
  L_F35_mt_out.initialize();
  #endif
  D_F35_knum_out.allocate("D_F35_knum_out");
  #ifndef NO_AD_INITIALIZE
  D_F35_knum_out.initialize();
  #endif
  D_F35_mt_out.allocate("D_F35_mt_out");
  #ifndef NO_AD_INITIALIZE
  D_F35_mt_out.initialize();
  #endif
  wgt_wgted_L_mt.allocate(1,nages,"wgt_wgted_L_mt");
  #ifndef NO_AD_INITIALIZE
    wgt_wgted_L_mt.initialize();
  #endif
  N_age_msy.allocate(1,nages,"N_age_msy");
  #ifndef NO_AD_INITIALIZE
    N_age_msy.initialize();
  #endif
  N_age_msy_spawn.allocate(1,nages,"N_age_msy_spawn");
  #ifndef NO_AD_INITIALIZE
    N_age_msy_spawn.initialize();
  #endif
  L_age_msy.allocate(1,nages,"L_age_msy");
  #ifndef NO_AD_INITIALIZE
    L_age_msy.initialize();
  #endif
  Z_age_msy.allocate(1,nages,"Z_age_msy");
  #ifndef NO_AD_INITIALIZE
    Z_age_msy.initialize();
  #endif
  F_L_age_msy.allocate(1,nages,"F_L_age_msy");
  #ifndef NO_AD_INITIALIZE
    F_L_age_msy.initialize();
  #endif
  F_msy.allocate(1,n_iter_msy,"F_msy");
  #ifndef NO_AD_INITIALIZE
    F_msy.initialize();
  #endif
  spr_msy.allocate(1,n_iter_msy,"spr_msy");
  #ifndef NO_AD_INITIALIZE
    spr_msy.initialize();
  #endif
  R_eq.allocate(1,n_iter_msy,"R_eq");
  #ifndef NO_AD_INITIALIZE
    R_eq.initialize();
  #endif
  L_eq_mt.allocate(1,n_iter_msy,"L_eq_mt");
  #ifndef NO_AD_INITIALIZE
    L_eq_mt.initialize();
  #endif
  L_eq_knum.allocate(1,n_iter_msy,"L_eq_knum");
  #ifndef NO_AD_INITIALIZE
    L_eq_knum.initialize();
  #endif
  SSB_eq.allocate(1,n_iter_msy,"SSB_eq");
  #ifndef NO_AD_INITIALIZE
    SSB_eq.initialize();
  #endif
  B_eq.allocate(1,n_iter_msy,"B_eq");
  #ifndef NO_AD_INITIALIZE
    B_eq.initialize();
  #endif
  FdF_msy.allocate(styr,endyr,"FdF_msy");
  #ifndef NO_AD_INITIALIZE
    FdF_msy.initialize();
  #endif
  FdF35.allocate(styr,endyr,"FdF35");
  #ifndef NO_AD_INITIALIZE
    FdF35.initialize();
  #endif
  SdSSB_msy.allocate(styr,endyr,"SdSSB_msy");
  #ifndef NO_AD_INITIALIZE
    SdSSB_msy.initialize();
  #endif
  SdSSB_msy_end.allocate("SdSSB_msy_end");
  #ifndef NO_AD_INITIALIZE
  SdSSB_msy_end.initialize();
  #endif
  FdF_msy_end.allocate("FdF_msy_end");
  #ifndef NO_AD_INITIALIZE
  FdF_msy_end.initialize();
  #endif
  FdF_msy_end_mean.allocate("FdF_msy_end_mean");
  #ifndef NO_AD_INITIALIZE
  FdF_msy_end_mean.initialize();
  #endif
  SdSSB_F35.allocate(styr,endyr,"SdSSB_F35");
  #ifndef NO_AD_INITIALIZE
    SdSSB_F35.initialize();
  #endif
  Sdmsst_F35.allocate(styr,endyr,"Sdmsst_F35");
  #ifndef NO_AD_INITIALIZE
    Sdmsst_F35.initialize();
  #endif
  SdSSB_F35_end.allocate("SdSSB_F35_end");
  #ifndef NO_AD_INITIALIZE
  SdSSB_F35_end.initialize();
  #endif
  Sdmsst_F35_end.allocate("Sdmsst_F35_end");
  #ifndef NO_AD_INITIALIZE
  Sdmsst_F35_end.initialize();
  #endif
  FdF35_end_mean.allocate("FdF35_end_mean");
  #ifndef NO_AD_INITIALIZE
  FdF35_end_mean.initialize();
  #endif
  Fend_mean_temp.allocate("Fend_mean_temp");
  #ifndef NO_AD_INITIALIZE
  Fend_mean_temp.initialize();
  #endif
  Fend_mean.allocate("Fend_mean");
  #ifndef NO_AD_INITIALIZE
  Fend_mean.initialize();
  #endif
  L_age_F35.allocate(1,nages,"L_age_F35");
  #ifndef NO_AD_INITIALIZE
    L_age_F35.initialize();
  #endif
  iter_inc_msy.allocate("iter_inc_msy");
  #ifndef NO_AD_INITIALIZE
  iter_inc_msy.initialize();
  #endif
  M.allocate(1,nages,"M");
  #ifndef NO_AD_INITIALIZE
    M.initialize();
  #endif
  smsy2msst.allocate("smsy2msst");
  #ifndef NO_AD_INITIALIZE
  smsy2msst.initialize();
  #endif
  smsy2msst75.allocate("smsy2msst75");
  #ifndef NO_AD_INITIALIZE
  smsy2msst75.initialize();
  #endif
  F.allocate(styr,endyr,1,nages,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  Fsum.allocate(styr,endyr,"Fsum");
  #ifndef NO_AD_INITIALIZE
    Fsum.initialize();
  #endif
  Fapex.allocate(styr,endyr,"Fapex");
  #ifndef NO_AD_INITIALIZE
    Fapex.initialize();
  #endif
  Z.allocate(styr,endyr,1,nages,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  log_avg_F_fleet1.allocate(log_avg_F_fleet1_LO,log_avg_F_fleet1_HI,log_avg_F_fleet1_PH,"log_avg_F_fleet1");
  log_avg_F_fleet1_out.allocate(1,8,"log_avg_F_fleet1_out");
  #ifndef NO_AD_INITIALIZE
    log_avg_F_fleet1_out.initialize();
  #endif
  log_F_dev_fleet1.allocate(styr_fleet1_L,endyr_fleet1_L,log_F_dev_fleet1_LO,log_F_dev_fleet1_HI,log_F_dev_fleet1_PH,"log_F_dev_fleet1");
  log_F_dev_fleet1_out.allocate(styr_fleet1_L,endyr_fleet1_L,"log_F_dev_fleet1_out");
  #ifndef NO_AD_INITIALIZE
    log_F_dev_fleet1_out.initialize();
  #endif
  F_fleet1.allocate(styr,endyr,1,nages,"F_fleet1");
  #ifndef NO_AD_INITIALIZE
    F_fleet1.initialize();
  #endif
  F_fleet1_out.allocate(styr,endyr,"F_fleet1_out");
  #ifndef NO_AD_INITIALIZE
    F_fleet1_out.initialize();
  #endif
  log_F_dev_init_fleet1.allocate("log_F_dev_init_fleet1");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_init_fleet1.initialize();
  #endif
  log_F_dev_end_fleet1.allocate("log_F_dev_end_fleet1");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_end_fleet1.initialize();
  #endif
  F_init.allocate(F_init_LO,F_init_HI,F_init_PH,"F_init");
  F_init_out.allocate(1,8,"F_init_out");
  #ifndef NO_AD_INITIALIZE
    F_init_out.initialize();
  #endif
  N_age_spr.allocate(1,nages,"N_age_spr");
  #ifndef NO_AD_INITIALIZE
    N_age_spr.initialize();
  #endif
  N_age_spr_spawn.allocate(1,nages,"N_age_spr_spawn");
  #ifndef NO_AD_INITIALIZE
    N_age_spr_spawn.initialize();
  #endif
  L_age_spr.allocate(1,nages,"L_age_spr");
  #ifndef NO_AD_INITIALIZE
    L_age_spr.initialize();
  #endif
  Z_age_spr.allocate(1,nages,"Z_age_spr");
  #ifndef NO_AD_INITIALIZE
    Z_age_spr.initialize();
  #endif
  spr_static.allocate(styr,endyr,"spr_static");
  #ifndef NO_AD_INITIALIZE
    spr_static.initialize();
  #endif
  F_L_age_spr.allocate(1,nages,"F_L_age_spr");
  #ifndef NO_AD_INITIALIZE
    F_L_age_spr.initialize();
  #endif
  F_spr.allocate(1,n_iter_spr,"F_spr");
  #ifndef NO_AD_INITIALIZE
    F_spr.initialize();
  #endif
  spr_spr.allocate(1,n_iter_spr,"spr_spr");
  #ifndef NO_AD_INITIALIZE
    spr_spr.initialize();
  #endif
  spr_ratio.allocate(1,n_iter_spr,"spr_ratio");
  #ifndef NO_AD_INITIALIZE
    spr_ratio.initialize();
  #endif
  L_spr.allocate(1,n_iter_spr,"L_spr");
  #ifndef NO_AD_INITIALIZE
    L_spr.initialize();
  #endif
  N_spr_F0.allocate(1,nages,"N_spr_F0");
  #ifndef NO_AD_INITIALIZE
    N_spr_F0.initialize();
  #endif
  N_bpr_F0.allocate(1,nages,"N_bpr_F0");
  #ifndef NO_AD_INITIALIZE
    N_bpr_F0.initialize();
  #endif
  N_spr_initial.allocate(1,nages,"N_spr_initial");
  #ifndef NO_AD_INITIALIZE
    N_spr_initial.initialize();
  #endif
  N_initial_eq.allocate(1,nages,"N_initial_eq");
  #ifndef NO_AD_INITIALIZE
    N_initial_eq.initialize();
  #endif
  F_initial.allocate(1,nages,"F_initial");
  #ifndef NO_AD_INITIALIZE
    F_initial.initialize();
  #endif
  Z_initial.allocate(1,nages,"Z_initial");
  #ifndef NO_AD_INITIALIZE
    Z_initial.initialize();
  #endif
  spr_initial.allocate("spr_initial");
  #ifndef NO_AD_INITIALIZE
  spr_initial.initialize();
  #endif
  spr_F0.allocate("spr_F0");
  #ifndef NO_AD_INITIALIZE
  spr_F0.initialize();
  #endif
  bpr_F0.allocate("bpr_F0");
  #ifndef NO_AD_INITIALIZE
  bpr_F0.initialize();
  #endif
  rec_mean.allocate("rec_mean");
  #ifndef NO_AD_INITIALIZE
  rec_mean.initialize();
  #endif
  iter_inc_spr.allocate("iter_inc_spr");
  #ifndef NO_AD_INITIALIZE
  iter_inc_spr.initialize();
  #endif
  sdnr_I_survey1.allocate("sdnr_I_survey1");
  #ifndef NO_AD_INITIALIZE
  sdnr_I_survey1.initialize();
  #endif
  sdnr_I_survey2.allocate("sdnr_I_survey2");
  #ifndef NO_AD_INITIALIZE
  sdnr_I_survey2.initialize();
  #endif
  sdnr_ac_survey1.allocate("sdnr_ac_survey1");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_survey1.initialize();
  #endif
  sdnr_ac_survey2.allocate("sdnr_ac_survey2");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_survey2.initialize();
  #endif
  sdnr_ac_fleet1.allocate("sdnr_ac_fleet1");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_fleet1.initialize();
  #endif
  w_L.allocate("w_L");
  #ifndef NO_AD_INITIALIZE
  w_L.initialize();
  #endif
  w_I_survey1.allocate("w_I_survey1");
  #ifndef NO_AD_INITIALIZE
  w_I_survey1.initialize();
  #endif
  w_I_survey2.allocate("w_I_survey2");
  #ifndef NO_AD_INITIALIZE
  w_I_survey2.initialize();
  #endif
  w_ac_fleet1.allocate("w_ac_fleet1");
  #ifndef NO_AD_INITIALIZE
  w_ac_fleet1.initialize();
  #endif
  w_ac_survey1.allocate("w_ac_survey1");
  #ifndef NO_AD_INITIALIZE
  w_ac_survey1.initialize();
  #endif
  w_ac_survey2.allocate("w_ac_survey2");
  #ifndef NO_AD_INITIALIZE
  w_ac_survey2.initialize();
  #endif
  w_Nage_init.allocate("w_Nage_init");
  #ifndef NO_AD_INITIALIZE
  w_Nage_init.initialize();
  #endif
  w_rec.allocate("w_rec");
  #ifndef NO_AD_INITIALIZE
  w_rec.initialize();
  #endif
  w_rec_early.allocate("w_rec_early");
  #ifndef NO_AD_INITIALIZE
  w_rec_early.initialize();
  #endif
  w_rec_end.allocate("w_rec_end");
  #ifndef NO_AD_INITIALIZE
  w_rec_end.initialize();
  #endif
  w_Ftune.allocate("w_Ftune");
  #ifndef NO_AD_INITIALIZE
  w_Ftune.initialize();
  #endif
  f_fleet1_L.allocate("f_fleet1_L");
  #ifndef NO_AD_INITIALIZE
  f_fleet1_L.initialize();
  #endif
  f_survey1_cpue.allocate("f_survey1_cpue");
  #ifndef NO_AD_INITIALIZE
  f_survey1_cpue.initialize();
  #endif
  f_survey2_cpue.allocate("f_survey2_cpue");
  #ifndef NO_AD_INITIALIZE
  f_survey2_cpue.initialize();
  #endif
  f_fleet1_agec.allocate("f_fleet1_agec");
  #ifndef NO_AD_INITIALIZE
  f_fleet1_agec.initialize();
  #endif
  f_survey1_agec.allocate("f_survey1_agec");
  #ifndef NO_AD_INITIALIZE
  f_survey1_agec.initialize();
  #endif
  f_survey2_agec.allocate("f_survey2_agec");
  #ifndef NO_AD_INITIALIZE
  f_survey2_agec.initialize();
  #endif
  f_Nage_init.allocate("f_Nage_init");
  #ifndef NO_AD_INITIALIZE
  f_Nage_init.initialize();
  #endif
  f_rec_dev.allocate("f_rec_dev");
  #ifndef NO_AD_INITIALIZE
  f_rec_dev.initialize();
  #endif
  f_rec_dev_early.allocate("f_rec_dev_early");
  #ifndef NO_AD_INITIALIZE
  f_rec_dev_early.initialize();
  #endif
  f_rec_dev_end.allocate("f_rec_dev_end");
  #ifndef NO_AD_INITIALIZE
  f_rec_dev_end.initialize();
  #endif
  f_Ftune.allocate("f_Ftune");
  #ifndef NO_AD_INITIALIZE
  f_Ftune.initialize();
  #endif
  f_priors.allocate("f_priors");
  #ifndef NO_AD_INITIALIZE
  f_priors.initialize();
  #endif
  fval.allocate("fval");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  fval_data.allocate("fval_data");
  #ifndef NO_AD_INITIALIZE
  fval_data.initialize();
  #endif
  grad_max.allocate("grad_max");
  #ifndef NO_AD_INITIALIZE
  grad_max.initialize();
  #endif
  denom.allocate("denom");
  #ifndef NO_AD_INITIALIZE
  denom.initialize();
  #endif
  numer.allocate("numer");
  #ifndef NO_AD_INITIALIZE
  numer.initialize();
  #endif
  F_reg_proj.allocate("F_reg_proj");
  #ifndef NO_AD_INITIALIZE
  F_reg_proj.initialize();
  #endif
  F_proj.allocate(styr_proj,endyr_proj,"F_proj");
  #ifndef NO_AD_INITIALIZE
    F_proj.initialize();
  #endif
  L_knum_proj.allocate(styr_proj,endyr_proj,"L_knum_proj");
  #ifndef NO_AD_INITIALIZE
    L_knum_proj.initialize();
  #endif
  L_mt_proj.allocate(styr_proj,endyr_proj,"L_mt_proj");
  #ifndef NO_AD_INITIALIZE
    L_mt_proj.initialize();
  #endif
  B_proj.allocate(styr_proj,endyr_proj,"B_proj");
  #ifndef NO_AD_INITIALIZE
    B_proj.initialize();
  #endif
  SSB_proj.allocate(styr_proj,endyr_proj,"SSB_proj");
  #ifndef NO_AD_INITIALIZE
    SSB_proj.initialize();
  #endif
  R_proj.allocate(styr_proj,endyr_proj,"R_proj");
  #ifndef NO_AD_INITIALIZE
    R_proj.initialize();
  #endif
  FL_age_proj.allocate(1,nages,"FL_age_proj");
  #ifndef NO_AD_INITIALIZE
    FL_age_proj.initialize();
  #endif
  N_proj.allocate(styr_proj,endyr_proj,1,nages,"N_proj");
  #ifndef NO_AD_INITIALIZE
    N_proj.initialize();
  #endif
  N_spawn_proj.allocate(styr_proj,endyr_proj,1,nages,"N_spawn_proj");
  #ifndef NO_AD_INITIALIZE
    N_spawn_proj.initialize();
  #endif
  Z_proj.allocate(styr_proj,endyr_proj,1,nages,"Z_proj");
  #ifndef NO_AD_INITIALIZE
    Z_proj.initialize();
  #endif
  L_age_proj.allocate(styr_proj,endyr_proj,1,nages,"L_age_proj");
  #ifndef NO_AD_INITIALIZE
    L_age_proj.initialize();
  #endif
}

void model_parameters::set_runtime(void)
{
  dvector temp("{1.0e-4;}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
  dvector temp1("{1000,1600,10000;}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
	
  Linf=set_Linf(1);
  K=set_K(1);
  t0=set_t0(1);
  len_cv_val=set_len_cv(1);
  M=set_M; 
  smsy2msst=0.5;
  smsy2msst75=0.75;
  
  log_R0=set_log_R0(1);
  steep=set_steep(1);
  R_autocorr=set_R_autocorr(1);
  rec_sigma=set_rec_sigma(1);
  log_dm_fleet1_ac=set_log_dm_fleet1_ac(1);
  log_dm_survey1_ac=set_log_dm_survey1_ac(1);
  log_dm_survey2_ac=set_log_dm_survey2_ac(1);
   
  log_q_survey1=set_log_q_survey1(1);
  log_q_survey2=set_log_q_survey2(1);
   
  w_L=set_w_L;
  w_I_survey1=set_w_I_survey1;
  w_I_survey2=set_w_I_survey2;
  
  w_ac_fleet1=set_w_ac_fleet1;
  w_ac_survey1=set_w_ac_survey1;
  w_ac_survey2=set_w_ac_survey2;
  w_Nage_init=set_w_Nage_init;
  w_rec=set_w_rec;
  w_rec_early=set_w_rec_early;
  w_rec_end=set_w_rec_end;
  w_Ftune=set_w_Ftune;
  F_init=set_F_init(1);
  log_avg_F_fleet1=set_log_avg_F_fleet1(1);
  
  log_F_dev_fleet1=set_log_F_dev_fleet1_vals;
  selpar_A50_survey1=set_selpar_A50_survey1(1);
  selpar_slope_survey1=set_selpar_slope_survey1(1);
 
  selpar_A50_survey2=set_selpar_A50_survey2(1);
  selpar_slope_survey2=set_selpar_slope_survey2(1);
 
  selpar_A50_fleet1_B1=set_selpar_A50_fleet1_B1(1);
  selpar_slope_fleet1_B1=set_selpar_slope_fleet1_B1(1);
  selpar_A50_fleet1_B2=set_selpar_A50_fleet1_B2(1);
  selpar_slope_fleet1_B2=set_selpar_slope_fleet1_B2(1);
 
 sqrt2pi=sqrt(2.*3.14159265);
 g2mt=0.000001;         //conversion of grams to metric tons
 g2kg=0.001;            //conversion of grams to kg 
 kg2mt=0.001;            //conversion of kilograms to metric tons  
 mt2kg=1000.0;            //conversion of metric tons to kilograms   
 mt2klb=2.20462;        //conversion of metric tons to 1000 lb 
 mt2lb=mt2klb*1000.0;   //conversion of metric tons to lb
 g2klb=g2mt*mt2klb;     //conversion of grams to 1000 lb 
 dzero=0.00001;         
 huge_number=1.0e+10;   
 
 SSB_msy_out=0.0;
 iter_inc_msy=max_F_spr_msy/(n_iter_msy-1);
 iter_inc_spr=max_F_spr_msy/(n_iter_spr-1); 
 prop_f=prop_f_obs;
 maturity_f=maturity_f_obs;
 
  nsamp_survey1_agec_allyr=missing;
  nsamp_survey2_agec_allyr=missing;
  nsamp_fleet1_agec_allyr=missing;
      						
  for (iyear=1; iyear<=nyr_survey1_agec; iyear++)
         {if (nsamp_survey1_agec(iyear)>=minSS_agec)
           {nsamp_survey1_agec_allyr(yrs_survey1_agec(iyear))=nsamp_survey1_agec(iyear);
            nfish_survey1_agec_allyr(yrs_survey1_agec(iyear))=nfish_survey1_agec(iyear);}} 
  for (iyear=1; iyear<=nyr_survey2_agec; iyear++)
         {if (nsamp_survey2_agec(iyear)>=minSS_agec)
           {nsamp_survey2_agec_allyr(yrs_survey2_agec(iyear))=nsamp_survey2_agec(iyear);
            nfish_survey2_agec_allyr(yrs_survey2_agec(iyear))=nfish_survey2_agec(iyear);}} 
  for (iyear=1; iyear<=nyr_fleet1_agec; iyear++)
         {if (nsamp_fleet1_agec(iyear)>=minSS_agec)
           {nsamp_fleet1_agec_allyr(yrs_fleet1_agec(iyear))=nsamp_fleet1_agec(iyear);
            nfish_fleet1_agec_allyr(yrs_fleet1_agec(iyear))=nfish_fleet1_agec(iyear);}}
			
  obs_survey1_cpue_allyr=missing; 
  pred_survey1_cpue_allyr=missing;
  survey1_cpue_cv_allyr=missing;	
  obs_survey2_cpue_allyr=missing; 
  pred_survey2_cpue_allyr=missing;
  survey2_cpue_cv_allyr=missing;  
	  
  F_msy(1)=0.0;  
  for (ff=2;ff<=n_iter_msy;ff++) {F_msy(ff)=F_msy(ff-1)+iter_inc_msy;}
  F_spr(1)=0.0;  
  for (ff=2;ff<=n_iter_spr;ff++) {F_spr(ff)=F_spr(ff-1)+iter_inc_spr;}
  F_fleet1.initialize();
  L_fleet1_num.initialize();
   
  F_fleet1_out.initialize();
  sel_survey1.initialize();
  sel_survey2.initialize();
  sel_fleet1.initialize();
  
  log_rec_dev_output.initialize();  
  log_rec_dev=set_log_rec_dev_vals;
  log_Nage_dev_output.initialize();
  log_Nage_dev=set_log_Nage_dev_vals;
 
}

void model_parameters::userfunction(void)
{
  fval =0.0;
 //cout<<"start"<<endl;
 get_length_weight_at_age(); 
 //cout << "got length, weight, fecundity transitions" <<endl;
 get_reprod();
 //cout << "got repro stuff" << endl;
 get_spr_F0();
 //cout << "got F0 spr" << endl;
 get_selectivity(); 
 //cout << "got selectivity" << endl;
 get_mortality(); 
 //cout << "got mortalities" << endl;
 get_bias_corr(); 
 //cout<< "got recruitment bias correction" << endl;
 get_numbers_at_age(); 
 //cout << "got numbers at age" << endl;
 get_landings_numbers();
 //cout << "got landings in numbers" << endl;
 get_landings_wgt();
 // cout << "got landings in wgt" << endl;
 get_indices();
 //cout << "got indices" << endl;
 get_age_comps();
 //cout<< "got age comps"<< endl;
 evaluate_objective_function();
 //cout << "objective function calculations complete" << endl;
}

void model_parameters::get_length_weight_at_age(void)
{
	//population total length in mm
    //compute mean length (mm TL) and weight (whole) at age
    meanlen_TL=Linf*(1.0-mfexp(-K*(agebins-t0)));     
    wgt_kg=wgtpar_a*pow(meanlen_TL,wgtpar_b);             //whole wgt in kg 
    wgt_g=wgt_kg/g2kg;                                    //convert wgt in kg to weight in g    
    wgt_mt=wgt_g*g2mt;                                    //convert weight in g to weight in mt
    wgt_klb=mt2klb*wgt_mt;                                //1000 lb of whole wgt
    wgt_lb=mt2lb*wgt_mt;                                  //lb of whole wgt
}

void model_parameters::get_reprod(void)
{
   //reprod is product of stuff going into reproductive capacity calcs
   reprod=elem_prod((elem_prod(prop_f,maturity_f)),wgt_mt);
}

void model_parameters::get_spr_F0(void)
{
  N_spr_F0(1)=1.0*mfexp(-1.0*M(1)*spawn_time_frac); //at peak spawning time
  N_bpr_F0(1)=1.0;      //at start of year
  for (iage=2; iage<=nages; iage++)
  {
    N_spr_F0(iage)=N_spr_F0(iage-1)*mfexp(-1.0*(M(iage-1)*(1.0-spawn_time_frac) + M(iage)*spawn_time_frac)); 
    N_bpr_F0(iage)=N_bpr_F0(iage-1)*mfexp(-1.0*(M(iage-1)));    
  }
  N_spr_F0(nages)=N_spr_F0(nages)/(1.0-mfexp(-1.0*M(nages))); //plus group (sum of geometric series)
  N_bpr_F0(nages)=N_bpr_F0(nages)/(1.0-mfexp(-1.0*M(nages)));
  spr_F0=sum(elem_prod(N_spr_F0,reprod)); 
  bpr_F0=sum(elem_prod(N_bpr_F0,wgt_mt));    
}

void model_parameters::get_selectivity(void)
{
  selvec_survey1=logistic(agebins, selpar_A50_survey1, selpar_slope_survey1);	
  selvec_survey2=logistic(agebins, selpar_A50_survey2, selpar_slope_survey2);	
  selvec_fleet1_B1=logistic(agebins, selpar_A50_fleet1_B1, selpar_slope_fleet1_B1);  
  selvec_fleet1_B2=selvec_fleet1_B1; //Block2 selex = Block1 selex; i.e., assuming no block structure here
  //selvec_fleet1_B2=logistic(agebins, selpar_A50_fleet1_B2, selpar_slope_fleet1_B2); //uncomment if using two selex blocks 	
   for (iyear=styr; iyear<=endyr_selex_phase1; iyear++)
   {     
 	   sel_survey1(iyear)(1,nages_agec)=selvec_survey1(1,nages_agec);
           sel_survey2(iyear)(1,nages_agec)=selvec_survey2(1,nages_agec);
	   sel_fleet1(iyear)(1,nages_agec)=selvec_fleet1_B1(1,nages_agec);
   }
 for (iyear=(endyr_selex_phase1+1); iyear<=endyr; iyear++)
   {
	   sel_survey1(iyear)(1,nages_agec)=selvec_survey1(1,nages_agec);
           sel_survey2(iyear)(1,nages_agec)=selvec_survey2(1,nages_agec);
	   sel_fleet1(iyear)(1,nages_agec)=selvec_fleet1_B2(1,nages_agec);
   }	   
}

void model_parameters::get_mortality(void)
{
  Fsum.initialize();
  Fapex.initialize();
  F.initialize();
  //log_F_dev_init_fleet1=sum(log_F_dev_fleet1(styr_fleet1_L,(styr_fleet1_L+2)))/3.0;  //initialization F is avg from first 3 yrs of observed landings                  
  log_F_dev_init_fleet1=log_F_dev_fleet1(styr_fleet1_L); //For sim study, only using first year to match operating model
  F_init_fleet1_prop=mfexp(log_avg_F_fleet1+log_F_dev_init_fleet1)/(mfexp(log_avg_F_fleet1+log_F_dev_init_fleet1));
  for (iyear=styr; iyear<=endyr; iyear++) 
  {
    if(iyear>=styr_fleet1_L & iyear<=endyr_fleet1_L)
    {  F_fleet1_out(iyear)=mfexp(log_avg_F_fleet1+log_F_dev_fleet1(iyear));     
       F_fleet1(iyear)=sel_fleet1(iyear)*F_fleet1_out(iyear);
       Fsum(iyear)+=F_fleet1_out(iyear);
    }
    //Total F at age
    F(iyear)=F_fleet1(iyear); //first in additive series (NO +=)
    Fapex(iyear)=max(F(iyear));
    Z(iyear)=M+F(iyear);  
  }  //end iyear 
}

void model_parameters::get_bias_corr(void)
{
  var_rec_dev=norm2(log_rec_dev(styr_rec_dev,endyr_rec_dev)-
              sum(log_rec_dev(styr_rec_dev,endyr_rec_dev))/nyrs_rec)
              /(nyrs_rec-1.0);                           
  //if (set_BiasCor <= 0.0) {BiasCor=mfexp(var_rec_dev/2.0);}   //bias correction based on empirical residuals
  rec_sigma_sq=square(rec_sigma);
  if (set_BiasCor <= 0.0) {BiasCor=mfexp(rec_sigma_sq/2.0);}   //bias correction based on Rsigma               
  else {BiasCor=set_BiasCor;}
}

void model_parameters::get_numbers_at_age(void)
{
  R0=mfexp(log_R0);
  S0=spr_F0*R0;
  R_virgin=SR_eq_func(R0, steep, spr_F0, spr_F0, BiasCor, SR_switch);
  B0=bpr_F0*R_virgin;   
  F_initial=F_init*F_fleet1(1);  //initialization F is a scalar multiple (F_init) of F in year one
  Z_initial=M+F_initial;
  N_spr_initial(1)=1.0*mfexp(-1.0*Z_initial(1)*spawn_time_frac); //at peak spawning time;
  for (iage=2; iage<=nages; iage++)
    {
      N_spr_initial(iage)=N_spr_initial(iage-1)*
                   mfexp(-1.0*(Z_initial(iage-1)*(1.0-spawn_time_frac) + Z_initial(iage)*spawn_time_frac)); 
    }
  N_spr_initial(nages)=N_spr_initial(nages)/(1.0-mfexp(-1.0*Z_initial(nages))); //plus group
  spr_initial=sum(elem_prod(N_spr_initial,reprod));
  if (styr==styr_rec_dev) {R1=SR_eq_func(R0, steep, spr_F0, spr_initial, 1.0, SR_switch);} //without bias correction (deviation added later)
  else {R1=SR_eq_func(R0, steep, spr_F0, spr_initial, BiasCor, SR_switch);} //with bias correction
  if(R1<10.0) {R1=10.0;} //Avoid unrealistically low popn sizes during search algorithm
  N_initial_eq(1)=R1;
  for (iage=2; iage<=nages; iage++)
  {
    N_initial_eq(iage)=N_initial_eq(iage-1)*
        mfexp(-1.0*(Z_initial(iage-1)));    
  }
  //plus group calculation
  N_initial_eq(nages)=N_initial_eq(nages)/(1.0-mfexp(-1.0*Z_initial(nages))); //plus group
  N(styr)(2,nages)=elem_prod(N_initial_eq(2,nages),mfexp(log_Nage_dev));
  if (styr==styr_rec_dev) {N(styr,1)=N_initial_eq(1)*mfexp(log_rec_dev(styr_rec_dev));}
  else {N(styr,1)=N_initial_eq(1);}
  N_mdyr(styr)(1,nages)=elem_prod(N(styr)(1,nages),(mfexp(-1.*(Z_initial(1,nages))*0.5))); //mid year 
  N_spawn(styr)(1,nages)=elem_prod(N(styr)(1,nages),(mfexp(-1.*(Z_initial(1,nages))*spawn_time_frac))); //peak spawning time 
  SSB(styr)=sum(elem_prod(N_spawn(styr),reprod));
  MatFemB(styr)=sum(elem_prod(N_spawn(styr),reprod));
  for (iyear=styr; iyear<endyr; iyear++)
  {
    if(iyear<(styr_rec_dev-1)||iyear>(endyr_rec_dev-1)) //recruitment follows S-R curve (with bias correction) exactly
    {
        N(iyear+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB(iyear),SR_switch);
        N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(mfexp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*mfexp(-1.*Z(iyear,nages)); //plus group
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year 
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); //peak spawning time 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod));   
        MatFemB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod));		
    }
    else   //recruitment follows S-R curve with lognormal deviation
    {
        N(iyear+1,1)=SR_func(R0, steep, spr_F0, SSB(iyear),SR_switch)*mfexp(log_rec_dev(iyear+1));
        N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(mfexp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*mfexp(-1.*Z(iyear,nages)); //plus group
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); //peak spawning time 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod));
		MatFemB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod));		
    }
  }
  //last year (projection) has no recruitment variability
  N(endyr+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB(endyr),SR_switch);
  N(endyr+1)(2,nages)=++elem_prod(N(endyr)(1,nages-1),(mfexp(-1.*Z(endyr)(1,nages-1))));
  N(endyr+1,nages)+=N(endyr,nages)*mfexp(-1.*Z(endyr,nages)); //plus group
  rec=column(N,1);
  SdS0=SSB/S0;
}

void model_parameters::get_landings_numbers(void)
{
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      L_fleet1_num(iyear,iage)=N(iyear,iage)*F_fleet1(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
    }          
    pred_fleet1_L_knum(iyear)=sum(L_fleet1_num(iyear))/1000.0;    
  }
}

void model_parameters::get_landings_wgt(void)
{
  for (iyear=styr; iyear<=endyr; iyear++)
  {    
    L_fleet1_mt(iyear)=elem_prod(L_fleet1_num(iyear),wgt_mt);     //in mt
    pred_fleet1_L_mt(iyear)=sum(L_fleet1_mt(iyear)); 
  }
}

void model_parameters::get_indices(void)
{
 //Survey 1: cpue
  for (iyear=1; iyear<=nyr_survey1_cpue; iyear++)
  {   //index in number units, start of yr
      N_survey1(iyear)=elem_prod(N(yrs_survey1_cpue(iyear)),sel_survey1(yrs_survey1_cpue(iyear)));
	  pred_survey1_cpue(iyear)=mfexp(log_q_survey1)*sum(N_survey1(iyear));
	  obs_survey1_cpue_allyr(yrs_survey1_cpue(iyear))=obs_survey1_cpue(iyear);
	  pred_survey1_cpue_allyr(yrs_survey1_cpue(iyear))=pred_survey1_cpue(iyear);
	  survey1_cpue_cv_allyr(yrs_survey1_cpue(iyear))=survey1_cpue_cv(iyear);
  }
  //Survey 2: cpue
  for (iyear=1; iyear<=nyr_survey2_cpue; iyear++)
  {   //index in number units, start of yr
      N_survey2(iyear)=elem_prod(N(yrs_survey2_cpue(iyear)),sel_survey2(yrs_survey2_cpue(iyear)));
	  pred_survey2_cpue(iyear)=mfexp(log_q_survey2)*sum(N_survey2(iyear));
	  obs_survey2_cpue_allyr(yrs_survey2_cpue(iyear))=obs_survey2_cpue(iyear);
	  pred_survey2_cpue_allyr(yrs_survey2_cpue(iyear))=pred_survey2_cpue(iyear);
	  survey2_cpue_cv_allyr(yrs_survey2_cpue(iyear))=survey2_cpue_cv(iyear);
  }
}

void model_parameters::get_age_comps(void)
{
 //survey1
 for (iyear=1;iyear<=nyr_survey1_agec;iyear++)
  {
    ErrorFree_survey1_agec(iyear)=N_survey1(iyear)/sum(N_survey1(iyear));
    pred_survey1_agec_allages(iyear)=age_error*ErrorFree_survey1_agec(iyear); 
    for (iage=1; iage<=nages_agec; iage++) {pred_survey1_agec(iyear,iage)=pred_survey1_agec_allages(iyear,iage);} 
    for (iage=(nages_agec+1); iage<=nages; iage++) {pred_survey1_agec(iyear,nages_agec)+=pred_survey1_agec_allages(iyear,iage);} //plus group                        
  }
  //survey2
 for (iyear=1;iyear<=nyr_survey2_agec;iyear++)
  {
    ErrorFree_survey2_agec(iyear)=N_survey2(iyear)/sum(N_survey2(iyear));
    pred_survey2_agec_allages(iyear)=age_error*ErrorFree_survey2_agec(iyear); 
    for (iage=1; iage<=nages_agec; iage++) {pred_survey2_agec(iyear,iage)=pred_survey2_agec_allages(iyear,iage);} 
    for (iage=(nages_agec+1); iage<=nages; iage++) {pred_survey2_agec(iyear,nages_agec)+=pred_survey2_agec_allages(iyear,iage);} //plus group                        
  }
  //fleet1
  for (iyear=1;iyear<=nyr_fleet1_agec;iyear++) 
  {
    ErrorFree_fleet1_agec(iyear)=L_fleet1_num(yrs_fleet1_agec(iyear))/sum(L_fleet1_num(yrs_fleet1_agec(iyear)));  
    pred_fleet1_agec_allages(iyear)=age_error*(ErrorFree_fleet1_agec(iyear)/sum(ErrorFree_fleet1_agec(iyear)));   
    for (iage=1; iage<=nages_agec; iage++) {pred_fleet1_agec(iyear,iage)=pred_fleet1_agec_allages(iyear,iage);} 
    for (iage=(nages_agec+1); iage<=nages; iage++) {pred_fleet1_agec(iyear,nages_agec)+=pred_fleet1_agec_allages(iyear,iage);} //plus group                             
  }
 //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
}

void model_parameters::get_weighted_current(void)
{
  F_temp_sum=0.0;
  F_fleet1_prop=1.0; //only used if multiple fleets; set to one otherwise
  sel_wgted_tot=selvec_fleet1_B2;
  sel_wgted_L=selvec_fleet1_B2;
  wgt_wgted_L_mt=wgt_mt; //weighted among fleets, only if multiple fleets are modeled                         
}

void model_parameters::get_msy(void)
{
  //compute values as functions of F
  for(ff=1; ff<=n_iter_msy; ff++)
  {
    //uses fleet-weighted F's; although only one fleet in this sim study
    Z_age_msy=0.0;
    F_L_age_msy=0.0;
    F_L_age_msy=F_msy(ff)*sel_wgted_L;
    Z_age_msy=M+F_L_age_msy;         
    N_age_msy(1)=1.0;
    for (iage=2; iage<=nages; iage++)
      {N_age_msy(iage)=N_age_msy(iage-1)*mfexp(-1.*Z_age_msy(iage-1));}
    N_age_msy(nages)=N_age_msy(nages)/(1.0-mfexp(-1.*Z_age_msy(nages)));
    N_age_msy_spawn(1,(nages-1))=elem_prod(N_age_msy(1,(nages-1)),
                                   mfexp((-1.*Z_age_msy(1,(nages-1)))*spawn_time_frac));                 
    N_age_msy_spawn(nages)=(N_age_msy_spawn(nages-1)*(mfexp(-1.*(Z_age_msy(nages-1)*(1.0-spawn_time_frac) + 
                            Z_age_msy(nages)*spawn_time_frac) )))/(1.0-mfexp(-1.*Z_age_msy(nages)));
    spr_msy(ff)=sum(elem_prod(N_age_msy_spawn,reprod));
    R_eq(ff)=SR_eq_func(R0, steep, spr_msy(1), spr_msy(ff), BiasCor, SR_switch);
    if (R_eq(ff)<dzero) {R_eq(ff)=dzero;}    
    N_age_msy*=R_eq(ff);
    N_age_msy_spawn*=R_eq(ff);
    for (iage=1; iage<=nages; iage++)
    {
      L_age_msy(iage)=N_age_msy(iage)*(F_L_age_msy(iage)/Z_age_msy(iage))*
                      (1.-mfexp(-1.*Z_age_msy(iage)));
    }
    SSB_eq(ff)=sum(elem_prod(N_age_msy_spawn,reprod));
	B_eq(ff)=sum(elem_prod(N_age_msy,wgt_mt));
    L_eq_mt(ff)=sum(elem_prod(L_age_msy,wgt_wgted_L_mt)); //in whole weight
    L_eq_knum(ff)=sum(L_age_msy)/1000.0;  
  }  
  msy_mt_out=max(L_eq_mt); //msy in whole weight 
  for(ff=1; ff<=n_iter_msy; ff++)
  {
   if(L_eq_mt(ff) == msy_mt_out) 
      {    
        SSB_msy_out=SSB_eq(ff);
        B_msy_out=B_eq(ff);
        R_msy_out=R_eq(ff);
        msy_knum_out=L_eq_knum(ff);
        F_msy_out=F_msy(ff);  
        spr_msy_out=spr_msy(ff);      
      }
  }
}

void model_parameters::get_per_recruit_stuff(void)
{
  //static per-recruit stuff 
  for(iyear=styr; iyear<=endyr; iyear++)
  {
    N_age_spr(1)=1.0;
    for(iage=2; iage<=nages; iage++)
     {N_age_spr(iage)=N_age_spr(iage-1)*mfexp(-1.*Z(iyear,iage-1));}
    N_age_spr(nages)=N_age_spr(nages)/(1.0-mfexp(-1.*Z(iyear,nages)));    
    N_age_spr_spawn(1,(nages-1))=elem_prod(N_age_spr(1,(nages-1)),
                                mfexp(-1.*Z(iyear)(1,(nages-1))*spawn_time_frac));
    N_age_spr_spawn(nages)=(N_age_spr_spawn(nages-1)*
                          (mfexp(-1.*(Z(iyear)(nages-1)*(1.0-spawn_time_frac) + Z(iyear)(nages)*spawn_time_frac) )))
                          /(1.0-mfexp(-1.*Z(iyear)(nages)));           
    spr_static(iyear)=sum(elem_prod(N_age_spr_spawn,reprod))/spr_F0;
  }
  //compute SSB/R and YPR as functions of F
  for(ff=1; ff<=n_iter_spr; ff++)
  {
    //uses fishery-weighted F's, same as in MSY calculations    
    F_L_age_spr=0.0;
    Z_age_spr=0.0;
	F_L_age_spr=F_spr(ff)*sel_wgted_L;
    Z_age_spr=M+F_L_age_spr; 
    N_age_spr(1)=1.0;
    for (iage=2; iage<=nages; iage++)
     {N_age_spr(iage)=N_age_spr(iage-1)*mfexp(-1.*Z_age_spr(iage-1));}
    N_age_spr(nages)=N_age_spr(nages)/(1-mfexp(-1.*Z_age_spr(nages)));
    N_age_spr_spawn(1,(nages-1))=elem_prod(N_age_spr(1,(nages-1)),
                                   mfexp((-1.*Z_age_spr(1,(nages-1)))*spawn_time_frac));                 
    N_age_spr_spawn(nages)=(N_age_spr_spawn(nages-1)*
                          (mfexp(-1.*(Z_age_spr(nages-1)*(1.0-spawn_time_frac) + Z_age_spr(nages)*spawn_time_frac) )))
                          /(1.0-mfexp(-1.*Z_age_spr(nages)));
    spr_spr(ff)=sum(elem_prod(N_age_spr_spawn,reprod));
    L_spr(ff)=0.0;
    for (iage=1; iage<=nages; iage++)
    {
      L_age_spr(iage)=N_age_spr(iage)*(F_L_age_spr(iage)/Z_age_spr(iage))*
                      (1.-mfexp(-1.*Z_age_spr(iage)));
      L_spr(ff)+=L_age_spr(iage)*wgt_wgted_L_mt(iage)*1000.0; //in whole wgt
    }   
  }
  spr_ratio=spr_spr/spr_F0;
  F30_dum=min(fabs(spr_ratio-0.3));
  F35_dum=min(fabs(spr_ratio-0.35));
  F40_dum=min(fabs(spr_ratio-0.4));
  for(ff=1; ff<=n_iter_spr; ff++)
  {   
      if (fabs(spr_ratio(ff)-0.3)==F30_dum) {F30_out=F_spr(ff);}	  
	  if (fabs(spr_ratio(ff)-0.35)==F35_dum) {F35_out=F_spr(ff);}	  
	  if (fabs(spr_ratio(ff)-0.4)==F40_dum) {F40_out=F_spr(ff);}
  }
  rec=column(N,1);
  rec_mean=sum(rec(styr_rec_dev, endyr_rec_dev))/nyrs_rec; //avg recruitment across all years
  R_F35_out=rec_mean;
  F_L_age_spr=F35_out*sel_wgted_L;
  Z_age_spr=M+F_L_age_spr; 
  N_age_spr(1)=R_F35_out;
  for (iage=2; iage<=nages; iage++)
     {N_age_spr(iage)=N_age_spr(iage-1)*mfexp(-1.*Z_age_spr(iage-1));}
  N_age_spr(nages)=N_age_spr(nages)/(1-mfexp(-1.*Z_age_spr(nages)));
  N_age_spr_spawn(1,(nages-1))=elem_prod(N_age_spr(1,(nages-1)),
                                   mfexp((-1.*Z_age_spr(1,(nages-1)))*spawn_time_frac));                 
  N_age_spr_spawn(nages)=(N_age_spr_spawn(nages-1)*
                          (mfexp(-1.*(Z_age_spr(nages-1)*(1.0-spawn_time_frac) + Z_age_spr(nages)*spawn_time_frac) )))
                          /(1.0-mfexp(-1.*Z_age_spr(nages)));
  for (iage=1; iage<=nages; iage++)
    {
      L_age_F35(iage)=N_age_spr(iage)*(F_L_age_spr(iage)/Z_age_spr(iage))*
                      (1.-mfexp(-1.*Z_age_spr(iage)));
    }
  SSB_F35_out=sum(elem_prod(N_age_spr_spawn,reprod));
  B_F35_out=sum(elem_prod(N_age_spr,wgt_mt));
  L_F35_mt_out=sum(elem_prod(L_age_F35,wgt_wgted_L_mt)); //in whole weight
  L_F35_knum_out=sum(L_age_F35)/1000.0;  
}

void model_parameters::get_miscellaneous_stuff(void)
{
  if(var_rec_dev>0.0)
   {sigma_rec_dev=sqrt(var_rec_dev);} //sample SD of predicted residuals (may not equal rec_sigma)  
   else{sigma_rec_dev=0.0;}
  //compute total landings- and discards-at-age in 1000 fish and mt whole weight
  L_total_num.initialize();
  L_total_mt.initialize();
  L_total_knum_yr.initialize();
  L_total_mt_yr.initialize();  
  for(iyear=styr; iyear<=endyr; iyear++)
  {
        L_total_mt_yr(iyear)=pred_fleet1_L_mt(iyear);
        L_total_knum_yr(iyear)=pred_fleet1_L_knum(iyear);
        B(iyear)=elem_prod(N(iyear),wgt_mt);
        totN(iyear)=sum(N(iyear));
        totB(iyear)=sum(B(iyear));   
  }
  L_total_num=L_fleet1_num;   //landings at age in number fish
  L_total_mt=L_fleet1_mt;     //landings at age in mt whole weight
  //Time series of interest  
  B(endyr+1)=elem_prod(N(endyr+1),wgt_mt);
  totN(endyr+1)=sum(N(endyr+1));
  totB(endyr+1)=sum(B(endyr+1));  
  SdS0=SSB/S0;
  Fend_mean_temp=1.0;
  for (iyear=1; iyear<=selpar_n_yrs_wgted; iyear++) {Fend_mean_temp*=Fapex(endyr-iyear+1);}
  Fend_mean=pow(Fend_mean_temp,(1.0/selpar_n_yrs_wgted)); //fishing status based on last selpar_n_yrs_wgted years of assessment	  
  if(F_msy_out>0)
    {
      FdF_msy=Fapex/F_msy_out;
      FdF_msy_end=FdF_msy(endyr);
      FdF_msy_end_mean=Fend_mean/F_msy_out;
    }
  if(SSB_msy_out>0)
    {
      SdSSB_msy=SSB/SSB_msy_out;
      SdSSB_msy_end=SdSSB_msy(endyr);
    }  
	if(F35_out>0)
    {
	  FdF35=Fapex/F35_out;
	  FdF35_end_mean=Fend_mean/F35_out;
	}
  if(SSB_F35_out>0)
    {
      SdSSB_F35=SSB/SSB_F35_out;
	  Sdmsst_F35=SSB/(smsy2msst75*SSB_F35_out);
      SdSSB_F35_end=SdSSB_F35(endyr);
	  Sdmsst_F35_end=Sdmsst_F35(endyr);
    }  	
   //fill in log recruitment deviations for yrs they are nonzero
   for(iyear=styr_rec_dev; iyear<=endyr_rec_dev; iyear++)
     {log_rec_dev_output(iyear)=log_rec_dev(iyear);}
   //fill in log Nage deviations for ages they are nonzero (ages2+)
   for(iage=2; iage<=nages; iage++)
     {log_Nage_dev_output(iage)=log_Nage_dev(iage);}
}

void model_parameters::get_projection(void)
{
    switch(Fproj_switch){
       case 1: //F=Fcurrent
          F_reg_proj=Fend_mean;
          break;
       case 2: //F=Fmsy
          F_reg_proj=F_msy_out;
          break;
       case 3: //F=F30
          F_reg_proj=F30_out;
          break;     
       case 4: //F=F40
          F_reg_proj=F40_out;
          break;          		  
       default: // no such switch available
          cout << "Error in input: Projection switch Fproj_switch must be set to 1, 2, 3, or 4." << endl;
          cout << "Presently it is set to " << Fproj_switch <<"."<< endl;
          exit(0);          
   }
  N_proj(styr_proj)=N(endyr+1); //initial conditions computed previously
  for (iyear=styr_proj; iyear<=endyr_proj; iyear++) //recruitment follows S-R curve (with bias correction) exactly
  {     
        if (iyear<styr_regs) {F_proj(iyear)=Fend_mean;}
		else {F_proj(iyear)=Fproj_mult*F_reg_proj;}
		FL_age_proj=sel_wgted_L*F_proj(iyear);
        Z_proj(iyear)=M+FL_age_proj; //+FD_age_proj;
        N_spawn_proj(iyear)(1,nages)=elem_prod(N_proj(iyear)(1,nages),(mfexp(-1.*(Z_proj(iyear)(1,nages))*spawn_time_frac))); //peak spawning time
		SSB_proj(iyear)= sum(elem_prod(N_spawn_proj(iyear),reprod));
        B_proj(iyear)=sum(elem_prod(N_proj(iyear),wgt_mt)); //uses spawning weight
		for (iage=1; iage<=nages; iage++)
			{L_age_proj(iyear,iage)=N_proj(iyear,iage)*FL_age_proj(iage)*(1.-mfexp(-1.*Z_proj(iyear,iage)))/Z_proj(iyear,iage);
			}          
        L_knum_proj(iyear)=sum(L_age_proj(iyear))/1000.0;
	    L_mt_proj(iyear)=sum(elem_prod(L_age_proj(iyear),wgt_wgted_L_mt));     //in mt
		if (iyear<endyr_proj) {
			N_proj(iyear+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB_proj(iyear),SR_switch);
			N_proj(iyear+1)(2,nages)=++elem_prod(N_proj(iyear)(1,nages-1),(mfexp(-1.*Z_proj(iyear)(1,nages-1))));
			N_proj(iyear+1,nages)+=N_proj(iyear,nages)*mfexp(-1.*Z_proj(iyear,nages)); //plus group		
		}
  }
   R_proj=column(N_proj,1);                          
      // neff_survey1_agec_allyr_out=missing;
      // neff_fleet1_agec_allyr_out=missing;
      // for (iyear=1; iyear<=nyr_survey1_agec; iyear++)
         // {if (nsamp_survey1_agec(iyear)>=minSS_survey1_agec)
            // {neff_survey1_agec_allyr_out(yrs_survey1_agec(iyear))=multinom_eff_N(pred_survey1_agec(iyear),obs_survey1_agec(iyear));}                            
          // else {neff_survey1_agec_allyr_out(yrs_survey1_agec(iyear))=-99;}
         // }
      // for (iyear=1; iyear<=nyr_fleet1_agec; iyear++)
         // {if (nsamp_fleet1_agec(iyear)>=minSS_fleet1_agec)
            // {neff_fleet1_agec_allyr_out(yrs_fleet1_agec(iyear))=multinom_eff_N(pred_fleet1_agec(iyear),obs_fleet1_agec(iyear));}                            
          // else {neff_fleet1_agec_allyr_out(yrs_fleet1_agec(iyear))=-99;}
         // }      
 //---------------------------------------------------------
}

void model_parameters::evaluate_objective_function(void)
{
  fval=0.0;
  fval_data=0.0;  
  //fval=square(xdum-2);
  f_survey1_cpue=lk_lognormal(pred_survey1_cpue, obs_survey1_cpue, survey1_cpue_cv, w_I_survey1);
  fval+=f_survey1_cpue;
  fval_data+=f_survey1_cpue;
  f_survey2_cpue=lk_lognormal(pred_survey2_cpue, obs_survey2_cpue, survey2_cpue_cv, w_I_survey2);
  fval+=f_survey2_cpue;
  fval_data+=f_survey2_cpue;
  //f_fleet1_L in mt
  f_fleet1_L=lk_lognormal(pred_fleet1_L_mt(styr_fleet1_L,endyr_fleet1_L), obs_fleet1_L(styr_fleet1_L,endyr_fleet1_L),
                      fleet1_L_cv(styr_fleet1_L,endyr_fleet1_L), w_L);
  fval+=f_fleet1_L;
  fval_data+=f_fleet1_L;
  //f_survey1_agec
  //f_survey1_agec=lk_robust_multinomial(nsamp_survey1_agec, pred_survey1_agec, obs_survey1_agec, nyr_survey1_agec, double(nages_agec), minSS_agec, w_ac_survey1);
  f_survey1_agec=lk_multinomial(nsamp_survey1_agec, pred_survey1_agec, obs_survey1_agec, nyr_survey1_agec, minSS_agec, w_ac_survey1);
  //f_survey1_agec=lk_dirichlet_multinomial(nsamp_survey1_agec, pred_survey1_agec, obs_survey1_agec, nyr_survey1_agec, double(nages_agec), minSS_agec, log_dm_survey1_ac);
  fval+=f_survey1_agec;
  fval_data+=f_survey1_agec;
  //f_survey2_agec
  //f_survey2_agec=lk_robust_multinomial(nsamp_survey2_agec, pred_survey2_agec, obs_survey2_agec, nyr_survey2_agec, double(nages_agec), minSS_agec, w_ac_survey2);
  f_survey2_agec=lk_multinomial(nsamp_survey2_agec, pred_survey2_agec, obs_survey2_agec, nyr_survey2_agec, minSS_agec, w_ac_survey2);
  //f_survey2_agec=lk_dirichlet_multinomial(nsamp_survey2_agec, pred_survey2_agec, obs_survey2_agec, nyr_survey2_agec, double(nages_agec), minSS_agec, log_dm_survey2_ac);
  fval+=f_survey2_agec;
  fval_data+=f_survey2_agec;
  //f_fleet1_agec
  //f_fleet1_agec=lk_robust_multinomial(nsamp_fleet1_agec, pred_fleet1_agec, obs_fleet1_agec, nyr_fleet1_agec, double(nages_agec), minSS_agec, w_ac_fleet1);
  f_fleet1_agec=lk_multinomial(nsamp_fleet1_agec, pred_fleet1_agec, obs_fleet1_agec, nyr_fleet1_agec, minSS_agec, w_ac_fleet1);
  //f_fleet1_agec=lk_dirichlet_multinomial(nsamp_fleet1_agec, pred_fleet1_agec, obs_fleet1_agec, nyr_fleet1_agec, double(nages_agec), minSS_agec, log_dm_fleet1_ac);
  fval+=f_fleet1_agec;
  fval_data+=f_fleet1_agec;
  //Light penalty applied to log_Nage_dev for deviation from zero. If not estimated, this penalty equals zero.
  f_Nage_init=norm2(log_Nage_dev);        
  fval+=w_Nage_init*f_Nage_init;
  f_rec_dev=0.0;
  rec_logL_add=nyrs_rec*log(rec_sigma);
  f_rec_dev=(square(log_rec_dev(styr_rec_dev) + rec_sigma_sq/2.0)/(2.0*rec_sigma_sq));
  for(iyear=(styr_rec_dev+1); iyear<=endyr; iyear++)
  {f_rec_dev+=(square(log_rec_dev(iyear)-R_autocorr*log_rec_dev(iyear-1) + rec_sigma_sq/2.0)/
               (2.0*rec_sigma_sq));}
  f_rec_dev+=rec_logL_add;            
  fval+=w_rec*f_rec_dev;
  f_rec_dev_early=0.0; //possible extra constraint on early rec deviations
  if (w_rec_early>0.0)
    { if (styr_rec_dev<endyr_rec_phase1)
        {  
          for(iyear=styr_rec_dev; iyear<=endyr_rec_phase1; iyear++)
          {f_rec_dev_early+=square(log_rec_dev(iyear));}
        }
  fval+=w_rec_early*f_rec_dev_early;
  }
  f_rec_dev_end=0.0; //possible extra constraint on ending rec deviations
  if (w_rec_end>0.0)
  { if (endyr_rec_phase2<endyr_rec_dev)
        {  
          for(iyear=(endyr_rec_phase2+1); iyear<=endyr_rec_dev; iyear++)
           {f_rec_dev_end+=square(log_rec_dev(iyear));}
        }
      fval+=w_rec_end*f_rec_dev_end;
   }  
  //Ftune penalty: does not apply in last phase
  f_Ftune=0.0; 
  if (w_Ftune>0.0)
  {if (set_Ftune>0.0 && !last_phase()) {f_Ftune=square(Fapex(set_Ftune_yr)-set_Ftune);}
   fval+=w_Ftune*f_Ftune;
  }
  f_priors=0.0; 
  f_priors+=neg_log_prior(len_cv_val,set_len_cv(5),set_len_cv(6),set_len_cv(7));
  f_priors+=neg_log_prior(steep,set_steep(5),set_steep(6),set_steep(7)); 
  f_priors+=neg_log_prior(log_R0,set_log_R0(5),set_log_R0(6),set_log_R0(7)); 
  f_priors+=neg_log_prior(R_autocorr,set_R_autocorr(5),set_R_autocorr(6),set_R_autocorr(7));
  f_priors+=neg_log_prior(rec_sigma,set_rec_sigma(5),set_rec_sigma(6),set_rec_sigma(7));
  f_priors+=neg_log_prior(selpar_A50_survey1,set_selpar_A50_survey1(5), set_selpar_A50_survey1(6), set_selpar_A50_survey1(7));
  f_priors+=neg_log_prior(selpar_slope_survey1,set_selpar_slope_survey1(5), set_selpar_slope_survey1(6), set_selpar_slope_survey1(7));
  f_priors+=neg_log_prior(selpar_A50_survey2,set_selpar_A50_survey2(5), set_selpar_A50_survey2(6), set_selpar_A50_survey2(7));
  f_priors+=neg_log_prior(selpar_slope_survey2,set_selpar_slope_survey2(5), set_selpar_slope_survey2(6), set_selpar_slope_survey2(7));
  f_priors+=neg_log_prior(selpar_A50_fleet1_B1,set_selpar_A50_fleet1_B1(5), set_selpar_A50_fleet1_B1(6), set_selpar_A50_fleet1_B1(7));
  f_priors+=neg_log_prior(selpar_slope_fleet1_B1,set_selpar_slope_fleet1_B1(5), set_selpar_slope_fleet1_B1(6), set_selpar_slope_fleet1_B1(7));
  f_priors+=neg_log_prior(selpar_A50_fleet1_B2,set_selpar_A50_fleet1_B2(5), set_selpar_A50_fleet1_B2(6), set_selpar_A50_fleet1_B2(7));
  f_priors+=neg_log_prior(selpar_slope_fleet1_B2,set_selpar_slope_fleet1_B2(5), set_selpar_slope_fleet1_B2(6), set_selpar_slope_fleet1_B2(7));
  f_priors+=neg_log_prior(log_q_survey1,set_log_q_survey1(5),set_log_q_survey1(6),set_log_q_survey1(7));
  f_priors+=neg_log_prior(log_q_survey2,set_log_q_survey2(5),set_log_q_survey2(6),set_log_q_survey2(7));
  f_priors+=neg_log_prior(F_init,set_F_init(5),set_F_init(6),set_F_init(7));
  f_priors+=neg_log_prior(log_avg_F_fleet1,set_log_avg_F_fleet1(5),set_log_avg_F_fleet1(6),set_log_avg_F_fleet1(7));
  f_priors+=neg_log_prior(log_dm_fleet1_ac,set_log_dm_fleet1_ac(5),set_log_dm_fleet1_ac(6),set_log_dm_fleet1_ac(7));
  f_priors+=neg_log_prior(log_dm_survey1_ac,set_log_dm_survey1_ac(5),set_log_dm_survey1_ac(6),set_log_dm_survey1_ac(7));
  f_priors+=neg_log_prior(log_dm_survey2_ac,set_log_dm_survey2_ac(5),set_log_dm_survey2_ac(6),set_log_dm_survey2_ac(7));
  fval+=f_priors;
  //cout << "fval = " << fval << "  fval_data = " << fval_data << endl;
}

dvar_vector model_parameters::logistic(const dvar_vector& ages, const dvariable& L50, const dvariable& slope)
{
  //ages=vector of ages, L50=age at 50% selectivity, slope=rate of increase
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=1./(1.+mfexp(-1.*slope*(ages-L50))); //logistic;  
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;
  //----------------------------------------------------------------------------------
}

dvar_vector model_parameters::logistic_peak(const dvar_vector& ages, const dvariable& L50, const dvariable& slope, const dvariable& peak)
{
  //ages=vector of ages, L50=age at 50% selectivity, slope=rate of increase, peak=asymptote
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=peak/(1.+mfexp(-1.*slope*(ages-L50))); //logistic;  
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;
}

dvar_vector model_parameters::logistic_exponential(const dvar_vector& ages, const dvariable& L50, const dvariable& slope, const dvariable& sigma, const dvariable& joint)
{
  //ages=vector of ages, L50=age at 50% sel (ascending limb), slope=rate of increase, sigma=controls rate of descent (descending)                               
  //joint=age to join curves                                                                                                                                    
  RETURN_ARRAYS_INCREMENT();                                                                                                                                    
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());                                                                                                         
  Sel_Tmp=1.0;                                                                                                                                                  
  for (iage=1; iage<=nages; iage++)                                                                                                                             
  {                                                                                                                                                             
   if (ages(iage)<joint) {Sel_Tmp(iage)=1./(1.+mfexp(-1.*slope*(ages(iage)-L50)));}                                                                             
   if (ages(iage)>joint){Sel_Tmp(iage)=mfexp(-1.*square((ages(iage)-joint)/sigma));}                                                                            
  }                                                                                                                                                             
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);                                                                                                                                 
  RETURN_ARRAYS_DECREMENT();                                                                                                                                    
  return Sel_Tmp;   
}

dvar_vector model_parameters::logistic_double(const dvar_vector& ages, const dvariable& L501, const dvariable& slope1, const dvariable& L502, const dvariable& slope2)
{
  //ages=vector of ages, L50=age at 50% selectivity, slope=rate of increase, L502=age at 50% decrease additive to L501, slope2=slope of decrease
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=elem_prod( (1./(1.+mfexp(-1.*slope1*(ages-L501)))),(1.-(1./(1.+mfexp(-1.*slope2*(ages-(L501+L502)))))) );     
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;
}

dvar_vector model_parameters::logistic_joint(const dvar_vector& ages, const dvariable& L501, const dvariable& slope1, const dvariable& L502, const dvariable& slope2, const dvariable& satval, const dvariable& joint)
{
  //ages=vector of ages, L501=age at 50% sel (ascending limb), slope1=rate of increase,L502=age at 50% sel (descending), slope1=rate of increase (ascending), 
  //satval=saturation value of descending limb, joint=location in age vector to join curves (may equal age or age + 1 if age-0 is infleet2uded)
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=1.0; 
  for (iage=1; iage<=nages; iage++)
  {
   if (double(iage)<joint) {Sel_Tmp(iage)=1./(1.+mfexp(-1.*slope1*(ages(iage)-L501)));}  
   if (double(iage)>joint){Sel_Tmp(iage)=1.0-(1.0-satval)/(1.+mfexp(-1.*slope2*(ages(iage)-L502)));}  
  }  
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;
}

dvar_vector model_parameters::gaussian_double(const dvar_vector& ages, const dvariable& peak, const dvariable& top, const dvariable& ascwid, const dvariable& deswid, const dvariable& init, const dvariable& final)
{
  //ages=vector of ages, peak=ascending inflection location (as logistic), top=width of plateau, ascwid=ascent width (as log(width))
  //deswid=descent width (as log(width))
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step1(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step2(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step3(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step4(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step5(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step6(ages.indexmin(),ages.indexmax());
  dvar_vector pars_tmp(1,6); dvar_vector sel_tmp_iq(1,2);
  pars_tmp(1)=peak;
  pars_tmp(2)=peak+1.0+(0.99*ages(nages)-peak-1.0)/(1.0+mfexp(-top));
  pars_tmp(3)=mfexp(ascwid);
  pars_tmp(4)=mfexp(deswid);
  pars_tmp(5)=1.0/(1.0+mfexp(-init));
  pars_tmp(6)=1.0/(1.0+mfexp(-final));
  sel_tmp_iq(1)=mfexp(-(square(ages(1)-pars_tmp(1))/pars_tmp(3)));
  sel_tmp_iq(2)=mfexp(-(square(ages(nages)-pars_tmp(2))/pars_tmp(4)));
  sel_step1=mfexp(-(square(ages-pars_tmp(1))/pars_tmp(3)));
  sel_step2=pars_tmp(5)+(1.0-pars_tmp(5))*(sel_step1-sel_tmp_iq(1))/(1.0-sel_tmp_iq(1));  
  sel_step3=mfexp(-(square(ages-pars_tmp(2))/pars_tmp(4)));
  sel_step4=1.0+(pars_tmp(6)-1.0)*(sel_step3-1.0)/(sel_tmp_iq(2)-1.0);
  sel_step5=1.0/ (1.0+mfexp(-(20.0* elem_div((ages-pars_tmp(1)), (1.0+sfabs(ages-pars_tmp(1)))) )));
  sel_step6=1.0/(1.0+mfexp(-(20.0*elem_div((ages-pars_tmp(2)),(1.0+sfabs(ages-pars_tmp(2)))) )));  
  Sel_Tmp=elem_prod(sel_step2,(1.0-sel_step5))+ 
          elem_prod(sel_step5,((1.0-sel_step6)+ elem_prod(sel_step4,sel_step6)) ); 
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;
}

dvariable model_parameters::SR_func(const dvariable& R0, const dvariable& h, const dvariable& spr_F0, const dvariable& SSB, int func)
{
  //R0=virgin recruitment, h=steepness, spr_F0=spawners per recruit @ F=0, SSB=spawning biomass
  //func=1 for Beverton-Holt, 2 for Ricker
  RETURN_ARRAYS_INCREMENT();
  dvariable Recruits_Tmp;
  switch(func) {
    case 1: //Beverton-Holt
      Recruits_Tmp=((0.8*R0*h*SSB)/(0.2*R0*spr_F0*(1.0-h)+(h-0.2)*SSB));       
    break;
    case 2: //Ricker
      Recruits_Tmp=((SSB/spr_F0)*mfexp(h*(1-SSB/(R0*spr_F0))));       
    break;
  }
  RETURN_ARRAYS_DECREMENT();
  return Recruits_Tmp;
}

dvariable model_parameters::SR_eq_func(const dvariable& R0, const dvariable& h, const dvariable& spr_F0, const dvariable& spr_F, const dvariable& BC, int func)
{
  //R0=virgin recruitment, h=steepness, spr_F0=spawners per recruit @ F=0, spr_F=spawners per recruit @ F, BC=bias correction
  //func=1 for Beverton-Holt, 2 for Ricker
  RETURN_ARRAYS_INCREMENT();
  dvariable Recruits_Tmp;
  switch(func) {
    case 1: //Beverton-Holt
      Recruits_Tmp=(R0/((5.0*h-1.0)*spr_F))*(BC*4.0*h*spr_F-spr_F0*(1.0-h));    
    break;
    case 2: //Ricker
      Recruits_Tmp=R0/(spr_F/spr_F0)*(1.0+log(BC*spr_F/spr_F0)/h);      
    break;
  }
  RETURN_ARRAYS_DECREMENT();
  return Recruits_Tmp;
}

dvariable model_parameters::multinom_eff_N(const dvar_vector& pred_comp, const dvar_vector& obs_comp)
{
  //pred_comp=vector of predicted comps, obscomp=vector of observed comps
  dvariable EffN_Tmp; dvariable numer; dvariable denom;
  RETURN_ARRAYS_INCREMENT();
  numer=sum( elem_prod(pred_comp,(1.0-pred_comp)) );
  denom=sum( square(obs_comp-pred_comp) );
  if (denom>0.0) {EffN_Tmp=numer/denom;}
  else {EffN_Tmp=-missing;}                            
  RETURN_ARRAYS_DECREMENT();
  return EffN_Tmp;
}

dvariable model_parameters::lk_lognormal(const dvar_vector& pred, const dvar_vector& obs, const dvar_vector& cv, const dvariable& wgt_dat)
{
  //pred=vector of predicted vals, obs=vector of observed vals, cv=vector of CVs in arithmetic space, wgt_dat=constant scaling of CVs
  //small_number is small value to avoid log(0) during searfleet1
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvariable small_number=0.00001;
  dvar_vector var(cv.indexmin(),cv.indexmax()); //variance in log space
  var=log(1.0+square(cv/wgt_dat));   // convert cv in arithmetic space to variance in log space
  LkvalTmp=sum(0.5*elem_div(square(log(elem_div((pred+small_number),(obs+small_number)))),var) );
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;
}

dvariable model_parameters::lk_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const double& minSS, const dvariable& wgt_dat)
{
  //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, minSS=min N threshold, wgt_dat=scaling of N's
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvariable small_number=0.00001;
  LkvalTmp=0.0;
  for (int ii=1; ii<=ncomp; ii++)
  {if (nsamp(ii)>=minSS)
    {LkvalTmp-=wgt_dat*nsamp(ii)*sum(elem_prod((obs_comp(ii)+small_number),
               log(elem_div((pred_comp(ii)+small_number), (obs_comp(ii)+small_number)))));
    }
  }  
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;
}

dvariable model_parameters::lk_robust_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const dvariable& mbin, const double& minSS, const dvariable& wgt_dat)
{
  //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, mbin=number of bins, minSS=min N threshold, wgt_dat=scaling of N's
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvariable small_number=0.00001;
  LkvalTmp=0.0;
  dvar_matrix Eprime=elem_prod((1.0-obs_comp), obs_comp)+0.1/mbin; //E' of Francis 2011, p.1131  
  dvar_vector nsamp_wgt=nsamp*wgt_dat;
  for (int ii=1; ii<=ncomp; ii++)
  {if (nsamp(ii)>=minSS)
    {LkvalTmp+= sum(0.5*log(Eprime(ii))-log(small_number+mfexp(elem_div((-square(obs_comp(ii)-pred_comp(ii))) , (Eprime(ii)*2.0/nsamp_wgt(ii)) ))) );
    }
  }  
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;
}

dvariable model_parameters::lk_dirichlet_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const dvariable& mbin, const double& minSS, const dvariable& log_dir_par)
{
  //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, mbin=number of bins, minSS=min N threshold, wgt_dat=scaling of N's
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  LkvalTmp=0.0; 
  dvar_vector nsamp_adjust=nsamp*mfexp(log_dir_par);
  for (int ii=1; ii<=ncomp; ii++)
  {
	if (nsamp(ii)>=minSS)
    {
		LkvalTmp-=gammln(nsamp_adjust(ii))-gammln(nsamp(ii)+nsamp_adjust(ii));
		LkvalTmp-=sum(gammln(nsamp(ii)*obs_comp(ii)+nsamp_adjust(ii)*pred_comp(ii)));
        LkvalTmp+=sum(gammln(nsamp_adjust(ii)*pred_comp(ii)));		
    }
  }  
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;
}

dvariable model_parameters::lk_logistic_normal(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const dvariable& mbin, const double& minSS)
{
  //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, mbin=number of bins, minSS=min N threshold
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvariable small_number=0.0001;
  LkvalTmp=0.0;
  dvar_matrix nu=pred_comp+0.0;
  dvar_matrix pred_plus=pred_comp+small_number;
  dvar_matrix obs_plus=obs_comp+small_number;
  dvariable nu_mean;
  dvariable nu_sum_sq;
  dvariable tau_hat_sq;
  dvariable year_count; //keeps track of years included in likelihood (i.e., that meet the sample size requirement)
  LkvalTmp=0.0;
  nu_sum_sq=0.0;
  year_count=0.0;
  for (int ii=1; ii<=ncomp; ii++)
  {if (nsamp(ii)>=minSS)
    {
		year_count+=1.0;
		nu_mean=sum( log(obs_plus(ii))-log(pred_plus(ii))  )/mbin;	//year-specific mean log residual
		for (int jj=1; jj<=mbin;jj++)
		{
			nu(ii,jj) = log(obs_plus(ii,jj)) - log(pred_plus(ii,jj)) - nu_mean;
			nu_sum_sq += square(nu(ii,jj));
		}
    }
  }  
  if (year_count>0.0)
  {
	  tau_hat_sq = nu_sum_sq/((mbin-1.0)*year_count);
	  LkvalTmp = (mbin-1.0)*year_count*log(tau_hat_sq);
  }
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;
}

dvariable model_parameters::neg_log_prior(dvariable pred, const double& prior, dvariable var, int pdf)
{
  //prior=prior point estimate, var=variance (if negative, treated as CV in arithmetic space), pred=predicted value, pdf=prior type (1=none, 2=lognormal, 3=normal, 4=beta)
    dvariable LkvalTmp;
    dvariable alpha, beta, ab_iq;
    dvariable big_number=1e10;
    LkvalTmp=0.0;
    // compute generic pdf's
    switch(pdf) {
        case 1: //option to turn off prior
          LkvalTmp=0.0;
          break;
        case 2: // lognormal
          if(prior<=0.0) cout << "YIKES: Don't use a lognormal distn for a negative prior" << endl;
          else if(pred<=0) LkvalTmp=big_number=1e10;
          else {
            if(var<0.0) var=log(1.0+var*var) ;      // convert cv to variance on log scale
            LkvalTmp= 0.5*( square(log(pred/prior))/var + log(var) );
          }
        break;
        case 3: // normal
          if(var<0.0 && prior!=0.0) var=square(var*prior);       // convert cv to variance on observation scale
          else if(var<0.0 && prior==0.0) var=-var;               // cv not really appropriate if prior value equals zero
          LkvalTmp= 0.5*( square(pred-prior)/var + log(var) );
          break;
        case 4: // beta
          if(var<0.0) var=square(var*prior);          // convert cv to variance on observation scale
          if(prior<=0.0 || prior>=1.0) cout << "YIKES: Don't use a beta distn for a prior outside (0,1)" << endl;
          ab_iq=prior*(1.0-prior)/var - 1.0; alpha=prior*ab_iq; beta=(1.0-prior)*ab_iq;
          if(ab_iq<=0) {
			cout << "Parameter input error: For beta priors, mu*(1-mu)/var must be greater than one. Try decreasing var." << endl;
            exit(0);
		  }
		  if(pred>=0 && pred<=1) LkvalTmp= (1.0-alpha)*log(pred)+(1.0-beta)*log(1.0-pred)-gammln(alpha+beta)+gammln(alpha)+gammln(beta);
          else LkvalTmp=big_number;
          break;
        default: // no such prior pdf currently available
          cout << "Parameter input error: Prior must be either 1(none), 2(lognormal), 3(normal), or 4(beta)." << endl;
          cout << "Presently at least one is " << pdf << endl;
          exit(0);
    }
    return LkvalTmp;
}

dvariable model_parameters::sdnr_multinomial(const double& ncomp, const dvar_vector& ages, const dvar_vector& nsamp, 
                                    const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const dvariable& wgt_dat)
{
  //ncomp=number of years of data, ages=vector of ages, nsamp=vector of N's, 
  //pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, wgt_dat=likelihood weight for data source
  RETURN_ARRAYS_INCREMENT();
  dvariable SdnrTmp;
  dvar_vector o(1,ncomp);  
  dvar_vector p(1,ncomp);  
  dvar_vector ose(1,ncomp);  
  dvar_vector res(1,ncomp);
  SdnrTmp=0.0;
  for (int ii=1; ii<=ncomp; ii++)
  {
    o(ii)=sum(elem_prod(ages,obs_comp(ii)));
    p(ii)=sum(elem_prod(ages,pred_comp(ii)));
    ose(ii)=sqrt((sum(elem_prod(square(ages),pred_comp(ii)))-square(p(ii)))/(nsamp(ii)*wgt_dat));
  }
  res=elem_div((o-p),ose); 
  SdnrTmp=sqrt(sum(square(res-(sum(res)/ncomp))/(ncomp-1.0))); 
  RETURN_ARRAYS_DECREMENT();
  return SdnrTmp;
}

dvariable model_parameters::sdnr_lognormal(const dvar_vector& pred, const dvar_vector& obs, const dvar_vector& cv, const dvariable& wgt_dat)
{
  //nyr=number of years of data, pred=vector of predicted data, obs=vector of observed data, cv=vector of cv's, wgt_dat=likelihood weight for data source
  RETURN_ARRAYS_INCREMENT();
  dvariable SdnrTmp;
  dvariable small_number=0.00001;
  dvariable n;
  dvar_vector res(cv.indexmin(),cv.indexmax());
  SdnrTmp=0.0;
  res=elem_div(log(elem_div(obs+small_number,pred+small_number)),sqrt(log(1+square(cv/wgt_dat))));
  n=cv.indexmax()-cv.indexmin()+1;
  SdnrTmp=sqrt(sum(square(res-(sum(res)/n))/(n-1.0))); 
  RETURN_ARRAYS_DECREMENT();
  return SdnrTmp; 
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
      cout << "TotalLikelihood " << fval << endl;
	  //cout << "xdum = " << xdum << endl;
    if (last_phase())  
    {
      //cout<<"start report"<<endl;
       get_weighted_current();
      //cout<<"got weighted"<<endl;
       get_msy();
      //cout<<"got msy"<<endl;
	  get_per_recruit_stuff();
      //cout<<"got per recruit"<<endl;  
      get_miscellaneous_stuff();
      //cout<<"got misc stuff"<<endl;
	  get_projection();
      //cout<<"got projection"<<endl;
      // get_effective_sample_sizes();
      //cout<<"got effective_sample_sizes"<<endl;
      grad_max=objective_function_value::pobjfun->gmax;
      time(&finish);
	  elapsed_time=difftime(finish,start);
	  hour=long(elapsed_time)/3600;
	  minute=long(elapsed_time)%3600/60;
	  second=(long(elapsed_time)%3600)%60;
	  cout<<endl<<endl<<"*******************************************"<<endl;
	  cout<<"--Start time: "<<ctime(&start)<<endl;
	  cout<<"--Finish time: "<<ctime(&finish)<<endl;
	  cout<<"--Runtime: ";
	  cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	  cout<<"*******************************************"<<endl;
      cout <<endl;     
      cout << "><>--><>--><>--><>--><>--><>--><>--><>--><>--><>"  <<endl;
      cout << "Fmsy=" << F_msy_out<< "   SSBmsy=" << SSB_msy_out <<endl;
      cout <<"F status="<<FdF_msy_end<<endl;
      cout <<"Pop status="<<SdSSB_msy_end<<endl;
      cout << "h="<<steep<<"   R0="<<R0<<endl;
      cout << "><>--><>--><>--><>--><>--><>--><>--><>--><>--><>"  <<endl;  
      //report << "TotalLikelihood " << fval << endl;
      //report << "N" << endl;
      //report << N<<endl;
      //report << "F" << endl;
      //report << F <<endl;   
      sdnr_ac_survey1=sdnr_multinomial(nyr_survey1_agec, agebins_agec, nsamp_survey1_agec, pred_survey1_agec, obs_survey1_agec, w_ac_survey1);
      sdnr_ac_survey2=sdnr_multinomial(nyr_survey2_agec, agebins_agec, nsamp_survey2_agec, pred_survey2_agec, obs_survey2_agec, w_ac_survey2);
      sdnr_ac_fleet1=sdnr_multinomial(nyr_fleet1_agec, agebins_agec, nsamp_fleet1_agec, pred_fleet1_agec, obs_fleet1_agec, w_ac_fleet1);
      sdnr_I_survey1=sdnr_lognormal(pred_survey1_cpue, obs_survey1_cpue, survey1_cpue_cv, w_I_survey1);
      sdnr_I_survey2=sdnr_lognormal(pred_survey2_cpue, obs_survey2_cpue, survey2_cpue_cv, w_I_survey2);
      //#################################################################################################
      //##  Pass parameters to vector for bounds check plotting
      //################################################################################################# 
       Linf_out(8)=Linf; Linf_out(1,7)=set_Linf; 
       K_out(8)=K; K_out(1,7)=set_K;
       t0_out(8)=t0; t0_out(1,7)=set_t0;
       len_cv_val_out(8)=len_cv_val; len_cv_val_out(1,7)=set_len_cv;
       log_R0_out(8)=log_R0; log_R0_out(1,7)=set_log_R0;
       steep_out(8)=steep; steep_out(1,7)=set_steep;
       rec_sigma_out(8)=rec_sigma; rec_sigma_out(1,7)=set_rec_sigma;
       R_autocorr_out(8)=R_autocorr; R_autocorr_out(1,7)=set_R_autocorr;
       log_dm_fleet1_ac_out(8)=log_dm_fleet1_ac; log_dm_fleet1_ac_out(1,7)=set_log_dm_fleet1_ac;
       log_dm_survey1_ac_out(8)=log_dm_survey1_ac; log_dm_survey1_ac_out(1,7)=set_log_dm_survey1_ac;
       log_dm_survey2_ac_out(8)=log_dm_survey2_ac; log_dm_survey2_ac_out(1,7)=set_log_dm_survey2_ac;
       selpar_A50_survey1_out(8)=selpar_A50_survey1; selpar_A50_survey1_out(1,7)=set_selpar_A50_survey1;
       selpar_slope_survey1_out(8)=selpar_slope_survey1; selpar_slope_survey1_out(1,7)=set_selpar_slope_survey1;
       selpar_A50_survey2_out(8)=selpar_A50_survey2; selpar_A50_survey2_out(1,7)=set_selpar_A50_survey2;
       selpar_slope_survey2_out(8)=selpar_slope_survey2; selpar_slope_survey2_out(1,7)=set_selpar_slope_survey2;
       selpar_A50_fleet1_B1_out(8)=selpar_A50_fleet1_B1; selpar_A50_fleet1_B1_out(1,7)=set_selpar_A50_fleet1_B1;
       selpar_slope_fleet1_B1_out(8)=selpar_slope_fleet1_B1; selpar_slope_fleet1_B1_out(1,7)=set_selpar_slope_fleet1_B1;
       selpar_A50_fleet1_B2_out(8)=selpar_A50_fleet1_B2; selpar_A50_fleet1_B2_out(1,7)=set_selpar_A50_fleet1_B2;
       selpar_slope_fleet1_B2_out(8)=selpar_slope_fleet1_B2; selpar_slope_fleet1_B2_out(1,7)=set_selpar_slope_fleet1_B2;
       log_q_survey1_out(8)=log_q_survey1; log_q_survey1_out(1,7)=set_log_q_survey1;
       log_q_survey2_out(8)=log_q_survey2; log_q_survey2_out(1,7)=set_log_q_survey2;
       log_avg_F_fleet1_out(8)=log_avg_F_fleet1; log_avg_F_fleet1_out(1,7)=set_log_avg_F_fleet1;
       F_init_out(8)=F_init; F_init_out(1,7)=set_F_init;
       log_rec_dev_out(styr_rec_dev, endyr_rec_dev)=log_rec_dev;
       log_F_dev_fleet1_out(styr_fleet1_L,endyr_fleet1_L)=log_F_dev_fleet1;
      #include "BAM-Sim.cxx"   // write the R-compatible report
    } //endl last phase loop     
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  time(&start);
  arrmblsize=20000000;
  gradient_structure::set_MAX_NVAR_OFFSET(1600);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(2000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(2000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(10000);
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
