// Create a file with an R object from AD Model Builder

// open the file using the default AD Model Builder file name, and 
// 6 digits of precision
open_r_file(adprogram_name + ".rdat", 6);

// Example of an INFO object
open_r_info_list("info", true);
	wrt_r_item("title", "Simulated data");
	wrt_r_item("species", "Sim");
	wrt_r_item("model", "Statistical Catch at Age");
	if(SR_switch==1)
	    {wrt_r_item("rec.model", "BH-steep");}
    if(SR_switch==2)
	    {wrt_r_item("rec.model", "Ricker-steep");}	
	wrt_r_item("base.run", "BAM-Sim1.tpl");  	
	wrt_r_item("units.length", "mm");	
	wrt_r_item("units.weight", "kg");		
	wrt_r_item("units.biomass", "mt");	
	wrt_r_item("units.ssb", "mt");		
	wrt_r_item("units.ypr", "kg"); 			
	wrt_r_item("units.landings", "mt");	
	//wrt_r_item("units.discards", "1000 dead fish");	
    wrt_r_item("units.numbers", "number fish");		
    wrt_r_item("units.naa", "number fish");
	wrt_r_item("units.rec", "number fish");	    	
close_r_info_list();

// VECTOR object of parameters and estimated quantities
open_r_info_list("parms", false); 
	wrt_r_item("styr", styr);
	wrt_r_item("endyr", endyr);
	wrt_r_item("styrR", styr_rec_dev);
	wrt_r_item("endyrR", endyr_rec_dev);
	wrt_r_item("endyrR.phase1", endyr_rec_phase1);
	wrt_r_item("endyrR.phase2", endyr_rec_phase2);
	wrt_r_item("endyr.sel.phase1", endyr_selex_phase1);
	wrt_r_item("max.F.spr.msy", max_F_spr_msy);
	wrt_r_item("selpar.n.yrs.wgted", selpar_n_yrs_wgted);
	
	wrt_r_item("minSS.agec", minSS_agec);
	
    wrt_r_item("spawn.time", spawn_time_frac);
	wrt_r_item("q.survey1", mfexp(log_q_survey1));     	       	    
	    
    wrt_r_item("F.prop.fleet1", F_fleet1_prop);			
	wrt_r_item("F.init.scale", F_init);
	wrt_r_item("log.F.avg.fleet1", log_avg_F_fleet1);
	wrt_r_item("F.tune", set_Ftune);
	wrt_r_item("F.tune.yr", set_Ftune_yr);
	
	wrt_r_item("B0", B0);
	wrt_r_item("Bstyr.B0", totB(styr)/B0);
	wrt_r_item("SSB0", S0);							
	wrt_r_item("SSBstyr.SSB0", SSB(styr)/S0);
	wrt_r_item("Rstyr.R0", rec(styr)/R0);
	wrt_r_item("SR.switch", SR_switch);
	if(SR_switch==1)
	{
      wrt_r_item("BH.biascorr",BiasCor);
	  wrt_r_item("BH.Phi0", spr_F0);	
	  wrt_r_item("BH.R0", R0);	
	  wrt_r_item("BH.steep", steep);
    }
    if(SR_switch==2)
	{
      wrt_r_item("Ricker.biascorr",BiasCor);
	  wrt_r_item("Ricker.Phi0", spr_F0);	
	  wrt_r_item("Ricker.R0", R0);	
	  wrt_r_item("Ricker.steep", steep);
    }	
	wrt_r_item("R.sigma.logdevs", sigma_rec_dev); 
	wrt_r_item("R.sigma.par",rec_sigma);
	wrt_r_item("R.autocorr",R_autocorr);
	wrt_r_item("R0", R0); //same as BH.R0, but used in BSR.time.plots	
	wrt_r_item("R.virgin.bc", R_virgin); //bias-corrected virgin recruitment			
	wrt_r_item("rec.lag", 1.0);
	wrt_r_item("set.BiasCor", set_BiasCor);		
	wrt_r_item("msy.mt", msy_mt_out);
	wrt_r_item("msy.knum", msy_knum_out);	
	wrt_r_item("Fmsy", F_msy_out);
	wrt_r_item("SSBmsy", SSB_msy_out);
	wrt_r_item("msst", smsy2msst*SSB_msy_out);	//For FishGraph; if using Fmsy benchmark
	wrt_r_item("Bmsy", B_msy_out);
	wrt_r_item("Rmsy", R_msy_out);
	wrt_r_item("sprmsy",spr_msy_out);
	wrt_r_item("Fend.Fmsy", FdF_msy_end);
	wrt_r_item("Fend.Fmsy.mean", FdF_msy_end_mean);	
	wrt_r_item("SSBend.SSBmsy", SdSSB_msy_end);
	wrt_r_item("SSBend.MSST", SdSSB_msy_end/smsy2msst);
	wrt_r_item("SSBend.MSST.75", SdSSB_msy_end/smsy2msst75);
	wrt_r_item("F35", F35_out); //For FishGraph
	wrt_r_item("SSB.F35", SSB_F35_out); //For FishGraph
	wrt_r_item("msst.F35", smsy2msst*SSB_F35_out);
	//wrt_r_item("msst", smsy2msst*SSB_F35_out); //For FishGraph; if using F35 benchmark
close_r_info_list();

 // VECTOR object of parameters and estimated quantities
open_r_info_list("spr.brps", false);
	wrt_r_item("F30", F30_out);
	wrt_r_item("F35", F35_out);
	wrt_r_item("F40", F40_out);
	wrt_r_item("SSB.F35", SSB_F35_out);
    wrt_r_item("msst.F35", smsy2msst*SSB_F35_out);   
    wrt_r_item("B.F35", B_F35_out);
	wrt_r_item("R.F35", R_F35_out);
    wrt_r_item("L.F35.knum", L_F35_knum_out);
    wrt_r_item("L.F35.mt", L_F35_mt_out);
    wrt_r_item("Fend.F35.mean", FdF35_end_mean);
    wrt_r_item("SSBend.SSBF35", SdSSB_F35_end);
	wrt_r_item("SSBend.MSSTF35", Sdmsst_F35_end);	
close_r_info_list();

// MATRIX object of parameter constraints
open_r_df("parm.cons",1,8,2);
    wrt_r_namevector(1,8);
 
    wrt_r_df_col("Linf",Linf_out);
	wrt_r_df_col("L",K_out);
	wrt_r_df_col("t0",t0_out);
	wrt_r_df_col("len.cv",len_cv_val_out);
	
    wrt_r_df_col("log_R0",log_R0_out);
    wrt_r_df_col("steep",steep_out);
    wrt_r_df_col("rec_sigma",rec_sigma_out);
    wrt_r_df_col("R_autocorr",R_autocorr_out);    

	wrt_r_df_col("log_dm_fleet1_ac",log_dm_fleet1_ac_out);
	wrt_r_df_col("log_dm_survey1_ac",log_dm_survey1_ac_out);
    
	wrt_r_df_col("selpar_A50_survey1",selpar_A50_survey1_out);
    wrt_r_df_col("selpar_slope_survey1",selpar_slope_survey1_out);
	
	wrt_r_df_col("selpar_A50_fleet1_B1",selpar_A50_fleet1_B1_out);
    wrt_r_df_col("selpar_slope_fleet1_B1",selpar_slope_fleet1_B1_out);
 	wrt_r_df_col("selpar_A50_fleet1_B2",selpar_A50_fleet1_B2_out);
    wrt_r_df_col("selpar_slope_fleet1_B2",selpar_slope_fleet1_B2_out);
   
    wrt_r_df_col("log_q_survey1",log_q_survey1_out);
	
    wrt_r_df_col("F_init.scale",F_init_out);    
    wrt_r_df_col("log_avg_F_fleet1",log_avg_F_fleet1_out);
       
close_r_df();

// DATA FRAME of time series deviation vector estimates
// names used in this object must match the names used in the "parm.tvec.cons" object
open_r_df("parm.tvec", styr, endyr, 2);
	wrt_r_namevector(styr,endyr);
	wrt_r_df_col("year", styr,endyr);
    wrt_r_df_col("log.rec.dev", log_rec_dev_out);
    wrt_r_df_col("log.F.dev.fleet1", log_F_dev_fleet1_out);
close_r_df();

// MATRIX object of deviation vector constraints
// names used in this object must match the names used in the "parm.vec" object
open_r_df("parm.tvec.cons",1,3,2);
    wrt_r_namevector(1,3);
    wrt_r_df_col("log.rec.dev",set_log_rec_dev);
    wrt_r_df_col("log.F.dev.fleet1",set_log_F_dev_fleet1);
close_r_df();

// DATA FRAME of age vector deviation estimates
// names used in this object must match the names used in the "parm.avec.cons" object
open_r_df("parm.avec", 1, nages, 2);
	wrt_r_namevector(1,nages);
	wrt_r_df_col("age", agebins); //deviations for first age not estimated
    wrt_r_df_col("log.Nage.dev", log_Nage_dev_output);
close_r_df();

// MATRIX object of age vector deviation constraints
// names used in this object must match the names used in the "parm.avec" object
open_r_df("parm.avec.cons",1,3,2);
    wrt_r_namevector(1,3);
    wrt_r_df_col("log.Nage.dev",set_log_Nage_dev);
close_r_df();

// VECTOR object of SDNR calculations
open_r_vector("sdnr"); 
    wrt_r_item("sdnr.U.survey1", sdnr_I_survey1);
    wrt_r_item("sdnr.agec.survey1", sdnr_ac_survey1);  
    wrt_r_item("sdnr.agec.fleet1", sdnr_ac_fleet1);    
close_r_vector();    


// VECTOR object of likelihood contributions
open_r_vector("like"); 
	wrt_r_item("lk.total", fval); //weighted likelihood
    wrt_r_item("lk.unwgt.data", fval_data);  //likelihood of just data components

    wrt_r_item("lk.L.fleet1", f_fleet1_L);  
    wrt_r_item("lk.U.survey1", f_survey1_cpue);    
    wrt_r_item("lk.agec.survey1", f_survey1_agec);	
    wrt_r_item("lk.agec.fleet1", f_fleet1_agec);	
    wrt_r_item("lk.Nage.init", f_Nage_init);
    wrt_r_item("lk.SRfit", f_rec_dev);
    wrt_r_item("lk.SRearly", f_rec_dev_early);
    wrt_r_item("lk.SRend", f_rec_dev_end);        
    wrt_r_item("lk.Ftune", f_Ftune);
    wrt_r_item("lk.priors",f_priors);
 	
 	wrt_r_item("w.L", w_L); 
    wrt_r_item("w.U.survey1", w_I_survey1);
    wrt_r_item("w.agec.survey1", w_ac_survey1);	
    wrt_r_item("w.agec.fleet1", w_ac_fleet1);	    
    wrt_r_item("w.Nage.init", w_Nage_init);
    wrt_r_item("w.SRfit", w_rec);
    wrt_r_item("w.SRearly", w_rec_early);
    wrt_r_item("w.SRend", w_rec_end);        
    wrt_r_item("w.Ftune", w_Ftune);
    
    wrt_r_item("grad.max",grad_max);      
       			
close_r_vector();
	
	open_r_matrix("N.age");
    wrt_r_matrix(N, 2, 2);
    wrt_r_namevector(styr, (endyr+1));
    wrt_r_namevector(agebins);    
    close_r_matrix();
    
    open_r_matrix("N.age.mdyr");
    wrt_r_matrix(N_mdyr, 2, 2);
    wrt_r_namevector(styr, endyr);
    wrt_r_namevector(agebins);   
    close_r_matrix();
	
    open_r_matrix("N.age.spawn");
    wrt_r_matrix(N_spawn, 2, 2);
    wrt_r_namevector(styr, endyr);
    wrt_r_namevector(agebins);   
    close_r_matrix();
    
    open_r_matrix("B.age");
    wrt_r_matrix(B, 2, 2);
    wrt_r_namevector(styr, (endyr+1));
    wrt_r_namevector(agebins);   
    close_r_matrix();

    open_r_matrix("F.age");
    wrt_r_matrix(F, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);   
    close_r_matrix();
        
    open_r_matrix("Z.age");
    wrt_r_matrix(Z, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);   
    close_r_matrix();
    
    open_r_matrix("L.age.pred.num");
    wrt_r_matrix(L_total_num, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);   
    close_r_matrix();
    
    open_r_matrix("L.age.pred.mt");
    wrt_r_matrix(L_total_mt, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);   
    close_r_matrix();

// LIST object with annual selectivity at age by fishery

open_r_list("sel.age");

    wrt_r_complete_vector("sel.v.wgted.L",sel_wgted_L, agebins);
    wrt_r_complete_vector("sel.v.wgted.tot",sel_wgted_tot, agebins);

    open_r_matrix("sel.m.survey1");
    wrt_r_matrix(sel_survey1, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);   
    close_r_matrix();  

    open_r_matrix("sel.m.fleet1");
    wrt_r_matrix(sel_fleet1, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);   
    close_r_matrix();  
	
close_r_list();


//LIST object with predicted and observed composition data
open_r_list("comp.mats");
	
    open_r_matrix("acomp.survey1.ob");		
    wrt_r_matrix(obs_survey1_agec, 2, 2);
    wrt_r_namevector(yrs_survey1_agec);
    wrt_r_namevector(agebins_agec);
    close_r_matrix();
    
    open_r_matrix("acomp.survey1.pr");		
    wrt_r_matrix(pred_survey1_agec, 2, 2);
    wrt_r_namevector(yrs_survey1_agec);
    wrt_r_namevector(agebins_agec);
    close_r_matrix();  
	
    open_r_matrix("acomp.fleet1.ob");		
    wrt_r_matrix(obs_fleet1_agec, 2, 2);
    wrt_r_namevector(yrs_fleet1_agec);
    wrt_r_namevector(agebins_agec);
    close_r_matrix();
    
    open_r_matrix("acomp.fleet1.pr");		
    wrt_r_matrix(pred_fleet1_agec, 2, 2);
    wrt_r_namevector(yrs_fleet1_agec);
    wrt_r_namevector(agebins_agec);
    close_r_matrix();  


close_r_list();

// DATA FRAME of time series
open_r_df("t.series", styr, (endyr+1), 2);
	wrt_r_namevector(styr,(endyr+1));
	wrt_r_df_col("year", styr,(endyr+1));
	wrt_r_df_col("F.Fmsy", FdF_msy);
	wrt_r_df_col("F.F35.ratio", FdF35);	//*.ratio extension is for FishGraph
	wrt_r_df_col("F.full", Fapex);		
	wrt_r_df_col("F.fleet1", F_fleet1_out);
	wrt_r_df_col("Fsum", Fsum);	
    wrt_r_df_col("N", totN); //abundance at start of year			    
    wrt_r_df_col("recruits", rec);
    wrt_r_df_col("logR.dev", log_rec_dev_output); //places zeros in yrs deviations not estimated  
    wrt_r_df_col("SSB", SSB);
    wrt_r_df_col("SSB.SSBmsy", SdSSB_msy);
    wrt_r_df_col("SSB.msst", SdSSB_msy/smsy2msst);
    wrt_r_df_col("SSB.msst.75", SdSSB_msy/smsy2msst75);
	wrt_r_df_col("SSB.SSBF35", SdSSB_F35);
    wrt_r_df_col("SSB.msstF35", Sdmsst_F35);
    wrt_r_df_col("B", totB);
    wrt_r_df_col("B.B0", totB/B0);
    wrt_r_df_col("MatFemB", MatFemB);    
    wrt_r_df_col("SPR.static", spr_static);
                                    
    wrt_r_df_col("total.L.mt", L_total_mt_yr);
    wrt_r_df_col("total.L.knum", L_total_knum_yr);
	
    wrt_r_df_col("U.survey1.ob", obs_survey1_cpue_allyr);
    wrt_r_df_col("U.survey1.pr", pred_survey1_cpue_allyr);
    wrt_r_df_col("cv.U.survey1", survey1_cpue_cv_allyr);    
 	
    wrt_r_df_col("L.fleet1.ob", obs_fleet1_L);			
    wrt_r_df_col("L.fleet1.pr", pred_fleet1_L_mt);  		
    wrt_r_df_col("cv.L.fleet1", fleet1_L_cv);		

    //comp sample sizes      
    wrt_r_df_col("acomp.survey1.n", nsamp_survey1_agec_allyr);   
    wrt_r_df_col("acomp.fleet1.n", nsamp_fleet1_agec_allyr);    
     
    wrt_r_df_col("acomp.survey1.nfish", nfish_survey1_agec_allyr); 
    wrt_r_df_col("acomp.fleet1.nfish", nfish_fleet1_agec_allyr);       //wrt_r_df_col("acomp.fleet1.D.nfish", nfish_fleet1_D_agec_allyr);    	

    wrt_r_df_col("acomp.survey1.neff",  w_ac_survey1*nsamp_survey1_agec_allyr);   
    wrt_r_df_col("acomp.fleet1.neff", w_ac_fleet1*nsamp_fleet1_agec_allyr); 
    //wrt_r_df_col("acomp.survey1.neff", (1+nsamp_survey1_agec_allyr*exp(log_dm_survey1_ac_out(8)))/(1+exp(log_dm_survey1_ac_out(8))) );   
    //wrt_r_df_col("acomp.fleet1.neff", (1+nsamp_fleet1_agec_allyr*exp(log_dm_fleet1_ac_out(8)))/(1+exp(log_dm_fleet1_ac_out(8))) );      
    
close_r_df();

// DATA FRAME of L and D time series by fishery
open_r_df("LD.pr.tseries", styr, endyr, 2);
	wrt_r_namevector(styr,endyr);
	wrt_r_df_col("year", styr, endyr);

    wrt_r_df_col("L.fleet1.mt", pred_fleet1_L_mt);   
    wrt_r_df_col("L.fleet1.knum", pred_fleet1_L_knum);
                  		
close_r_df();	


open_r_df("a.series", 1, nages, 2);
	wrt_r_namevector(1,nages);
	wrt_r_df_col("age", agebins);
	wrt_r_df_col("weight", wgt_kg);     //for FishGraph
    wrt_r_df_col("length", meanlen_TL);
	wrt_r_df_col("wgt.wgted.L.mt", wgt_wgted_L_mt);	
	wrt_r_df_col("wgt.mt.spawn.eq", wgt_mt);		
	wrt_r_df_col("reprod", reprod);	
	wrt_r_df_col("M", M);	
	wrt_r_df_col("mat.female", maturity_f);
	wrt_r_df_col("prop.female", prop_f);
	wrt_r_df_col("F.initial", F_initial);
	wrt_r_df_col("Z.initial", Z_initial);
	wrt_r_df_col("Nage.eq.init",N_initial_eq);
	wrt_r_df_col("log.Nage.init.dev",log_Nage_dev_output);	
close_r_df();	

open_r_df("eq.series", 1, n_iter_msy, 2);
	wrt_r_namevector(1,n_iter_msy);
	wrt_r_df_col("F.eq", F_msy);
	wrt_r_df_col("spr.eq", spr_msy);
	wrt_r_df_col("R.eq", R_eq);
	wrt_r_df_col("SSB.eq", SSB_eq);
	wrt_r_df_col("B.eq", B_eq);
	wrt_r_df_col("L.eq.mt", L_eq_mt);	
    wrt_r_df_col("L.eq.knum", L_eq_knum);
close_r_df();

open_r_df("pr.series", 1, n_iter_spr, 2);
	wrt_r_namevector(1,n_iter_spr);
	wrt_r_df_col("F.spr", F_spr);
	wrt_r_df_col("spr", spr_spr);
	wrt_r_df_col("SPR", spr_ratio);
	wrt_r_df_col("ypr.kg", L_spr);
close_r_df();


open_r_list("CLD.est.mats");

    open_r_matrix("Lw.fleet1");			
       wrt_r_matrix(L_fleet1_mt, 1,1);		
    close_r_matrix();
    
    open_r_matrix("Ln.fleet1");			
       wrt_r_matrix(L_fleet1_num, 1,1);
    close_r_matrix();

	open_r_matrix("Lw.total");
    wrt_r_matrix(L_total_mt, 1,1);
    close_r_matrix();

    open_r_matrix("Ln.total");
        wrt_r_matrix(L_total_num, 1,1);
    close_r_matrix();
		
close_r_list();


//LIST object with predicted and observed composition data
open_r_list("age.error");

    open_r_matrix("error.mat");		
    wrt_r_matrix(age_error, 2, 2);
    wrt_r_namevector(agebins);
    wrt_r_namevector(agebins);
    close_r_matrix();
  
close_r_list();

// DATA FRAME of time series
open_r_df("projection", styr_proj, endyr_proj, 2);
    wrt_r_namevector(styr_proj,endyr_proj);
    wrt_r_df_col("year", styr_proj,endyr_proj);
    wrt_r_df_col("F.proj", F_proj);
    wrt_r_df_col("SSB.proj", SSB_proj);
    wrt_r_df_col("B.proj", B_proj);
    wrt_r_df_col("R.proj", R_proj);
	wrt_r_df_col("L.knum.proj", L_knum_proj);
	wrt_r_df_col("L.mt.proj", L_mt_proj);
 close_r_df();

close_r_file();
