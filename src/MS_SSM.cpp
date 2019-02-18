#include <TMB.hpp> 
#include <iostream> // Header that defines the standard input/output stream objects (cin, cout, etc.)
#include <cmath> /// for functions such as power pow()  
#include "MS_SSM_functions.h" // File that contains the functions used within this script

#define see(object) std::cout << #object ":\n" << object << "\n"; // To print things of interest

  
  
 
///////////////////////////////////////// START OF OBJECTIVE FUNCTION ///////////////////////////////////////////////////////

  
template <class Type>   
Type objective_function<Type>::operator() () 
{
  // Input data
  DATA_INTEGER(sp); // Number of prey species  
  DATA_MATRIX(interaction); // Interactions matrix (sp,sp) (row=prey,col=predator)
  DATA_INTEGER(Y); // Last year position for optimization
  //DATA_IVECTOR(Y1); // First year position per species
  DATA_IVECTOR(Aplus); // Age max per prey species (age plus group)
  DATA_INTEGER(max_A); // Age max for all prey species, max(A)
  DATA_IVECTOR(max_A_catch); // Age max for catch
  DATA_IVECTOR(min_A_catch); // Age min for catch
  DATA_ARRAY(M); // Natural mortality
  DATA_ARRAY(wc); // Weight at age in catch (kg)
  DATA_IVECTOR(n_surv); // Number of surveys per species
  DATA_INTEGER(n_surv_max); // Total number of different surveys
  DATA_MATRIX(prop_y_elapsed_surv); // proportion of year elapsed for surveys (Y,n_surv_max)
  DATA_MATRIX(obs_aggr_Cw); // Input aggregated catch in tons
  DATA_MATRIX(sd_obs_aggr_Cw); // Standard deviation for aggregated catch
  DATA_MATRIX(flag_aggr_Cw); // flags which say if annual aggregated catch estimates should be used or not for specific species (Y,sp)
  DATA_ARRAY(obs_aggr_I); // Input aggregated indices (Y,sp,n_surv_max)
  DATA_ARRAY(sd_obs_aggr_I); // Standard deviation of aggregated indices (Y,sp,n_surv_max)
  DATA_ARRAY(flag_aggr_I); // flags which say if annual aggregated index estimates should be used or not for specific species (Y,sp,n_surv_max)
  DATA_IMATRIX(flag_Caa); // flags which say if annual catch at age should be used or not for specific species (Y,sp)
  DATA_ARRAY(obs_prop_Caa); // Observed age composition (prop) in catch (Y,max_A,sp)
  DATA_ARRAY(flag_Iaa); // flags which say if annual index at age estimates should be used or not for specific species (Y,sp,n_surv_max)
  DATA_MATRIX(Neff_C); // Effective sample size for age composition catch (Y,sp)
  DATA_IVECTOR(age_comp_model_catch); // Index for choice likelihood distribution for age composition catch (sp)
  DATA_IMATRIX(catch_aref); // last positive age class (Y,sp)
  DATA_ARRAY(obs_prop_Iaa); // Observed age composition (prop) in surveys (Y,max_A,sp,n_surv_max)
  DATA_ARRAY(Neff_surv); // Effective sample size for age composition abundance indices (Y,sp,n_surv_max)
  DATA_IVECTOR(age_comp_model_indices); // Index for choice likelihood distribution for age composition catch (sp)
  DATA_IARRAY(index_aref); // last positive age class (Y,sp,n_surv_max)
  DATA_INTEGER(process_rec); // 0=no process error on recruitment, 1=process error on recruitment
  DATA_INTEGER(process_survival); // 0=no process error on recruitment, 1=process error on recruitment
  DATA_INTEGER(recruit_model); // if process error on, 1=random walk, 2=random about the mean
  DATA_IVECTOR(M_model);
  DATA_INTEGER(process_M);
  DATA_VECTOR(scale_M_upper); // max scaling factor on M
  DATA_IMATRIX(sel_model_surv); // selectivity option survey 1=selectivity at age, 2=logistic
  DATA_IMATRIX(min_A_surv); // Min age to consider for survey selectivity
  DATA_IMATRIX(max_A_surv); // Max age to consider for survey selectivity
  DATA_INTEGER(data_simulate); // 1=data simulated
  DATA_INTEGER(error_simulate); // 1=process errors simulated
  DATA_ARRAY(wSSB); // weight at age in SSB in kg
  DATA_ARRAY(mature); // proportion of mature at age
  DATA_VECTOR(prop_y_elapsed_SSB); // proportion of year elapsed before spawning
  DATA_IVECTOR(sel_model_F); //selectivity option catch 1=selectivity at age, 2=logistic
  // For trophic interactions
  DATA_INTEGER(n_prey); // sp+1 for Other food
  DATA_INTEGER(n_pred);
  DATA_IVECTOR(Bplus); // Age max per predator species (age plus group)
  DATA_INTEGER(max_B); // Age max for all predator species, max(B)
  DATA_VECTOR(log_shape_gamma_dat); // Log of shape parameter for gamma distribution size preference (n_pred)
  DATA_VECTOR(log_scale_gamma_dat); // Log of scale parameter for gamma distribution size preference (n_pred)
  DATA_INTEGER(gamma_pref_estim); // Option for estimation gamma distribution for size pref, 0=no estimation, 1=estimation
  DATA_INTEGER(n_ratio); // Size of log_ratio_w_diet 
  DATA_MATRIX(log_ratio_w_diet); // Log ratio w predator over w prey from food habits data (n_ratio,n_pred)
  //DATA_SCALAR(B_ecosystem); // Total biomass of ecosystem
  DATA_SCALAR(B_other); // Biomass other food in tons
  DATA_ARRAY(ratio_diet); // Proportion of prey i in diet of predator j of age b (n_stom_max,Y,sp,n_pred)
  DATA_INTEGER(cons_rate_estim); // Option for estimation normal distribution for log consumption rates, 0=no estimation, 1=estimation
  //DATA_INTEGER(n_cons_rate); // number of observed consumption estimates per predator (n_pred)
  //DATA_MATRIX(log_obs_cons_rate); // Log of consumption rates empirically estimated 
  //DATA_IMATRIX(flag_cons_rate); // Flag to avoid the NAs in observed consumption rates (n_cons_rate,n_pred)
  DATA_INTEGER(Y_cons); // Last year position for optimization
  DATA_ARRAY(spring_cons); // array (Y,2,n_pred), 1st column=meansw (sum meansw for 10 sp), 2nd column=bottemp
  DATA_ARRAY(fall_cons); // array (Y,2,n_pred), 1st column=meansw (sum meansw for 10 sp), 2nd column=bottemp
  DATA_ARRAY(spring_cons_other); // array (Y,2,n_pred), 1st column=meansw (sum meansw for 10 sp), 2nd column=bottemp
  DATA_ARRAY(fall_cons_other); // array (Y,2,n_pred), 1st column=meansw (sum meansw for 10 sp), 2nd column=bottemp
  DATA_IVECTOR(flag_cons_rate); // 0 if NA in meansw or bottemp, 1 otherwise
  DATA_INTEGER(predation_on); // Are trophic interactions on? 0=no, 1=yes
  DATA_IMATRIX(n_stom); // Number of stomach in diet data (Y,n_pred)
  DATA_INTEGER(n_stom_max); // Max number of stomach for all predators
  DATA_ARRAY(prob_l_given_b); // Probabilty of predator length l given its age b (max(n_stom),Y,max_B,n_pred)
  DATA_INTEGER(biomass_other_option); // Option to model biomasss_other, 1= cst, 2=surplus production
  DATA_INTEGER(diet_model); // ddeltadir=1, ddirichet=2
  DATA_IARRAY(age_pred); // Age of predator for each stomach in diet data (n_stom_max,Y,n_pred)
  DATA_INTEGER(flag_nll_diet); // 1=initial model with probs, 2=new simplified model per stomach, 3=highly simplified model aggregated over stomach
  DATA_ARRAY(ratio_diet3); // Proportion of prey i in diet of predator j of age b for flag_nll_diet=3 (Y,sp,b,n_pred)
  DATA_IVECTOR(t1);
  DATA_IVECTOR(t2);
  DATA_VECTOR(length_t);
  DATA_INTEGER(xmax);
  DATA_INTEGER(process_F);
  DATA_IVECTOR(min_Fbar_idx); // minimum age use for calculation mean mortality rates
  DATA_IVECTOR(max_Fbar_idx); // maximum age use for calculation mean mortality rates
  DATA_VECTOR(Fbar_range); // length age range use for calculation mean mortality rates
  
  
  
  // Parameters to estimate
  PARAMETER_VECTOR(logit_gamma_F); // log slope param for s_F 
  PARAMETER_VECTOR(logit_A50_F); // logit age at 50% selectivity for calculation of s_F
  PARAMETER_MATRIX(log_N1); // Log initial numbers at age in year 1 (numbers)
  PARAMETER_MATRIX(log_rec); // Log annual recruitment year 2 to Y (numbers) (Y-1,sp)
  PARAMETER_MATRIX(log_E); // Log year component of F
  PARAMETER_ARRAY(log_NAA); // Log fish numbers for process error on N (Y-1,max_A-1,sp) because no ages at Y1 and rec at age 1
  PARAMETER_ARRAY(log_sd_log_NAA); // standard deviation log random effect on survival (max_A-1,sp). Needs to be an array or doesn't work!
  PARAMETER_VECTOR(log_sd_log_rec); // standard deviation random effect on log recruitment (sp)
  vector<Type> sd_log_rec=exp(log_sd_log_rec);
  PARAMETER_MATRIX(logit_q); // Logit catchability surveys so q between 0 and 1 (sp,n_surv_max)
  PARAMETER_MATRIX(logit_gamma_surv); // log slope param for s_surv (sp,n_surv_max)
  PARAMETER_MATRIX(logit_A50_surv); // logit age at 50% selectivity for calculation of s_surv (sp,n_surv_max)
  PARAMETER_MATRIX(acomp_pars_temp_catch) // (sp,3) To estimate for age composition distribution
  PARAMETER_ARRAY(acomp_pars_temp_index); // (sp,3) To estimate for age composition distribution
  PARAMETER_VECTOR(mean_log_rec); // Average log recruitment over time in numbers (sp)
  PARAMETER_ARRAY(log_M); // (Y,max_A,sp)
  PARAMETER_ARRAY(log_M1); // (max_A,sp)
  PARAMETER_ARRAY(log_MAA); // (Y-1,max_A,sp)
  PARAMETER_VECTOR(log_lorenzen1); // (sp)
  PARAMETER_VECTOR(lorenzen2); // (sp)
  PARAMETER_VECTOR(log_sd_log_MAA); // log sd random walk M for M_model 2 (sp)
  vector<Type> sd_log_MAA=exp(log_sd_log_MAA);
  PARAMETER_MATRIX(logit_scale_M); // scale on M from assessment (max_A,sp)
  PARAMETER_ARRAY(logit_s_surv); // survey selectivity when sel_model_surv=1 (max_A,sp,n_surv_max)
  PARAMETER_MATRIX(logit_s_F); // fishing selectivity when sel_model_F=1 (max_A,sp)
  // For trophic inetractions
  PARAMETER_MATRIX(vuln_par); // General vulnerability for a pair prey-predator (sp,n_pred)
  PARAMETER_VECTOR(log_shape_gamma_par); // Log of shape parameter for gamma distribution size preference (n_pred)
  PARAMETER_VECTOR(log_scale_gamma_par); // Log of scale parameter for gamma distribution size preference (n_pred)
  PARAMETER_VECTOR(log_power_typeIII); // Exponent of type III functional response (n_pred)
  PARAMETER_ARRAY(par_deltadir); // phi, slope and intercept parameters to estimate in dirichlet (Y,3,n_pred)
  PARAMETER_ARRAY(log_cons_rate); // log of per capita consumption rate (food consumed per ind per year) (max_B,n_pred). Need to be array!!!
  PARAMETER_ARRAY(log_sd_cons_rate); // log of sd of consumption rate distribution (max_B,n_pred). Need to be array!!!
  PARAMETER_ARRAY(log_alpha_cons); // alpha of evacuation rate equation (sp,n_pred). Needs to be array!!!!!
  PARAMETER_ARRAY(log_beta_cons); // beta of evacuation rate equation (sp,n_pred) Needs to be array!!!!!
  array<Type> alpha_cons(sp,n_pred);
  alpha_cons = exp(log_alpha_cons);
  array<Type> beta_cons(sp,n_pred);
  beta_cons = exp(log_beta_cons);
  // PARAMETER_MATRIX(log_cons_rate_other); // log of per capita consumption rate of other food (food consumed per ind per year) (max_B,n_pred)
  // PARAMETER_MATRIX(log_sd_cons_rate_other); // log of sd of consumption rate other food distribution (max_B,n_pred)
  PARAMETER(log_biomass_other_y1); // Initial biomass other food
  PARAMETER(log_growth_other);
  PARAMETER(log_K_other);
  
  // Ar1 on log_E
  PARAMETER_VECTOR(log_sd_process_F); // log sd on 1st pred_log_E when process_F=1 (sp) 
  vector<Type> sd_process_F=exp(log_sd_process_F);
  PARAMETER_VECTOR(logit_phi_process_F); // parameter of AR1 when process_F=1 (sp) 
  vector<Type> phi_process_F(sp);
  for(int i = 0; i < sp; i++) phi_process_F(i) = inverse_logit(logit_phi_process_F(i), Type(-1), Type(1)); // logit transform so -1 < phi_process_F < 1
  PARAMETER_VECTOR(cst_process_F); // cst to calculate mean of AR1 when process_F=1 (sp)
  vector<Type> mean_log_E = cst_process_F/(1-phi_process_F);
    
  // // Test estimate CV for aggregated catch
  // PARAMETER_ARRAY(cv_obs_aggr_Cw); // Coeff deviation for aggregated catch, needs to be an array or log()+1 below doesn't work
  // array<Type> sd_obs_aggr_Cw(Y,sp);
  // sd_obs_aggr_Cw=sqrt(log(square(cv_obs_aggr_Cw)+1));
  // //
  
  // Define dimensions for model objects
  vector<Type> A50_F(sp);
  vector<Type> gamma_F(sp);
  matrix<Type> s_F(max_A,sp); // fishing selectivity
  array<Type> F(Y,max_A,sp);
  array<Type> Z(Y,max_A,sp);
  array<Type> NAA(Y,max_A,sp); // true numbers at age 
  array<Type> pred_NAA(Y,max_A,sp); // matrix to fill for fish numbers
  array<Type> C(Y,max_A,sp);
  array<Type> Cw(Y,max_A,sp);
  array<Type> aggr_Cw(Y,sp); // Needs to be an array even if 2D otherwise it doesn't work since coming from operation on 3D array!
  array<Type> aggr_C(Y,sp); // Needs to be an array even if 2D otherwise it doesn't work since coming from operation on 3D array!
  matrix<Type> q_surv(sp,n_surv_max); // survey catchability
  matrix<Type> A50_surv(sp,n_surv_max);
  matrix<Type> gamma_surv(sp,n_surv_max);
  array<Type> s_surv(max_A,sp,n_surv_max); // Survey selectivity
  array<Type> I(Y,max_A,sp,n_surv_max); // Survey indices of abundance
  array<Type> aggr_I(Y,sp,n_surv_max); // Aggregated indices over ages
  array<Type> prop_Caa(Y,max_A,sp); // Age composition in catch
  array<Type> prop_Iaa(Y,max_A,sp,n_surv_max); // Age composition in surveys
  vector<Type> acomp_pars(3), t_pred_paa(max_A), t_paa(max_A); // 3=max number of param for age comp distribution
  array<Type> sd_log_NAA(max_A-1,sp); // Needs to be an array or doesn't work! No age 1 because no need since already r.e with log_rec
  sd_log_NAA=exp(log_sd_log_NAA);
  Type NLL= 0; // total negative log-likelihood
  array<Type> pred_MAA(Y,max_A,sp);
  array<Type> MAA(Y,max_A,sp);
  matrix<Type> scale_M(max_A,sp);
  matrix<Type> SSB(Y,sp);
  // For trophic interactions
  array<Type> log_ratio_w(Y,max_A,sp,max_B,n_pred); // Log ratio of weight predator against weight prey for size preference
  vector<Type> shape_gamma(n_pred);
  vector<Type> scale_gamma(n_pred);
  vector<Type> mode_gamma(n_pred); // Mode of size pref distribution
  vector<Type> density_mode_gamma(n_pred); // Density of size pref distribution at the mode (max density)
  vector<Type> sum_vuln_par(n_pred); // Sum of vulnerability parameters
  matrix<Type> vuln(sp,n_pred); // General vulnerability for modelled prey
  vector<Type> sum_vuln(n_pred);
  vector<Type> vuln_other(n_pred); // General vulnerability for other food = 1-sum_vuln
  array<Type> size_pref(Y,max_A,sp,max_B,n_pred); // Preference of predator j for prey i
  array<Type> suit(Y,max_A,sp,max_B,n_pred); // Suitability of prey to pred
  array<Type> suit_other(Y,max_B,n_pred); // Suitability of other food to pred
  //array<Type> scaled_suit(Y,max_A,sp,max_B,n_pred); // Scaled suitability of prey to pred relative to all food
  //array<Type> scaled_suit_other(Y,max_B,n_pred); // Scaled suitability of other food relative to all food
  //array<Type> sum_suit(Y,max_B,n_pred);
  array<Type> biomass_prey_avail(Y,max_A,sp,max_B,n_pred); // Biomass of prey available to predator
  array<Type> biomass_prey_avail_no_age(Y,sp,max_B,n_pred); // Biomass of other prey available summed over age of prey
  vector<Type> biomass_other(Y); // Biomass of other food whcih depends on ecosystem biomass
  array<Type> biomass_other_avail(Y,max_B,n_pred); // Biomass of other food available to j
  array<Type> total_biomass_prey_avail(Y,max_B,n_pred); // Total biomass of modelled prey vaialble to j
  array<Type> ratio_biomass_all_prey_avail_no_age(Y,n_prey,max_B,n_pred); // merged ratio for modelled species and other food
  vector<Type> obs_ratio(n_prey); // ratio of observed stomach content per predator, predator age and year
  vector<Type> pred_ratio(n_prey); // ratio of predicted stomach content per predator, predator age and year
  vector<Type> par_ratio(3); // log(phi), log(slope) and intercept parameters for delta dirichlet to estimates for each predator and year
  array<Type> sd_cons_rate(max_B,n_pred); // sd of consumption rate distribution. Need to be array!!!
  array<Type> cons_rate(max_B,n_pred); //per capita consumption rate (food consumed per ind per year). Need to be array!!!
  sd_cons_rate = exp(log_sd_cons_rate);
  cons_rate = exp(log_cons_rate); 
  array<Type> obs_cons_rate(Y_cons,n_pred);
  array<Type> PAA(Y,max_A,sp); // Predation mortality (predation rate)
  array<Type> Cfish(Y,max_A,sp); // Catch of prey species from total predators
  array<Type> Cfishw(Y,max_A,sp); // Same than Cfish but in tonnes
  matrix<Type> N_pred(Y,n_pred); // Total number of predators, needs to be array
  array<Type> prob_b(Y,max_B,n_pred); // Probability of age B
  array<Type> prob_l_b(n_stom_max,Y,max_B,n_pred); // Probability of l and b
  array<Type> prob_l(n_stom_max,Y,n_pred); // Probability of l
  array<Type> prob_b_given_l(n_stom_max,Y,max_B,n_pred); 
  array<Type> pred_ratio_l(n_stom_max,Y,n_prey,max_B,n_pred); // Predictive biomass ratio of i in each stomach of j of age b
  array<Type> pred_ratio_l_no_age(n_stom_max,Y,n_prey,n_pred); // Predictive biomass ratio summed over predator age
  vector<Type> power_typeIII = exp(log_power_typeIII);
  array<Type> mean_ratio(xmax,n_prey,max_B,n_pred);
  array<Type> mean_pred_ratio(xmax,n_prey,max_B,n_pred);
  
  array<Type> total_predation_other(Y,max_B,n_pred);
  // matrix<Type> sd_cons_rate_other(max_B,n_pred); // sd of consumption rate distribution for other food
  // matrix<Type> cons_rate_other(max_B,n_pred); //per capita consumption rate of other food (food consumed per ind per year)
  // for(int j = 0; j < n_pred; j++){
  //   for(int b = 0; b < Bplus(j); b++){
  //     sd_cons_rate_other(b,j) = exp(log_sd_cons_rate_other(b,j));
  //     cons_rate_other(b,j) = exp(log_cons_rate_other(b,j));
  //   }
  // }
  Type growth_other = exp(log_growth_other);
  Type K_other = exp(log_K_other);
  
  matrix<Type> pred_log_E((Y-1),sp); // year component of F when process_F=1
  vector<Type> sd_increment_process_F(sp); // sd for increment noice for F process
  sd_increment_process_F=exp(log_sd_process_F)*sqrt(1-square(phi_process_F)); 
  
  // Objects used to ADREPORT to avoid adreporting large matrices
  matrix<Type> sumFAA(max_A,sp);
  matrix<Type> sumFy(Y,sp);
  matrix<Type> sumMAA(max_A,sp);
  matrix<Type> sumMy(Y,sp);
  matrix<Type> sumPAA(max_A,sp);
  matrix<Type> sumPy(Y,sp);
  matrix<Type> sumZAA(max_A,sp);
  matrix<Type> sumZy(Y,sp);
  matrix<Type> mean_FAA(max_A,sp);
  matrix<Type> mean_Fy(Y,sp);
  matrix<Type> mean_MAA(max_A,sp);
  matrix<Type> mean_My(Y,sp);
  matrix<Type> mean_PAA(max_A,sp);
  matrix<Type> mean_Py(Y,sp);
  matrix<Type> mean_ZAA(max_A,sp);
  matrix<Type> mean_Zy(Y,sp);
  matrix<Type> recruits(Y,sp);
  

  // Model
  
  
  ////////////////////////////////////// 1. Calculate F //////////////////////////////////////////////////////////
  
  //fishing selectivity
  A50_F.setZero();
  gamma_F.setZero();
  s_F.setZero();
  for(int i = 0; i < sp; i++){
    if (sel_model_F(i)==1){ // estimated selectivity at age
      for(int a = (min_A_catch(i)-1); a < max_A_catch(i); a++){
        s_F(a,i)=1/(1+exp(-logit_s_F(a,i)));
      }
    } else {
      if (sel_model_F(i)==2){ //logistic selectivity
        //fishing parameters
        A50_F(i)=max_A_catch(i)/(1+exp(-logit_A50_F(i))); // logit scale for A50 with lower bound =0 and upper band =Aplus
        gamma_F(i)=max_A_catch(i)/(1+exp(-logit_gamma_F(i)));
        //fishing selectivity
        for(int a = (min_A_catch(i)-1); a < max_A_catch(i); a++){
          s_F(a,i)=1/(1+exp(-((a+1)-A50_F(i))/gamma_F(i)));
        }
        for(int a = (min_A_catch(i)-1); a < max_A_catch(i); a++){
          s_F(a,i) = s_F(a,i)/s_F((max_A_catch(i)-1),i); // so selectivity is forced to be 1 for last age
        }
      } else {
        Rf_error("sel_model_F option does not exist");
      }
    }
  }
  
  // log_E as a AR1 random effect
  pred_log_E.setZero();
  Type nll_process_F = 0;
  using namespace density; // if use the option SCALE and AR1
  if (process_F==1){
    for (int i = 0; i < sp; i++){
      nll_process_F -= dnorm(log_E(0,i),mean_log_E(i),sd_process_F(i),1); // nll for pred_log_E value in 1st year (assume mean of 0)
      if (error_simulate==1){
        SIMULATE {
          log_E(0,i) = rnorm(mean_log_E(i),sd_process_F(i));
          REPORT(log_E);
        }
      }
      for (int t = 1; t < (Y); t++){
        pred_log_E(t-1,i) = cst_process_F(i)+phi_process_F(i)*log_E(t-1,i);
        nll_process_F -= dnorm(log_E(t,i),pred_log_E(t-1,i),sd_increment_process_F(i),1);
        if (error_simulate==1){
          SIMULATE {
            log_E(t,i) = rnorm(pred_log_E(t-1,i),sd_increment_process_F(i));
            REPORT(log_E);
          }
        }
      }
    }
    // // or if want to use build-in AR1
    // for (int i = 0; i < sp; i++){
    //     vector<Type> temp_log_E = log_E.col(i).array()-mean_log_E(i);
    //     nll_process_F += SCALE(AR1(phi_process_F(i)),sd_process_F(i))(temp_log_E);
    //     SIMULATE {
    //       SCALE(AR1(phi_process_F(i)),sd_process_F(i)).simulate(temp_log_E);
    //       log_E.col(i) = temp_log_E+mean_log_E(i);
    //       REPORT(log_E);
    //     }
    // }
  }
  
  
  //F
  F.setZero();
  for(int i = 0; i < sp; i++){
    for(int t = 0; t < Y; t++){
      for(int a = (min_A_catch(i)-1); a < max_A_catch(i); a++){
        F(t,a,i)=exp(log_E(t,i))*s_F(a,i);
      }
    }
  }
  
  
  ////////////////////////////////////// 2. Calculate MAA //////////////////////////////////////////////////////////
  
  //M
  pred_MAA.setZero();
  MAA.setZero();
  scale_M.setZero();
  Type nll_log_MAA = 0;
  for(int i = 0; i < sp; i++){
    if (M_model(i)==4){
      for(int a = 0; a < Aplus(i); a++){
        scale_M(a,i)=scale_M_upper(i)/(1+exp(-logit_scale_M(a,i)));
      }
    }
  }
  for(int i = 0; i < sp; i++){
    for(int a = 0; a < Aplus(i); a++){
      if (M_model(i)==2){
        if (process_M==1){
          MAA(0,a,i)=exp(log_M1(a,i));
          pred_MAA(0,a,i)=MAA(0,a,i);
        }
      }
      for(int t = 0; t < Y; t++){
        if (M_model(i)==1){
          MAA(t,a,i)=M(t,a,i); // M is given as input data
        } else {
          if (M_model(i)==2){
            if (process_M==1){
              if (t > 0){
                MAA(t,a,i)=exp(log_MAA(t-1,a,i)); // random walk M
                pred_MAA(t,a,i)=MAA(t-1,a,i);
                //////// PROCESS ERRORS AND PRIORS M ////////
                // NLL for process error M
                nll_log_MAA -= dnorm(log_MAA(t-1,a,i),log(pred_MAA(t,a,i)),sd_log_MAA(i),1); // Random walk M
                // OR nll_log_MAA -= dnorm(log(MAA(t,a,i)),log(pred_MAA(t,a,i)),sd_log_MAA(i),1); // Random walk M
                if (error_simulate==1){
                  SIMULATE {
                    log_MAA(t-1,a,i) = rnorm(log(pred_MAA(t,a,i)), sd_log_MAA(i));
                    REPORT(log_MAA);
                    MAA(t,a,i)=exp(log_MAA(t-1,a,i));
                    REPORT(MAA);
                  }
                }
              }
            } else {
              MAA(t,a,i)=exp(log_M(t,a,i)); // M is an estimated matrix
            }
          } else {
            if (M_model(i)==3){
              MAA(t,a,i)=exp(log_lorenzen1(i)+lorenzen2(i)*log(wc(t,a,i)*1000)); // Lorenzen M with priors on parameters, w must be in grams
            } else {
              if (M_model(i)==4){
                MAA(t,a,i)=scale_M(a,i)*M(t,a,i); // rescaled M from input data
              } else {
                if (M_model(i)!=5){
                  Rf_error("M_model option does not exist");
                }
              }
            }
          }
        }
      }
    }
  }
  
  
  // //////// PROCESS ERRORS AND PRIORS M ////////
  // 
  // // NLL for process error M
  // Type nll_log_MAA = 0;
  // if (process_M==1){
  //   for(int i = 0; i < sp; i++){
  //     if (M_model[i]==2){
  //       for(int a = 0; a < (Aplus[i]); a++){
  //         for(int t = 1; t < (Y); t++){
  //           nll_log_MAA -= dnorm(log_MAA(t-1,a,i),log(pred_MAA(t,a,i)),sd_log_MAA(i),1); // Random walk M
  //           // OR nll_log_MAA -= dnorm(log(MAA(t,a,i)),log(pred_MAA(t,a,i)),sd_log_MAA(i),1); // Random walk M
  //           if (error_simulate==1){
  //             SIMULATE {
  //               log_MAA(t-1,a,i) = rnorm(log(pred_MAA(t,a,i)), sd_log_MAA(i));
  //               REPORT(log_MAA);
  //               MAA(t,a,i)=exp(log_MAA(t-1,a,i));
  //               REPORT(MAA);
  //             }
  //           }
  //         }
  //       }
  //     }
  //   }
  // }
  
  // Priors Lorenzen parameters
  Type nll_log_lorenzen1 = 0;
  Type nll_lorenzen2 = 0;
  Type mean_lorenzen1=3.69;
  Type sd_lorenzen1=0.5;
  Type mean_lorenzen2=-0.305;
  Type sd_lorenzen2=0.028;
  Type log_sd_lorenzen1=sqrt(log(square(sd_lorenzen1/mean_lorenzen1)+1));
  for(int i = 0; i < sp; i++){
    if (M_model(i)==3){
      // nll_log_lorenzen1 -= dnorm(exp(log_lorenzen1(i)),mean_lorenzen1,sd_lorenzen1,1);
      nll_log_lorenzen1 -= dnorm(log_lorenzen1(i),log(mean_lorenzen1),log_sd_lorenzen1,1);
      nll_lorenzen2 -= dnorm(lorenzen2(i),mean_lorenzen2,sd_lorenzen2,1);
      // if (error_simulate==1){
      //   SIMULATE {
      //     log_lorenzen1(i) = rnorm(log(mean_lorenzen1),log_sd_lorenzen1);
      //     //log_lorenzen1(i) = log(rnorm(mean_lorenzen1,sd_lorenzen1));
      //     lorenzen2(i) = rnorm(mean_lorenzen2,sd_lorenzen2);
      //     REPORT(log_lorenzen1);
      //     REPORT(lorenzen2);
      //     for(int t = 0; t < Y; t++){
      //       for (int a = 0; a < Aplus(i); a++){
      //         MAA(t,a,i)=exp(log_lorenzen1(i)+lorenzen2(i)*log(wc(t,a,i)*1000));
      //       }
      //     }
      //     REPORT(MAA);
      //   }
      // }
    }
  }
  /////////////////////////////////////////////
  
  
  ////////////////////////////////////// 3. Calculate Recruitment and NAA in first year //////////////////////////////////////////////////////////
  
  //NAA true number at age
  // pred_NAA, predicted number at age for random effects
  NAA.setZero();
  pred_NAA.setZero();
  Type nll_log_rec = 0;
  // fill N at age in first year Y1
  for(int i = 0; i < sp; i++){
    for(int a = 0; a < Aplus(i); a++){
      NAA(0,a,i)=exp(log_N1(a,i)); //Give N at age in first year (numbers)
      pred_NAA(0,a,i)=NAA(0,a,i); //fill ages 1st year for predicted N
    }
    // fill recruitment column for NAA
    for(int t = 1; t < Y; t++){
      NAA(t,0,i)=exp(log_rec(t-1,i)); //Give recruitment at age 1 every year (numbers) with possibility of random effect
      // fill recruitment column for pred_NAA
      if (process_rec==1){ // random recruitment
        if (recruit_model==1){
          pred_NAA(t,0,i)=NAA(t-1,0,i); //random walk recruitment
        } else {
          if (recruit_model==2){ 
            pred_NAA(t,0,i)=exp(mean_log_rec(i)); //random about the mean
          } else {
            Rf_error("recruit_model does not exist");
          }
        }
        //////// PROCESS ERRORS RECRUITMENT ////////
        // NLL for process error recruitment
        nll_log_rec -= dnorm(log(NAA(t,0,i)),log(pred_NAA(t,0,i)),sd_log_rec(i),1); // Random recruitment
        // nll_log_rec -= dnorm(log_rec(t-1,i),log(pred_NAA(t,0,i)),sd_log_rec(i),1); // Random walk recruitment
        if (error_simulate==1){
          SIMULATE {
            log_rec(t-1,i) = rnorm(log(pred_NAA(t,0,i)), sd_log_rec(i));
            REPORT(log_rec);
            NAA(t,0,i)=exp(log_rec(t-1,i));
            REPORT(NAA);
          }
        }
      }
      if (process_rec==0){ // deterministic recruitment
        pred_NAA(t,0,i)=NAA(t,0,i); //same than pred_NAA(t,0,i)=exp(log_rec(t-1,i)) and log_rec is estimated
      } 
    }
  }

  //////// PROCESS ERRORS RECRUITEMENT ////////
  
  // // NLL for process error recruitment
  // //parallel_accumulator<Type> nll_log_rec(this);
  // Type nll_log_rec = 0;
  // if (process_rec==1){
  //   for(int i = 0; i < sp; i++){
  //     for(int t = 1; t < (Y); t++){ //Recruitment start t=2 so max t = Y-1
  //       nll_log_rec -= dnorm(log(NAA(t,0,i)),log(pred_NAA(t,0,i)),sd_log_rec(i),1); // Random recruitment
  //       // nll_log_rec -= dnorm(log_rec(t-1,i),log(pred_NAA(t,0,i)),sd_log_rec(i),1); // Random walk recruitment
  //       if (error_simulate==1){
  //         SIMULATE {
  //           log_rec(t-1,i) = rnorm(log(pred_NAA(t,0,i)), sd_log_rec(i));
  //           REPORT(log_rec);
  //           NAA(t,0,i)=exp(log_rec(t-1,i));
  //           REPORT(NAA);
  //         }
  //       }
  //     }
  //   }
  // }
  
  
  
  ///////////////////////////////// 4. Calculate predation variables if predation is on //////////////////////////////////////////////////////////
  
  
  Type nll_gamma_pref = 0;
  Type nll_cons_rate = 0;
  Type nll_alpha_cons = 0;
  Type nll_beta_cons = 0;
  
  if (predation_on==1){ // Only if predation is on

    /////////// Size preference ///////////////////////
    
    // Gamma parameters for calculation of size pref distribution given or estimated
    // If given shape and scale empirically size_pref calculated given the observed ratio of weight predator over prey
    // If estimated, same thing but shape and scale estimated within the model to get s.e.
    shape_gamma.setZero();
    scale_gamma.setZero();
    mode_gamma.setZero();
    density_mode_gamma.setZero();
    for(int j = 0; j < n_pred; j++){
      shape_gamma(j)=exp(log_shape_gamma_dat(j)); // >0
      scale_gamma(j)=exp(log_scale_gamma_dat(j)); // >0
      mode_gamma(j)=(shape_gamma(j)-1)*scale_gamma(j);
      density_mode_gamma(j) = dgamma(mode_gamma(j),shape_gamma(j),scale_gamma(j));
      if (gamma_pref_estim==1){ // estimated given input log_ratio_w_diet
        // NLL gamma distribution size preference when parameters estimated
        for (int k = 0; k < n_ratio; k++){
          nll_gamma_pref -= dgamma(log_ratio_w_diet(k,j),exp(log_shape_gamma_par(j)),exp(log_scale_gamma_par(j)),1);
          if (data_simulate==1){
            SIMULATE {
              log_ratio_w_diet(k,j) = rgamma(exp(log_shape_gamma_par(j)),exp(log_scale_gamma_par(j)));
              REPORT(log_ratio_w_diet);
            }
          }
        }
      }
    }


    // Size preference of predator j for prey i calculated given the log ratio of weights in the model but shape and scale previously estimated from food habits data
    log_ratio_w.setZero();
    size_pref.setZero();
    for(int j = 0; j < n_pred; j++){
      for(int i = 0; i < sp; i++){
        for(int b = 0; b < Bplus(j); b++){
          for(int a = 0; a < Aplus(i); a++){
            for(int t = 0; t < Y; t++){
              log_ratio_w(t,a,i,b,j) = log(wc(t,b,j)/wc(t,a,i));
              if (log_ratio_w(t,a,i,b,j)>0){
                size_pref(t,a,i,b,j) = dgamma(log_ratio_w(t,a,i,b,j),shape_gamma(j),scale_gamma(j))/density_mode_gamma(j); // divided by max so between 0 and 1
              } else {
                size_pref(t,a,i,b,j) = 0;
              }
            }
          }
        }
      }
    }
    
    /////////// Prey suitability ///////////////////////
    
    // Vulnerability of prey i for predator j
    // vuln estimated for sp species and other food vuln = 1-sum(vuln)
    // Criteria for vuln: 
    // 1. vuln > 0 (hence logs)
    // 2. sum_vuln + vuln_other = 1
    // 3. 0 <= sum_vuln over modelled prey <= 1 (hence logit transformation for sum_vuln)
    
    sum_vuln_par.setZero();
    vuln.setZero();
    vuln_other.setZero();
    suit.setZero();
    sum_vuln.setZero();
    suit_other.setZero();
    for(int j = 0; j < n_pred; j++){
      for(int i = 0; i < sp; i++){
        if (interaction(i,j)==1) sum_vuln_par(j) += exp(vuln_par(i,j));
      }
      for(int i = 0; i < sp; i++){
        if (interaction(i,j)==1){
          vuln(i,j) = exp(vuln_par(i,j))/(1+sum_vuln_par(j)); // multinomial logistic transformation
          sum_vuln(j) += vuln(i,j); // sum vuln over modelled prey
        }
      }
      vuln_other(j) = 1-sum_vuln(j); // vuln-other=1-sum-vuln but transform so sum_vuln+vuln_other=1
    }
    // Suitability of prey i for predator j
    for(int j = 0; j < n_pred; j++){
      for(int b = 0; b < Bplus(j); b++){
        for(int t = 0; t < Y; t++){
          suit_other(t,b,j) = vuln_other(j); // size_pref_other=1 so suit_other=vuln_other*1 (may be change in future if calculate a size_pref_other)
          for(int i = 0; i < sp; i++){
            for(int a = 0; a < Aplus(i); a++){
              suit(t,a,i,b,j) = vuln(i,j)*size_pref(t,a,i,b,j);
            }
          }
        }
      }
    }

    // // Scaled suitability relative to all prey so sum to 1
    // sum_suit.setZero();
    // scaled_suit.setZero();
    // scaled_suit_other.setZero();
    // for(int t = 0; t < Y; t++){
    //   for(int j = 0; j < n_pred; j++){
    //     for(int b = 0; b < Bplus(j); b++){
    //       for(int i = 0; i < sp; i++){
    //         for(int a = 0; a < Aplus(i); a++){
    //           sum_suit(t,b,j) += suit(t,a,i,b,j);
    //         }
    //       }
    //       for(int i = 0; i < sp; i++){
    //         for(int a = 0; a < Aplus(i); a++){
    //           scaled_suit(t,a,i,b,j) = suit(t,a,i,b,j)/(sum_suit(t,b,j)+suit_other(t,b,j)); // scaled suit for prey i
    //         }
    //       }
    //       scaled_suit_other(t,b,j) = suit_other(t,b,j)/(sum_suit(t,b,j)+suit_other(t,b,j)); // scaled suit for other food
    //     }
    //   }
    // }
    
    
    /////////// Consumption rate of predator ///////////////////////
    
    // Estimate consumption rate parameters within the model with a penalty function, otherwise alpha and beta = mean in NLL
    obs_cons_rate.setZero();
    
    if (cons_rate_estim==1){ 
      // NLL penalties on alpha and beta to help estimation given meta-analysis on 13 papers for 18 values of alpha and beta
      for (int j = 0; j < n_pred; j++){
        for(int i = 0; i < sp; i++){
          nll_alpha_cons -= dnorm(alpha_cons(i,j),Type(0.02835569),Type(0.01642670),1);
          nll_beta_cons -= dnorm(beta_cons(i,j),Type(0.11982778),Type(0.02921756),1);
        }
      }
    //} 
    
    // Estimate annual consumption rate for each predator in kg of prey eaten per predator given food habit data (temperatures and mean w in stomach)
    
    for (int j = 0; j < n_pred; j++){
        for (int t = 0; t < Y_cons; t++){ 
          if (flag_cons_rate(t)!=0)
            //consumption rate spring+fall in kg/predator/year
            obs_cons_rate(t,j) = 182.5/1000*24*alpha_cons(j)*exp(beta_cons(j)*spring_cons(t,1,j))*spring_cons(t,0,j) + 182.5/1000*24*alpha_cons(j)*exp(beta_cons(j)*fall_cons(t,1,j))*fall_cons(t,0,j);
        }
    }
    
    // NLL per-capita consumption rates in log(kg consumed/fish/year)
    // The observed data is the one calculated from spring and fall food habit data
    /////////// Need to do the same for cons_rate_other if interested in estimating it
    for (int j = 0; j < n_pred; j++){
      for (int b = 0; b < max_B; b++){
        for(int i = 0; i < sp; i++){
          for (int t = 0; t < Y_cons; t++){
            if (obs_cons_rate(t,j)!=0)
              nll_cons_rate -= dnorm(log(obs_cons_rate(t,j)),log_cons_rate(b,j),sd_cons_rate(b,j),1);
            if (data_simulate==1){
              SIMULATE {
                obs_cons_rate(t,j) = exp(rnorm(log_cons_rate(b,j),sd_cons_rate(b,j))); 
                REPORT(obs_cons_rate);
              }
            }
          }
        }
      }
    }
    
    }
    
    
  } // end of predation model for now

  
  
  ///////////////////////////////// 5. Calculate PAA and fill pred_NAA with Z=F+MAA+PAA //////////////////////////////////////////////////////////
  
  ///////////////// NAA needed from now on ////////////////
    
    
    biomass_other.setZero();
    biomass_prey_avail.setZero();
    biomass_prey_avail_no_age.setZero();
    biomass_other_avail.setZero();
    total_biomass_prey_avail.setZero();
    ratio_biomass_all_prey_avail_no_age.setZero();
    total_predation_other.setZero();
    PAA.setZero();
    Z.setZero();
    Type nll_log_NAA = 0;
    
    if (predation_on==1){
      if (biomass_other_option==2){ // surplus production model
        biomass_other(0) = exp(log_biomass_other_y1); // Fill 1st year biomass other food
      }
    }
    for(int t = 0; t < Y; t++){
      if (predation_on==1){ // Only if predation is on otherwise PAA=0
        for(int j = 0; j < n_pred; j++){
          for(int b = 0; b < Bplus(j); b++){
            for(int i = 0; i < sp; i++){
              for(int a = 0; a < (Aplus(i)); a++){
                // Biomass of modelled prey i available to predator j in tons
                biomass_prey_avail(t,a,i,b,j) = pow((NAA(t,a,i)*wc(t,a,i)/1000),power_typeIII(j))*suit(t,a,i,b,j);// Type II if power_typeIII not estimated (NAs in map argument) otherwise type III
                biomass_prey_avail_no_age(t,i,b,j) += biomass_prey_avail(t,a,i,b,j); // Biomass of available prey summed over prey ages
              }
              total_biomass_prey_avail(t,b,j) += biomass_prey_avail_no_age(t,i,b,j); // Total available biomass of modelled prey
            }
            // Biomass of other food available to j in tons
            if (biomass_other_option==1){ // biomass other food is cst
              biomass_other_avail(t,b,j) = pow(B_other,power_typeIII(j))*suit_other(t,b,j);// Type II if power_typeIII not estimated (NAs in map argument) otherwise type III, assume suit_other=1
            } else {
              if (biomass_other_option==2){
                total_predation_other(t) += cons_rate(b,j)*NAA(t,b,j)*(biomass_other_avail(t,b,j) / (total_biomass_prey_avail(t,b,j)+biomass_other_avail(t,b,j))); // kg/y
                Type eps = Type(1.0);
                Type pen = Type(0.0);
                if (t<(Y-1))
                  biomass_other(t+1) = posfun(biomass_other(t)+growth_other*biomass_other(t)*(1-biomass_other(t)/K_other)-(total_predation_other(t)/1000),eps,pen);
                biomass_other_avail(t,b,j) = biomass_other(t);
                NLL += 1.0E+20*pen;
              } else {
                Rf_error("biomass_other_option does not exist");
              }
            }
            // Ratio for nll_ratio_diet for modelled prey
            for(int i = 0; i < sp; i++){
              ratio_biomass_all_prey_avail_no_age(t,i,b,j) = biomass_prey_avail_no_age(t,i,b,j)/(total_biomass_prey_avail(t,b,j)+biomass_other_avail(t,b,j)); // size ratio_biomass_all_prey_avail_no_age(Y,n_prey,max_B,n_pred)
            }
            // Fill last column by other food ratio
            ratio_biomass_all_prey_avail_no_age(t,(n_prey-1),b,j) = biomass_other_avail(t,b,j)/(total_biomass_prey_avail(t,b,j)+biomass_other_avail(t,b,j)); 
            for(int i = 0; i < sp; i++){
              for(int a = 0; a < (Aplus(i)); a++){
                // Predation mortality (PAA) on modelled species in year^-1
                PAA(t,a,i) += (NAA(t,b,j)*cons_rate(b,j)/1000*suit(t,a,i,b,j))/(total_biomass_prey_avail(t,b,j)+biomass_other_avail(t,b,j));
              }
            }
          }
        }
      }
      for(int i = 0; i < sp; i++){
        for(int a = 0; a < (Aplus(i)-1); a++){
          if (M_model(i)==5){ // MAA+PAA=input M
            MAA(t,a,i) = M(t,a,i)-PAA(t,a,i);
          }
          // Z
          Z(t,a,i)=F(t,a,i)+MAA(t,a,i)+PAA(t,a,i);
          // Fill pred_NAA now that we have Z
          if (t < (Y-1)){ // because t+1
            pred_NAA(t+1,a+1,i)=NAA(t,a,i)*exp(-Z(t,a,i)); //in numbers
            if (process_survival==1){
              if (error_simulate==1){
                SIMULATE {
                  log_NAA(t,a,i) = rnorm(log(pred_NAA(t+1,a+1,i)), sd_log_NAA(a,i));
                  REPORT(log_NAA);
                }
              }
              NAA(t+1,a+1,i)=exp(log_NAA(t,a,i)); // Random effect for N at age
            } else {
              NAA(t+1,a+1,i)=pred_NAA(t+1,a+1,i); // deterministic survival
            }
          }
        }
        //// Age plus group
        if (M_model(i)==5) MAA(t,Aplus(i)-1,i) = M(t,Aplus(i)-1,i)-PAA(t,Aplus(i)-1,i);
        Z(t,Aplus(i)-1,i)=F(t,Aplus(i)-1,i)+MAA(t,Aplus(i)-1,i)+PAA(t,Aplus(i)-1,i);
        if (t < (Y-1)){ // because t+1
          pred_NAA(t+1,Aplus(i)-1,i)+=NAA(t,Aplus(i)-1,i)*exp(-Z(t,Aplus(i)-1,i)); // age plus group
          if (process_survival==1){
            if (error_simulate==1){
              SIMULATE {
                log_NAA(t,Aplus(i)-2,i) = rnorm(log(pred_NAA(t+1,Aplus(i)-1,i)), sd_log_NAA(Aplus(i)-2,i));
                REPORT(log_NAA);
              }
            }
            // NLL for process error fish survival
            NAA(t+1,Aplus(i)-1,i)=exp(log_NAA(t,Aplus(i)-2,i)); // Random effect for N at age
          } else {
            NAA(t+1,Aplus(i)-1,i)=pred_NAA(t+1,Aplus(i)-1,i); // deterministic survival
          }
        }
        //////// PROCESS ERRORS SURVIVAL ////////
        // NLL for process error fish survival
        if (process_survival==1){
          for(int a = 1; a < Aplus(i); a++){
            if (t > 0){ 
              nll_log_NAA -= dnorm(log_NAA(t-1,a-1,i),log(pred_NAA(t,a,i)),sd_log_NAA(a-1,i),1); // Random walk on fish survival
              // OR nll_log_NAA -= dnorm(log(NAA(t,a,i)),log(pred_NAA(t,a,i)),sd_log_NAA(a-1,i),1); // Random walk on fish survival
            }
          }
        }
      }
    }
    


  // //////// PROCESS ERRORS SURVIVAL ////////
  // 
  // // NLL for process error fish survival
  // //parallel_accumulator<Type> nll_log_NAA(this);
  // Type nll_log_NAA = 0;
  // if (process_survival==1){
  //   for(int i = 0; i < sp; i++){
  //     for(int a = 1; a < (Aplus[i]); a++){
  //       for(int t = 1; t < (Y); t++){
  //         nll_log_NAA -= dnorm(log_NAA(t-1,a-1,i),log(pred_NAA(t,a,i)),sd_log_NAA(a-1,i),1); // Random walk on fish survival
  //         // OR nll_log_NAA -= dnorm(log(NAA(t,a,i)),log(pred_NAA(t,a,i)),sd_log_NAA(a-1,i),1); // Random walk on fish survival
  //         if (error_simulate==1){
  //           SIMULATE {
  //             log_NAA(t-1,a-1,i) = rnorm(log(pred_NAA(t,a,i)), sd_log_NAA(a-1,i));
  //             REPORT(log_NAA);
  //             NAA(t,a,i)=exp(log_NAA(t-1,a-1,i));
  //             REPORT(NAA);
  //           }
  //         }
  //       }
  //     }
  //   }
  // }


  //////// NLLs for interaction model ////////
  
  Type nll_ratio_diet = 0;
  
  N_pred.setZero();
  prob_b.setZero();
  prob_l_b.setZero();
  prob_l.setZero();
  prob_b_given_l.setZero();
  pred_ratio_l.setZero();
  pred_ratio_l_no_age.setZero();
  mean_ratio.setZero();
  mean_pred_ratio.setZero();
  obs_ratio.setZero();
  pred_ratio.setZero();
  
  if (predation_on==1){ // Only if predation is on


    // Prepare predictive ratio for NLL
    if (flag_nll_diet==1){

      for(int j = 0; j < n_pred; j++){
        for(int t = 0; t < Y; t++){
          for(int b = 0; b < Bplus(j); b++){
            N_pred(t,j) += NAA(t,b,j); // total number of predators
          }
          for(int b = 0; b < Bplus(j); b++){
            prob_b(t,b,j) = NAA(t,b,j)/N_pred(t,j); // prob of j to be of age b in year t (depends on abundance in system)
          }
        }
      }

    }

    for(int j = 0; j < n_pred; j++){
      if ( (flag_nll_diet==1) || (flag_nll_diet==2) || (flag_nll_diet==3) ){
      for(int t = 0; t < Y; t++){
        // For NLL ratio diet
        for (int p = 0; p < 3; p++){ // log(phi)=precision, log(slope) and intercept for delta dirichlet
          par_ratio(p) = par_deltadir(t,p,j);
        }

        if (flag_nll_diet!=3){ // flag_nll_diet = 1 or 2

        for (int l = 0; l < n_stom(t,j); l++){

          if (flag_nll_diet==1){

          for(int b = 0; b < Bplus(j); b++){
            prob_l_b(l,t,b,j) = prob_l_given_b(l,t,b,j)*prob_b(t,b,j); // joint probability of b and l
            prob_l(l,t,j) += prob_l_b(l,t,b,j); // prob of length l
          }
          for(int i = 0; i < n_prey; i++){ // n_prey because modelled prey + other food
            for(int b = 0; b < Bplus(j); b++){
              prob_b_given_l(l,t,b,j) = prob_l_b(l,t,b,j)/prob_l(l,t,j); // prob of age b given length l
              pred_ratio_l(l,t,i,b,j) = ratio_biomass_all_prey_avail_no_age(t,i,b,j)*prob_b_given_l(l,t,b,j);
              pred_ratio_l_no_age(l,t,i,j) += pred_ratio_l(l,t,i,b,j); // predictive biomass ratio for a given predator length
            }
          }
        //}

          } else { // flag_nll_diet=2, new simplified model per stomach

            for(int i = 0; i < n_prey; i++){ // n_prey because modelled prey + other food
              pred_ratio_l_no_age(l,t,i,j) = ratio_biomass_all_prey_avail_no_age(t,i,(age_pred(l,t,j)-1),j);
            }

          }
        // for (int l = 0; l < n_stom(t,j); l++){
          int n_size = 0;
          for(int i = 0; i < n_prey; i++){ // n_prey because modelled prey + other food
            obs_ratio(i) = ratio_diet(l,t,i,j);
            pred_ratio(i) = pred_ratio_l_no_age(l,t,i,j);
            if (pred_ratio(i)!=0 && obs_ratio(i)!=0){
              n_size++;
            }
          }
          //see(n_size);
          vector<Type> pred_ratio2(n_size);
          vector<Type> obs_ratio2(n_size);
          int n = 0;
          for(int i = 0; i < n_prey; i++){
            if (pred_ratio(i)!=0 && obs_ratio(i)!=0){
              pred_ratio2(n) = pred_ratio(i);
              obs_ratio2(n) = obs_ratio(i);
              n++;
            }
          }
          pred_ratio2 = pred_ratio2/pred_ratio2.sum();
          obs_ratio2 = obs_ratio2/obs_ratio2.sum();
          //see(pred_ratio2.sum());
          //see(pred_ratio2);
          //if (t==0 && j==0 && l==4) see(obs_ratio2);
           if (pred_ratio2.size()!=0 && obs_ratio2.sum()!=0){
             if (diet_model==1) nll_ratio_diet -= ddeltadir(obs_ratio2,pred_ratio2,par_ratio,1);
             if (diet_model==2) nll_ratio_diet -= ddirichlet(obs_ratio2,pred_ratio2,exp(par_ratio(0)),1);
           }
           //if (l==1) see(nll_ratio_diet);
           if (data_simulate==1){
             SIMULATE {
               //if (pred_ratio2.size()!=0){
               //vector<Type> ratio(pred_ratio2.size());
               obs_ratio2.setZero();
               if (diet_model==1) obs_ratio2 = rdeltadir(pred_ratio2,par_ratio);
               if (diet_model==2) obs_ratio2 = rdirichlet(pred_ratio2,exp(par_ratio(0)));
               int k = 0;
               for(int i = 0; i < n_prey; i++){
                 if (pred_ratio(i)!=0 && obs_ratio(i)!=0){
                   ratio_diet(l,t,i,j) = obs_ratio2(k);
                   k++;
                  } else {
                    ratio_diet(l,t,i,j) = 0;
                  }
                }
               //}
                REPORT(ratio_diet);
              }
            }
        }

        } else { // flag_nll_diet=3, new highly simplified model where stomachs aggregated

          for(int b = 0; b < Bplus(j); b++){
            int n_size = 0;
            for(int i = 0; i < n_prey; i++){ // n_prey because modelled prey + other food
              obs_ratio(i) = ratio_diet3(t,i,b,j);
              pred_ratio(i) = ratio_biomass_all_prey_avail_no_age(t,i,b,j);
              if (pred_ratio(i)!=0 && obs_ratio(i)!=0){
                n_size++;
              }
            }
            //see(n_size);
            vector<Type> pred_ratio2(n_size);
            vector<Type> obs_ratio2(n_size);
            int n = 0;
            for(int i = 0; i < n_prey; i++){
              if (pred_ratio(i)!=0 && obs_ratio(i)!=0){
                pred_ratio2(n) = pred_ratio(i);
                obs_ratio2(n) = obs_ratio(i);
                n++;
              }
            }
            pred_ratio2 = pred_ratio2/pred_ratio2.sum();
            obs_ratio2 = obs_ratio2/obs_ratio2.sum();
            // see(pred_ratio2);
            // see(obs_ratio2);
            if (pred_ratio2.size()!=0 && obs_ratio2.sum()!=0){
              if (diet_model==1) nll_ratio_diet -= ddeltadir(obs_ratio2,pred_ratio2,par_ratio,1);
              if (diet_model==2) nll_ratio_diet -= ddirichlet(obs_ratio2,pred_ratio2,exp(par_ratio(0)),1);
            }
            //see(nll_ratio_diet);
            //if (l==1) see(nll_ratio_diet);
            if (data_simulate==1){
              SIMULATE {
                //if (pred_ratio2.size()!=0){
                //vector<Type> ratio(pred_ratio2.size());
                obs_ratio2.setZero();
                if (diet_model==1) obs_ratio2 = rdeltadir(pred_ratio2,par_ratio);
                if (diet_model==2) obs_ratio2 = rdirichlet(pred_ratio2,exp(par_ratio(0)));
                int k = 0;
                for(int i = 0; i < n_prey; i++){
                  if (pred_ratio(i)!=0 && obs_ratio(i)!=0){
                    ratio_diet3(t,i,b,j) = obs_ratio2(k);
                    k++;
                  } else {
                    ratio_diet3(t,i,b,j) = 0;
                  }
                }
                //}
                REPORT(ratio_diet3);
              }
            }
          }

        }
      } // close t
      } else { // flag_nll_diet != 1 or 2 or 3 so flag_nll_diet = 4 or 5

        if (flag_nll_diet==4){ // xmax years average

          for(int b = 0; b < Bplus(j); b++){
            for(int i = 0; i < n_prey; i++){ // n_prey because modelled prey + other food
              for (int x = 0; x < xmax; x++){
                for (int t = t1(x); t < t2(x); t++){
                  mean_ratio(x,i,b,j) += ratio_diet3(t,i,b,j)/length_t(x); // diet prop averaged over length_t years
                  mean_pred_ratio(x,i,b,j) += ratio_biomass_all_prey_avail_no_age(t,i,b,j)/length_t(x);
                }
              }
            }
          }
          //see(mean_ratio);

          for (int x = 0; x < xmax; x++){

            par_ratio.setZero();
            for (int p = 0; p < 3; p++){ // log(phi)=precision, log(slope) and intercept for delta dirichlet
              for (int t = t1(x); t < t2(x); t++){
                par_ratio(p) += par_deltadir(t,p,j)/length_t(x);
              }
            }
            //see(par_ratio);

            for(int b = 0; b < Bplus(j); b++){
              int n_size = 0;
              for(int i = 0; i < n_prey; i++){ // n_prey because modelled prey + other food
                obs_ratio(i) = mean_ratio(x,i,b,j);
                pred_ratio(i) = mean_pred_ratio(x,i,b,j);
                if (pred_ratio(i)!=0 && obs_ratio(i)!=0){
                  n_size++;
                }
              }
              //see(obs_ratio);
              vector<Type> pred_ratio2(n_size);
              vector<Type> obs_ratio2(n_size);
              int n = 0;
              for(int i = 0; i < n_prey; i++){
                if (pred_ratio(i)!=0 && obs_ratio(i)!=0){
                  pred_ratio2(n) = pred_ratio(i);
                  obs_ratio2(n) = obs_ratio(i);
                  n++;
                }
              }
              pred_ratio2 = pred_ratio2/pred_ratio2.sum();
              obs_ratio2 = obs_ratio2/obs_ratio2.sum();
              //see(pred_ratio2);
              if (pred_ratio2.size()!=0 && obs_ratio2.sum()!=0){
                if (diet_model==1) nll_ratio_diet -= ddeltadir(obs_ratio2,pred_ratio2,par_ratio,1);
                if (diet_model==2) nll_ratio_diet -= ddirichlet(obs_ratio2,pred_ratio2,exp(par_ratio(0)),1);
              }
              //see(nll_ratio_diet);
            }
          }
        } else {
          if (flag_nll_diet==5){ // flag_nll_diet==5, ratio averaged over whole time series
            for (int p = 0; p < 3; p++){ // log(phi)=precision, log(slope) and intercept for delta dirichlet
              for (int t = 0; t < Y; t++){
                par_ratio(p) += par_deltadir(t,p,j)/Y;
              }
            }
            //see(par_ratio);
  
            for(int b = 0; b < Bplus(j); b++){
              int n_size = 0;
              obs_ratio.setZero();
              pred_ratio.setZero();
              for(int i = 0; i < n_prey; i++){ // n_prey because modelled prey + other food
                for (int t = 0; t < Y; t++){
                  obs_ratio(i) += ratio_diet3(t,i,b,j)/Y;
                  pred_ratio(i) += ratio_biomass_all_prey_avail_no_age(t,i,b,j)/Y;
                }
                if (pred_ratio(i)!=0 && obs_ratio(i)!=0){
                  n_size++;
                }
              }
              //see(obs_ratio);
              //see(n_size);
              vector<Type> pred_ratio2(n_size);
              vector<Type> obs_ratio2(n_size);
              int n = 0;
              for(int i = 0; i < n_prey; i++){
                if (pred_ratio(i)!=0 && obs_ratio(i)!=0){
                  pred_ratio2(n) = pred_ratio(i);
                  obs_ratio2(n) = obs_ratio(i);
                  n++;
                }
              }
              pred_ratio2 = pred_ratio2/pred_ratio2.sum();
              obs_ratio2 = obs_ratio2/obs_ratio2.sum();
              //see(pred_ratio2);
              //see(obs_ratio2);
              //see(pred_ratio2.size());
              if (pred_ratio2.size()!=0 && obs_ratio2.sum()!=0){
                if (diet_model==1) nll_ratio_diet -= ddeltadir(obs_ratio2,pred_ratio2,par_ratio,1);
                if (diet_model==2) nll_ratio_diet -= ddirichlet(obs_ratio2,pred_ratio2,exp(par_ratio(0)),1);
              }
              //see(nll_ratio_diet);
            }
          } else {
            Rf_error("flag_nll_diet option does not exist");
          }
        }
      }
    } // close j



  } // final end of the interaction model

 //////////////////////////////////////////////////////// 
 
 
 
  
  /////////////////////// 6. Calculate SSB //////////////////////////////////
  
  SSB.setZero(); //Really important or SSB not well calculated
  for(int i = 0; i < sp; i++){
    for(int t = 0; t < Y; t++){
      for(int a = 0; a < Aplus(i); a++){
        SSB(t,i) += NAA(t,a,i)*(wSSB(t,a,i)/1000)*mature(t,a,i)*exp(-Z(t,a,i)*prop_y_elapsed_SSB(i)); //in tons
      }
    }
  }
  
  
  
  /////////////////////// 7. Calculate Catch //////////////////////////////////
  
  C.setZero();
  Cw.setZero();
  Cfish.setZero();
  Cfishw.setZero();
  aggr_C.setZero();
  aggr_Cw.setZero();
  Type nll_aggr_Cw = 0;
  prop_Caa.setZero();
  Type nll_prop_Caa = 0;
  acomp_pars.setZero();
  t_pred_paa.setZero();
  t_paa.setZero();
  //C.fill() = 0.0
  
  for(int i = 0; i < sp; i++){
    // For nll_prop_Caa
    for (int j = 0; j < 3; j++){ // 3=max number of param for age comp distribution
      acomp_pars(j) = acomp_pars_temp_catch(i,j); 
    }
    for(int t = 0; t < Y; t++){
      for(int a = (min_A_catch(i)-1); a < max_A_catch(i); a++){
        
        // Catches in numbers and weight
        C(t,a,i)=F(t,a,i)/Z(t,a,i)*NAA(t,a,i)*(1-exp(-Z(t,a,i))); // in numbers
        Cw(t,a,i)=C(t,a,i)*wc(t,a,i)/1000; //in tons
        if (predation_on==1){
          Cfish(t,a,i)=PAA(t,a,i)/Z(t,a,i)*NAA(t,a,i)*(1-exp(-Z(t,a,i))); // in numbers
          Cfishw(t,a,i)=Cfish(t,a,i)*wc(t,a,i)/1000; //in tons
        }
        
        // Aggregated catch        
        aggr_C(t,i) += C(t,a,i); // in numbers
        aggr_Cw(t,i) += Cw(t,a,i); // in tons
      }
      
      // NLL total aggregated catch
      if (flag_aggr_Cw(t,i)==1){
        nll_aggr_Cw -= dnorm(log(obs_aggr_Cw(t,i)),log(aggr_Cw(t,i)),sd_obs_aggr_Cw(t,i),1);
        if (data_simulate==1){
          SIMULATE {
            obs_aggr_Cw(t,i) = exp(rnorm(log(aggr_Cw(t,i)), sd_obs_aggr_Cw(t,i)));
            REPORT(obs_aggr_Cw);
          }
        }
      }
      
      // Age composition catch 
      for(int a = (min_A_catch(i)-1); a < max_A_catch(i); a++){
        prop_Caa(t,a,i) = C(t,a,i)/aggr_C(t,i);
      }
      
      //NLL age composition catch
      if (flag_Caa(t,i)==1){
        for(int a = (min_A_catch(i)-1); a < max_A_catch(i); a++){
          t_pred_paa(a) = prop_Caa(t,a,i);
          t_paa(a) = obs_prop_Caa(t,a,i);
        } 
        nll_prop_Caa -= get_acomp_ll(Y, min_A_catch(i)-1, max_A_catch(i), Neff_C(t,i), age_comp_model_catch(i), t_paa, t_pred_paa, acomp_pars, catch_aref(t,i));
        //std::cout << "nll_prop_Caa=" << nll_prop_Caa << std::endl;
        if (data_simulate==1){
          SIMULATE {
            vector<Type> paa = sim_acomp(Y, min_A_catch(i)-1, max_A_catch(i), Neff_C(t,i), age_comp_model_catch(i), t_paa, t_pred_paa, acomp_pars, catch_aref(t,i));
            for(int a = (min_A_catch(i)-1); a < Aplus(i); a++){
              obs_prop_Caa(t,a,i) = paa(a);
            }
            if (max_A_catch(i) < max_A){
              for(int a = max_A_catch(i); a < max_A; a++){
                obs_prop_Caa(t,a,i) = 0;
              }
            }
            REPORT(obs_prop_Caa);
          }
        }
      }
    }
  }
  


  /////////////////////// 8. Calculate survey indices //////////////////////////////////
  

  q_surv.setZero();
  A50_surv.setZero();
  gamma_surv.setZero();
  s_surv.setZero();
  I.setZero();
  aggr_I.setZero();
  prop_Iaa.setZero();
  acomp_pars.setZero();
  t_pred_paa.setZero();
  t_paa.setZero();
  Type nll_aggr_I = 0;
  Type nll_prop_Iaa = 0;

  for(int i = 0; i < sp; i++){
    for(int k = 0; k < n_surv(i); k++){
      
      // For nll_prop_Iaa
      for (int j = 0; j < 3; j++){ // 3=max number of param for age comp distribution
        acomp_pars(j) = acomp_pars_temp_index(i,j,k);
      }
      // Survey catchability
      q_surv(i,k)=1/(1+exp(-logit_q(i,k))); // q between 0 and 1
      
      // Survey selectivity
      if (sel_model_surv(i,k)==1){ // estimated selectivity at age
        for(int a = (min_A_surv(i,k)-1); a < max_A_surv(i,k); a++){
          s_surv(a,i,k)=1/(1+exp(-logit_s_surv(a,i,k)));
        }
      }
      if (sel_model_surv(i,k)==2){ //logistic selectivity
        //Survey parameters
        A50_surv(i,k)=max_A_surv(i,k)/(1+exp(-logit_A50_surv(i,k))); // logit scale for A50 with lower bound =0 and upper band =Aplus
        gamma_surv(i,k)=max_A_surv(i,k)/(1+exp(-logit_gamma_surv(i,k))); 
        //Survey selectivity
        for(int a = (min_A_surv(i,k)-1); a < max_A_surv(i,k); a++){
          s_surv(a,i,k)=1/(1+exp(-((a+1)-A50_surv(i,k))/gamma_surv(i,k)));
        }
        for(int a = (min_A_surv(i,k)-1); a < max_A_surv(i,k); a++){
          s_surv(a,i,k) = s_surv(a,i,k)/s_surv((max_A_surv(i,k)-1),i,k); // so selectivity is forced to be 1 for last age
        }
      }
      
      //Abundance indices from surveys
      for(int t = 0; t < Y; t++){
        for(int a = (min_A_surv(i,k)-1); a < max_A_surv(i,k); a++){
          I(t,a,i,k) = q_surv(i,k)*s_surv(a,i,k)*NAA(t,a,i)*exp(-Z(t,a,i)*prop_y_elapsed_surv(t,k));
          
          //Aggregated indices
          aggr_I(t,i,k) += I(t,a,i,k);
        }
        
        // NLL aggregated survey indices of abundance 
        if (flag_aggr_I(t,i,k)==1){ // flags which say if annual index estimates should be used or not for specific species
          nll_aggr_I -= dnorm(log(obs_aggr_I(t,i,k)),log(aggr_I(t,i,k)),sd_obs_aggr_I(t,i,k),1);
          if (data_simulate==1){
            SIMULATE {
              obs_aggr_I(t,i,k) = exp(rnorm(log(aggr_I(t,i,k)), sd_obs_aggr_I(t,i,k)));
              REPORT(obs_aggr_I);
            }
          } 
        }
        
        // Age composition surveys
        for(int a = (min_A_surv(i,k)-1); a < max_A_surv(i,k); a++){
          prop_Iaa(t,a,i,k) = I(t,a,i,k)/aggr_I(t,i,k);

          // NLL age composition surveys
          t_pred_paa(a) = prop_Iaa(t,a,i,k);
          t_paa(a) = obs_prop_Iaa(t,a,i,k);
        }
        // flags which say if annual index estimates should be used or not for specific species
        if (flag_Iaa(t,i,k)==1){
          nll_prop_Iaa -= get_acomp_ll(Y, min_A_surv(i,k)-1, max_A_surv(i,k), Neff_surv(t,i,k), age_comp_model_indices(i), t_paa, t_pred_paa, acomp_pars, index_aref(t,i,k));
          if (data_simulate==1){
            SIMULATE {
              vector<Type> iaa = sim_acomp(Y, min_A_surv(i,k)-1, max_A_surv(i,k), Neff_surv(t,i,k), age_comp_model_indices(i), t_paa, t_pred_paa, acomp_pars, index_aref(t,i,k));
              for(int a = (min_A_surv(i,k)-1); a < max_A_surv(i,k); a++){
                obs_prop_Iaa(t,a,i,k) = iaa(a);
                if (max_A_surv(i,k) < max_A){
                  for(int a = max_A_surv(i,k); a < max_A; a++){
                    obs_prop_Iaa(t,a,i,k) = 0;
                  }
                }
              }
              REPORT(obs_prop_Iaa);
            }
          }
        }
      }
    }
  }

  
  // Things to calculate for ADREPORT
  sumFAA.setZero();
  sumFy.setZero();
  sumMAA.setZero();
  sumMy.setZero();
  sumPAA.setZero();
  sumPy.setZero();
  sumZAA.setZero();
  sumZy.setZero();
  mean_FAA.setZero();
  mean_Fy.setZero();
  mean_MAA.setZero();
  mean_My.setZero();
  mean_PAA.setZero();
  mean_Py.setZero();
  mean_ZAA.setZero();
  mean_Zy.setZero();
  // To avoid adreporting large matrices
  for(int i = 0; i < sp; i++){
    for(int t = 0; t < Y; t++){
      for(int a = min_Fbar_idx(i)-1; a < max_Fbar_idx(i); a++){
        sumFy(t,i) += F(t,a,i);
        sumMy(t,i) += MAA(t,a,i);
        sumPy(t,i) += PAA(t,a,i);
        sumZy(t,i) += Z(t,a,i);
      }
    }
    for(int a = min_Fbar_idx(i)-1; a < max_Fbar_idx(i); a++){
      for(int t = 0; t < Y; t++){
        sumFAA(a,i) += F(t,a,i);
        sumMAA(a,i) += MAA(t,a,i);
        sumPAA(a,i) += PAA(t,a,i);
        sumZAA(a,i) += Z(t,a,i);
      }
    }
    for(int t = 0; t < Y; t++){
      for(int a = min_Fbar_idx(i)-1; a < max_Fbar_idx(i); a++){
        mean_FAA(a,i) = sumFAA(a,i)/Y;
        mean_Fy(t,i) = sumFy(t,i)/Fbar_range(i);
        mean_MAA(a,i) = sumMAA(a,i)/Y;
        mean_My(t,i) = sumMy(t,i)/Fbar_range(i);
        mean_PAA(a,i) = sumPAA(a,i)/Y;
        mean_Py(t,i) = sumPy(t,i)/Fbar_range(i);
        mean_ZAA(a,i) = sumZAA(a,i)/Y;
        mean_Zy(t,i) = sumZy(t,i)/Fbar_range(i);
      }
    }
    for(int t = 0; t < Y; t++){
      recruits(t,i) = NAA(t,0,i);
    }
  }

  
  
  NLL += nll_process_F;
  NLL += nll_aggr_Cw;
  NLL += nll_aggr_I;
  NLL += nll_prop_Caa;
  NLL += nll_prop_Iaa;
  NLL += nll_log_rec;
  NLL += nll_log_NAA;
  NLL += nll_log_MAA;
  NLL += nll_log_lorenzen1;
  NLL += nll_lorenzen2;
  // For trophic interactions
  NLL += nll_gamma_pref;
  NLL += nll_ratio_diet;
  NLL += nll_cons_rate;
  NLL += nll_alpha_cons;
  NLL += nll_beta_cons;
 
  
  
  
  
  //To report back to R
  REPORT(A50_F);
  REPORT(gamma_F);
  REPORT(s_F);
  REPORT(log_E);
  REPORT(pred_log_E);
  REPORT(sd_increment_process_F);
  REPORT(sd_process_F);
  REPORT(phi_process_F);
  REPORT(mean_log_E);
  REPORT(F);
  REPORT(MAA);
  // REPORT(pred_MAA);
  REPORT(Z);
  REPORT(NAA);
  REPORT(pred_NAA);
  REPORT(C);
  REPORT(Cw);
  REPORT(aggr_C);
  REPORT(aggr_Cw);
  REPORT(q_surv);
  REPORT(A50_surv);
  REPORT(gamma_surv);
  REPORT(s_surv);
  REPORT(I);
  REPORT(aggr_I);
  REPORT(prop_Caa);
  REPORT(prop_Iaa);
  REPORT(nll_aggr_Cw);
  REPORT(nll_aggr_I);
  REPORT(nll_prop_Caa);
  REPORT(nll_prop_Iaa);
  REPORT(nll_log_rec);
  REPORT(nll_log_NAA);
  REPORT(nll_log_MAA);
  REPORT(nll_log_lorenzen1);
  REPORT(nll_lorenzen2);
  // REPORT(exp(log_lorenzen1));
  // REPORT(lorenzen2);
  // REPORT(scale_M);
  REPORT(SSB);
  REPORT(nll_gamma_pref);
  REPORT(nll_cons_rate);
  REPORT(nll_alpha_cons);
  REPORT(nll_beta_cons);
  REPORT(nll_ratio_diet);
  REPORT(PAA);
  REPORT(NLL);
  REPORT(nll_process_F);
  REPORT(mean_Fy);
  REPORT(mean_FAA);
  REPORT(mean_My);
  REPORT(mean_MAA);
  REPORT(mean_Zy);
  REPORT(mean_ZAA);
  REPORT(recruits);
  if (predation_on==1){
  // For trophic interactions
    REPORT(mean_Py);
    REPORT(mean_PAA);
    // REPORT(Cfish);
    // REPORT(Cfishw);
    // REPORT(shape_gamma);
    // REPORT(scale_gamma);
    // REPORT(mode_gamma);
    // REPORT(density_mode_gamma);
    // REPORT(size_pref);
    // REPORT(log_ratio_w);
    // REPORT(sum_vuln_par);
    REPORT(vuln);
    REPORT(vuln_other);
    // REPORT(sum_vuln);
    // REPORT(suit);
    // REPORT(suit_other);
    // REPORT(sum_suit);
    // REPORT(scaled_suit);
    // REPORT(scaled_suit_other);
    // REPORT(biomass_prey_avail);
    // REPORT(biomass_prey_avail_no_age);
    // REPORT(biomass_other);
    // REPORT(biomass_other_avail);
    // REPORT(total_biomass_prey_avail);
    // REPORT(ratio_biomass_all_prey_avail_no_age);
    // REPORT(obs_cons_rate);
    // REPORT(cons_rate);
    // REPORT(alpha_cons);
    // REPORT(beta_cons);
    // REPORT(Cfish);
    // REPORT(Cfishw);
    // REPORT(N_pred);
    // REPORT(prob_b);
    // REPORT(prob_l_b);
    // REPORT(prob_l);
    // REPORT(prob_b_given_l);
    // REPORT(pred_ratio_l);
    // REPORT(pred_ratio_l_no_age);
    // REPORT(mean_ratio);
    // REPORT(mean_pred_ratio);
    // if (biomass_other_option==2){
    //   REPORT(total_predation_other);
    //   REPORT(growth_other);
    //   REPORT(K_other);
    // }
  }



  //For sdreport 
  ADREPORT(mean_Fy);
  ADREPORT(mean_FAA);
  ADREPORT(mean_My);
  ADREPORT(mean_MAA);
  ADREPORT(mean_Zy);
  ADREPORT(mean_ZAA);
  ADREPORT(recruits);
  ADREPORT(SSB);
  //ADREPORT(NAA);
  ADREPORT(aggr_Cw);
  ADREPORT(aggr_I);
  //ADREPORT(MAA);
  //ADREPORT(F);
  //ADREPORT(Z);
  if (predation_on==1){
    ADREPORT(mean_Py);
    ADREPORT(mean_PAA);
  //   ADREPORT(PAA);
  //   ADREPORT(biomass_other);
  //   ADREPORT(size_pref);
  //   ADREPORT(cons_rate);
  //   ADREPORT(alpha_cons);
  //   ADREPORT(beta_cons);
  }


  return(NLL);
  
  
}
