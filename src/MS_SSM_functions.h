/////////// Functions needed in the cpp code MS_SM.cpp ///////////////////

// Function that does the inverse logit transformation
template <class Type> 
Type inverse_logit(Type logit_x, Type lower, Type upper){
  Type x = (upper-lower)/(1+exp(-logit_x))+lower;
  return x;
}



// Function to force an argument to be positive
template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x, eps, Type(0.01) * pow(x-eps,2), Type(0));
  return CppAD::CondExpGe(x, eps, x, eps/(Type(2)-x/eps));
}

// Function that square x
template <class Type> 
Type square(Type x){return x*x;}


// Function with options for age composition data
template<class Type>
Type get_acomp_ll(int year, int min_age, int n_ages, Type Neff, int age_comp_model, vector<Type> paa_obs, vector<Type> paa_pred, vector<Type> age_comp_pars, int aref)
{
  Type zero = Type(0);
  Type one = Type(1);
  Type half = Type(0.5);
  Type two = Type(2);
  vector<Type> temp_n(n_ages);
  Type temp_Neff = zero, ll = zero;
  if(age_comp_model == 1) //multinomial
  {
    temp_Neff = Neff * exp(age_comp_pars(0));
    temp_n = temp_Neff * paa_obs;
    ll = lgamma(temp_Neff + one);
    for(int a = min_age; a < n_ages; a++) ll += -lgamma(temp_n(a) + one) + temp_n(a) * log(paa_pred(a) + Type(1e-15));
  } else {
    if(age_comp_model == 2) //dirichlet-multinomial
    {
      temp_Neff = Neff;
      temp_n = temp_Neff * paa_obs;
      ll = lgamma(temp_Neff + one) + lgamma(exp(age_comp_pars(0))) - lgamma(temp_Neff + exp(age_comp_pars(0)));
      for(int a = min_age; a < n_ages; a++) ll += -lgamma(temp_n(a) + one) + lgamma(temp_n(a) + exp(age_comp_pars(0)) * paa_pred(a)) -
        lgamma(exp(age_comp_pars(0)) * paa_pred(a));
    } else {
      if(age_comp_model == 3) //dirichlet
      {
        Type obs = zero, pred = zero, obs_2 = zero, pred_2 = zero;
        for(int a = aref-1; a < n_ages; a++)
        {
          obs_2 += paa_obs(a);
          pred_2 += paa_pred(a);
        }
        ll = lgamma(exp(age_comp_pars(0)));
        for(int a = min_age; a < aref-1; a++)
        {
          pred += paa_pred(a);
          obs += paa_obs(a);
          if(paa_obs(a) > Type(1.0e-15))
          {
            ll +=  -lgamma(exp(age_comp_pars(0)) * pred) + (exp(age_comp_pars(0)) * pred - one) * log(obs);
            pred = zero;
            obs = zero;
          }
          //else pooling with next age
        }
        //add in the last age class(es).
        ll += -lgamma(exp(age_comp_pars(0)) * pred_2) + (exp(age_comp_pars(0)) * pred_2 - one) * log(obs_2);
      } else {
        if(age_comp_model == 4) //zero-one inflated logistic normal. Inspired by zero-one inflated beta in Ospina and Ferrari (2012).
        {
          vector<Type> X(n_ages), p0(n_ages);
          Type mu = zero, sd = zero, pos_obs = zero, pos_pred = zero, pos_obs_l = zero, pos_pred_l = zero, pos_obs_sum = zero;
          Type pos_pred_sum = zero, y = zero;
          X = log(paa_pred + Type(1.0e-15)) - log(one - paa_pred + Type(1.0e-15));
          p0 = one/(one + exp(exp(age_comp_pars(1))*(X - age_comp_pars(0)))); //prob of zero declines with proportion caught
          sd = exp(age_comp_pars(2));
          int last_pos = 0;
          pos_obs_sum = sum(paa_obs);
          for(int a = min_age; a < n_ages; a++) if(paa_obs(a) > Type(1.0e-15))
          {
            pos_pred_sum += paa_pred(a);
            last_pos = a;
          }
          //logistic applies only to proportions of non-zero observations
          pos_obs_l = paa_obs(last_pos)/pos_obs_sum;
          pos_pred_l = paa_pred(last_pos)/pos_pred_sum;
          for(int a = min_age; a < n_ages; a++)
          {
            if(paa_obs(a) < Type(1.0e-15)) ll += log(p0(a) + Type(1.0e-15));
            else
            {
              ll += log(one - p0(a) + Type(1.0e-15));
              if(a < last_pos) //add in logistic-normal for positive observations less than last observed age class
              {
                pos_pred = paa_pred(a)/pos_pred_sum;
                pos_obs = paa_obs(a)/pos_obs_sum;
                y = log(pos_obs) - log(pos_obs_l);
                mu = log(pos_pred + Type(1.0e-15)) - log(pos_pred_l + Type(1.0e-15));
                ll += -half * (log(two * M_PI) + square((y - mu)/sd)) - log(sd) - log(pos_obs);
              }
            }
          }
          ll -= log(pos_obs_l); //add in the last observed age class(es).
        } else {
          if(age_comp_model == 5) //logistic normal. Pool zero observations with adjacent age classes.
          {
            Type mu = zero, sd = zero, obs = zero, pred = zero, obs_2 = zero, pred_2 = zero, y = zero;
            for(int a = aref-1; a < n_ages; a++)
            {
              obs_2 += paa_obs(a);
              pred_2 += paa_pred(a);
            }
            for(int a = min_age; a < aref-1; a++)
            {
              pred += paa_pred(a);
              obs += paa_obs(a);
              if(paa_obs(a) > Type(1.0e-15))
              {
                sd = exp(age_comp_pars(0)-half*log(Neff));
                y = log(obs) - log(obs_2);
                mu = log(pred + Type(1.0e-15)) - log(pred_2 + Type(1.0e-15));
                ll += -half * (log(two * M_PI) + square((y - mu)/sd)) - log(sd) - log(obs);
                pred = zero;
                obs = zero;
              }
              //else pooling with next age
            }
            ll -= log(obs_2); //add in the last age class(es).
          } else {
            if(age_comp_model == 6) //zero-one inflated logistic normal where p0 is a function of binomial sample size. 2 parameters
            {
              vector<Type> p0(n_ages);
              Type n_e = zero, mu = zero, sd = zero, pos_obs = zero, pos_pred = zero, pos_obs_l = zero, pos_pred_l = zero, pos_obs_sum = zero;
              Type pos_pred_sum = zero, y = zero;
              n_e = exp(age_comp_pars(0));
              p0 = exp(n_e * log(one-paa_pred + Type(1.0e-15))); //prob of zero declines with proportion caught
              sd = exp(age_comp_pars(1));
              int last_pos = 0;
              pos_obs_sum = sum(paa_obs);
              for(int a = min_age; a < n_ages; a++) if(paa_obs(a) > Type(1.0e-15))
              {
                pos_pred_sum += paa_pred(a);
                last_pos = a;
              }
              //logistic applies only to proportions of non-zero observations
              pos_obs_l = paa_obs(last_pos)/pos_obs_sum;
              pos_pred_l = paa_pred(last_pos)/pos_pred_sum;
              for(int a = min_age; a < n_ages; a++)
              {
                if(paa_obs(a) < Type(1.0e-15)) ll += log(p0(a) + Type(1.0e-15));
                else
                {
                  ll += log(one - p0(a) + Type(1.0e-15));
                  if(a < last_pos) //add in logistic-normal for positive observations less than last observed age class
                  {
                    pos_pred = paa_pred(a)/pos_pred_sum;
                    pos_obs = paa_obs(a)/pos_obs_sum;
                    y = log(pos_obs) - log(pos_obs_l);
                    mu = log(pos_pred + Type(1.0e-15)) - log(pos_pred_l + Type(1.0e-15));
                    ll += -half * (log(two * M_PI) + square((y - mu)/sd)) - log(sd) - log(pos_obs);
                  }
                }
              }
              ll -= log(pos_obs_l); //add in the last observed age class(es).
            } else {
              Rf_error("age_comp_model option for catch and/or survey does not exist, estimation impossible");
            }
          }
        }
      }
    }
  }
  return ll;
}


// Function to simulate age composition data
template<class Type>
vector<Type> sim_acomp(int year, int min_age, int n_ages, Type Neff, int age_comp_model, vector<Type> paa_obs, vector<Type> paa_pred, vector<Type> age_comp_pars, int aref)
{
  Type zero = Type(0);
  Type one = Type(1);
  Type half = Type(0.5);
  //Type two = Type(2);
  vector<Type> obs(n_ages);
  obs.setZero();
  if(age_comp_model == 1) //multinomial generated from condition binomials
  { 
    Type N = Neff * exp(age_comp_pars(0));
    obs(min_age) = rbinom(N, paa_pred(min_age));
    if (n_ages>1) {
      for(int a = (min_age+1); a < n_ages-1; a++) 
      {
        Type denom = one - paa_pred.head(a).sum();
        Type cond_N = N-obs.head(a).sum();
        if(denom > Type(1.0e-15)) if(cond_N > Type(1.0e-15))
        {
          Type cond_p = paa_pred(a)/denom; //.head first a components
          if(one - cond_p > Type(1.0e-15)) obs(a) = rbinom(cond_N,cond_p);
          else obs(a) = cond_N; //p pretty much 1
        }
      }
      obs(n_ages-1) = N - obs.sum();
    }
    obs = obs/obs.sum();// proportions
  } else {
    if(age_comp_model == 2) //dirichlet-multinomial. dirichlet generated from iid gammas and multinomial from conditional binomials
    { 
      Type N = Neff;
      vector<Type> ps(n_ages);
      for(int a = min_age; a < n_ages; a++) ps(a) = rgamma(paa_pred(a)*exp(age_comp_pars(0)),one);
      ps = ps/ps.sum();
      obs(min_age) = rbinom(N, ps(min_age));
      for(int a = (min_age+1); a < n_ages-1; a++) 
      {
        Type denom = one-ps.head(a).sum();
        Type cond_N = N-obs.head(a).sum(); //.head(a) first a components
        if(denom > Type(1.0e-15)) if(cond_N > Type(1.0e-15))
        {
          Type cond_p = ps(a)/denom; 
          if(one - cond_p > Type(1.0e-15)) obs(a) = rbinom(cond_N,cond_p);
          else obs(a) = cond_N; //p pretty much 1
        }
      }
      obs(n_ages-1) = N - obs.sum();
      obs = obs/obs.sum();// proportions
    } else {
      if(age_comp_model == 3) //dirichlet generated from iid gammas
      { 
        Type obs_2 = 0.0;
        vector<Type> best_obs(n_ages);
        for(int a = min_age; a < n_ages; a++) best_obs(a) = rgamma(paa_pred(a)*exp(age_comp_pars(0)),one);
        best_obs = best_obs/best_obs.sum();
        obs_2 = best_obs.tail(n_ages-aref+1).sum(); // .tail last n_ages-aref+1 components
        //obs_2 = best_obs.segment(aref-1,n_ages-1).sum(); 
        for(int a = aref-1; a < n_ages; a++) if(paa_obs(a) > Type(1.0e-15)) obs(a) = obs_2;
        obs_2 = zero;
        for(int a = min_age; a < aref-1; a++) 
        {
          obs_2 += best_obs(a);
          if(paa_obs(a) > Type(1.0e-15))
          {
            obs(a) = obs_2;
            obs_2 = zero;
          }
          else obs(a) = zero;
          //else pooling with next age
        }
      } else {
        if(age_comp_model == 4) //zero-one inflated logistic normal. Inspired by zero-one inflated beta in Ospina and Ferrari (2012).
        {
          vector<Type> X = log(paa_pred + Type(1.0e-15)) - log(one - paa_pred + Type(1.0e-15));
          vector<Type> p0 = one/(one + exp(exp(age_comp_pars(1))*(X - age_comp_pars(0)))); //prob of zero declines with proportion caught
          Type sd = exp(age_comp_pars(2));
          for(int a = min_age; a < n_ages; a++) obs(a) = rbinom(one, one - p0(a)); // generate instances of positive observations
          int n_pos = 0;
          for(int a = min_age; a < n_ages; a++) if(obs(a) > 0.5) n_pos++;
          if(n_pos>0)
          {
            vector<Type> pos_pred(n_pos);
            int k = 0;
            for(int a = min_age; a < n_ages; a++) if(obs(a) > 0.5)
            {
              pos_pred(k) = paa_pred(a);
              k++;
            }
            vector<Type> pos_obs(n_pos);
            pos_obs.setZero();
            for(int a = min_age; a < n_pos-1; a++) 
            {
              pos_obs(a) = exp(rnorm(log(pos_pred(a)) - log(pos_pred(n_pos-1)), sd));
            }
            pos_obs = pos_obs/(one + pos_obs.sum());
            pos_obs(n_pos-1) = one - pos_obs.sum();
            k = 0;
            for(int a = min_age; a < n_ages; a++) if(obs(a) > 0.5)
            {
              obs(a) = pos_obs(k);
              k++;
            }
          }
        } else {
          if(age_comp_model == 5) //logistic normal. Pool zero observations with adjacent age classes.
          {
            vector<Type> best_obs(n_ages);
            best_obs.setZero();
            Type sd = exp(age_comp_pars(0)-half*log(Neff));
            for(int a = min_age; a < n_ages-1; a++) best_obs(a) = exp(rnorm(log(paa_pred(a)) - log(paa_pred(n_ages-1)), sd));
            best_obs = best_obs/(one + best_obs.sum());
            best_obs(n_ages-1) = one - best_obs.sum();
            // for(int a = 0; a < n_ages; a++) best_obs(a) = exp(rnorm(log(paa_pred(a)) - log(paa_pred(n_ages-1)), sd));
            // best_obs = best_obs/(best_obs.sum());
            
            
            Type obs_2 = best_obs.tail(n_ages-aref+1).sum(); // .tail last n_ages-aref+1 components
            for(int a = aref-1; a < n_ages; a++) if(paa_obs(a) > Type(1.0e-15)) obs(a) = obs_2;
            obs_2 = zero;
            for(int a = min_age; a < aref-1; a++) 
            {
              obs_2 += best_obs(a);
              if(paa_obs(a) > Type(1.0e-15))
              {
                obs(a) = obs_2;
                obs_2 = zero;
              }
              else obs(a) = zero;
              //else pooling with next age
            }
          } else {
            if(age_comp_model == 6) //zero-one inflated logistic normal where p0 is a function of binomial sample size. 2 parameters
            {
              Type n_e = exp(age_comp_pars(0));
              vector<Type> p0 = exp(n_e * log(one-paa_pred + Type(1.0e-15))); //prob of zero declines with proportion caught
              Type sd = exp(age_comp_pars(1));
              
              for(int a = min_age; a < n_ages; a++) obs(a) = rbinom(one, one - p0(a)); // generate instances of positive observations
              int n_pos = min_age;
              for(int a = min_age; a < n_ages; a++) if(obs(a) > 0.5) n_pos++;
              if(n_pos>0)
              {
                vector<Type> pos_pred(n_pos);
                int k = min_age;
                for(int a = min_age; a < n_ages; a++) if(obs(a) > 0.5)
                {
                  pos_pred(k) = paa_pred(a);
                  k++;
                }
                vector<Type> pos_obs(n_pos);
                pos_obs.setZero();
                for(int a = min_age; a < n_pos-1; a++) 
                {
                  pos_obs(a) = exp(rnorm(log(pos_pred(a)) - log(pos_pred(n_pos-1)), sd));
                }
                pos_obs = pos_obs/(one + pos_obs.sum());
                pos_obs(n_pos-1) = one - pos_obs.sum();
                k = 0;
                for(int a = min_age; a < n_ages; a++) if(obs(a) > 0.5)
                {
                  obs(a) = pos_obs(k);
                  k++;
                }
              }
            } else {
              Rf_error("age_comp_model option for catch and/or survey does not exist, simulation impossible");
            }
          }
        }
      }
    }
  }
  return obs;
}



// Function delta-dirichlet for proportion of prey in diet
template<class Type>
Type ddeltadir(vector<Type> prop_obs, vector<Type> prop_pred, vector<Type> pars, int do_log)
{
  //application of the delta-dirichlet from Lewy 1996, Biometrics 52: 1394-1409
  //it is a a special case where we view observations as imperfectly observed Dirichlet variable: 
  //Categories are detected with probability p which is a declining function of a expected true Dirichlet proportions
  //Also assume that occurences of zeros are independent. Not sure how or why to model dependence.
  int n = prop_obs.size();
  vector<int> isdir(n);
  vector<Type> ll_1(n);
  ll_1.setZero();
  Type d = 0.0;
  Type phi = exp(pars(0)); //sum of dirichlet alphas
  Type beta1 = exp(pars(1)); //increasing
  Type beta0 = pars(2);
  vector<Type> pi(n), alphas = prop_pred*phi;
  
  for(int i = 0; i < n; i++) 
  {
    pi(i) = 1.0/(1.0 + exp(-(beta0 + beta1 * prop_pred(i)))); //probability of detection increases with true proportion
    if(prop_obs(i) > 1e-10) //greater than zero
    { 
      isdir(i) = 1;
      d += prop_pred(i);
      ll_1(i) += log(pi(i));
    }
    else 
    {
      isdir(i) = 0;
      ll_1(i) += log(1.0 - pi(i));
    }
  }
  //see(pi);
  
  int ndir = isdir.sum();
  Type ll = ll_1.sum(); //add up all the log-probabilities of zeros and non-zeros.
  if(ndir>0)
  {
    vector<Type> dir_alphas(ndir);
    int k = 0;
    for(int i = 0; i < n; i++) if(isdir(i) == 1) 
    {
      dir_alphas(k) = alphas(i);
      ll += (alphas(i)-1.0) * log(prop_obs(i)) -lgamma(alphas(i));
      k++;
    }
    ll += lgamma(dir_alphas.sum());
  }
  if(do_log) return(ll);
  else return(exp(ll));
}


// Function to simulate delta-dirichlet for proportion of prey in diet
template<class Type>
vector<Type> rdeltadir(vector<Type> prop_pred, vector<Type> pars)
{
  //application of the delta-dirichlet from Lewy 1996, Biometrics 52: 1394-1409
  //it is a a special case where we view observations as imperfectly observed Dirichlet variable: 
  //Categories are detected with probability p which is a declining function of a expected true Dirichlet proportions
  //Also assume that occurences of zeros are independent. Not sure how or why to model dependence.
  int n = prop_pred.size();
  Type phi = exp(pars(0)); //sum of dirichlet alphas
  Type beta1 = exp(pars(1)); //increasing
  Type beta0 = pars(2);
  vector<Type> alphas = prop_pred*phi;
  
  vector<Type> pi = Type(1)/(Type(1) + exp(-(beta0 + beta1 * prop_pred))); //probability of detection increases with true proportion
  vector<Type> obs = rbinom(Type(1), pi);
  // while (obs.sum()==0){
  //   obs = rbinom(Type(1), pi);
  // }
  vector<int> isdir(n);
  for(int i = 0; i < n; i++) 
  {
    if(obs(i) >0.9) isdir(i) = 1;
    else isdir(i) = 0;
  }
  int ndir = isdir.sum();
  //see(ndir);
  //see(isdir);
  //see(rbinom(Type(1),pi));
  //see(pi);
  if(ndir>0)
  {
    //vector<Type> dir_alphas(ndir);
    //int k = 0;
    for(int i = 0; i < n; i++) if(isdir(i) == 1) 
    {
      //dir_alphas(k) = alphas(i);
      obs(i) = rgamma(alphas(i),Type(1));
      //k++;
    }
    if (obs.sum()!=0)
      obs = obs/obs.sum(); //dirichlet using gamma random variables.
  }
  //see(obs);
  //see(obs.sum());
  return(obs);
}



template <class Type>
Type ddirichlet(vector<Type> obs, vector<Type>p, Type phi, int do_log) 
{
  int n = obs.size();
  Type ll = lgamma(phi);
  vector<Type> alphas = (p)*phi;
  for(int i = 0; i < n; i++) ll +=  -lgamma(alphas(i)) + (alphas(i) - Type(1.0)) * log(obs(i));
  if(do_log == 1) return(ll);
  else return(exp(ll));
}



template<class Type>
vector<Type> rdirichlet(vector<Type> p, Type phi)
{
  int n = p.size();
  vector<Type> alpha = p * phi;
  vector<Type> obs(n);
  for(int i = 0; i < n; i++) obs(i) = rgamma(alpha(i),Type(1.0));
  if (obs.sum()!=0) obs = obs/obs.sum();
  return (obs);
}
