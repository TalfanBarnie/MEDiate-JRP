
functions {

    //***************************************************************************************
    //
    // Aki Vehtari's Generalised Pareto functions from 
    // https://mc-stan.org/users/documentation/case-studies/gpareto_functions.html
    //
    //***************************************************************************************
    

    real gpareto_lpdf(vector y, real ymin, real k, real sigma) {
        // generalised Pareto log pdf 
        int N = rows(y);
        real inv_k = inv(k);
        if (k<0 && max(y-ymin)/sigma > -inv_k)
          reject("k<0 and max(y-ymin)/sigma > -1/k; found k, sigma =", k, sigma);
        if (sigma<=0)
          reject("sigma<=0; found sigma =", sigma);
        if (fabs(k) > 1e-15)
          return -(1+inv_k)*sum(log1p((y-ymin) * (k/sigma))) -N*log(sigma);
        else
          return -sum(y-ymin)/sigma -N*log(sigma); // limit k->0
    }
  
  
    real gpareto_cdf(vector y, real ymin, real k, real sigma) {
        // generalised Pareto cdf
        real inv_k = inv(k);
        if (k<0 && max(y-ymin)/sigma > -inv_k)
          reject("k<0 and max(y-ymin)/sigma > -1/k; found k, sigma =", k, sigma);
        if (sigma<=0)
          reject("sigma<=0; found sigma =", sigma);
        if (fabs(k) > 1e-15)
          return exp(sum(log1m_exp((-inv_k)*(log1p((y-ymin) * (k/sigma))))));
        else
          return exp(sum(log1m_exp(-(y-ymin)/sigma))); // limit k->0
    }
  
    real gpareto_lcdf(vector y, real ymin, real k, real sigma) {
        // generalised Pareto log cdf
        real inv_k = inv(k);
        if (k<0 && max(y-ymin)/sigma > -inv_k)
          reject("k<0 and max(y-ymin)/sigma > -1/k; found k, sigma =", k, sigma);
        if (sigma<=0)
          reject("sigma<=0; found sigma =", sigma);
        if (fabs(k) > 1e-15)
          return sum(log1m_exp((-inv_k)*(log1p((y-ymin) * (k/sigma)))));
        else
          return sum(log1m_exp(-(y-ymin)/sigma)); // limit k->0
    }
  
    real gpareto_lccdf(vector y, real ymin, real k, real sigma) {
        // generalised Pareto log ccdf
        real inv_k = inv(k);
        if (k<0 && max(y-ymin)/sigma > -inv_k)
          reject("k<0 and max(y-ymin)/sigma > -1/k; found k, sigma =", k, sigma);
        if (sigma<=0)
          reject("sigma<=0; found sigma =", sigma);
        if (fabs(k) > 1e-15)
          return (-inv_k)*sum(log1p((y-ymin) * (k/sigma)));
        else
          return -sum(y-ymin)/sigma; // limit k->0
    }
  
    real gpareto_rng(real ymin, real k, real sigma) {
        // generalised Pareto rng
        if (sigma<=0)
          reject("sigma<=0; found sigma =", sigma);
        if (fabs(k) > 1e-15)
          return ymin + (uniform_rng(0,1)^-k -1) * sigma / k;
        else
          return ymin - sigma*log(uniform_rng(0,1)); // limit k->0
    }
    



    //***************************************************************************************
    //
    // Ben Goodrich's Copula functions from 
    // https://spinkney.github.io/helpful_stan_functions/group__copula.html
    //
    //***************************************************************************************



    real gumbel_copula_lpdf(real u, real v, real theta) {
      real neg_log_u = -log(u);
      real log_neg_log_u = log(neg_log_u);
      real neg_log_v = -log(v);
      real log_neg_log_v = log(neg_log_v);
      real log_temp = log_sum_exp(theta * log_neg_log_u, theta * log_neg_log_v);
      real theta_m1 = theta - 1;

      if (theta < 1) 
        reject("theta must be >= 1");

      if (is_inf(theta)) {
        if (u == v) 
          return 0;
        else 
          return negative_infinity();
      }

      return theta_m1 * log_neg_log_u + theta_m1 * log_neg_log_v + neg_log_u + neg_log_v
             - exp(log_temp / theta)
             + log_sum_exp(2 * theta_m1 / -theta * log_temp,
                           log(theta_m1) + (1 - 2 * theta) / theta * log_temp);
    }

    //***************************************************************************************
    //
    // Other needed functions 
    //
    //***************************************************************************************
    
    real gpareto_ppf(real F, real mu, real xi, real sigma) {

        

        if (fabs(xi) > 1e-15)
          return mu + (sigma/xi)*((1-F)^(-xi)-1);
          
        else
          
          return mu - sigma*log(1-F); // limit k->0
    }



    real gumbel_generator(real t, real theta){
    
        return (-log(t))^theta;
    }
    
    real gumbel_inverse_generator(real t, real theta){
    
        return exp(-t^(1/theta));
    }
    
    real gumbel_K_function(real w, real theta ){
        
        return w*(1- log(w)/theta);
    
    }
    
    vector gumbel_invK_system(
              vector y,        // unknowns
              vector theta,    // parameters
              real[] x_r,      // data (real)
              int[] x_i) {     // data (integer)
              
              vector[1] z;
              real v;
              real theta0;
              real w;
              
              w = y[1];
              theta0 = theta[1];
              v = x_r[1];
              z[1] = v - gumbel_K_function(w, theta0);
              
              
              
              return z;
              
              }

    real gumbel_copula_cdf(real u, real v, real theta) {
    
        real cdf;
        
        cdf = exp(
                    -(
                        (-log(u))^theta +
                        (-log(v))^theta
                    )^(1/theta)
                );
                
        return cdf;
    
    }
    
    real P_and(real tide, real flow, 
                real mut, real muf, 
                real xit, real xif, 
                real sigmat, real sigmaf,
                real theta){
    
        real p;
        real ctide;
        real cflow;
        real ccop;
        vector[1] t;
        vector[1] f;
        
        t[1] = tide;
        f[1] = flow;
    
        ctide =  gpareto_cdf(t, mut, xit, sigmat);
        cflow =  gpareto_cdf(f, muf, xif, sigmaf);
        ccop = gumbel_copula_cdf(ctide, cflow, theta);

        p = 1 -ctide -cflow + ccop;

        return p;
    
    }
    
    real JRP(real tide, real flow, 
                real mut, real muf, 
                real xit, real xif, 
                real sigmat, real sigmaf,
                real theta, 
                real lam){

        real pand;
        real jrp;

        pand = P_and( tide,  flow, 
                 mut,  muf, 
                 xit,  xif, 
                 sigmat,  sigmaf,
                 theta);

        jrp = lam/pand;

        return jrp;
        }
   

    //***************************************************************************************
    //
    // Joint PDF functions 
    //
    //***************************************************************************************

    real joint_pdf_lpdf(
        row_vector tf, 
        real mut, real xit, real sigmat, 
        real muf, real xif, real sigmaf, 
        real theta){
    
        vector[1] t;
        vector[1] f;
        real lpt;
        real lpf;
        real ct;
        real cf;
        real lpcopula;
        
        t[1] = tf[1];
        f[1] = tf[2];
        
        lpt = gpareto_lpdf(t | mut, xit, sigmat); 
        lpf = gpareto_lpdf(f | muf, xif, sigmaf);
        ct = gpareto_cdf(t, mut, xit, sigmat); 
        cf = gpareto_cdf(f, muf, xif, sigmaf); 

        lpcopula = gumbel_copula_lpdf(ct | cf, theta);
        
        
        return lpcopula + lpt + lpf;
        
    }


}

data {

    int N; // number of events
    int M; // number of return probabilities
    int L; // number of evaluations of gpdf
    int O; // number of tide x flow pairs
    int P; // number of ppc's
    int Q; // number of combinations
    
    matrix[N,2] events; // N x [tide, flow]
    real mut; // tide threshold used for mean of tide marginal GPD
    real muf; // flow threshold used for mean of flow marginal GPD
    real mut100; // tide threshold used for mean of tide marginal GPD
    real muf100; // flow threshold used for mean of flow marginal GPD
    
    vector[M] returnProb;
    vector[L] tides_for_pdf;
    vector[L] flows_for_pdf;
    real tmax;
    real fmax;
    
    vector[O] tide_grid;
    vector[O] flow_grid;

    real lam;
    
    //vector[P] v1;
    //vector[P] v2;
    
    vector[N] t_quantiles;
    vector[N] f_quantiles;
    
    vector[Q] t_RP_for_JRPS;
    vector[Q] f_RP_for_JRPS;
}

transformed data {

  int x_i[0];
  //real x_r[P];
  
  
  //for (i in 1:P){
  //    x_r[i] = v2[i];
  //}
  

}


parameters {


    // generalised pareto distribution parameters
    // for tides ...
    real<lower=0> sigmat;
    real<lower=-sigmat/(tmax-mut)> xit; 
    
    // ... and for flow:
    real<lower=0> sigmaf;
    real<lower=-sigmaf/(fmax-muf)> xif; 

    
    // copula parameters
    real<lower=1> theta;

}



model {


    for (i in 1:N){
        events[i] ~ joint_pdf(  mut,  xit,  sigmat, 
                                muf,  xif,  sigmaf, 
                                theta 
                                );

    };

}

generated quantities {

    vector[M] return_t;
    vector[M] return_f;
    vector[L] tides_gpdf;
    vector[L] flows_gpdf;
    vector[O] jrp_grid;
    vector[N] jrp_events;
    vector[N] rp_tides;
    vector[N] rp_flows;
    vector[1] t;
    vector[1] f;
    matrix[N,2] events_hat;
    vector[N] t_for_quantile;
    vector[N] f_for_quantile;
    
    vector[Q] t_for_JRPS;
    vector[Q] f_for_JRPS; 
    vector[Q] t_for_JRPS100;
    vector[Q] f_for_JRPS100; 
    
    vector[Q] JRP_for_RP_pair;
    vector[Q] JRP_for_RP_pair100;
    

    vector[1] thetas;
    vector[1] w_guess;
    vector[1] w;
    vector[1] u1;
    vector[1] u2;
    real t_hat;
    real f_hat;
    real v1;
    real v2[1];
    
    
        

    for (i in 1:N){
        t_for_quantile[i] = gpareto_ppf(t_quantiles[i], mut, xit, sigmat);
        f_for_quantile[i] = gpareto_ppf(f_quantiles[i], muf, xif, sigmaf);
    
    }

    
    for (i in 1:M){
        return_t[i] = gpareto_ppf(returnProb[i], mut, xit, sigmat);
        return_f[i] = gpareto_ppf(returnProb[i], muf, xif, sigmaf);
    
    }
    
    
    for (i in 1:L){
        tides_gpdf[i] = exp(gpareto_lpdf(tides_for_pdf[i:i]| mut, xit, sigmat));
        flows_gpdf[i] = exp(gpareto_lpdf(flows_for_pdf[i:i]| muf, xif, sigmaf));
    }
    
    for (i in 1:O){
        jrp_grid[i] = JRP(tide_grid[i], flow_grid[i], 
                mut, muf, 
                xit, xif, 
                sigmat, sigmaf,
                theta, 
                lam);
    }
    
    for (i in 1:N){
        jrp_events[i] = JRP(events[i,1], events[i,2], 
                mut, muf, 
                xit, xif, 
                sigmat, sigmaf,
                theta, 
                lam);
                
        t[1] = events[i,1];
        f[1] = events[i,2];
                
        rp_tides[i] = lam/(1-gpareto_cdf(t| mut, xit, sigmat));
        
        rp_flows[i] = lam/(1-gpareto_cdf(f| muf, xif, sigmaf));
    }
    
    
    for (i in 1:Q){
    
        t_for_JRPS[i] = gpareto_ppf(1 - lam/t_RP_for_JRPS[i], mut, xit, sigmat);
        f_for_JRPS[i] = gpareto_ppf(1 - lam/f_RP_for_JRPS[i], muf, xif, sigmaf);
        t_for_JRPS100[i] = gpareto_ppf(1 - lam/t_RP_for_JRPS[i], mut100, xit, sigmat);
        f_for_JRPS100[i] = gpareto_ppf(1 - lam/f_RP_for_JRPS[i], muf100, xif, sigmaf);
        
        JRP_for_RP_pair[i] = JRP(t_for_JRPS[i], f_for_JRPS[i], 
                mut, muf, 
                xit, xif, 
                sigmat, sigmaf,
                theta, 
                lam);
        JRP_for_RP_pair100[i] = JRP(t_for_JRPS100[i], f_for_JRPS100[i], 
                mut100, muf100, 
                xit, xif, 
                sigmat, sigmaf,
                theta, 
                lam);        
    }
  
    for (i in 1:N){

       
        v1  = uniform_rng(0,1);
        v2[1]  = uniform_rng(0,1);
        
        thetas[1] = theta;
        //w_guess[1] =  x_r[i]; //0.5;
        w_guess[1] =  v2[1]; //0.5;
        
        //w = algebra_solver(gumbel_invK_system, w_guess, thetas,  x_r[i:i], x_i, 1e-10, 1e-6,1e4);
        w = algebra_solver(gumbel_invK_system, w_guess, thetas,  v2, x_i, 1e-10, 1e-6,1e4);

        
        //u1 = exp(v1[i]^(1/theta) * log(w));
        
        //u2 = exp((1-v1[i])^(1/theta)*log(w));
        
        u1 = exp(v1^(1/theta) * log(w));
        
        u2 = exp((1-v1)^(1/theta)*log(w));
    
        t_hat = gpareto_ppf(u1[1], mut, xit, sigmat); 
        
        f_hat = gpareto_ppf(u2[1], muf, xif, sigmaf);
        
        events_hat[i,1] = t_hat;
        events_hat[i,2] = f_hat;
    
    }    
   




}

