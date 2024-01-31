
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
                real mu1, real mu2, 
                real xi1, real xi2, 
                real sigma1, real sigma2,
                real theta){
    
        real p;
        real ctide;
        real cflow;
        real ccop;
        vector[1] t;
        vector[1] f;
        
        t[1] = tide;
        f[1] = flow;
    
        ctide =  gpareto_cdf(t, mu1, xi1, sigma1);
        cflow =  gpareto_cdf(f, mu2, xi2, sigma2);
        ccop = gumbel_copula_cdf(ctide, cflow, theta);

        p = 1 -ctide -cflow + ccop;

        return p;
    
    }
    
    real JRP(real tide, real flow, 
                real mu1, real mu2, 
                real xi1, real xi2, 
                real sigma1, real sigma2,
                real theta, 
                real lam){

        real pand;
        real jrp;

        pand = P_and( tide,  flow, 
                 mu1,  mu2, 
                 xi1,  xi2, 
                 sigma1,  sigma2,
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
        real mu1, real xi1, real sigma1, 
        real mu2, real xi2, real sigma2, 
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
        
        lpt = gpareto_lpdf(t | mu1, xi1, sigma1); 
        lpf = gpareto_lpdf(f | mu2, xi2, sigma2);
        ct = gpareto_cdf(t, mu1, xi1, sigma1); 
        cf = gpareto_cdf(f, mu2, xi2, sigma2); 

        lpcopula = gumbel_copula_lpdf(ct | cf, theta);
        
        
        return lpcopula + lpt + lpf;
        
    }


}

data {


    int N_evnt; //  number of extreme events prodicuced by Eleonora's code
    int N_grid; //  total number of points in grid for JRP contour plot
    int N_mrp; // number of points to calculate marginal return periods

    matrix[N_evnt,2] events; // selected extreme events
    vector[N_evnt] quantiles_1; // marginal empirical quantiles for driver 1 to calculate for events for Q-Q plot
    vector[N_evnt] quantiles_2; // marginal empirical quantiles for driver 2 to calculate for events for Q-Q plot

    vector[N_grid] grid1; // points for JRP grid for driver 1
    vector[N_grid] grid2; // points for JRP grid for driver 2


    vector[N_mrp] mrp_prob; // marginal return probabilities for calculating driver magnitudes for a given return period (note these are probabilities!)



    real mu1; // driver 1 threshold used for mean of driver 1 marginal GPD
    real mu2; // driver 2 threshold used for mean of driver 2 marginal GPD
    real max1; // driver 1 max value for placing limits on xi1
    real max2; // driver 2 max value for placing limits on xi2
    real lam; // mean time interval between events in days


    
    
}

transformed data {
    int x_i[0];

}


parameters {


    // generalised pareto distribution parameters
    // for tides ...
    real<lower=0> sigma1;
    real<lower=-sigma1/(max1-mu1)> xi1; 
    
    // ... and for flow:
    real<lower=0> sigma2;
    real<lower=-sigma2/(max2-mu2)> xi2; 

    
    // copula parameters
    real<lower=1> theta;

}



model {


    for (i in 1:N_evnt){
        events[i] ~ joint_pdf(  mu1,  xi1,  sigma1, 
                                mu2,  xi2,  sigma2, 
                                theta 
                                );

    };

}

generated quantities {

    vector[N_grid] jrp_grid;
    vector[N_mrp]  mrp_1;
    vector[N_mrp]  mrp_2;

    vector[N_evnt] jrp_for_events; // JRP values for each event events
    vector[N_evnt] rp1_for_events; // marginal RP for driver 1 for each event
    vector[N_evnt] rp2_for_events; // marginal RP for driver 2 for each event

    vector[N_evnt] d1_for_quantile; // values for quantiles for QQ plot
    vector[N_evnt] d2_for_quantile; // values for quantiles for QQ plot

    matrix[N_evnt,2] events_hat; // posterior predictive for pppc plot



    // calculate JRP grid for making contour plot
   for (i in 1:N_grid){
        jrp_grid[i] = JRP(grid1[i], grid2[i], 
                mu1, mu2, 
                xi1, xi2, 
                sigma1, sigma2,
                theta, 
                lam);
    }

    // calculate values of the drivers for different returnperiods
    // (here return period is expressed as return probability for the ppf function)
    for (i in 1:N_mrp){
        mrp_1[i] = gpareto_ppf(mrp_prob[i], mu1, xi1, sigma1);
        mrp_2[i] = gpareto_ppf(mrp_prob[i], mu2, xi2, sigma2);
    
    }

    // calculate the joint return period and marginal return periods for each event
    for (i in 1:N_evnt){

        vector[1] d1;
        vector[1] d2;

        jrp_for_events[i] = JRP(events[i,1], events[i,2], 
                mu1, mu2, 
                xi1, xi2, 
                sigma1, sigma2,
                theta, 
                lam);

        d1[1] = events[i,1];
        d2[1] = events[i,2];

        rp1_for_events[i] = lam/(1-gpareto_cdf(d1| mu1, xi1, sigma1));
        
        rp2_for_events[i] = lam/(1-gpareto_cdf(d2| mu2, xi2, sigma2));    
    }

    // calculate the values for the Q-Q plot
    for (i in 1:N_evnt){
        d1_for_quantile[i] = gpareto_ppf(quantiles_1[i], mu1, xi1, sigma1);
        d2_for_quantile[i] = gpareto_ppf(quantiles_2[i], mu2, xi2, sigma2);
    
    }


    // ppc plot
    for (i in 1:N_evnt){


      vector[1] thetas;
      vector[1] w_guess;
      vector[1] w;
      vector[1] u1;
      vector[1] u2;
      real v1;
      real v2[1];  
      real d1_hat;
      real d2_hat;  

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
    
        d1_hat = gpareto_ppf(u1[1], mu1, xi1, sigma1); 
        
        d2_hat = gpareto_ppf(u2[1], mu2, xi2, sigma2);
        
        events_hat[i,1] = d1_hat;
        events_hat[i,2] = d2_hat;
    
    }    

}

