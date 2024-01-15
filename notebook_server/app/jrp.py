import nest_asyncio
# hack to get Stan to work
nest_asyncio.apply()
import stan
import pandas as pd
import datetime
import arviz as az
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import ScalarFormatter
from matplotlib.patches import ConnectionPatch
from statsmodels.distributions.empirical_distribution import ECDF



class JRP:
    
    def __init__(self,file_data="/app/DATA/data.csv", file_thresholds="/app/DATA/thresholds.csv"):
        print("Loading event data ...")
        self.df_data = pd.read_csv(file_data,parse_dates=['date'])
        self.df_thresholds  = pd.read_csv(file_thresholds)

    def get_stan_data(self, 
                      marginal_years = np.geomspace(0.1,500,100),
                      L=20, 
                      P = 1000, 
                      n_grid=20,
                      t_RP_for_JRPS =np.array([20,20,200,200])*365, 
                      f_RP_for_JRPS = np.array([5,100,5,100])*365, 
                      random_seed=1
                      ):

        print("Preprocessing data for Stan ...")
        # random seed for Stan
        self.random_seed = random_seed

        # format data for Stan
        data = {
                    'events':self.df_data[['surge','flow']].values,
                }
        data['N'] = len(data['events'])

        # contemporray thresholds
        data['muf'] = self.df_thresholds[self.df_thresholds['year']==2023]['flow'].item()
        data['mut'] = self.df_thresholds[self.df_thresholds['year']==2023]['surge'].item()

        # future projections due to climate change
        data['muf100'] = self.df_thresholds[self.df_thresholds['year']==2100]['flow'].item()
        data['mut100'] = self.df_thresholds[self.df_thresholds['year']==2100]['surge'].item()

        # mean inter event period
        lam = self.df_data['date'].sort_values().diff().dropna().mean().total_seconds()/(60*60*24)
        data['lam'] = lam


        # get the return probability for the marginal return years we are interested in
        marginal_days = marginal_years*365

        returnProb = 1-lam/marginal_days

        self.marginal_days = marginal_days

        self.marginal_years = marginal_years

        data['returnProb'] = returnProb
        data['M'] = len(data['returnProb'])

        # get the tides and flows for plotting pdfs

        
        tides_for_pdf = np.linspace(
            self.df_data['surge'].min(),
            self.df_data['surge'].max(),
            L
        )
        flows_for_pdf = np.linspace(
            self.df_data['flow'].min(),
            self.df_data['flow'].max(),
            L
        )
        data['tides_for_pdf'] = tides_for_pdf
        data['flows_for_pdf'] = flows_for_pdf

        data['L'] = len(data['tides_for_pdf'])


        data['tmax'] = self.df_data['surge'].max()
        data['fmax'] = self.df_data['flow'].max()

        _grid_t = np.linspace(self.df_data['surge'].min(), self.df_data['surge'].max(), n_grid)
        _grid_f = np.linspace(self.df_data['flow'].min(), self.df_data['flow'].max(), n_grid)
        grid_t, grid_f = np.meshgrid(_grid_t, _grid_f)
        grid_t, grid_f = grid_t.flatten(), grid_f.flatten()
        grid_t, grid_f 

        data['tide_grid'] = grid_t
        data['flow_grid'] = grid_f
        data['O'] = len(data['flow_grid'])
        data['lam'] = lam


        
        data['v1'] = np.random.uniform(0,1,P)
        data['v2'] = np.random.uniform(0,1,P)
        data['P'] = P



        tides_ecdf = ECDF(self.df_data['surge'])
        flow_ecdf = ECDF(self.df_data['flow'])


        data['t_quantiles'] = tides_ecdf(self.df_data['surge'])
        data['f_quantiles'] = flow_ecdf(self.df_data['flow'])


        data['t_RP_for_JRPS'] = t_RP_for_JRPS 
        data['f_RP_for_JRPS'] = f_RP_for_JRPS 
        data['Q'] = 4

        self.data = data

    def build_stan_model(self, file_model="/app/model.stan"):

        print("Building Stan model ...")
        with open(file_model) as f:
            model = f.read()
    

        self.model = stan.build(model, data=self.data, random_seed=1)

    def run_stan_model(self, num_chains=4, num_samples=2000):
        print("Running Stan model ...")

        self.result = self.model.sample(num_chains=num_chains, num_samples=num_samples)

        self.summary = az.summary(self.result)

    def generate_plots(self):
        print("Generating plots ...")
        
        self.generate_trace_plot()

        self.generate_corner_plot()

        self.generate_QQ_plot()

        self.generate_jrp_plot()



    def generate_trace_plot(self,filename='/app/DATA/trace_plot.png'):
        print("Generating trace plot ...")
        az.plot_trace(self.result)
        plt.tight_layout()
        plt.savefig(filename)
        return()


    def generate_corner_plot(self,filename='/app/DATA/corner_plot.png'):
        print("Generating corner plot ...")
        az.plot_pair(self.result,
            var_names = ['sigmat','sigmaf','xit','xif','theta'],
            marginals=True,  kind='kde', divergences=True,figsize=(10,10))
        plt.savefig(filename)


    def generate_jrp_plot(self,filename='app/DATA/jrp_plot.png'):
        print("Generating JRP plot ...")
	
        N_samples =20

        every_n_samples = int(8000/N_samples)

        alpha =0.5

        temp = (az
			.from_pystan(self.result)
			.posterior
			.stack(dimensions={'sample':['chain','draw']})
		)

	

        lam = self.df_data['date'].sort_values().diff().dropna().mean().total_seconds()/(60*60*24)


        fig, axs = plt.subplots(2,2,figsize=(10,10))

        self.df_data.plot.scatter(
			    x='surge',
			    y='flow',
			    ax=axs[0,1],
			    color='lightgrey',
			    marker='.'
			)

        year_contours = np.array([10,50,100,200,500])

        day_contours =  year_contours*365

        marginal_years = self.marginal_years

        jrp_50 = np.percentile(temp['jrp_grid'].values, 50,axis=1)

        for jrp_i in temp['jrp_grid'].values[:,::every_n_samples].T:
                axs[0,1].tricontour(
		            self.data['tide_grid'],
		            self.data['flow_grid'],
		            jrp_i,
		            levels=day_contours,
		            cmap='hsv',
		            alpha=0.1,
		            norm=colors.LogNorm(
		                        vmin=jrp_50.min(), 
		                        vmax=jrp_50.max()
		                                )  
	    )



        CS = axs[0,1].tricontour(
		            self.data['tide_grid'],
		            self.data['flow_grid'],
		            jrp_50,
		            levels=day_contours,
		            alpha=1,
		            cmap='hsv',
		            norm=colors.LogNorm(
		                        vmin=jrp_50.min(), 
		                        vmax=jrp_50.max()
		                                )
		            )    
	    
	 
        axs[0,1].clabel(
	    CS, 
	    CS.levels,
	    fmt={ key:item for key, item in zip(day_contours, year_contours)},
	    inline=True, 
	    fontsize=10
	)


        ################################### TIDES #########################################

        axs[1,1].fill_betweenx(
	    marginal_years,
	    np.percentile(temp['return_t'].values,20,axis=1),
	    np.percentile(temp['return_t'].values,80,axis=1),
	    alpha=0.1      
	)
        axs[1,1].plot(
		np.percentile(temp['return_t'].values,50,axis=1),
	    marginal_years
	)
        axs[1,1].plot(
		np.percentile(temp['return_t'].values,20,axis=1),
	    marginal_years,
	    linestyle='dotted',
	     c='#1f77b4'
	)
        axs[1,1].plot(
		np.percentile(temp['return_t'].values,80,axis=1),
		marginal_years,
		linestyle='dotted',
	     c='#1f77b4'

	)

        axs[1,1].grid(True, which="both", ls="-", color='0.65')
        axs[1,1].set_yscale('log')
        axs[1,1].yaxis.set_major_formatter(ScalarFormatter())

        ################################### Flow #########################################

        axs[0,0].fill_between(
	    marginal_years,
	    np.percentile(temp['return_f'].values,20,axis=1),
	    np.percentile(temp['return_f'].values,80,axis=1),
	    alpha=0.1      
	)
        axs[0,0].plot(
	    marginal_years,
		np.percentile(temp['return_f'].values,50,axis=1),
	    
	)
        axs[0,0].plot(
		marginal_years,

		np.percentile(temp['return_f'].values,20,axis=1),
	    linestyle='dotted',
	     c='#1f77b4'
	)
        axs[0,0].plot(
		    marginal_years,

		np.percentile(temp['return_f'].values,80,axis=1),
		linestyle='dotted',
	     c='#1f77b4'

	)

        axs[0,0].grid(True, which="both", ls="-", color='0.65')
        axs[0,0].set_xscale('log')
        axs[0,0].xaxis.set_major_formatter(ScalarFormatter())



        tide_lim = [0.4,2.0]
        axs[0,1].set_xlim(tide_lim)
        axs[1,1].set_xlim(tide_lim)

        flow_lim = [0.5,7]
        axs[0,1].set_ylim(flow_lim)
        axs[0,0].set_ylim(flow_lim)



        axs[1,0].axis('off')



        # Hide X and Y axes label marks
        axs[0,1].xaxis.set_tick_params(
	    labelbottom=False,
	    labeltop=True,
	    bottom=False, 
	    top=True, 

	)
        axs[0,1].yaxis.set_tick_params(
	    labelleft=False,
	    labelright=True,
	    left=False, 
	    right=True
	)
        # Hide X and Y axes tick marks
        #axs[0,1].set_xticks([])
        #axs[0,1].set_yticks([])
        axs[0,1].set_xlabel('')
        axs[0,1].set_ylabel('')


        axs[1,1].set_xlabel("Tide m")
        axs[1,1].set_ylabel("Return period yr")


        axs[0,0].set_ylabel("Discharge m3 s-1")
        axs[0,0].set_xlabel("Return period yr")


        axs[0,1].set_xlabel("Tide m")
        axs[0,1].set_ylabel("Discharge m3 s-1")
        axs[0,1].yaxis.set_label_position("right")
        axs[0,1].xaxis.set_label_position("top")



        axs[0,0].xaxis.set_tick_params(
	    labeltop=True,
	    top=True
	)

        axs[1,1].yaxis.set_tick_params(
	    labelright=True,
	    right=True
	)




        get_t = lambda rp :  np.interp(

		                rp,
		                marginal_years,
		                np.percentile(temp['return_t'].values,50,axis=1),


		                )
        get_f = lambda rp :  np.interp(

		                rp,
		                marginal_years,
		                np.percentile(temp['return_f'].values,50,axis=1),


		                )


        plt.subplots_adjust(hspace=0.02,wspace=0.02)

        plt.savefig(filename)


    def generate_QQ_plot(self,filename='/app/DATA/qq_plot.png'):
        print("Generating QQ plot ...")
        pass


    def infer(self):

        self.get_stan_data()

        self.build_stan_model()

        self.run_stan_model()

        self.generate_plots()



