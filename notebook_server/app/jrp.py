from statsmodels.distributions.empirical_distribution import ECDF
from dataclasses import dataclass
import pandas as pd
import numpy as np
import nest_asyncio
# hack to get Stan to work
nest_asyncio.apply()
import stan
import arviz as az
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import ScalarFormatter
from statsmodels.distributions.empirical_distribution import ECDF
from matplotlib.patches import ConnectionPatch
import datetime
import json



def convert_detectExtremeEvents_xls_to_json(file):
    """Helper function to convert the output of detectExtremeEvents to json for 
    testing the api
    """
    data = pd.read_excel(file,sheet_name='data')

    info = pd.read_excel(file,sheet_name='info')

    inpt = pd.read_excel(file,sheet_name='input')
    
    data_json = data.to_json()
    
    info_json = info.to_json()
    
    inpt_json = inpt.to_json()
    
    data = {}

    data['inpt'] = json.loads(inpt_json)
    
    data['info'] = json.loads(info_json)
    
    data['data'] = json.loads(data_json)
        
    with open("data.json","w") as f:
        json.dump(data, f)



@dataclass
class Driver:
    name : str
    unit : str
    threshold_value: float
    threshold_value_unit : str
    threshold_percent : float
    threshold_percent_unit : str
    min_peaks_distance : float
    min_peaks_distance_unit : str
    POT : float
    minimum : float
    maximum : float

        
@dataclass
class PotsInfo:
    correlation_coefficient : float
    search_window_value : float
    search_window_units : str
    number_of_joint_events : int
    dt_value : str
    dt_units : str
    

class JRP:
    
    def __init__(self, data, info, inpt,
                num_chains=4, num_samples=2000, 
                 year_contours = np.array([10,50,100,200,500]),
                 marginal_years = np.geomspace(0.1,500,100)
                
                ):
        

        # get values we need later
        min_1 = data['Driver1'].min()
        max_1 = data['Driver1'].max()
        min_2 = data['Driver2'].min()
        max_2 = data['Driver2'].max()
        
        driver1 = {
         'name':info.loc[info['info']=='Driver1','value'].astype(str).item(),
         'unit':info.loc[info['info']=='Driver1','unit'].astype(str).item(),
         'threshold_value':info.loc[info['info']=='Threshold Driver1','value'].astype(float).item(),
         'threshold_value_unit':info.loc[info['info']=='Threshold Driver1','unit'].astype(str).item(),
         'threshold_percent':inpt.loc[inpt['INPUT']=='threshold 1','value'].astype(float).item(),
         'threshold_percent_unit':inpt.loc[inpt['INPUT']=='threshold 1','unit'].astype(str).item(),
         'min_peaks_distance':inpt.loc[inpt['INPUT']=='min peaks distance 1','value'].astype(int).item(),
         'min_peaks_distance_unit':inpt.loc[inpt['INPUT']=='min peaks distance 1','unit'].astype(str).item(),
         'POT':inpt.loc[inpt['INPUT']=='POT 1','value'].astype(int).item(),
         'minimum':min_1,
         'maximum':max_1
        }

        driver2 = {
         'name':info.loc[info['info']=='Driver2','value'].astype(str).item(),
         'unit':info.loc[info['info']=='Driver2','unit'].astype(str).item(),
         'threshold_value':info.loc[info['info']=='Threshold Driver 2','value'].astype(float).item(),
         'threshold_value_unit':info.loc[info['info']=='Threshold Driver 2','unit'].astype(str).item(),
         'threshold_percent':inpt.loc[inpt['INPUT']=='threshold 2','value'].astype(float).item(),
         'threshold_percent_unit':inpt.loc[inpt['INPUT']=='threshold 2','unit'].astype(str).item(),
         'min_peaks_distance':inpt.loc[inpt['INPUT']=='min peaks distance 2','value'].astype(int).item(),
         'min_peaks_distance_unit':inpt.loc[inpt['INPUT']=='min peaks distance 2','unit'].astype(str).item(),
         'POT':inpt.loc[inpt['INPUT']=='POT 2','value'].astype(int).item(),
         'minimum':min_2,
         'maximum':max_2
        }
        
        potsinfo = {
        'correlation_coefficient' : info.loc[info['info']=='Correlation coefficient','value'].astype(float).item(),
        'search_window_value' : info.loc[info['info']=='Search Window','value'].astype(int).item(),
        'search_window_units' : info.loc[info['info']=='Search Window','unit'].astype(str).item(),
        'number_of_joint_events' : info.loc[info['info']=='n° of joint events','value'].astype(int).item(),
        'dt_value' : inpt.loc[inpt['INPUT']=='dt','value'].astype(str).item(),
        'dt_units' : inpt.loc[inpt['INPUT']=='dt','unit'].astype(str).item()
        }
        
        driver1 = Driver(**driver1)

        driver2 = Driver(**driver2)
        
        potsinfo = PotsInfo(**potsinfo)

        self.data = data
        self.driver1 = driver1
        self.driver2 = driver2
        self.potsinfo = potsinfo
        self.num_chains=num_chains
        self.num_samples=num_samples
        self.year_contours = year_contours
        self.marginal_years = marginal_years
        
       
    @classmethod
    def from_dict(cls, dct,num_chains=4, num_samples=2000, 
                 year_contours = np.array([10,50,100,200,500]),
                 marginal_years = np.geomspace(0.1,500,100)):
        """returns a JRP object built from json
        """
        data = pd.DataFrame(dct['data'])

        data['Time'] = data['Time'].apply(lambda r: datetime.datetime.fromtimestamp(r/1e3))

        info = pd.DataFrame(dct['info'])

        inpt = pd.DataFrame(dct['inpt'])

        jrp = cls(data, info, inpt, num_chains=num_chains, num_samples=num_samples, year_contours=year_contours, marginal_years=marginal_years)
        
        return(jrp)



    @classmethod
    def from_detectExtremeEvents(cls, file,num_chains=4, num_samples=2000, 
                 year_contours = np.array([10,50,100,200,500]),
                 marginal_years = np.geomspace(0.1,500,100)):
        """Returns a JRP object built from the output
        of Eleonora Perugini's detectExtremeEvents 
        Matlab code
        """
        
        data = pd.read_excel(file,sheet_name='data')
        
        info = pd.read_excel(file,sheet_name='info')
        
        inpt = pd.read_excel(file,sheet_name='input')
        
        jrp = cls(data, info, inpt, num_chains=num_chains, num_samples=num_samples, year_contours=year_contours, marginal_years=marginal_years)
        
        return(jrp)
    
    def format_data_for_stan(self,
                            marginal_years = np.geomspace(0.1,500,100),
                            #driver1_RP_for_JRPS =np.array([20,20,200,200])*365, 
                            #driver2_RP_for_JRPS = np.array([5,100,5,100])*365, 
                            n_grid=20,
                            #P = 1000, 
                            #L=20
                            u_rp_d1_years=np.array([20,200]),
                            u_rp_d2_years=np.array([20,200]),
                            u_val_d1 = np.array([300,400]),
                            u_val_d2 = np.array([60,80]),
                            grid_d1_min = None,
                            grid_d1_max = None,
                            grid_d2_min = None,
                            grid_d2_max = None
                            ):
        """Builds the dict containing all the information Stan needs to fit 
        the JRP distribution and calculate all the necessary values
        for plotting
        """

        
        # initialise the stan data dict with the event values
        data = {
            'events':self.data[['Driver1','Driver2']].values
        }
        
        # the number of events
        data['N_evnt'] = len(data['events'])
        
        # thresholds - these are the mu parameters of the fit GPD
        data['mu1'] = self.driver1.threshold_value
        data['mu2'] = self.driver2.threshold_value
        
        # max values  - for setting boundaries on the xi parameter
        data['max1'] = self.driver1.maximum
        data['max2'] = self.driver2.maximum
        
        # mean inter event period in days
        lam = self.data['Time'].sort_values().diff().dropna().mean().total_seconds()/(60*60*24)
        data['lam'] = lam
        
        # get points for sampling JRP grid as a function of driver value
        # we want the minimum and maximum JRP grid values to span the range 
        # of the data as well as the desired marginal return periods
        if grid_d1_min is None:
            grid_d1_min = self.driver1.minimum

        if grid_d1_max is None:
            grid_d1_max = self.driver1.maximum

        if grid_d2_min is None:
            grid_d2_min = self.driver2.minimum

        if grid_d2_max is None:
            grid_d2_max = self.driver2.maximum


        _grid_1 = np.linspace(grid_d1_min, grid_d1_max, n_grid)
        _grid_2 = np.linspace(grid_d2_min, grid_d2_max, n_grid)


        grid_1, grid_2 = np.meshgrid(_grid_1, _grid_2)
        grid_1, grid_2 = grid_1.flatten(), grid_2.flatten()
        grid_1, grid_2 

        data['grid1'] = grid_1
        data['grid2'] = grid_2
        data['N_grid'] = n_grid*n_grid

        # get points for marginal return period plots (plotted next to the JRP grid)
        marginal_days = self.marginal_years*365
        marginal_returnprob = 1-lam/marginal_days
        data['mrp_prob'] = marginal_returnprob
        data['N_mrp'] = len(data['mrp_prob'])
        
        
        # get points for sampling JRP grid as a function of return period
        grid_rp_days = np.geomspace(lam, 500*365, n_grid)
    
        grid_rp_days_returnprob = 1-lam/grid_rp_days
        grid_1b, grid_2b = np.meshgrid(grid_rp_days_returnprob, grid_rp_days_returnprob)
        grid_1b, grid_2b = grid_1b.flatten(), grid_2b.flatten()
        data['grid1b'] = grid_1b
        data['grid2b'] = grid_2b
        data['N_gridb'] = len(data['grid2b'])

       
        
        
        # quantiles for Q-Q plot
        driver1_ecdf = ECDF(self.data['Driver1'])
        driver2_ecdf = ECDF(self.data['Driver2'])


        data['quantiles_1'] = driver1_ecdf(self.data['Driver1'])
        data['quantiles_2'] = driver2_ecdf(self.data['Driver2'])

        # user defined marginal return periods for generating driver values and JRP in years
        data['u_rp_d1'] = u_rp_d1_years*365
        data['u_rp_d2'] = u_rp_d2_years*365
        data['N_umrp'] = len(data['u_rp_d2']) 


        # user defined driver values for calculating MRPs and JRPs
        data['u_val_d1'] = u_val_d1
        data['u_val_d2'] = u_val_d2
        data['N_uval'] = len(data['u_val_d2'])
        
        
        self.stan_data = data
        
        print("-------------------------------------------------------------------------")        
        print("DATA SUMMARY")
        print("-------------------------------------------------------------------------")
        print("Driver 1:", self.driver1.name,
              ", threshold", self.driver1.threshold_percent, self.driver1.threshold_percent_unit,
              ",", self.driver1.threshold_value, self.driver1.threshold_value_unit,
             )
        print("Driver 2:", self.driver2.name,
              ", threshold", self.driver2.threshold_percent, self.driver2.threshold_percent_unit,
              ",", self.driver2.threshold_value, self.driver2.threshold_value_unit,
             )        
        print(len(self.data),"extreme events selected using:\n",
              "(1) A search window of",
             self.potsinfo.search_window_value,self.potsinfo.search_window_units,
              "\n","(2) minimum peaks distance of", 
              self.driver1.min_peaks_distance,self.driver1.min_peaks_distance_unit,
              "for driver 1\n","(3) minimum peaks distance of",
              self.driver2.min_peaks_distance,self.driver2.min_peaks_distance_unit,
              "for driver 2."
             )
        print("Giving a correlation coefficient of ", self.potsinfo.correlation_coefficient)
        
        print("-------------------------------------------------------------------------")
        print("STAN CONFIGURATION")
        print("-------------------------------------------------------------------------")
        print("Stan will be run with",self.num_chains, "chains for",self.num_samples,"samples each.")
        print("All other settings set to defaults")
        print("-------------------------------------------------------------------------")
        print("JRP CONTOUR PLOT CONFIGURATION")
        print("-------------------------------------------------------------------------")
        print("JRP grid for the contour plot will be calculated at the following points")
        print(n_grid, "points between",grid_d1_min,self.driver1.unit, 
              "and",grid_d1_max,self.driver1.unit,  "for driver 1")
        print(n_grid, "points between", grid_d2_min,self.driver2.unit, 
              "and", grid_d2_max,self.driver2.unit,  "for driver 2")
        print("giving a grid with a total of",n_grid*n_grid, "points, contoured")
        print("for the following years:",self.year_contours)
        print("The marginal return period will be calculated at",len(marginal_years),
              "points between ",min(marginal_years),"and",max(marginal_years),"years.")
        print("-------------------------------------------------------------------------")
        print("USER DEFINED DRIVER RETURN PERIODS")
        print("-------------------------------------------------------------------------")
        print("Distributions over driver values and JRP will be calculated for user defined")
        print("return periods of",u_rp_d1_years,"years for",self.driver1.name,"and")
        print(u_rp_d2_years, "years for",self.driver2.name)
        print("-------------------------------------------------------------------------")
        print("USER DEFINED DRIVER VALUES")
        print("-------------------------------------------------------------------------")
        print("Distributions over driver marginal return periods and JRP will be calculated for user defined")
        print("values of",u_val_d1, self.driver1.unit,"for",self.driver1.name,"and")
        print(u_val_d2, self.driver2.unit, "for",self.driver2.name)
        print("-------------------------------------------------------------------------")

        return()
 
    def build_stan_model(self, file_model="app/model.stan"):

        print("Building Stan model ...")
        with open(file_model) as f:
            model = f.read()
    

        self.model = stan.build(model, data=self.stan_data, random_seed=1)   
        

    def run_stan_model(self):
        print("Running Stan model ...")

        self.result = self.model.sample(
                            num_chains=self.num_chains, 
                            num_samples=self.num_samples
                            )

        #self.summary = az.summary(self.result)
        
        
    def generate_trace_plot(self,filename=None):
        print("Generating trace plot ...")
        az.plot_trace(self.result,var_names=['sigma1','xi1','sigma2','xi2','theta'])
        plt.tight_layout()
        if filename is not None:
            plt.savefig(filename)
        return()
    
    
    def generate_QQ_plot(self,filename=None):
        posterior = az.from_pystan(self.result).posterior.stack({'sample':['chain','draw']})

        fig, axs = plt.subplots(1,2,figsize=(10,5))

        lims1 = [self.driver1.minimum, self.driver1.maximum]

        width1 = (lims1[1]-lims1[0])/20

        data1 =  posterior['d1_for_quantile'].values

        for i, absissa in enumerate(self.data['Driver1'].values):

            d1 = data1[i]
            absissa = np.array([absissa])

            d1 = d1[~np.isnan(d1)]


            violin_parts = axs[0].violinplot(
                d1,
                absissa,
                #positions = self.data['Driver1'].values,
                #manage_ticks = False,
                widths=width1,
                );
        
            for pc in violin_parts['bodies']:
                pc.set_facecolor('#1f77b4')
                pc.set_edgecolor('#1f77b4')

            
            for partname in ('cbars', 'cmins', 'cmaxes'):#, 'cmeans', 'cmedians'):
                vp = violin_parts[partname]
                vp.set_edgecolor('#1f77b4')
                #vp.set_linewidth(1)


        axs[0].plot(lims1, lims1, c=u'#ff7f0e')
        axs[0].set_xlim(lims1)
        axs[0].set_ylim(lims1)
        axs[0].set_title("Q-Q plot for marginal distribution of\n"+ self.driver1.name+" "+ self.driver1.unit)
        axs[0].set_xlabel("Observed "+self.driver1.name+" "+self.driver1.unit)
        axs[0].set_ylabel("Equivalent quantile from distribution\n"+self.driver1.name+" "+self.driver1.unit)  


        lims2 = [self.driver2.minimum, self.driver2.maximum]

        width2 = (lims2[1]-lims2[0])/20

        data2 =  posterior['d2_for_quantile'].values

        for i, absissa in enumerate(self.data['Driver2'].values):

            d2 = data2[i]
            absissa = np.array([absissa])

            d2 = d2[~np.isnan(d2)]

            violin_parts = axs[1].violinplot(
                d2,
                absissa,
                #positions = self.data['Driver1'].values,
                #manage_ticks = False,
                widths=width2,
            );


            for pc in violin_parts['bodies']:
                pc.set_facecolor('#1f77b4')
                pc.set_edgecolor('#1f77b4')

            
            for partname in ('cbars', 'cmins', 'cmaxes'):#, 'cmeans', 'cmedians'):
                vp = violin_parts[partname]
                vp.set_edgecolor('#1f77b4')
                #vp.set_linewidth(1)



        axs[1].set_xlim(lims2)
        axs[1].set_ylim(lims2)

        axs[1].plot(lims2, lims2, c=u'#ff7f0e')



        axs[1].set_title("Q-Q plot for marginal distribution of\n"+ self.driver2.name+" "+ self.driver1.unit)
        axs[1].set_xlabel("Observed "+self.driver2.name+" "+self.driver2.unit)
        axs[1].set_ylabel("Equivalent quantile from distribution\n"+self.driver2.name+" "+self.driver2.unit)        

        plt.tight_layout()

        if filename is not None:
            plt.savefig(filename)
        
        
    def generate_table(self):
        
        # get the posterior as an xarray dataset and stack the chains
        posterior = az.from_pystan(self.result).posterior.stack({'sample':['chain','draw']})
        
        # get the return periods and joint return periods, convert to years
        jrp_for_events_years = posterior['jrp_for_events'].values/365
        rp1_for_events_years = posterior['rp1_for_events'].values/365
        rp2_for_events_years = posterior['rp2_for_events'].values/365


        # the statistics we want to calculate
        # new stats can be appended here
        funcs = {
            "mean":lambda x: np.nanmean(x,axis=1), 
            "median":lambda x: np.nanmedian(x,axis=1), 
            "std":lambda x: np.nanstd(x, axis=1),
            "pc10":lambda x: np.nanpercentile(x,10, axis=1),
            "pc90":lambda x: np.nanpercentile(x,90, axis=1)
        }


        # get JRP stats
        df_jrp= pd.DataFrame({name:func(jrp_for_events_years) for name, func in funcs.items()})
        df_jrp.columns = pd.MultiIndex.from_tuples([
            ("Joint Return Period (JRP, years)",x) for x in  df_jrp.columns ])

        # get RP stats for driver 1
        df_rp1= pd.DataFrame({name:func(rp1_for_events_years) for name, func in funcs.items()})
        df_rp1.columns = pd.MultiIndex.from_tuples([
            ("Marginal Return period for "+self.driver1.name +" (RP, years)",x) for x in  df_rp1.columns ])

        # get RP stats for driver 2
        df_rp2= pd.DataFrame({name:func(rp2_for_events_years) for name, func in funcs.items()})
        df_rp2.columns = pd.MultiIndex.from_tuples([
            ("Marginal Return period for "+self.driver2.name +" (RP, years)",x) for x in  df_rp2.columns ])

        # get the original data, fix the column names
        data = self.data.copy()

        data = data.rename(
            columns={
                'Driver1':self.driver1.name + " "+ self.driver1.unit,
                'Driver2':self.driver2.name + " "+ self.driver2.unit})

        data.columns = pd.MultiIndex.from_tuples([("data",x) for x in data.columns])

        # combine all into a new dataset and return
        df = pd.concat([data,df_jrp, df_rp1,df_rp2],axis=1)
        
        return(df)
    
    def generate_ppc_plot(self, filename=None):
        idata = az.from_pystan(
            posterior=self.result, 
            posterior_predictive=["events_hat"], 
            observed_data=["events"],
            posterior_model= self.model
        )

        fig, ax = plt.subplots(1,1)
        az.plot_ppc(idata, data_pairs={"events": "events_hat"}, ax=ax,mean=False)

        lims= [min(self.driver1.minimum,self.driver2.minimum),
               max(self.driver1.maximum,self.driver2.maximum)]

        lims = [lims[0]*0.5,lims[1]*1.5]

        ax.set_xlim(lims)

        ax.set_xlabel("Raw data for both drivers ( "+self.driver1.name+", "+ self.driver2.name+" )")
        ax.set_ylabel("Probability density")
        
        if filename is not None:
            plt.savefig(filename)

    def generate_jrp_plot(self,filename=None,driver1_coords=None,
                                driver2_coords=None,jrp_coords=None):
        print("Generating JRP plot ...")

        

        alpha =0.5

        temp = (az
            .from_pystan(self.result)
            .posterior
            .stack(dimensions={'sample':['chain','draw']})
        )



        lam = self.stan_data['lam']


        fig, axs = plt.subplots(2,2,figsize=(10,10))

        self.data.plot.scatter(
                x='Driver1',
                y='Driver2',
                ax=axs[0,1],
                color='lightgrey',
                marker='.'
            )

        

        day_contours =  self.year_contours*365

        print("NOTE: np.nanperctile used to calculate jrp_50")
        jrp_50 = np.nanpercentile(temp['jrp_grid'].values, 50,axis=1)

        # We plot every_n_samples
        N_samples =20
        every_n_samples = int(8000/N_samples)
        for jrp_i in temp['jrp_grid'].values[:,::every_n_samples].T:
                axs[0,1].tricontour(
                    self.stan_data['grid1'],
                    self.stan_data['grid2'],
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
                    self.stan_data['grid1'],
                    self.stan_data['grid2'],
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
            fmt={ key:item for key, item in zip(day_contours, self.year_contours)},
            inline=True, 
            fontsize=10
        )
        


        
        ################################### Driver1 #########################################

        axs[1,1].fill_betweenx(
            self.marginal_years,
            np.nanpercentile(temp['mrp_1'].values,20,axis=1),
            np.nanpercentile(temp['mrp_1'].values,80,axis=1),
            alpha=0.1      
        )
        axs[1,1].plot(
            np.nanpercentile(temp['mrp_1'].values,50,axis=1),
            self.marginal_years
        )
        axs[1,1].plot(
            np.nanpercentile(temp['mrp_1'].values,20,axis=1),
            self.marginal_years,
            linestyle='dotted',
             c='#1f77b4'
        )
        axs[1,1].plot(
            np.nanpercentile(temp['mrp_1'].values,80,axis=1),
            self.marginal_years,
            linestyle='dotted',
             c='#1f77b4'

        )

        axs[1,1].grid(True, which="both", ls="-", color='0.65')
        axs[1,1].set_yscale('log')
        axs[1,1].yaxis.set_major_formatter(ScalarFormatter())
        
        axs[1,1].set_xlim([
            self.driver1.minimum,
            self.driver1.maximum
        ])

        
        ################################### Driver 2 #########################################

        axs[0,0].fill_between(
            self.marginal_years,
            np.nanpercentile(temp['mrp_2'].values,20,axis=1),
            np.nanpercentile(temp['mrp_2'].values,80,axis=1),
            alpha=0.1      
        )
        axs[0,0].plot(
            self.marginal_years,
            np.nanpercentile(temp['mrp_2'].values,50,axis=1),

        )
        axs[0,0].plot(
            self.marginal_years,
            np.nanpercentile(temp['mrp_2'].values,20,axis=1),
            linestyle='dotted',
             c='#1f77b4'
        )
        axs[0,0].plot(
            self.marginal_years,
            np.nanpercentile(temp['mrp_2'].values,80,axis=1),
            linestyle='dotted',
             c='#1f77b4'

        )

        axs[0,0].grid(True, which="both", ls="-", color='0.65')
        axs[0,0].set_xscale('log')
        axs[0,0].xaxis.set_major_formatter(ScalarFormatter())

        axs[0,0].set_ylim([
            self.driver2.minimum,
            self.driver2.maximum
        ])
        
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

        axs[0,1].set_xlabel('')
        axs[0,1].set_ylabel('')

    

        axs[1,1].set_xlabel(self.driver1.name+" "+self.driver1.unit)
        axs[1,1].set_ylabel("Return period yr")


        axs[0,0].set_ylabel(self.driver2.name+" "+self.driver2.unit)
        axs[0,0].set_xlabel("Return period yr")


        axs[0,1].set_xlabel(self.driver1.name+" "+self.driver1.unit)
        axs[0,1].set_ylabel(self.driver2.name+" "+self.driver2.unit)
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
        
        
            

        plt.subplots_adjust(hspace=0.02,wspace=0.02)
        
        
        
        if driver1_coords is not None:
            
            for d1, d2, jrp in zip(driver1_coords,
                                  driver2_coords,
                                  jrp_coords):

                con = ConnectionPatch(
                                xyA=d2, xyB=jrp, 
                                coordsA="data", coordsB="data",
                                axesA=axs[0,0], axesB=axs[0,1],
                                color="red",linestyle='dashed')
                fig.add_artist(con)


                con = ConnectionPatch(
                                xyB=d1, xyA=jrp, 
                                coordsB="data", coordsA="data",
                                axesB=axs[1,1], axesA=axs[0,1],
                                color="red",linestyle='dashed')

                fig.add_artist(con)


        
        if filename is not None:
            plt.savefig(filename)


        return(fig)
    
    def get_summary(self):
        return(az.summary(self.result,skipna=True))
    
    
    def generate_corner_plot(self,filename=None):
        print("Generating corner plot ...")
        az.plot_pair(self.result,
            var_names = ['sigma1','sigma2','xi1','xi2','theta'],
            marginals=True,  kind='kde', divergences=True,figsize=(10,10))
        if filename is not None:
            plt.savefig(filename)
            



    def generate_jrp_plot2(self,filename=None,driver1_coords=None,
                                driver2_coords=None,jrp_coords=None,
                                PLOT_ALL_DATA = False,  
                                PLOT_USER_D1D2 = False,  
                                PLOT_USER_MRPS = False):
                

        temp = (az
            .from_pystan(self.result)
            .posterior
            .stack(dimensions={'sample':['chain','draw']})
        )


        posterior = az.from_pystan(self.result).posterior.stack({'sample':['chain','draw']})

        lam = self.stan_data['lam']


        fig, axs = plt.subplots(2,2,figsize=(10,10))





        day_contours =  self.year_contours*365

        jrp_50 = np.nanpercentile(temp['jrp_grid'].values, 50,axis=1)

        # We plot every_n_samples
        N_samples =20
        every_n_samples = int(8000/N_samples)
        for jrp_i in temp['jrp_grid'].values[:,::every_n_samples].T:
                axs[0,1].tricontour(
                    self.stan_data['grid1'],
                    self.stan_data['grid2'],
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
                    self.stan_data['grid1'],
                    self.stan_data['grid2'],
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
            fmt={ key:item for key, item in zip(day_contours, self.year_contours)},
            inline=True, 
            fontsize=10
        )


        ################################### JRP contours function of RP #####################

        jrp_50b = np.nanpercentile(temp['jrp_gridb'].values, 50,axis=1)

        # We plot every_n_samples
        N_samples =20
        every_n_samples = int(8000/N_samples)
        for jrp_i in temp['jrp_gridb'].values[:,::every_n_samples].T:
                axs[1,0].tricontour(
                    (1/365)*lam/(1-self.stan_data['grid1b']),
                    (1/365)*lam/(1-self.stan_data['grid2b']),
                    jrp_i,
                    levels=day_contours,
                    cmap='hsv',
                    alpha=0.1,
                    norm=colors.LogNorm(
                                vmin=jrp_50.min(), 
                                vmax=jrp_50.max()
                                        )  
        )
                
                
        
        CSb = axs[1,0].tricontour(
                    (1/365)*lam/(1-self.stan_data['grid1b']),
                    (1/365)*lam/(1-self.stan_data['grid2b']),
                    jrp_50b,
                    levels=day_contours,
                    alpha=1,
                    cmap='hsv',
                    norm=colors.LogNorm(
                                vmin=jrp_50.min(), 
                                vmax=jrp_50.max()
                                        )
                    )    

        axs[1,0].clabel(
            CSb, 
            CSb.levels,
            fmt={ key:item for key, item in zip(day_contours, self.year_contours)},
            inline=False, 
            fontsize=10,
            inline_spacing=1
        )





        ################################### Driver1 #########################################

        axs[1,1].fill_betweenx(
            self.marginal_years,
            np.nanpercentile(temp['mrp_1'].values,20,axis=1),
            np.nanpercentile(temp['mrp_1'].values,80,axis=1),
            alpha=0.1      
        )
        axs[1,1].plot(
            np.nanpercentile(temp['mrp_1'].values,50,axis=1),
            self.marginal_years
        )
        axs[1,1].plot(
            np.nanpercentile(temp['mrp_1'].values,20,axis=1),
            self.marginal_years,
            linestyle='dotted',
            c='#1f77b4'
        )
        axs[1,1].plot(
            np.nanpercentile(temp['mrp_1'].values,80,axis=1),
            self.marginal_years,
            linestyle='dotted',
            c='#1f77b4'

        )



        ################################### Driver 2 #########################################

        axs[0,0].fill_between(
            self.marginal_years,
            np.nanpercentile(temp['mrp_2'].values,20,axis=1),
            np.nanpercentile(temp['mrp_2'].values,80,axis=1),
            alpha=0.1      
        )
        axs[0,0].plot(
            self.marginal_years,
            np.nanpercentile(temp['mrp_2'].values,50,axis=1),

        )
        axs[0,0].plot(
            self.marginal_years,
            np.nanpercentile(temp['mrp_2'].values,20,axis=1),
            linestyle='dotted',
            c='#1f77b4'
        )
        axs[0,0].plot(
            self.marginal_years,
            np.nanpercentile(temp['mrp_2'].values,80,axis=1),
            linestyle='dotted',
            c='#1f77b4'

        )





        ###### Plot all data ###############################################

        if PLOT_ALL_DATA:
            
            # Plot extreme events in JRP as function of values
            self.data.plot.scatter(
                x='Driver1',
                y='Driver2',
                ax=axs[0,1],
                #color='lightgrey',
                #marker='.'
                color='#1f77b4',
                marker='o'
            )

            # plot the median and distribution over MRP space also
            df=self.generate_table()
            axs[1,0].scatter(
                df['Marginal Return period for r24h (RP, years)']['median'],
                df['Marginal Return period for API (RP, years)']['median'],
                marker='.'
            )
            
            # plot distributions 
            for i in range(self.stan_data['N_evnt']):

                az.plot_kde(
                    posterior['rp2_for_events'].isel(rp2_for_events_dim_0=i).values/365,
                    posterior['rp1_for_events'].isel(rp1_for_events_dim_0=i).values/365,
                    ax=axs[1,0],
                    hdi_probs=[0.5],
                    contourf_kwargs={"colors": ['#1f77b4'],'alpha':0.1},
                    contour_kwargs={"colors": ['#1f77b4'], "linestyles":["dotted"]}
                )


            
        
                

        ###### Axes tick params ###############################################

        axs[0,0].xaxis.set_tick_params(
            labeltop=True,
            top=True,
            labelbottom=False,
            bottom=False
            )

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

        axs[1,1].yaxis.set_tick_params(
            labelleft=False,
            labelright=True,
            left=False, 
            right=True 
        )

        ###### Axes tick label formatters ###############################################

        axs[0,0].xaxis.set_major_formatter(ScalarFormatter())
        axs[1,0].xaxis.set_major_formatter(ScalarFormatter())
        axs[0,1].yaxis.set_major_formatter(ScalarFormatter())
        axs[1,1].yaxis.set_major_formatter(ScalarFormatter())



        ###### Axes label positions ###############################################
        axs[0,0].xaxis.set_label_position("top")
        axs[0,1].xaxis.set_label_position("top")
        axs[0,1].yaxis.set_label_position("right")
        axs[1,1].yaxis.set_label_position("right")
        axs[1,1].xaxis.set_label_position("bottom")



        ###### All axes scales ###############################################

        axs[0,0].set_xscale('log')
        axs[1,0].set_xscale('log')
        axs[1,0].set_yscale('log')
        axs[1,1].set_yscale('log')


        ###### All axes labels ###############################################

        axs[0,0].set_ylabel(self.driver2.name + " values ("+self.driver2.unit+")")
        axs[0,0].set_xlabel(self.driver2.name + " marginal return priod (yr)")
        axs[0,1].set_xlabel(self.driver1.name + " values ("+self.driver1.unit+")")
        axs[0,1].set_ylabel(self.driver2.name + " values ("+self.driver2.unit+")")
        axs[1,0].set_xlabel(self.driver2.name + " marginal return priod (yr)")
        axs[1,0].set_ylabel(self.driver1.name + " marginal return priod (yr)")
        axs[1,1].set_xlabel(self.driver1.name + " values ("+self.driver1.unit+")")
        axs[1,1].set_ylabel(self.driver1.name + " marginal return priod (yr)")


        ###### Limits for all axes ###############################################

        axs[0,0].set_ylim([
            self.driver2.minimum,
            self.driver2.maximum
        ])

        axs[0,1].set_ylim([
            self.driver2.minimum,
            self.driver2.maximum
        ])

        axs[0,1].set_xlim([
            self.driver1.minimum,
            self.driver1.maximum
        ])

        xlim = axs[0,0].get_xlim()
        ylim = axs[1,1].get_ylim()

        axs[0,0].set_xlim((lam/365,100))
        axs[1,0].set_xlim((lam/365,100))
        axs[1,0].set_ylim((lam/365,100))
        axs[1,1].set_ylim((lam/365,100))
        axs[1,1].set_xlim([
            self.driver1.minimum,
            self.driver1.maximum
        ])




        ###### Plot distribution over (mrp1, mrp2) for user specified (d1,d2)  ##
        if PLOT_USER_MRPS:

            for i in range(self.stan_data['N_uval']):
                
                # get the coordinates we need
                u_r1 = posterior['u_rp_for_d1'].isel(u_rp_for_d1_dim_0=i).median().item()/365
                u_r2 = posterior['u_rp_for_d2'].isel(u_rp_for_d2_dim_0=i).median().item()/365
                u_d1 = self.stan_data['u_val_d1'][i]
                u_d2 = self.stan_data['u_val_d2'][i]
                
                
                # plot location on (mrp1, mrp2)
                axs[0,1].scatter(u_d1, u_d2, color='#1f77b4')
                
                con = ConnectionPatch(
                                xyA=[u_r2,u_r1], xyB=[u_r2,u_d2], 
                                coordsA="data", coordsB="data",
                                axesA=axs[1,0], axesB=axs[0,0],
                                color='#1f77b4',linestyle='dotted')
                
                fig.add_artist(con)
                
                
                
                con = ConnectionPatch(
                                xyA=[u_r2,u_d2], xyB=[u_d1,u_d2], 
                                coordsA="data", coordsB="data",
                                axesA=axs[0,0], axesB=axs[0,1],
                                color='#1f77b4',linestyle='dotted')
                
                fig.add_artist(con)
                
                
                con = ConnectionPatch(
                                xyA=[u_r2,u_r1], xyB=[u_d1,u_r1], 
                                coordsA="data", coordsB="data",
                                axesA=axs[1,0], axesB=axs[1,1],
                                color='#1f77b4',linestyle='dotted')
                
                fig.add_artist(con)
                
                
                con = ConnectionPatch(
                                xyA=[u_d1,u_r1], xyB=[u_d1,u_d2], 
                                coordsA="data", coordsB="data",
                                axesA=axs[1,1], axesB=axs[0,1],
                                color='#1f77b4',linestyle='dotted')
                
                fig.add_artist(con)
                
                
                # plot the distribution over (d1, d2)
                az.plot_kde(
                        posterior['u_rp_for_d2'].isel(u_rp_for_d2_dim_0=i).values/365,
                        posterior['u_rp_for_d1'].isel(u_rp_for_d1_dim_0=i).values/365,
                        ax=axs[1,0],
                        hdi_probs=[0.5],
                        contourf_kwargs={"colors": ['#1f77b4'],'alpha':0.1},
                        contour_kwargs={"colors": ['#1f77b4'], "linestyles":["dotted"]}
                    )
                
                
                u_d1_ax, u_d2_ax  = (axs[0,1].transData + axs[0,1].transAxes.inverted()).transform([u_d1,u_d2])
                
                ax_inset = axs[0,1].inset_axes([u_d1_ax, u_d2_ax,0.2,0.2]) #transform=axs[1,0].transData
                
                (posterior['u_JRP_for_vals'].isel(u_JRP_for_vals_dim_0=i)/365).plot.hist(ax=ax_inset,alpha=0.5)

                med = (posterior['u_JRP_for_vals'].isel(u_JRP_for_vals_dim_0=i)/365).median()
                ymax = ax_inset.get_ylim()[1]
                ax_inset.plot(
                            [med,med],
                            [0,ymax],
                    color='#1f77b4'
                )
                ax_inset.text(med, ymax, str(round(med.item(),3)),horizontalalignment='center',color='#1f77b4')
                
                
                ax_inset.xaxis.set_major_locator(plt.MaxNLocator(2))
                ax_inset.spines['top'].set_visible(False)
                ax_inset.spines['right'].set_visible(False)
                ax_inset.spines['left'].set_visible(False)
                ax_inset.yaxis.set_tick_params(labelleft=False)
                ax_inset.set_yticks([])
                ax_inset.set_title("");
                ax_inset.set_xlabel("JRP years")
                ax_inset.set_facecolor = None
                ax_inset.patch.set_alpha(0)
                ax_inset.patch.set_alpha(0)
                ax_inset.patch.set_facecolor('#909090')
                ax_inset.spines['bottom'].set_color('#1f77b4')
                ax_inset.tick_params(axis='x', colors='#1f77b4')
                ax_inset.xaxis.label.set_color('#1f77b4')
                

        ###### Plot distribution over (d1,d2) for user specified (mrp1, mrp2) ##

        if PLOT_USER_D1D2:

            for i in range(self.stan_data['N_umrp']):
                
                # get the coordinates we need
                u_r1 = self.stan_data['u_rp_d1'][i]/365
                u_r2 = self.stan_data['u_rp_d2'][i]/365
                u_d1 = posterior['u_d1_for_rp'].isel(u_d1_for_rp_dim_0=i).median().item()
                u_d2 = posterior['u_d2_for_rp'].isel(u_d2_for_rp_dim_0=i).median().item()
            
            
                # plot location on (mrp1, mrp2)
                axs[1,0].scatter(u_r2, u_r1, color='#1f77b4')
                
                con = ConnectionPatch(
                                xyA=[u_r2,u_r1], xyB=[u_r2,u_d2], 
                                coordsA="data", coordsB="data",
                                axesA=axs[1,0], axesB=axs[0,0],
                                color='#1f77b4',linestyle='dotted')
                
                fig.add_artist(con)
                
                
                
                con = ConnectionPatch(
                                xyA=[u_r2,u_d2], xyB=[u_d1,u_d2], 
                                coordsA="data", coordsB="data",
                                axesA=axs[0,0], axesB=axs[0,1],
                                color='#1f77b4',linestyle='dotted')
                
                fig.add_artist(con)
                
                
                con = ConnectionPatch(
                                xyA=[u_r2,u_r1], xyB=[u_d1,u_r1], 
                                coordsA="data", coordsB="data",
                                axesA=axs[1,0], axesB=axs[1,1],
                                color='#1f77b4',linestyle='dotted')
                
                fig.add_artist(con)
                
                
                con = ConnectionPatch(
                                xyA=[u_d1,u_r1], xyB=[u_d1,u_d2], 
                                coordsA="data", coordsB="data",
                                axesA=axs[1,1], axesB=axs[0,1],
                                color='#1f77b4',linestyle='dotted')
                
                fig.add_artist(con)
                
            
                # plot the distribution over (d1, d2)
                az.plot_kde(
                        posterior['u_d1_for_rp'].isel(u_d1_for_rp_dim_0=i).values,
                        posterior['u_d2_for_rp'].isel(u_d2_for_rp_dim_0=i).values,
                        ax=axs[0,1],
                        hdi_probs=[0.5],
                        contourf_kwargs={"colors": ['#1f77b4'],'alpha':0.1},
                        contour_kwargs={"colors": ['#1f77b4'], "linestyles":["dotted"]}
                    )
                
                
                u_r2_ax, u_r1_ax  = (axs[1,0].transData + axs[1,0].transAxes.inverted()).transform([u_r2,u_r1])
                #u_r2_ax, u_r1_ax = axs[1,0].transData.inverted().transform([u_r2,u_r1])
                
                ax_inset = axs[1,0].inset_axes([u_r2_ax, u_r1_ax,0.2,0.2]) #transform=axs[1,0].transData
                
                (posterior['u_JRP'].isel(u_JRP_dim_0=i)/365).plot.hist(ax=ax_inset,alpha=0.5)

                med = (posterior['u_JRP'].isel(u_JRP_dim_0=i)/365).median()
                ymax = ax_inset.get_ylim()[1]
                ax_inset.plot(
                            [med,med],
                            [0,ymax],
                    color='#1f77b4'
                )
                ax_inset.text(med, ymax, str(round(med.item(),3)),horizontalalignment='center',color='#1f77b4')
                
                
                ax_inset.xaxis.set_major_locator(plt.MaxNLocator(2))
                ax_inset.spines['top'].set_visible(False)
                ax_inset.spines['right'].set_visible(False)
                ax_inset.spines['left'].set_visible(False)
                ax_inset.yaxis.set_tick_params(labelleft=False)
                ax_inset.set_yticks([])
                ax_inset.set_title("");
                ax_inset.set_xlabel("JRP years")
                ax_inset.set_facecolor = None
                ax_inset.patch.set_alpha(0)
                ax_inset.patch.set_alpha(0)
                ax_inset.patch.set_facecolor('#909090')
                ax_inset.spines['bottom'].set_color('#1f77b4')
                ax_inset.tick_params(axis='x', colors='#1f77b4')
                ax_inset.xaxis.label.set_color('#1f77b4')

                
        ###### Grid settings for all axes ########################################

        axs[0,0].grid(True, which="both", ls="-", color='0.9')
        axs[1,0].grid(True, which="both", ls="-", color='0.9')
        axs[1,1].grid(True, which="both", ls="-", color='0.9')
        axs[0,1].grid(True, which="both", ls="-", color='0.9')


        ###### Figure labels ########################################

            
        axs[0,0].text(0.05, 0.95, '(a)',
            horizontalalignment='center',
            verticalalignment='center',
                    fontsize=12,
            transform = axs[0,0].transAxes) 

        axs[0,1].text(0.05, 0.95, '(b)',
            horizontalalignment='center',
            verticalalignment='center',
                    fontsize=12,
            transform = axs[0,1].transAxes) 

        axs[1,0].text(0.05, 0.95, '(c)',
            horizontalalignment='center',
            verticalalignment='center',
                    fontsize=12,
            transform = axs[1,0].transAxes) 

        axs[1,1].text(0.05, 0.95, '(d)',
            horizontalalignment='center',
            verticalalignment='center',
                    fontsize=12,
            transform = axs[1,1].transAxes) 

        ###### Final adjustments ########################################


        plt.subplots_adjust(hspace=0.02,wspace=0.02)

        if filename is not None:
            plt.savefig(filename)
