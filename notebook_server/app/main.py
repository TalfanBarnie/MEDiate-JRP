from typing import Union
from fastapi import FastAPI, Request
from fastapi.staticfiles import StaticFiles
from . import jrp as JRP
import pandas as pd
import json

app = FastAPI()

# Finished figure .pngs are stored as static files and relative urls are sent to the user 
# in response to the POST method
app.mount("/static", StaticFiles(directory="/app/DATA"), name="static")


# The user sends the data to perform inference on as json via a POST method
@app.post("/")
async def get_body(request: Request):

    # we wait for the json to be decoded
    user_data =  await request.json()

    # json is then converted to csv and stored locally for use by the JRP module ..

    df_data = pd.DataFrame({
                    'date':user_data['data']['date'],
                    'surge':user_data['data']['surge'],
                    'flow':user_data['data']['flow']
                })
    df_data.to_csv("/app/DATA/temporary_data.csv")

    # ... the user data is stored as two separate .csvs, one conctaining the 
    # bivariate data, once containing threshold used for contemporray
    # JRPs, one containing thershold projected into the future
    # to account for climate change

    df_threshold = pd.DataFrame({
                    'year':user_data['threshold']['year'],
                    'surge':user_data['threshold']['surge'],
                    'flow':user_data['threshold']['flow']
                })
    df_threshold.to_csv("/app/DATA/temporary_threshold.csv")

    # initialise the JRP object with the user data
    jrp = JRP.JRP(file_data="/app/DATA/temporary_data.csv", file_thresholds="/app/DATA/temporary_threshold.csv")

    # perform inference on the user data using stan
    # this produces a summary pandas dataframe containing 
    # the statistics of the parameters and produces figures
    # as .pngs which are written to the static directory
    jrp.infer()

    summary = jrp.summary.to_json()

    summary = json.loads(summary)

    # we append the relative url of the figures
    # so the user can just append them to the 
    # base url and fetch them
    summary['corner_plot']='static/corner_plot.png'
    summary['trace_plot']='static/trace_plot.png'


    return json.dumps(summary)

