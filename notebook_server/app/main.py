from typing import Union
from fastapi import FastAPI, Request
from fastapi.staticfiles import StaticFiles
from . import jrp as JRP
import pandas as pd
import json

app = FastAPI(debug=True)

# Finished figure .pngs are stored as static files and relative urls are sent to the user 
# in response to the POST method
app.mount("/static", StaticFiles(directory="/app/DATA"), name="static")


# The user sends the data to perform inference on as json via a POST method
@app.post("/")
async def get_body(request: Request):

    # we wait for the json to be decoded
    user_data =  await request.json()

    #user_dict = json.loads(user_data)


    # initialise the JRP object with the user data
    jrp = JRP.JRP.from_dict(user_data)
    # perform inference on the user data using stan
    # this produces a summary pandas dataframe containing 
    # the statistics of the parameters and produces figures
    # as .pngs which are written to the static directory
    jrp.format_data_for_stan()
    jrp.build_stan_model()
    jrp.run_stan_model()
    jrp.generate_trace_plot("/app/DATA/trace_plot.png")
    jrp.generate_jrp_plot("/app/DATA/jrp_plot.png")
    jrp.generate_QQ_plot("/app/DATA/qq_plot.png")
    jrp.generate_ppc_plot("/app/DATA/ppc_plot.png")
    jrp.generate_corner_plot("/app/DATA/corner_plot.png")
    summary = jrp.get_summary()
    table = jrp.generate_table()

    out_data = {}
    out_data["summary"] = json.loads(summary.to_json())
    out_data["table"] = json.loads(table.to_json())
    
    # we append the relative url of the figures
    # so the user can just append them to the 
    # base url and fetch them
    out_data['corner_plot']='static/corner_plot.png'
    out_data['trace_plot']='static/trace_plot.png'
    out_data['ppc_plot']='static/ppc_plot.png'
    out_data['qq_plot']='static/qq_plot.png'
    out_data['jrp_plot']='static/jrp_plot.png'


#json.dumps({"hello":"hello"})
    return json.dumps(out_data)

