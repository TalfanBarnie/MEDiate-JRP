# MEDiate-JRP
Package for generating Joint Return Probabilities for paired extreme events data. Can be run as an API from a container or in a regular Jupyter Notebook. The package consists of the following files:


    MEDiate-JRP/
    ├── api_server
    │   ├── app
    │   │   ├── __init__.py                                     - initialises folder as a python package
    │   │   ├── jrp.py                                          - preprocesses, runs model, make figures
    │   │   ├── main.py                                         - handles API
    │   │   └── model.stan                                      - specifies the statistical model
    │   ├── DATA
    │   │   └── README.md
    │   ├── docker-compose.yaml                                 - orchestrates containers, ports, volumes
    │   ├── Dockerfile                                          - specifies the container that runs the API
    │   ├── environment.yaml                                    - the required python packages
    │   └── example.json                                        - example data  for testing API
    ├── documentation
    │   └── Overview_of_JRP_estimation_and_Oslo_example.pdf
    ├── notebook_server
    │   ├── app
    │   │   ├── __init__.py                                     - initialises folder as a python package
    │   │   ├── jrp.py                                          - preprocesses, runs model, make figures
    │   │   ├── main.py                                         - handles API
    │   │   └── model.stan                                      - specifies the statistical model
    │   ├── DATA
    │   │   ├── jextr_api9824h98pot5_daily_7dw.xlsx             - example output from detectExtremeEvents
    │   │   └── README.md
    │   ├── docker-compose.yaml                                 - orchestrates containers, ports, volumes
    │   ├── Dockerfile                                          - specifies the container that runs the notebook server
    │   ├── environment.yaml                                    - the required python packages
    │   └── notebook.ipynb                                      - notebook for using jrp locally
    └── README.md


## To download the repository
Clone the repository to your machine

	git clone https://github.com/TalfanBarnie/MEDiate-JRP.git

You will be asked for your username and password. For the password github now requires a token:
- on github, click yur user icon in the top right corner
- settings -> developer settings -> personal access tokens -> Tokens (classic) -> Generate new token -> Generate new token (classic) 
- enter you authentifcation code
- under note give it a name, click "repo" to select al check boxes, then click generate token
- copy result enter it as password


## To run the API server
Navigate into the api server directory

	cd MEDiate-JRP/api_server

Build the container, start the server, run in background

    sudo docker compose up --build -d
    
Test with the example json data file

	curl -X POST -H "Content-Type: application/json" -d @example.json http://127.0.0.1/
    
NOTE! This can take a few minutes to run. Returns JSON containing summary statistics and relative urls to figures.

	sudo docker compose down
    
Shuts down the container

## To run the notebook server
Navigate into the notebook server directory

	cd MEDiate-JRP/notebook_server

Build the container, start the server, run in background

	sudo docker compose up --build -d

Open a browser go to 

	http://localhost/notebooks/notebook.ipynb

double click the notebook.ipynb to open. Run the contents of the notebook, figures are saved in shared volumed DATA on the host.

	sudo docker compose down
    
Shuts down the container

## Todo list
- app folder is the same in both directories, must be a way to have it in build context for both Dockerfiles 
- unit tests for Stan code
- integration tests - synthetic data drawn from a known distribution
- add missing figures
- save compiled Stan model between runs
- extend api to give user more control (number of chains, warmup, sampling, etc.)
- make sure all necessary diagnostics are returned
