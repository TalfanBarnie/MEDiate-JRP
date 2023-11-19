# MEDiate-JRP
Package for generating Joint Return Probabilities for paired extreme events data. Can be run as an API from a container or in a regular Jupyter Notebook. The package consists of the following files:

    MEDiate-JRP/
    ├── app
    │   ├── __init__.py                                     - initialises folder as a python package
    │   ├── jrp.py                                          - preprocesses, runs model, make figures
    │   ├── main.py                                         - handles API
    │   └── model.stan                                      - specifies the statistical model
    ├── DATA                                                - folder for output, static files shared by API
    │   └── README.md
    ├── docker-compose.yaml                                 - orchestrates containers, ports, volumes 
    ├── Dockerfile                                          - specifies the container that runs the API
    ├── documentation
    │   └── Overview_of_JRP_estimation_and_Oslo_example.pdf
    ├── environment.yaml                                    - the required python packages
    ├── example_data.csv                                    - example data (Oslo) for testing notebook
    ├── example.json                                        - example data (Oslo) for testing API
    ├── example_thresholds.csv                              - example data (Oslo) for testing notebook
    ├── notebook.ipynb                                      - notebook for using jrp locally
    └── README.md




## To run as an API
Clone the repository to your machine

	git clone https://github.com/TalfanBarnie/MEDiate-JRP.git
Build the container, start the server, run in background

    sudo docker compose up --build -d
    
Test with the example json data file

	curl -X POST -H "Content-Type: application/json" -d @example.json http://127.0.0.1/
    
NOTE! This can take a few minutes to run. Returns JSON containing summary statistics and relative urls to figures.

	sudo docker compose down
    
Shuts down the container

## To run as a notebook
Set up a python environment with your preferred environment manager (conda, mamba, pyenv ...) and the packages in environment.yaml.
(Note Jupyter is not included in the environment.yaml to avoid cluttering up the docker container with
packages we don't need - to run the notebook you have to install that too). Then

	jupyter run notebook
    
In the window that opens select notebook.ipynb, proceed as usual. 
## Todo list
- unit tests for Stan code
- integration tests - synthetic data drawn from a known distribution
- add missing figures
- save compiled Stan model between runs
- extend api to give user more control (number of chains, warmup, sampling, etc.)
- make sure all necessary diagnostics are returned
