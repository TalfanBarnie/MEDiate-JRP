FROM continuumio/miniconda3
RUN apt-get --allow-releaseinfo-change update
RUN apt-get install -y libgl1-mesa-glx ffmpeg libsm6 libxext6 libarchive-dev
COPY environment.yaml .
RUN conda install -c conda-forge mamba
RUN /bin/bash -c "mamba env create -f environment.yaml"
RUN ["conda", "run", "-n", "myenv", "python", "-m", "pip", "install", "pystan"] 
COPY ./app ./app
ENTRYPOINT ["conda", "run", "-n", "myenv","uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "80"]

