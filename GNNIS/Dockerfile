FROM dockerhub.aganitha.ai:4443/chem/pyscf_molconsul:ovhv1

# Set the working directory in the container
WORKDIR /app

# Install mamba
RUN conda install -y mamba -c conda-forge

# Clone the GitHub repository
RUN git clone https://github.com/rinikerlab/GNNImplicitSolvent.git

# Change directory to the cloned repo
WORKDIR /app/GNNImplicitSolvent

# Create the environment using mamba
RUN mamba env create -f environment.yml


# Use explicit source of conda.sh and activate env to install the package
RUN /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && conda activate GNNImplicitSolvent && pip install ."


# Default to an interactive shell
CMD ["/bin/bash"]

