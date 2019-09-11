FROM digitalproteomes/gosu

LABEL maintainer="Patrick Pedrioli"
LABEL Description="A container for the TP data vizualization tools" Version="1.0"

## Let apt-get know we are running in noninteractive mode
ENV DEBIAN_FRONTEND noninteractive

## Install python and R
RUN apt-get update
RUN apt-get install -y python3 r-base

# Install R dependencies
COPY r_deps.r /tmp
RUN Rscript /tmp/r_deps.r

## Install any additional package specified in the requirements.txt
# Anaconda installing
RUN wget https://repo.continuum.io/archive/Anaconda3-5.0.1-Linux-x86_64.sh \
    && bash Anaconda3-5.0.1-Linux-x86_64.sh -b \
    && rm Anaconda3-5.0.1-Linux-x86_64.sh

# Set path to conda
#ENV PATH /root/anaconda3/bin:$PATH
ENV PATH /home/user/anaconda3/bin:$PATH

# Updating Anaconda packages
#RUN conda update anaconda
RUN conda update --all

COPY requirements.txt /tmp/
RUN while read requirement; do conda install --yes $requirement; done < /tmp/requirements.txt

## Copy scripts and data into the container
COPY src /usr/local/src
COPY data /usr/local/data

RUN ln -s /usr/local/src/plot_proteo_fmi_table/plot_proteo_fmi_table.py /usr/local/bin \
    && ln -s /usr/local/src/plot_proteo_fmi_table/preprocessor2fmi.py /usr/local/bin \
    && ln -s /usr/local/src/plot_bar/plot_bar.r /usr/local/bin \
    && ln -s /usr/local/src/plot_word_cloud/plot_word_cloud.r /usr/local/bin \
    && ln -s /usr/local/src/preprocessor/preprocessor.r /usr/local/bin \
    && ln -s /usr/local/src/preprocessor/preprocessor.py /usr/local/bin \
    && ln -s /usr/local/src/qc/tpviz_input_qc.py /usr/local/bin
