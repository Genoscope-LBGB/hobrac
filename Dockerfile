FROM centos/python-38-centos7

USER root
	
# python setup
RUN pip install --upgrade pip
RUN pip install --upgrade setuptools
RUN pip install wheel

# copy the git repo to the container
COPY . /gitlab/repo

# change directory for the remaining run commands
WORKDIR /gitlab/repo

# install the project module
RUN pip install .

