FROM ubuntu:22.04
LABEL authors="Mila van Rooij, Dept. of Rheumatology, Leiden University Medical Centre"
LABEL description="BpresenT framework for BCR analysis. For academic use only."
LABEL version="1.1"
LABEL source="https://github.com/LUMC-Rheumatology/BpresenT"

# set workdir
WORKDIR /usr/src/bpresent

# set $PATH to include netmhcpan
ENV PATH="$PATH:/usr/src/bpresent/netMHCpan-4.1"

# set environment variables (disallow .pyc, display python stdout in term)
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

# expose TCP/UDP port 8000
EXPOSE 8000

# install dependencies
RUN apt-get update -y && apt-get upgrade -y
RUN apt-get install python3.10 python3-pip tcsh -y
COPY ./requirements.txt /usr/src/bpresent
RUN pip3 install -r requirements.txt
RUN pip3 install https://github.com/openvax/mhctools/tarball/441a650

# copy and initialize BpresenT
COPY . /usr/src/bpresent
RUN find . -type f,d -exec touch {} +
RUN python3 manage.py migrate
RUN python3 manage.py updategenedb && python3 manage.py initadmin && python3 manage.py populatedb

# on run, start server at port 8000
CMD ["python3", "manage.py", "runserver", "0.0.0.0:8000"]
