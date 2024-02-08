# Load Base Image
FROM ubuntu:20.04

# User definitions
ARG USER_NAME="lmmp"
ARG PASS="puc-rio"

# Update apt and install python3.7
RUN apt-get update && apt-get install -y \
        software-properties-common
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt-get update && apt-get install -y \
    python3.9 \
    python3-pip
RUN python3.9 -m pip install pip
RUN apt-get update && apt-get install -y \
    python3-distutils \
    python3-setuptools \
    iproute2  \
    dhcpcd5 \
    iputils-ping
ENV PYTHONUTF8=1

# Install OpenSSH
RUN apt install openssh-server sudo -y

# Install tmux 
RUN apt-get install tmux -y

# Install FEniCs
RUN apt-get install software-properties-common -y
RUN add-apt-repository ppa:fenics-packages/fenics -y
RUN apt-get update -y
RUN apt-get install --no-install-recommends fenics -y

# Create User (System account, HomeDirectory, bash as default, SUDO Group)
RUN useradd -rm -d /home/${USER_NAME} -s /bin/bash -g root -G sudo ${USER_NAME}
RUN echo "$USER_NAME:$PASS" | chpasswd

# Define container Volume and Working Directory
VOLUME /home/${USER_NAME}
WORKDIR /home/${USER_NAME}

# Install final dependencies
## Make requirements
COPY ./requirements.txt /home/${USER_NAME}
# RUN sudo chown -R ${USER_NAME} ..
RUN pip3 install -r requirements.txt

# Set User
USER ${USER_NAME} 
ENTRYPOINT sudo service ssh restart && bash

