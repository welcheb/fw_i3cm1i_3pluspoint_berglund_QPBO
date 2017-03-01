# Dockerfile for https://index.docker.io/u/welcheb/fw_i3cm1i_3pluspoint_berglund_qpbo/
FROM welcheb/ubuntu:8.04

MAINTAINER E. Brian Welch <brian.welch@vanderbilt.edu>

# update and install dependencies
#RUN apt-get -y update && apt-get install -y \
RUN apt-get install -y \
	libgl1-mesa-glx \
	libglu1-mesa \
	libxt6 \
	unzip

# copy DixonApp_LINUX.exe.zip
COPY ./berglund/QPBO/DixonApp/LINUX/DixonApp_LINUX.exe.zip /usr/local/bin/DixonApp_LINUX.exe.zip

# copy lib files
COPY ./berglund/QPBO/DixonApp/LINUX/lib/*.* /usr/local/lib/

# unzip DixonApp_LINUX.exe.zip
RUN cd /usr/local/bin && unzip DixonApp_LINUX.exe.zip

# scratch directory for volume mount
RUN mkdir /scratch && echo "This is the scratch folder" > /scratch/scratch.txt

# array form of entrypoint to support additional argument from docker run command
ENTRYPOINT ["/bin/bash"]
