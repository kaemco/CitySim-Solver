FROM ubuntu:latest AS builder

# get the necessary libraries
RUN apt-get update &&\
	apt-get install --yes git make g++ gfortran libglu1-mesa-dev

# prepare the directory for the executables
RUN mkdir ~/bin

# prepare the directory for the library
RUN mkdir ~/lib
RUN mkdir ~/lib/CitySim

# prepare the directory for the sources
WORKDIR /src

# get the open-source code
# RUN git clone https://github.com/kaemco/CitySim-Solver.git
# or copy it from the current directory
COPY . /src/CitySim-Solver

# compile the code
WORKDIR CitySim-Solver
RUN make -j $(nproc)

# make a final image with the executables
FROM alpine:latest
RUN apk add --no-cache --update gcompat gfortran glu

# prepare the environnement for the executables
ENV PATH "$PATH:/root/bin"

# copy the build solver from the builder image
COPY --from=builder /root/bin /root/bin

# execute the Solver by default in the data directory
# the data directory can point to the location where you store your CitySim XML files
# for instance:
# docker build . -t citysim
# docker run -v /mnt/c/myXMLFiles:/data -it citysim myFile.xml
WORKDIR /data
ENTRYPOINT ["CitySim"]
