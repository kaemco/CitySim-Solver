FROM alpine:latest AS builder

# get the necessary libraries
RUN apk add --no-cache \
	--update make g++ gfortran mesa-dev glu-dev

# prepare the directory for the executables
RUN mkdir ~/bin

# prepare the directory for the library
RUN mkdir ~/lib
RUN mkdir ~/lib/CitySim

# prepare the directory for the sources
WORKDIR /src

# get the open-source code
# RUN apk add --no-cache --update git
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

# execute the Solver with the xml file in the mounted directory
# it will also save a log.txt file
# for instance:
# docker build . -f Dockerfile.ubuntu -t citysim
# docker run --rm -v /mnt/c/mySimulationDir/:/data citysim
# As an alternative, it is possible to override the CMD:
# docker run --rm -v /mnt/c/mySimulationDir/:/data citysim CitySim myxml.xml
WORKDIR /data
CMD ["sh", "-c", "CitySim *.xml | tee log.txt"]
