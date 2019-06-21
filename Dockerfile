FROM ubuntu:latest as builder
#RUN apt-get update && apt-get install -y wget \
#    gpg
RUN apt-get update && apt-get install -y   git make gcc wget gpg vim python3-dev libopenblas-dev libfreetype6-dev pkgconf gfortran
RUN wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB 
RUN apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
RUN wget https://apt.repos.intel.com/setup/intelproducts.list -O /etc/apt/sources.list.d/intelproducts.list
RUN apt-get update &&  apt-get install -y intel-mkl-2018.2-046 libxml2-dev
COPY genomesim /genomesim
COPY mkl_funcs_list /genomesim
RUN cd /opt/intel/mkl/tools/builder && make libintel64 export=/genomesim/mkl_funcs_list threading=sequential && mv mkl_custom.so /genomesim/lib/mkl_small.so && cd /genomesim
WORKDIR /genomesim/
RUN /bin/bash -c "source /opt/intel/bin/compilervars.sh -arch intel64 -platform linux;make"
RUN make
CMD ["/bin/bash"]
