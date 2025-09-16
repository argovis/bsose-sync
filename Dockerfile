FROM rust:1.82.0

RUN apt-get update -y && apt-get install -y nano curl wget cmake libhdf5-serial-dev libnetcdff-dev netcdf-bin
WORKDIR /app
COPY . .
RUN cargo build --release
#RUN chown -R 1000660000 /app
#CMD bash run.sh
