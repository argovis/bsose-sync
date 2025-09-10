FROM python:3.9
RUN apt-get update -y && apt-get install -y nano
RUN pip install pymongo pandas numpy h5py netCDF4 xarray

WORKDIR /app
COPY . .