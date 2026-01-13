#!/usr/bin/env bash
set -euo pipefail

#for year in {2013..2024}; do
for year in {2021..2024}; do
#for year in {2020..2020}; do
    echo "=== Year $year ==="
    for start in $(seq 0 10 588); do
    #for start in $(seq 440 10 588); do
        if [ "$start" -eq 580 ]; then
            end=588
        else
            end=$((start + 10))
        fi

        echo "Running container for year=$year, range=$start,$end"
        until docker container run --network bsose --env MONGODB_URI=mongodb://database/argo -v /data/SOSE/SOSE/SO6/ITER156/5DAY/${year}:/tmp argovis/bsose-sync:hax target/release/bsose-sync /tmp/Theta_bsoseI156_${year}_5day.nc THETA $start $end 0 2160; do
            echo "❌ Failed for year=$year, range=$start,$end — retrying in 300s..."
            sleep 300
        done
        echo "✅ Completed for year=$year, range=$start,$end"
    done
done
