#!/usr/bin/env bash

createTables ()
{
  if [[ $1 -eq 600 ]]; then range=100; else range=99; fi

  start_depth=$1
  end_depth=$(($1 + $range))

  model=$2

  echo "Creating tables from $start_depth to $end_depth km depth..."

  for depth in $(seq $start_depth $end_depth);
  do
    ./bin/tables $depth $model
  done
  
  sleep 1
}

if [[ $# -eq 0 ]]; then
  echo "No models provided!"
  echo "Usage: $0 model1 model2 model3..."
  echo "Example: $0 iasp91 ak135f prem stw105 rem1d"
  
  exit 1
fi

mkdir -p tables

for model in $@;
do
  echo "Creating tables for model $model..."

  for depth in {0..600..100};
  do
    createTables $depth $model &
  done
  
  wait
  mkdir ./tables/$model
  mv ./tables/*.list ./tables/*.bin ./tables/$model/
  echo "Saving tables for model $model at pytracer/tables/$model..."
done

wait
echo "Done!"