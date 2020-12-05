#!/bin/bash

for inst in ../instances/*.txt; do
  $1 $(basename "$inst")
done;