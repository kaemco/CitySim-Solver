#!/bin/bash
FILES="./*/*.xml"

for f in $FILES
do
  echo "Simulation $f..."
  # take action on each file. $f store current file name
  ~/bin/CitySim $f
done

echo "Done."

