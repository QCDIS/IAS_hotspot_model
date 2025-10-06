#! /bin/bash


for f in /mnt/inputs/*.zip; do
  [ -e "$f" ] && unzip -o "$f" -d /mnt/inputs
done

cd /wrp

# Call main.r script with input arguments
Rscript --vanilla main.r "$@"

cd /mnt/outputs
for f in *; do
  [ -f "$f" ] &&  zip -q "${f%.*}.zip" "$f"
  # If we have a folder zip it
  [ -d "$f" ] &&  zip -rq "${f%.*}.zip" "$f"
done

cd /wrp