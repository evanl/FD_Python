#!/bin/bash

for nx in 128 256 512 
do
  python fd2_particle.py $nx >> gridprofile.out
  echo " completed $nx cells " 
done
