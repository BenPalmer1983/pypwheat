#!/bin/bash
gnuplot force-heating.plot
convert -density 300 force-heating.eps force-heating.png 
