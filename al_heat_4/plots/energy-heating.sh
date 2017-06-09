#!/bin/bash
gnuplot energy-heating.plot
convert -density 300 energy-heating.eps energy-heating.png 
