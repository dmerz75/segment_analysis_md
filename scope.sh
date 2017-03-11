#!/bin/bash

# find . -name *.cpp
# find . -name *.h
# comment src/deprecated/*

# cscope -b -q -R
# cscope -b -q -f cscope.files
# cscope -b -q -k -R -i cscope.files
cscope -b -i cscope.files
