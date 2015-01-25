#!/bin/sh

cat timestamp.hpp json.hpp seq.hpp ABCSpeedup.hpp | grep -v "#include \"" > submission.cpp
g++ -std=c++11 -c submission.cpp
gvim submission.cpp &
