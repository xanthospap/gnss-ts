#! /bin/bash

FLAGS="-std=c++11 -g -Wall -Wextra -pedantic -W -Wshadow -Winline -Wdisabled-optimization"

g++ -c $FLAGS top2daz.cpp -o top2daz.o

g++ $FLAGS main.cpp -o test.e
