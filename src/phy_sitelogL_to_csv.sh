#! /bin/bash

tail -n1 $1 | cut -f2- -d" " | sed 's/^   //g' | tr " " "\n"