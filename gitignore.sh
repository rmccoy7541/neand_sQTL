#!/bin/sh
# This is a shell script that updates the gitignore file with all of the files in the repo tthat exceed 100MB

find ./* -size +100M | cut -c 2-|  cat >> .gitignore
