#!/bin/sh
# This is a shell script that updates the gitignore file with all of the files in the repo tthat exceed 100MB

#find ./* -size +100M | cut -c 2-|  cat >> .gitignore

# Here's another script that can do the same thing but avoid re-adding the same files
find . -name .git -prune -o -size +100M -print | cut -c3-|
while read -r name; do
     git check-ignore -q "$name" || printf '%s\n' $name >> .gitignore
done
