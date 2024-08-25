#!/bin/bash

# Directory to search
dir_to_search="/path/to/your/directory"

# Find files larger than 1M, calculate their MD5 hash, sort by hash, find duplicates based on hash, print directory
find "$dir_to_search" -type f -size +1M -exec md5sum {} \; | sort | uniq -w32 -dD | awk '{ print $2 }' | xargs -I {} dirname {}