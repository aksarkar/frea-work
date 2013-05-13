#!/bin/bash
tr -d "\000" | sed 's/^.*hb3uucp00000000000000//' | awk 'NF > 0'
