#!/bin/bash
fgrep -H ">= 0.0" *_summary | awk '{print $1, $NF}'
