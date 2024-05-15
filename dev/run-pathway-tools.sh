#!/bin/bash

Xvfb $DISPLAY &
/opt/pathway-tools/pathway-tools "$@"
