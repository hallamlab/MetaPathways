#!/bin/bash

(
  echo "/opt/pathway-tools"
  echo "/data"
  echo ""
  echo "n"
  echo "Y"
  echo "n"
  echo ""
) | $INSTALLER_DIR/$INSTALL_EXE

## Questions:
## 1. Destination folder: /opt/pathway-tools (default)
## 2. Data folder: "/opt/data"
## 3. Enter ""
## 4. Shortcuts? "n"
## 5. Continue? "Y"
## 6. Launch PTools now? "n"
## 7. Enter ""
