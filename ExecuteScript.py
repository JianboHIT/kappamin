#!/usr/bin/env python3
# 
# Calculate the Cahill minimum limit to thermal conductivity
# 
# Usage:
#   python3 ExecuteScript.py
#


from kappamin import execute


fileinput = None
# Option Values:
#     'None': auto detect filename of configuration
#   filename: assign the filename of configuration

execute(filename=fileinput)
