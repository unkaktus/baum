#!/usr/bin/env python

import subprocess
import sys

def execute(cmd):
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)

# Allow for 5 lines after the nan
grace = 5

for path in execute(sys.argv[1:]):
    # Print all the lines
    print(path, end="")
    if path.__contains__('found NAN/INF'):
        grace -= 1
    if grace == 0:
        print("Found NAN in the output, terminating.")
        sys.exit(1)
