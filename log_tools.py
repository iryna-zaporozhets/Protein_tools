"""
Script contains functions required for logging important
info about the script execution.

Capabilities:
    Logging execution timestamp, location of the original script
    (including files pointed to with symbolic links),


"""

import sys
import subprocess
import time
import os 
import argparse
import inspect

def log_execution_info(output=sys.stdout,
                       log_python_info=True,
                       log_repo_version=True, 
                       log_doc_string=True,
                       ):

    output.write(70*"="+"\n")
    # log the execution command and execution time
    output.write(' '.join(sys.argv) + '\n')
    timestamp =  time.strftime("%a, %d %b %Y %H:%M:%S  %z", time.localtime())
    output.write(f"Execution began on {timestamp} \n")
    output.write("\n")

    # log location of the script that calls the current function
    stack = inspect.stack()
    path_parent = os.path.abspath(stack[1][1])
    output.write(f"Script source code located at {os.path.dirname(path_parent)}\n")
    output.write(f"Working directory: {os.getcwd()} \n")
    output.write("\n")

    # log the doc string of the script
    if log_doc_string:
        output.write("\n")
        output.write("DESCRIPTION\n")
        parentframe = stack[1][0]
        module = inspect.getmodule(parentframe)
        doc = inspect.getdoc(module)
        if doc is None:
            output.write("No docstring found")
        else:
            output.write(doc)
        output.write("\n")

    if log_python_info:
        output.write('\n')
        output.write("PYTHON INFO \n")
        output.write(f"Version : {sys.version} \n")
        output.write(f"Executable location : {sys.executable} \n")
        output.write('\n')

    if log_repo_version:
        output.write("VERSION CONTROL INFORMATION \n")
        try:
            label = subprocess.run(["git", "describe", "--always", "--first-parent", "--long", '--abbrev=14'], capture_output=True).stdout.decode(sys.stdout.encoding)
            output.write("Current version of the git repository: ")
            output.write(label)
        except: 
            output.write("Was not able to get git repository info.\n")
        try:
            result = subprocess.run(["git", "status", "-s", "--porcelain"], capture_output=True).stdout
            if len(result) > 0:
                if result.split()[0]  == 'fatal:':
                    output.write("No information about git tree state was found. \n")
                else:
                    output.write("Repository tree state \n")
                    output.write(result.decode(sys.stdout.encoding))
            else:
                output.write("Clean repository tree")
        except:
            output.write("No information about git tree state was found. \n")
    
    output.write(70*"="+"\n")
    return

if __name__ == '__main__':
    log_execution_info()
