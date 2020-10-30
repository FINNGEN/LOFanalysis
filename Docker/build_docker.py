#!/usr/bin/env python3

import shlex
from subprocess import Popen, PIPE,call,check_output
import argparse


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Build Docker file for variant filtering")

    parser.add_argument("--image", type= str,
                        help="name of image",default = 'lof')
    parser.add_argument("--version", type= str,
                        help="version value, e.g.0.001",required = True)
    parser.add_argument("--push",action = 'store_true')
    parser.add_argument("--args",type = str,default = '')
    args = parser.parse_args()

    docker_path = "eu.gcr.io/finngen-refinery-dev/"
    cmd = f"docker build -t {docker_path}{args.image}:{args.version} -f Dockerfile .. {args.args}"    
    print(cmd)
    call(shlex.split(cmd))

    if args.push:
        cmd = f"docker -- push {docker_path}{args.image}:{args.version}"       
        call(shlex.split(cmd))

