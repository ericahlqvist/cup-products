#!/bin/sh
screen -d -m salloc -C amd --time=3-00:00:0 --nodelist diablo06 --exclusive --mail-type=BEGIN,END --mail-user=eric.ahlqvist@ed.ac.uk -J diablo06 --no-shell
