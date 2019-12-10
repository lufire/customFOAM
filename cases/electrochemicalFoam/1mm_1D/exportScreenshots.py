#!/usr/bin/python
import os
import shutil
import sys
from paraview.simple import *
working_dir = os.path.dirname(os.path.abspath(__file__))
result_dir = os.path.join(working_dir,'results')
state_file = os.path.join(working_dir,'postprocess.pvsm')

LoadState(state_file)
layout1 = GetLayoutByName('Layout #1')
SaveScreenshot(os.path.join(result_dir, 'CurrentDensity.png'), viewOrLayout=layout1)
