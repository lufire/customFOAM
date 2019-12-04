#!/usr/bin/python
import os
import shutil
import csv
import sys
import subprocess

result_file_name = 'result_data.csv'
working_dir = '.'
result_dir = '.' 
patch_name = 'currentCollector'
value_name = 'currentDensity'

dict_dir = os.path.join(working_dir,'constant/electrochemicalProperties')
#control_dir = os.path.join(working_dir,'system/controlDict')

abs_working_dir = os.path.abspath(working_dir)
path_array = abs_working_dir.split(os.path.sep)
path_length = len(path_array)
simulation_name = path_array[path_length-1]
log_file_name = 'CurrentDensity-Potential.txt'
print(simulation_name)

# Initialization run: building mesh 
subprocess.call(os.path.join(working_dir,"./Allclean"))
subprocess.call(os.path.join(working_dir,"./initRun"))

result_list = []
with open(os.path.join(result_dir,result_file_name)) as result_file:
    for line in result_file:
        line = line.strip()
        line_list = line.split(',')
        result_list.append(line_list)
        
experiment_potential = result_list[0]
experiment_current = result_list[1]

n_operating_points = len(experiment_potential) - 3

simulation_current = []
simulation_potential = []
simulation_current.append(simulation_name)
simulation_potential.append(simulation_name)
simulation_current.append('Current Density')
simulation_current.append('A/m^2')
simulation_potential.append('Potential')
simulation_potential.append('V')

find_string1 = patch_name 
find_string2 = value_name
#find_string3 = end_name

dict_position = 0
control_position = 0

linecounter = 0
with open(dict_dir, 'r') as dict_file:
    dict_list = dict_file.readlines()
#    for line in dict_file:
#        dict_list.append(line)

with open(dict_dir, 'r') as dict_file:
    for line in dict_file:
        linecounter += 1
        if find_string1 in line:
#            print linecounter
#            print line
            for line in dict_file:
                linecounter += 1
                if find_string2 in line:
                    dict_position = linecounter  
                    break
              
#linecounter = 0
#with open(control_dir, 'r') as control_file:
#    control_list = control_file.readlines()
##    for line in dict_file:
##        dict_list.append(line)
#
#with open(control_dir, 'r') as control_file:
#    for line in control_file:
#        linecounter += 1
#        if find_string1 in line:
#            control_position = linecounter  
#            break

for i in range(0,n_operating_points):
   
    new_dict_list = dict_list
    #replace_string = '    potential       '+(experiment_potential[i+3])+';\n'
    replace_string = '    currentDensity  '+str(experiment_current[i+3])+';\n'

    new_dict_list[dict_position-1] = replace_string
#    print new_dict_list[valueposition-1]
    #new_dict_list[valueposition-1] = re.sub(str(BaseCaseCellPotential),str(CellPotential[i+3]),dict_list[valueposition-1])

    #createCase( projDir, mycasename, rowdata[ 'FLMFilename' ], dict_dir )
    
    with open(dict_dir, 'w') as dict_file:
        for item in new_dict_list:
            dict_file.write(item)

    #new_control_list = control_list
    ##replace_string = '    potential       '+(experiment_potential[i+3])+';\n'
    #if experiment_current[i+3] > 3000:
    #    end_time = 1.0
    #else:
    #    end_time = 0.5

    #replace_string = 'endTime         '+str(end_time)+';\n'

    #new_control_list[control_position-1] = replace_string
#    print new_dict_list[valueposition-1]
    #new_dict_list[valueposition-1] = re.sub(str(BaseCaseCellPotential),str(CellPotential[i+3]),dict_list[valueposition-1])

    #createCase( projDir, mycasename, rowdata[ 'FLMFilename' ], dict_dir )
    
    #with open(control_dir, 'w') as control_file:
    #    for item in new_control_list:
    #        control_file.write(item)
    # Prepare operating point
    #os.system("pyFoamWriteDictionary.py os.path.join(cfd_dir,'0','phiEs') 
    #    "boundaryField['currentCollector']['value']" 
    #    str("uniform "+str(experiment_potential[i+3]))")

    # Simulation for operating point 
    subprocess.call(os.path.join(working_dir,"./runDataPoint"))
    
    # Read average current density from log file 
    with open(os.path.join(working_dir,log_file_name)) as log_file:
        lines = log_file.readlines()
        currD = lines[0].split()[4]
        simulation_current.append(currD)
        pot = lines[1].split()[3]
        simulation_potential.append(pot)

result_list.append(simulation_current)
result_list.append(simulation_potential)

with open(os.path.join(result_dir,result_file_name),"a") as result_file:
    writeString = ','.join([str(item) for item in simulation_current])
    result_file.write(writeString + "\n")
    writeString = ','.join([str(item) for item in simulation_potential])
    result_file.write(writeString + "\n")
