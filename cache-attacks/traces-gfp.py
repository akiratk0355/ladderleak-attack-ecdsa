import sys
import time
import datetime
import re
import glob
import os.path

PROBE_THR = int(100)

#Parse traces
traces_begin = []
traces_end  = []
dir = sys.argv[1] + "/"

max = 0
while(True):
    file_name = "out."+str(max).zfill(5)+".dat"
    if os.path.exists(dir+file_name) == False:
         break
    max += 1

for ids in range(max):
    file_name = "out."+str(ids).zfill(5)+".dat"
    file = open(dir+file_name, "r")
        
    first_line = file.readline()
    list_line = first_line.split("./FR-trace starting at ")
    str_time = list_line[-1].rstrip()
    str_time = time.mktime(datetime.datetime.strptime(str_time, "%a %b %d %H:%M:%S %Y").timetuple())
    aux_1 = []
    aux_1.append(int(str_time))
    aux_2 = []
    aux_2.append(int(str_time))

    i = 0 ;
    for line in file:
        i = i + 1
        if (i < 40):
            continue
        if(i > 60):
            break
        if line[0] == "#":
            continue

        times = line.split(" ")
        if times[0] == "-":
            aux_1.append(-1)
        else:
            time_cycle_1 = int(times[0].rstrip())
            #Probe hit
            if(time_cycle_1 < PROBE_THR):
                aux_1.append(1)
            #Probe miss
            else:
                aux_1.append(0)

        if times[1] == "-":
            aux_2.append(-1)
        else:
            time_cycle_2 = int(times[1].rstrip())
            #Probe hit
            if(time_cycle_2 < PROBE_THR):
                aux_2.append(1)
            #Probe miss
            else:
                aux_2.append(0)
    
    #Add trace to dictionary
    traces_begin.append(aux_1)
    traces_end.append(aux_2)
    file.close()
    
#Parse log file
log = []
aux_time = ""
for line in open(glob.glob(dir+"secp192*1.log")[0]):
    if "WARNING" in line:
        continue
    if "cycles" in line:
        continue
    if  "Bit" in line:
        bit = int(line[-2])
        aux = []
        aux.append(aux_time)
        aux.append(bit)
        log.append(aux)
        aux_time = ""
    else:
        aux_time = line.rstrip()
        aux_time = time.mktime(datetime.datetime.strptime(aux_time, "%a %b %d %H:%M:%S %Y").timetuple())


for idx in range(len(traces_begin)):
    trace_b = traces_begin[idx]
    trace_e = traces_end[idx]
    if (idx > 1 and trace_b[0] != traces_begin[idx-1][0]):
        print("\n")
    else:
        print("")
    sys.stdout.write(time.asctime(time.localtime(trace_b[0])) + ", trace " + str(idx))
    for t in range(len(log)):
        if (log[t][0] == trace_b[0]):
            sys.stdout.write(", sig " + str(t) + ", " + str(log[t][1]))
    sys.stdout.write(": ")
    for idx2 in range(len(trace_b)):
        if (trace_b[idx2] == 1 and trace_e[idx2] == 1):
            sys.stdout.write("(BE)")
        elif (trace_b[idx2] == -1 or trace_e[idx2] == -1):
            sys.stdout.write("-")
        elif (trace_b[idx2] == 1):
            sys.stdout.write("B")
        elif (trace_e[idx2] == 1):
            sys.stdout.write("E")
        else:
            sys.stdout.write("|")
