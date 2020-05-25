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
    #print(ids)
    #print(aux_1)
    #print(aux_2)
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

guessing = []
guess_one = 0
guess_zero = 0

# Using guessing 
for idx in range(len(traces_begin)):
        trace_b = traces_begin[idx]
        trace_e = traces_end[idx]
        aux = []
        aux.append(trace_b[0])
        trace_b.remove(trace_b[0])
        trace_e.remove(trace_e[0])
        pos = 1
        second = 2
        start = 0
        v = -1

        idx2 = 0
        while idx2 < len(trace_e)-5:
                # Look for the second hit (on the end monitor)
                if (trace_e[idx2] == 1):
                        if (trace_b[idx2+2] == 1) and (trace_e[idx2+1] == 1) and (trace_e[idx2+2] == 1) and (trace_e[idx2+5] == 0):
                                v = 0
                                guess_zero += 0
                        else:
                                v = 1
                                guess_one += 1
                        break
                idx2 += 1
        aux.append(v)
        guessing.append(aux)

one_right = 0
one_x = 0
one_0 = 0
one_err = 0
zero_right = 0
zero_x = 0
zero_1 = 0
zero_err = 0

for i in range(len(log)):
    found = False
    for j in range(len(guessing)):
        if log[i][0] == guessing[j][0]:
            found = True
            if log[i][1] == 1:
                    if log[i][1] == guessing[j][1]:
                            one_right += 1
                    elif guessing[j][1] == 0:
                            one_0 += 1
                    elif guessing[j][1] == -1:
                            one_x += 1
                    else:
                            print("Invalid", guessing[j][1])

            elif log[i][1] == 0:
                    if log[i][1] == guessing[j][1]:
                            zero_right += 1
                    elif guessing[j][1] == 1:
                            zero_1 += 1
                    elif guessing[j][1] == -1:
                            zero_x += 1
                    else:
                            print("Invalid", guessing[j][1])
            else:
                print("Invalid log: ", log[i][1])
            break
    if found == False:
        if log[i][0] == 0:
            zero_err += 1
        else:
            one_err += 1

print("### Results ###")
print("Correct        :           1("+str(one_right)+")\t  0("+ str(zero_right) +")")
print("Unknown        :           1("+str(one_x)+")\t  0("+ str(zero_x)+")")
print("False Positives:           1("+str(one_0)+")\t  0("+ str(zero_1)+")")
print("Mismatched time:           1("+str(one_err)+")\t  0("+ str(zero_err)+")")
print("")
