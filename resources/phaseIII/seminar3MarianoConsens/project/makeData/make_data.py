# this initially just made admissions but it makes sense to generate 
# ICU stays from this same file

import time
from datetime import *
import random
import csv

# Initialize before generating a particular admission
# use Jan 1, 2000 as start date for generated dates
first = datetime(2000, 1, 1, 0, 0)
# at first just make a list of los for 4999 real visits 
f = open('los.data','r')
los_list = f.readlines()
# now los is list of 4999 strings 'ddd:hh:mm:00'
admission_types = ['EMERGENCY','URGENT','NEWBORN','ELECTIVE']

##
admissions = []
header = ['row_id', 'subject_id', 'hadm_id', 'admittime', 'dischtime', 'deathtime', 'admission_type']
admissions.append(header)

icustays = []
header = ['row_id','subject_id','hadm_id','icustay_id','intime','outtime']
icustays.append(header)

icu_counter = 1

for i in range(1, 2000):
    
    # randomly add between 0 and 355 days
    admittime = first + timedelta(random.randint(0,355))
    # now add randomly between 0 and 24*60-1 minutes
    admittime += timedelta(minutes=random.randint(0,24*60-1))
    
    # pick randomly from real los and create timedelta object
    choice = random.randint(0,4998) # this is inclusive of endpoints
    try: 
       day,hour,min,sec = los_list[choice].strip().split(':')
    except ValueError:
       print(choice, los_list[choice])
    los = timedelta(days=int(day),hours=int(hour),minutes=int(min))
    dischtime = admittime + los
    
    # pick a random subject from our subject list of 1000
    # not great since distribution over subjects won't be uniform in real data
    subject_id = random.randint(1,1000)
    
    # row_id and hadm_id will be the same generated in sequence
    # warn students not to count on this
    hadm_id = i
    
    # distribution from real data
    # ELECTIVE	7706 URGENT	1336 EMERGENCY	42071 NEWBORN	7863  
    # ELECTIVE	.14 URGENT	.02 EMERGENCY	.71 NEWBORN	.13
    choice = random.randint(0,99)
    if choice < 71:
        admission_type = 'EMERGENCY'
    elif choice < 84:
        admission_type = 'ELECTIVE'
    elif choice < 98:
        admission_type = 'NEWBORN'
    else:
        admission_type = 'URGENT'
        

    row = [i, subject_id, i, admittime.isoformat(), dischtime.isoformat(),"",admission_type]
    admissions.append(row)

    # We are done the admission. Now generate the corresponding ICU stays (if any)
    choice = random.randint(1, 1000)
    if choice < 20:
        pass    # 2% don't have an ICU stay
    #elif choice < 952:  # get exactly one for now don't make patients with 2+
    else:
        # generate the start time of the icu visit 
        # from 1/10th, 1/9th, ... 1/2 was through the hadm visit
        denom = random.randint(2,10)
        offset = los / denom;
        intime = admittime + offset;
        # generate end time randomly as 1/1, 1/2, 1/3, ... 1/10 of visit
        denom = random.randint(1,10)
        length = los / denom;
        outtime = intime + length;
        # but if outtime > dischtime then leave it as "" to become Null
        if outtime > dischtime:
            outtime = ""
        else:
            outtime = outtime.isoformat()
  
        # outtime already converted to string intime is not  
        icu_row = [icu_counter, subject_id, hadm_id, icu_counter, 
            intime.isoformat(), outtime]
        icustays.append(icu_row)
        icu_counter += 1

# create diagnoses data
# read in csv file of codes and counts
# create a list to pick from randomly so that we get some sort of distribution
f = open('real.diagnosis.csv','r')
diag_list = f.readlines()
diagnosis_codes_dist = []
for line in diag_list[1:]:
    count = int(line.strip().split(',')[0]) // 1000
    icd_code = line.strip().split(',')[1]
    for j in range(count):
        diagnosis_codes_dist.append(icd_code)

# on average 11 dignosis codes per admission, some have 0 but very few 47
# select how many diagnosis codes to generate by adding random number from -11 to +11# to 11

diagnoses_table = []
header = ['row_id', 'subject_id', 'hadm_id', 'seq_num', 'icd9_code']
diagnoses_table.append(header)
row_id = 1
for row in admissions[1:]:    # don't create diagnoses for the header row
    subject_id = row[1];
    hadm_id = row[2];
    # randomly decide how many diagnoses
    num_diagnoses = 11 + random.randint(-11,11)
    # now pick a diagnoses from real distribution
    for seq_num in range(1, num_diagnoses + 1):
        icd9_code = random.choice(diagnosis_codes_dist)
        diagnoses_table.append([row_id, subject_id, hadm_id, seq_num, icd9_code])
        row_id += 1
         
with open('diagnoses_icd.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(diagnoses_table)

with open('admissions.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(admissions)

with open('icustays.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(icustays)
