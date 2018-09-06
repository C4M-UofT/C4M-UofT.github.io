# row_id     | integer                        | not null
# subject_id | integer                        | not null
# gender     | character varying(5)           | not null
# dob        | timestamp(0) without time zone | not null
# dod        | timestamp(0) without time zone | 

import time
from datetime import *
import random
import csv

# Initialize before generating any particular patient
# use Jan 1, 1960 as start date for generated dates
initial = datetime(1960, 1, 1, 0, 0)

##
patients = []
header = ['row_id', 'subject_id', 'gender', 'dob', 'dod']
patients.append(header)

# use i for patient_id number which is 1 to 1000 inclusive
for i in range(1, 1001):

    # in the real data it is 44% female
    choice = random.randint(1,100)
    if choice <= 44:
        gender = 'F'
    else:
        gender = 'M'
 
    # randomly add somewhere between +/- 40 years (resolution is days)
    # now add randomly between 0 and 24*60-1 minutes
    dob = initial + timedelta(random.randint(-40*365,40*365))
    
    # row_id and subject_id will be the same generated in sequence
    # warn students not to count on this
    
    # for now everyone had dod of Null
    # actually Null for 30761 of X -- others range from over 300 to 0 or 1 day
    row = [i, i, gender, dob.isoformat(),""]
    patients.append(row)


with open('patients.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(patients)
