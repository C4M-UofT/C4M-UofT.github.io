
f = open('../january06.txt')

f.readline()
f.readline()

# line = f.readline()
# measurements = line.split()
#temperature = float(measurements[3])
temperature =  4000

while temperature != -12.14:
    line = f.readline()
    measurements = line.split()
    temperature = float(measurements[3])

# when we get here we know that temperature == -12.14
print(measurements[0], measurements[1])

#Task: print the date when we get -12.14C for the first time