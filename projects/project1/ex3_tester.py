# A program to test the read_in_symptoms function

# The dictionary we want to use for this test
symp_dict = {"Pat_1": ["symptom1", "symptom2"], "Pat_2": ["symptom2"]}

# Create a file based on this test dictionary
filename = "test.txt"
f = open(filename, "w")

# Iterate over each item in the dictionary
for name, symps in symp_dict.items():
    # build the appropriate line for this patient and write it out
    f.write(" ".join([name] + symps) + "\n")
f.close()

# Next, use the function being tested to read this new file
# back into a dictionary and compare that to the expected dictionary
if read_in_symptoms1(filename) == symp_dict:
    print("TEST 1 PASSED")
else:
    print("TEST 1 FAILED")


# Test 2 uses a different dictionary for the test
symp_dict = {"Patient1": ["Sym1", "Sym2"], "Patient2": [], "Patient3": ["Sym3"]}

# The remainder of test 2 is essentially the same
filename = "test.txt"
f = open(filename, "w")

for name, symps in symp_dict.items():
    f.write(" ".join([name] + symps) + "\n")
f.close()

if read_in_symptoms1(filename) == symp_dict:
    print("TEST 2 PASSED")
else:
    print("TEST 2 FAILED")
