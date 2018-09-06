HR = 80
QT = 300
sex = "male"
RR = 

QTc = QT / (RR ** (1 / 3))   # Compute the corrected QT interval

if sex == "female":
    if QTc > 450:
        print("Abnormal QTc")
    elif QTc >= 431:
        print("Borderline QTc")
    else:
        print("Normal QTc")
elif sex == "male":
    if QTc > 470:
        print("Abnormal QTc")
    elif QTc >= 451:
        print("Borderline QTc")
    else:
        print("Normal QTc")

else:
    print("Enter 'male' or 'female' as the sex")

