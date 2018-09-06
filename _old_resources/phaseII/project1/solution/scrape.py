import urllib
import os

url = "https://en.wikipedia.org/wiki/List_of_infectious_diseases"

page = urllib.request.urlopen(url).read().decode()

os.chdir("/home/guerzhoy/Desktop/diagnosis")

end_ind = page.find("<tr>")-1

while page.find("<tr", end_ind) >= 0:
    start_ind = page.find("<a href=\"", end_ind) + len("<td><b><a href=")
    end_ind = page.find("\" title", start_ind)
    disease = page[start_ind:end_ind]
    url = "https://en.wikipedia.org/wiki/" + disease
    try:
        dpage = urllib.request.urlopen(url).read().decode()
        f = open("wikipages/" + disease + ".html", "w")
        f.write(dpage)
        f.close()
    except:
        print("Skipping", disease)
    
    end_ind = page.find("<tr", end_ind)-1