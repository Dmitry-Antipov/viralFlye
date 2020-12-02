import sys

len_min = 5000
len_max = 1000000
cov_thr = 10

contigs_list = {}
links = set()


cov = 0
for line in open(sys.argv[1]):
    t = line.split()
    if t[0] == "S":
      #print (t[3].split(":"))
      if t[3].split(":")[0] == "dp":
        cov = float(t[3].split(":")[2])
      elif t[3].split(":")[0] == "KC":
        cov = float(t[3].split(":")[2])/len(t[2])
      else:
        print ("coverage info error")
        quit()
    if cov >= cov_thr and len(t[2])>= len_min and len(t[2])<=len_max:
        contigs_list[t[1]] = [t[2],len(t[2]),cov]

    elif t[0] == "L":
        links.add (t[1])

for i in contigs_list:
    if i not in links:
        print (f">{i}_length_{contigs_list[i][1]}_cov_{contigs_list[i][2]}")
        print (contigs_list[i][0])
