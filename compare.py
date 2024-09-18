def compare_strings(file1,file2):
    # print("equality of:" , file1 , " and " , file2 , ":")
    # s1 = open(file1,"r")
    # s2 = open(file2,"r")
    # return s1 == s2'
    print("differences between", file1, " and ", file2 , ":")
    from difflib import Differ
    with open(file1) as f1, open(file2) as f2:
        differ = Differ()
        diff = differ.compare(f1.readlines(),f2.readlines())
        c = 0               #counts lines   
        e = -1               #marker for first error
        for line in diff:
            if not line.startswith(" "):
                # print(line)
                e = c
                print("first difference in line:", e)
                break
            c += 1
    # print("(end list of different lines)")
    # if d != 0:
    #     print("first difference in line:", e)
    # else:
    if e < 0:
        print("no difference; linecount:" , c)

def uniform_zero(s):        # takes string of number, writes zeros in same way(mainly: remove "-")
    nice_s = ""
    for i in range(len(s)):
        if s[i] == "-":
            if s[i+1] == "0" and s[i+2] == ".":
                if s[i+3] == " " or s[i+3] == "+":
                    nice_s += " "                   # -0. -> 0. 
                elif s[i+3] == "j":
                    nice_s += "+"                   # -0.j -> +0.j
                else:
                    nice_s += "-"
            else:
                nice_s += "-"
        else:
            nice_s += s[i]                          # add next character
    return nice_s

def nice_string(file):
    with open(file,"r") as f:
        s = f.readlines()
    new_f = open(file[:-4]+"_nice.txt","w")
    for lines in s:
        nice_s = uniform_zero(lines)
        new_f.write(nice_s)
    new_f.close()