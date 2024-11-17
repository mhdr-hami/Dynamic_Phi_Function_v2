import numpy as np 
import sys
import matplotlib.pyplot as plt
import scipy.stats as st 
from collections import OrderedDict


plt.rcParams["figure.figsize"] = [7.00, 7.00]
plt.rcParams["legend.framealpha"] = 1.00
## Args: PythonAdrress Domain #Experiment #Policies #Weights DataAdrress

weight_to_num = {"1.50":0, "2.00":1, '3.00':2, '4.00':3, '5.00':4, '6.00':5, '7.00':6, '8.00':7, '9.00':8, '10.00':9}
# weight_to_num = {"1.50":0, "2.00":1, "2.50":2, "3.00":3, "3.50":4, "4.00":5, "4.50":6, "5.00":7, "5.50":8, "6.00":9, "6.50":10, "7.00":11, "7.50":12, "8.00":13, "8.50":14, "9.00":15, "9.50":16, "10.00":17}
# weight_to_num = {"1.20":0, "1.40":1, "1.60":2, "1.80":3, "2.00":4, "2.20":5, "2.40":6, "2.60":7, "2.80":8, "3.00":9, "3.20":10, "3.40":11, "3.60":12, "3.80":13, "4.00":14, "4.20":15, "4.40":16, "4.60":17, "4.80":18, "5.00":19, "5.20":20, "5.40":21, "5.60":22, "5.80":23, "6.00":24, "6.20":25, "6.40":26, "6.60":27, "6.80":28, "7.00":29, "7.20":30, "7.40":31, "7.60":32, "7.80":33, "8.00":34, "8.20":35, "8.40":36, "8.60":37, "8.80":38, "9.00":39, "9.20":40, "9.40":41, "9.60":42, "9.80":43, "10.00":44}

cost='5.0'
mapType = {0:'Obstacle Square', 1:'Swamped Square Cost='+cost, 2:'Obstacle Diamond', 3:'Swamped Diamond Cost='+cost, 4:'Obstacle Circle', 5:'Swamped Circle Cost='+cost, 6: sys.argv[5]+" Cost="+sys.argv[8]}

# markers = [',', '^', 'X', '8', 's', 'o'] ##MAP=5
# int_to_alg = {0:'WA*', 1:'PWXD', 2:'PWXU', 3:'XDP', 4:'XUP', 5:'MAP'}
# colours = ['tab:gray', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:cyan']

# markers = [',', '^', 'X', '8', 's', '*'] ##DWP=5
# int_to_alg = {0:'WA*', 1:'PWXD', 2:'PWXU', 3:'XDP', 4:'XUP', 5:'DWP'}
# colours = ['tab:gray', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:blue']

# markers = [',', '^', 'X', '8', 's', 'o', 'd', '*'] ##DWP=7, MAP=5
# int_to_alg = {0:'WA*', 1:'PWXD', 2:'PWXU', 3:'XDP', 4:'XUP', 5:'MAP', 6:'DPS', 7:'DWP'}
# colours = ['tab:gray', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:cyan', 'tab:olive', 'tab:blue']

markers = [',', '^', 'X', '8', 's', '*', 'd', 'o'] ##DWP=5, MAP=7
int_to_alg = {0:'WA*', 1:'PWXD', 2:'PWXU', 3:'XDP', 4:'XUP', 5:'DWP', 6:'DPS', 7:'MAP'}
colours = ['tab:gray', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:blue', 'tab:olive', 'tab:cyan']

# markers = [',', '^', 'X', '8', 's', '*', 'o'] ##DWP=5, MAP=6
# int_to_alg = {0:'WA*', 1:'PWXD', 2:'PWXU', 3:'XDP', 4:'XUP', 5:'DWP', 6:'MAP'}
# colours = ['tab:gray', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:blue', 'tab:cyan']

# markers = [',', '^', 'X', '8', 's', 'o', '*'] ##DWP=6, MAP=5
# int_to_alg = {0:'WA*', 1:'PWXD', 2:'PWXU', 3:'XDP', 4:'XUP', 5:'MAP', 6:'DWP'}
# colours = ['tab:gray', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:cyan', 'tab:blue']


linestyles = OrderedDict(
    [('dashed',              (0, (5, 5))),
     ('dashdotted',          (0, (3, 5, 1, 5))),
     ('dotted',              (0, (1, 5))),
     ('long dash with offset', (5, (10, 3))),
     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('solid',               (0, ())),
     ('loosely dotted',      (0, (1, 10))),
     ('densely dashed',      (0, (5, 1))),
     ('densely dotted',      (0, (1, 1))),
     ('loosely dashed',      (0, (5, 10))),
     ('loosely dashdotted',  (0, (3, 10, 1, 10))),
     ('densely dashdotted',  (0, (3, 1, 1, 1))),
     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))])
showErrorBar = False
tableType = 0

if sys.argv[1] == '-stp':
    if sys.argv[2] == '1':

        ##Experiment 1: Creates a table from average of runs of different problems over different weights
        table = np.zeros((int(sys.argv[3]), int(sys.argv[4])))
        count_table = np.zeros((int(sys.argv[3]), int(sys.argv[4])))
        TheDataSet = [[[] for _ in range(int(sys.argv[4]))] for _ in range(int(sys.argv[3]))]

        with open("./papers/DSDWA/results/"+sys.argv[5]+".txt", "r") as f:
            for line in f:
                data = line.split()
                if(len(data)): ## To check for empy lines
                    if data[0] == "STP" and (data[5] in list(weight_to_num.keys())) and (int(data[3]) in list(int_to_alg.keys())) : #and data[3]!='0' and data[3]!='1' and data[3]!='6' and data[3]!='8': # and data[9]!='0':
                        table[int(data[3])][weight_to_num[data[5]]] += float(data[7])
                        count_table[int(data[3])][weight_to_num[data[5]]] += 1
                        TheDataSet[int(data[3])][weight_to_num[data[5]]].append(int(data[7]))
                        # print(float(data[7]), count_table[int(data[3])][weight_to_num[data[5]]], sep=" ")

        result = np.divide(table, count_table)

        print()
        print('====================================================  Average Expansions Table  =====================================================')
        print('Alg / Weight|   1.50    |   2.00    |   3.00    |   4.00    |   5.00    |   6.00    |   7.00    |   8.00    |   9.00    |   10.0    |')
        print('_________________' * 7)
        for i in range(len(table)):
            if i!=-1:
                print(int_to_alg[i],end="")
                for k in range(12-len(int_to_alg[i])):
                    print(' ',end="")
                print('|', end="")
                for j in range(len(table[i])):
                    print(round(result[i][j], 2), end="")
                    for k in range(11-len(str(round(result[i][j], 2)))):
                        print(' ',end="")
                    print('|', end="")
                print()
        print('__________________' * 7)
        for cnt in range(len(table)):
            print(str(cnt+1)+' place Alg |', end="")
            for i in range(len(table[0])):
                col = table[:,i]
                print(int_to_alg[np.argsort(col)[cnt]], end="")
                # print(int_to_alg[np.argmin(col)], end="")
                for k in range(11-len(int_to_alg[np.argsort(col)[cnt]])):
                    print(' ',end="")
                print('|', end="")
            print()
        
        print('=====================================================',end='')
        for _ in range((28 - len(sys.argv[5]))//2):
            print(' ', end='')
        print(sys.argv[5], end='')
        for _ in range((28 - len(sys.argv[5]))//2):
            print(' ', end='')
        print('====================================================')
        print()

        if tableType==0:
            with open("./papers/DSDWA/results/"+sys.argv[5]+"_table0.txt", "w") as f:
                for i in range(len(table)):
                    if i!=-1:
                        f.write("$\Phi_{\\text{"+int_to_alg[i]+"}}$ ")
                        for j in range(len(table[i])):
                            if  j!=3 and j!=5 and j!=6 and j!=7 and j!=9:
                                tmp = str(round(result[i][j], 2))
                                firstPart = tmp.split(".")[0]
                                secondPart = tmp.split(".")[1]
                                f.write(" & ")
                                for k in range(len(firstPart)):
                                    if (len(firstPart)-k)%3 == 0 and k!=0:
                                        f.write(",")
                                    f.write(firstPart[k])
                                if len(secondPart):
                                    f.write(".")
                                f.write(secondPart)
                                f.write(" ")
                                # f.write("& " + str(round(result[i][j], 2))+" ")
                        f.write("\\\\ \n \n")
        
        elif tableType ==1:
            Weights = [1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
            xpoints = np.array(Weights)
            works = []
            for policy in int_to_alg:
                work_i = []
                for w in list(weight_to_num.keys()):
                    work_i.append(result[policy][weight_to_num[w]])
                works.append(work_i)
            
            with open("./papers/DSDWA/results/"+sys.argv[5]+"_table1.txt", "w") as f:
                for policy in int_to_alg:
                    f.write("$\Phi_{\\text{"+int_to_alg[policy]+"}}$ ")
                    ypoints = []
                    for i in works[policy]:
                        ypoints.append(i)
                    ypoints = np.array(ypoints)
                    for w in range(len(xpoints)):
                        if w!=3 and w!=5 and w!=6 and w!=7 and w!=9: ##This is index of w, not w itself
                            d = np.array(TheDataSet[policy][w])
                            avg = np.mean(d)
                            l, u = st.norm.interval(confidence=0.95, loc=np.mean(d), scale=st.sem(d))
                            confInt = u-l
                            tmp = str(round(avg, 2))
                            firstPart = tmp.split(".")[0]
                            secondPart = tmp.split(".")[1]
                            f.write(" & ")
                            for k in range(len(firstPart)):
                                if (len(firstPart)-k)%3 == 0 and k!=0:
                                    f.write(",")
                                f.write(firstPart[k])
                            if len(secondPart):
                                f.write(".")
                            f.write(secondPart)

                            tmp = str(round(confInt, 2))
                            firstPart = tmp.split(".")[0]
                            secondPart = tmp.split(".")[1]
                            f.write(" & ")
                            for k in range(len(firstPart)):
                                if (len(firstPart)-k)%3 == 0 and k!=0:
                                    f.write(",")
                                f.write(firstPart[k])
                            if len(secondPart):
                                f.write(".")
                            f.write(secondPart)

                        f.write(" ")
                    f.write("\\\\ \n \n")
                            
    ##############################################
    elif sys.argv[2] == '2':
        ##Experiment 2: creates a plot for one specific weight sorting the hardness of problems
        # x-axis is all the problems in decending order and y-axix is the umber of node expansions

        dataset = [{} for _ in range(int(sys.argv[3]))]
        problems_hardness = {}
        for i in range(1,101):
            problems_hardness[i] = 0

        with open("./papers/DSDWA/results.txt", "r") as f:
            for line in f:
                data = line.split()
                if data[0] == 'STP' and data[5] == '1.20':# and data[9]!='0':
                    problems_hardness[int(data[1])] += int(data[7])
                    dataset[int(data[3])][int(data[1])] = int(data[7])

        sorted_hardness = sorted(problems_hardness.items(), key=lambda x:x[1], reverse=True)
        sorted_hardness_dict = {}
        for i in (sorted_hardness):
            sorted_hardness_dict[i[0]] = i[1]

        print("HARDEST PROBLEM IS: ", sorted_hardness_dict)

        x_ax = list(sorted_hardness_dict.keys())
        x_axis = [str(i) for i in x_ax]

        y_axis = [[] for _ in range(int(sys.argv[3]))]

        for i in range(int(sys.argv[3])):
            for p in x_axis:
                if int(p) in dataset[i]:
                    y_axis[i].append(dataset[i][int(p)])
                else:
                    y_axis[i].append(0)

        # plt.plot(x_axis, y_axis[0], 'ko-', label='WA*')
        # plt.plot(x_axis, y_axis[1], 'g*-', label='XDP')
        # plt.plot(x_axis, y_axis[2], 'bs-', label='XUP')
        # plt.plot(x_axis, y_axis[3], 'rv-', label='HalfEdgeDrop')
        # plt.legend(loc='upper right')
        # plt.show()

    ##############################################
    elif sys.argv[2] == '3':
        ##Experiment 3: For each problem, creates a plot of applying different algorithms using different weights

        for problem in range(81, 83):
            dataset = [{} for _ in range(int(sys.argv[3]))]

            with open("./papers/DSDWA/results/STP-results.txt", "r") as f:
                for line in f:
                    data = line.split()
                    if data[0] == 'STP' and data[1] ==str(problem):# and data[9]!='0':
                        dataset[int(data[3])][float(data[5])] = int(data[7]) #dataset[1(xdp)][1.20] = 23455

            x_axis = list(dataset[0].keys())
            y_axis = [list(i.values()) for i in dataset]

            # for i in range(int(sys.argv[3])):
                # plt.plot(x_axis, y_axis[i], markers[i], label=int_to_alg[i])
                

            # plt.legend(loc='upper right')
            # plt.yscale('log')
            # plt.xlabel('weight')
            # plt.ylabel('Node Expansions')
            # plt.title('STP Problem '+str(problem))
            # plt.show()

    ##############################################
    elif sys.argv[2] == '4':

        ##Experiment 4: Creates the work/weight plot
        ##arg[3] is the number of algs and arg[4] is the number of weights
        table = np.zeros((int(sys.argv[3]), int(sys.argv[4])))
        count_table = np.zeros((int(sys.argv[3]), int(sys.argv[4])))
        TheDataSet = [[[] for _ in range(int(sys.argv[4]))] for _ in range(int(sys.argv[3]))]

        with open("./papers/DSDWA/results/"+sys.argv[5]+".txt", "r") as f:
            for line in f:
                data = line.split()
                if(len(data)): ## To check for empy lines
                    if data[0] == "STP" and (data[5] in list(weight_to_num.keys())) and (int(data[3]) in list(int_to_alg.keys())) : #and data[3]!='0' and data[3]!='1' and data[3]!='6' and data[3]!='8': # and data[9]!='0':
                        table[int(data[3])][weight_to_num[data[5]]] += float(data[7])
                        count_table[int(data[3])][weight_to_num[data[5]]] += 1
                        TheDataSet[int(data[3])][weight_to_num[data[5]]].append(int(data[7]))

        result = np.divide(table, count_table)

        # Weights = [1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        Weights = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]
        # Weights = [1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0]
        
        xpoints = np.array(Weights)
        works = []
        for policy in int_to_alg:
            work_i = []
            for w in list(weight_to_num.keys()):
                work_i.append(result[policy][weight_to_num[w]])
            works.append(work_i)
        for policy in int_to_alg:
            ypoints = []
            for i in works[policy]:
                ypoints.append(i)
            ypoints = np.array(ypoints)

            if showErrorBar:
                for w in range(len(xpoints)):
                    d = np.array(TheDataSet[policy][w])
                    l, u = st.norm.interval(confidence=0.95, loc=np.mean(d), scale=st.sem(d))
                    plt.errorbar(x=xpoints[w], y=ypoints[w], yerr=[[np.mean(d)-l], [u-np.mean(d)]], color=colours[policy], capsize=8)

            if int_to_alg[policy]!='DWP' and int_to_alg[policy]!='MAP':
                plt.plot(xpoints, ypoints, linestyle=linestyles[list(linestyles.keys())[5]], linewidth = '2.5', marker=markers[policy], ms='10', markeredgecolor="k", label=int_to_alg[policy], color=colours[policy], alpha=0.5)
            else:
                plt.plot(xpoints, ypoints, linestyle=linestyles[list(linestyles.keys())[5]], linewidth = '2.5', marker=markers[policy], ms='10', markeredgecolor="k", label=int_to_alg[policy], color=colours[policy])

        font = {'family':'serif','color':'darkred','size':12}
        plt.ylabel("Work", fontdict=font)
        plt.xlabel("Weights", fontdict=font)
        plt.title(mapType[int(sys.argv[7])]+", Size="+str(sys.argv[6]))
        # plt.legend(fontsize="26")
        plt.xticks(xpoints) 
        plt.yscale('log')
        plt.grid(axis='y', color='0.80', which='major')
        # plt.savefig(sys.argv[5]+"E4.pdf", format="pdf", bbox_inches="tight")
        plt.show()
    
    ##############################################
    elif sys.argv[2] == '5':

        ##Experiment 5: Creates the solutionQuality/weight plot
        ##arg[3] is the number of algs and arg[4] is the number of weights
        table = np.zeros((int(sys.argv[3]), int(sys.argv[4])))
        count_table = np.zeros((int(sys.argv[3]), int(sys.argv[4])))
        TheDataSet = [[[] for _ in range(int(sys.argv[4]))] for _ in range(int(sys.argv[3]))]

        with open("./papers/DSDWA/results/"+sys.argv[5]+".txt", "r") as f:
            for line in f:
                data = line.split()
                if(len(data)): ## To check for empy lines
                    if data[0] == "STP" and (data[5] in list(weight_to_num.keys())) and (int(data[3]) in list(int_to_alg.keys())) : #and data[3]!='0' and data[3]!='1' and data[3]!='6' and data[3]!='8': # and data[9]!='0':
                        table[int(data[3])][weight_to_num[data[5]]] += float(data[9])
                        count_table[int(data[3])][weight_to_num[data[5]]] += 1
                        TheDataSet[int(data[3])][weight_to_num[data[5]]].append(int(data[9]))

        result = np.divide(table, count_table)

        # Weights = [1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
        Weights = [1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        
        xpoints = np.array(Weights)
        works = []
        for policy in int_to_alg:
            work_i = []
            for w in list(weight_to_num.keys()):
                work_i.append(result[policy][weight_to_num[w]])
            works.append(work_i)
        for policy in int_to_alg:
            ypoints = []
            for i in works[policy]:
                ypoints.append(i)
            ypoints = np.array(ypoints)

            if showErrorBar:
                for w in range(len(xpoints)):
                    d = np.array(TheDataSet[policy][w])
                    l, u = st.norm.interval(confidence=0.95, loc=np.mean(d), scale=st.sem(d))
                    # plt.errorbar(x=xpoints[w], y=ypoints[w], yerr=(u-l)/2, color=colours[policy])
                    plt.errorbar(x=xpoints[w], y=ypoints[w], yerr=[[np.mean(d)-l], [u-np.mean(d)]], color=colours[policy], capsize=8)

            if policy!=5:
                plt.plot(xpoints, ypoints, linestyle=linestyles[list(linestyles.keys())[5]], linewidth = '2.5', marker=markers[policy], ms='10', markeredgecolor="k", label=int_to_alg[policy], color=colours[policy], alpha=0.5)
            else:
                plt.plot(xpoints, ypoints, linestyle=linestyles[list(linestyles.keys())[5]], linewidth = '2.5', marker=markers[policy], ms='10', markeredgecolor="k", label=int_to_alg[policy], color=colours[policy])

        font = {'family':'serif','color':'darkred','size':12}
        plt.ylabel("Work", fontdict=font)
        plt.xlabel("Weights", fontdict=font)
        plt.title(mapType[int(sys.argv[7])]+", Size="+str(sys.argv[6]))
        plt.legend(fontsize="26")
        plt.xticks(xpoints) 
        # plt.yscale('log')
        plt.savefig(sys.argv[5]+"E5.pdf", format="pdf", bbox_inches="tight")
        plt.show()
    
############################################################################################
elif sys.argv[1] == '-map':
    if sys.argv[2] == '1':

        ##Experiment 1: Creates a table from average of runs of different problems over different weights
        table = np.zeros((int(sys.argv[3]), int(sys.argv[4])))
        count_table = np.zeros((int(sys.argv[3]), int(sys.argv[4])))
        TheDataSet = [[[] for _ in range(int(sys.argv[4]))] for _ in range(int(sys.argv[3]))]

        with open("./papers/DSDWA/results/"+sys.argv[5]+".txt", "r") as f:
            for line in f:
                data = line.split()
                if(len(data)): ## To check for empy lines
                    if data[0] == "MAP" and (data[7] in list(weight_to_num.keys())) and (int(data[5]) in list(int_to_alg.keys())): #and data[5]!='0' and data[5]!='1' and data[5]!='6' and data[5]!='8':# and data[9]!='0':
                        table[int(data[5])][weight_to_num[data[7]]] += int(data[9])
                        count_table[int(data[5])][weight_to_num[data[7]]] += 1
                        TheDataSet[int(data[5])][weight_to_num[data[7]]].append(int(data[9]))

        result = np.divide(table, count_table)
        print()
        print('==========================================  Average Expansions Table  ===========================================')
        print('Alg / Weight|  1.50   |  2.00   |  3.00   |  4.00   |  5.00   |  6.00   |  7.00   |  8.00   |  9.00   |  10.0   |')
        print('_________________' * 7)
        for i in range(len(table)):
            if i!=-1:
                print(int_to_alg[i],end="")
                for k in range(12-len(int_to_alg[i])):
                    print(' ',end="")
                print('|', end="")
                for j in range(len(table[i])):
                    print(round(result[i][j], 2), end="")
                    for k in range(9-len(str(round(result[i][j], 2)))):
                        print(' ',end="")
                    print('|', end="")
                print()
        print('__________________' * 7)
                    
        for cnt in range(len(table)):
            print(str(cnt+1)+' place Alg |', end="")
            for i in range(len(table[0])):
                col = table[:,i]
                print(int_to_alg[np.argsort(col)[cnt]], end="")
                # print(int_to_alg[np.argmin(col)], end="")
                for k in range(9-len(int_to_alg[np.argsort(col)[cnt]])):
                    print(' ',end="")
                print('|', end="")
            print()
        
        print('=====================================================',end='')
        for _ in range((26 - len(sys.argv[5]))//2):
            print(' ', end='')
        print(sys.argv[5], end='')
        for _ in range((26 - len(sys.argv[5]))//2):
            print(' ', end='')
        print('====================================================')
        print()

        if tableType==0:
            with open("./papers/DSDWA/results/"+sys.argv[5]+"_table0.txt", "w") as f:
                for i in range(len(table)):
                    if i!=-1:
                        f.write("$\Phi_{\\text{"+int_to_alg[i]+"}}$ ")
                        for j in range(len(table[i])):
                            if  j!=3 and j!=5 and j!=6 and j!=7 and j!=9:
                                tmp = str(round(result[i][j], 2))
                                firstPart = tmp.split(".")[0]
                                secondPart = tmp.split(".")[1]
                                f.write(" & ")
                                for k in range(len(firstPart)):
                                    if (len(firstPart)-k)%3 == 0 and k!=0:
                                        f.write(",")
                                    f.write(firstPart[k])
                                if len(secondPart):
                                    f.write(".")
                                f.write(secondPart)
                                f.write(" ")
                                # f.write("& " + str(round(result[i][j], 2))+" ")
                        f.write("\\\\ \n \n")
        
        elif tableType ==1:
            Weights = [1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
            xpoints = np.array(Weights)
            works = []
            for policy in int_to_alg:
                work_i = []
                for w in list(weight_to_num.keys()):
                    work_i.append(result[policy][weight_to_num[w]])
                works.append(work_i)
            
            with open("./papers/DSDWA/results/"+sys.argv[5]+"_table1.txt", "w") as f:
                for policy in int_to_alg:
                    f.write("$\Phi_{\\text{"+int_to_alg[policy]+"}}$ ")
                    ypoints = []
                    for i in works[policy]:
                        ypoints.append(i)
                    ypoints = np.array(ypoints)
                    for w in range(len(xpoints)):
                        if w!=3 and w!=5 and w!=6 and w!=7 and w!=9: ##This is index of w, not w itself
                            d = np.array(TheDataSet[policy][w])
                            avg = np.mean(d)
                            l, u = st.norm.interval(confidence=0.95, loc=np.mean(d), scale=st.sem(d))
                            confInt = u-l
                            tmp = str(round(avg, 2))
                            firstPart = tmp.split(".")[0]
                            secondPart = tmp.split(".")[1]
                            f.write(" & ")
                            for k in range(len(firstPart)):
                                if (len(firstPart)-k)%3 == 0 and k!=0:
                                    f.write(",")
                                f.write(firstPart[k])
                            if len(secondPart):
                                f.write(".")
                            f.write(secondPart)

                            tmp = str(round(confInt, 2))
                            firstPart = tmp.split(".")[0]
                            secondPart = tmp.split(".")[1]
                            f.write(" & ")
                            for k in range(len(firstPart)):
                                if (len(firstPart)-k)%3 == 0 and k!=0:
                                    f.write(",")
                                f.write(firstPart[k])
                            if len(secondPart):
                                f.write(".")
                            f.write(secondPart)

                        f.write(" ")
                    f.write("\\\\ \n \n")

    ##############################################
    elif sys.argv[2] == '2':
        # Second idea: create plots, each plot for one weight, 
        # x-axis is all the problems in decending order and y-axix is the umber of node expansions

        dataset = [{} for _ in range(int(sys.argv[3]))]
        problems_hardness = {}
        for i in range(1,101):
            problems_hardness[i] = 0

        with open("./papers/DSDWA/results.txt", "r") as f:
            for line in f:
                data = line.split()
                if data[0] == 'MAP' and data[5] == '1.20':# and data[9]!='0':
                    problems_hardness[int(data[1])] += int(data[7])
                    dataset[int(data[3])][int(data[1])] = int(data[7])

        sorted_hardness = sorted(problems_hardness.items(), key=lambda x:x[1], reverse=True)
        sorted_hardness_dict = {}
        for i in (sorted_hardness):
            sorted_hardness_dict[i[0]] = i[1]

        print("HARDEST PROBLEM IS: ", sorted_hardness_dict)

        x_ax = list(sorted_hardness_dict.keys())
        x_axis = [str(i) for i in x_ax]

        y_axis = [[] for _ in range(int(sys.argv[3]))]

        # print(dataset[0])
        for i in range(int(sys.argv[3])):
            for p in x_axis:
                if int(p) in dataset[i]:
                    # print(dataset[i][p])
                    y_axis[i].append(dataset[i][int(p)])
                else:
                    y_axis[i].append(0)

        # print(y_axis[0])
        # print(x_axis)
        # plt.plot(x_axis, y_axis[0], 'ko-', label='WA*')
        # plt.plot(x_axis, y_axis[1], 'g*-', label='XDP')
        # plt.plot(x_axis, y_axis[2], 'bs-', label='XUP')
        # plt.plot(x_axis, y_axis[3], 'rv-', label='HalfEdgeDrop')
        # plt.legend(loc='upper right')
        # plt.show()

    ##############################################
    elif sys.argv[2] == '3':
        ##Experiment 3: For each problem, creates a plot of applying different algorithms using different weights

        problemsList = []
        with open("./papers/DSDWA/ALL-random-40-results.txt", "r") as f:
            for line in f:
                data = line.split()
                if data[0] == 'MAP' and (data[1] not in problemsList):
                    problemsList.append(data[1])
        f.close()
        
        for problem in problemsList:
            numberOfScenarios = 0
            with open("./papers/DSDWA/ALL-random-40-results.txt", "r") as f:
                dataset = [{} for _ in range(int(sys.argv[3]))]

                for line in f:
                    data = line.split()
                    if data[0] == 'MAP' and data[1] == problem :# and data[11]!='0':
                        if float(data[7]) not in list(dataset[int(data[5])].keys()):
                            dataset[int(data[5])][float(data[7])] = 0

                        dataset[int(data[5])][float(data[7])] += int(data[9]) #dataset[1(xdp)][1.20] = 23455

                        numberOfScenarios += 1
                        

                x_axis = list(dataset[0].keys())
                y_axis = [list(alg.values()) for alg in dataset]

                # print(y_axis)
                # print(y_axis)
                numberOfScenarios /= int(sys.argv[3])
                numberOfScenarios /= int(sys.argv[4])
                print(numberOfScenarios)
                for alg in range(len(y_axis)):
                    for i in range(len(y_axis[alg])):
                        y_axis[alg][i] = y_axis[alg][i]/numberOfScenarios
                # print(y_axis)

                # for i in range(int(sys.argv[3])):
                    # plt.plot(x_axis, y_axis[i], markers[i], label=int_to_alg[i])
                    

                # plt.legend(loc='upper right')
                # plt.yscale('log')
                # plt.xlabel('weight')
                # plt.ylabel('Node Expansions')
                # plt.title('MAP Problem '+str(problem))
                # plt.show()
    
    ##############################################
    elif sys.argv[2] == '4':

        ##Experiment 4: Creates the work/weight plot
        ##arg[3] is the number of algs and arg[4] is the number of weights
        table = np.zeros((int(sys.argv[3]), int(sys.argv[4])))
        count_table = np.zeros((int(sys.argv[3]), int(sys.argv[4])))
        TheDataSet = [[[] for _ in range(int(sys.argv[4]))] for _ in range(int(sys.argv[3]))]

        with open("./papers/DSDWA/results/"+sys.argv[5]+".txt", "r") as f:
            for line in f:
                data = line.split()
                if(len(data)): ## To check for empy lines
                    if data[0] == "MAP" and (data[7] in list(weight_to_num.keys())) and (int(data[5]) in list(int_to_alg.keys())) : #and data[5]!='0' and data[5]!='1' and data[5]!='6' and data[5]!='8':# and data[9]!='0':
                        table[int(data[5])][weight_to_num[data[7]]] += int(data[9])
                        count_table[int(data[5])][weight_to_num[data[7]]] += 1
                        TheDataSet[int(data[5])][weight_to_num[data[7]]].append(int(data[9]))

        result = np.divide(table, count_table)

        Weights = [1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        # Weights = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]
        # Weights = [1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0]
        
        xpoints = np.array(Weights)
        works = []
        for policy in int_to_alg:
            work_i = []
            for w in list(weight_to_num.keys()):
                work_i.append(result[policy][weight_to_num[w]])
            works.append(work_i)
        for policy in int_to_alg:
            ypoints = []
            for i in works[policy]:
                ypoints.append(i)
            ypoints = np.array(ypoints)

            if int_to_alg[policy]!='DWP' and int_to_alg[policy]!='MAP':
                plt.plot(xpoints, ypoints, linestyle=linestyles[list(linestyles.keys())[5]], linewidth = '2.5', marker=markers[policy], ms='10', markeredgecolor="k", label=int_to_alg[policy], color=colours[policy], alpha=0.5)
            else:
                plt.plot(xpoints, ypoints, linestyle=linestyles[list(linestyles.keys())[5]], linewidth = '2.5', marker=markers[policy], ms='10', markeredgecolor="k", label=int_to_alg[policy], color=colours[policy])
            
            if showErrorBar:
                for w in range(len(xpoints)):
                    d = np.array(TheDataSet[policy][w])
                    l, u = st.norm.interval(confidence=0.95, loc=np.mean(d), scale=st.sem(d))
                    plt.errorbar(x=xpoints[w], y=ypoints[w], yerr=[[np.mean(d)-l], [u-np.mean(d)]], color=colours[policy], capsize=8)
            
        font = {'family':'serif','color':'darkred','size':12}
        plt.ylabel("Work", fontdict=font)
        plt.xlabel("Weights", fontdict=font)
        plt.title(sys.argv[5])
        # plt.legend(fontsize="25")
        plt.xticks(xpoints) 
        plt.yscale('log')
        # plt.xscale('log')
        plt.grid(axis='y', color='0.80', which='major')
        # plt.savefig(sys.argv[5]+"E4.pdf", format="pdf", bbox_inches="tight")
        plt.show()
    
    ##############################################
    elif sys.argv[2] == '5':
        ##Experiment 5: Creates the solutionQuality/weight plot
        ##arg[3] is the number of algs and arg[4] is the number of weights
        table = np.zeros((int(sys.argv[3]), int(sys.argv[4])))
        count_table = np.zeros((int(sys.argv[3]), int(sys.argv[4])))
        TheDataSet = [[[] for _ in range(int(sys.argv[4]))] for _ in range(int(sys.argv[3]))]

        with open("./papers/DSDWA/results/"+sys.argv[5]+".txt", "r") as f:
            for line in f:
                data = line.split()
                if(len(data)): ## To check for empy lines
                    if data[0] == "MAP" and (data[7] in list(weight_to_num.keys())) and (int(data[5]) in list(int_to_alg.keys())) : #and data[5]!='0' and data[5]!='1' and data[5]!='6' and data[5]!='8':# and data[9]!='0':
                        table[int(data[5])][weight_to_num[data[7]]] += float(data[11])
                        count_table[int(data[5])][weight_to_num[data[7]]] += 1
                        TheDataSet[int(data[5])][weight_to_num[data[7]]].append(float(data[11]))

        result = np.divide(table, count_table)
        # Weights = [1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        Weights = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]

        xpoints = np.array(Weights)
        works = []
        for policy in int_to_alg:
            work_i = []
            for w in list(weight_to_num.keys()):
                work_i.append(result[policy][weight_to_num[w]])
            works.append(work_i)
        for policy in int_to_alg:
            ypoints = []
            for i in works[policy]:
                ypoints.append(i)
            ypoints = np.array(ypoints)

            if policy!=5:
                plt.plot(xpoints, ypoints, linestyle=linestyles[list(linestyles.keys())[5]], linewidth = '2.5', marker=markers[policy], ms='10', markeredgecolor="k", label=int_to_alg[policy], color=colours[policy], alpha=0.5)
            else:
                plt.plot(xpoints, ypoints, linestyle=linestyles[list(linestyles.keys())[5]], linewidth = '2.5', marker=markers[policy], ms='10', markeredgecolor="k", label=int_to_alg[policy], color=colours[policy])
            
            if showErrorBar:
                for w in range(len(xpoints)):
                    d = np.array(TheDataSet[policy][w])
                    l, u = st.norm.interval(confidence=0.95, loc=np.mean(d), scale=st.sem(d))
                    # plt.errorbar(x=xpoints[w], y=ypoints[w], yerr=(u-l)/2, color=colours[policy])
                    plt.errorbar(x=xpoints[w], y=ypoints[w], yerr=[[np.mean(d)-l], [u-np.mean(d)]], color=colours[7])
            
        font = {'family':'serif','color':'darkred','size':12}
        plt.ylabel("Work", fontdict=font)
        plt.xlabel("Weights", fontdict=font)
        plt.title(mapType[int(sys.argv[7])]+", Size="+str(sys.argv[6]))
        plt.legend(fontsize="25")
        plt.xticks(xpoints) 
        # plt.yscale('log')
        # plt.xscale('log')
        plt.savefig(sys.argv[5]+"E5.pdf", format="pdf", bbox_inches="tight")
        plt.show()
    ##############################################
    elif sys.argv[2] == '6':
        
        ##Experiment 6.0: Creates the work/weight plot
        ##arg[3] is the number of algs and arg[4] is the number of weights
        table = np.zeros((int(sys.argv[3]), int(sys.argv[4])))
        count_table = np.zeros((int(sys.argv[3]), int(sys.argv[4])))
        TheDataSet = [[[] for _ in range(int(sys.argv[4]))] for _ in range(int(sys.argv[3]))]

        with open("./papers/DSDWA/results/"+sys.argv[5]+".txt", "r") as f:
            for line in f:
                data = line.split()
                if(len(data)): ## To check for empy lines
                    if data[0] == "Time" and (data[9] in list(weight_to_num.keys())) and (int(data[7]) in list(int_to_alg.keys())) : #and data[5]!='0' and data[5]!='1' and data[5]!='6' and data[5]!='8':# and data[9]!='0':
                        table[int(data[7])][weight_to_num[data[9]]] += int(data[11])
                        count_table[int(data[7])][weight_to_num[data[9]]] += 1
                        TheDataSet[int(data[7])][weight_to_num[data[9]]].append(int(data[11]))

        result = np.divide(table, count_table)

        Weights = [1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        # Weights = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]
        # Weights = [1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0]
        
        xpoints = np.array(Weights)
        works = []
        for policy in int_to_alg:
            work_i = []
            for w in list(weight_to_num.keys()):
                work_i.append(result[policy][weight_to_num[w]])
            works.append(work_i)
        for policy in int_to_alg:
            ypoints = []
            for i in works[policy]:
                ypoints.append(i)
            ypoints = np.array(ypoints)

            if int_to_alg[policy]!='DWP' and int_to_alg[policy]!='MAP':
                plt.plot(xpoints, ypoints, linestyle=linestyles[list(linestyles.keys())[5]], linewidth = '2.5', marker=markers[policy], ms='10', markeredgecolor="k", label=int_to_alg[policy], color=colours[policy], alpha=0.5)
            else:
                plt.plot(xpoints, ypoints, linestyle=linestyles[list(linestyles.keys())[5]], linewidth = '2.5', marker=markers[policy], ms='10', markeredgecolor="k", label=int_to_alg[policy], color=colours[policy])
            
            if showErrorBar:
                for w in range(len(xpoints)):
                    d = np.array(TheDataSet[policy][w])
                    l, u = st.norm.interval(confidence=0.95, loc=np.mean(d), scale=st.sem(d))
                    plt.errorbar(x=xpoints[w], y=ypoints[w], yerr=[[np.mean(d)-l], [u-np.mean(d)]], color=colours[policy], capsize=8)
            
        font = {'family':'serif','color':'darkred','size':12}
        plt.ylabel("Work", fontdict=font)
        plt.xlabel("Weights", fontdict=font)
        plt.title(sys.argv[5])
        # plt.legend(fontsize="25")
        plt.xticks(xpoints) 
        plt.yscale('log')
        # plt.xscale('log')
        plt.grid(axis='y', color='0.80', which='major')
        # plt.savefig(sys.argv[5]+"E4.pdf", format="pdf", bbox_inches="tight")
        plt.show()


        ##Experiment 6.1: Creates the RunTimePerNode/weight plot
        ##arg[3] is the number of algs and arg[4] is the number of weights
        table = np.zeros((int(sys.argv[3]), int(sys.argv[4])))
        count_table = np.zeros((int(sys.argv[3]), int(sys.argv[4])))
        TheDataSet = [[[] for _ in range(int(sys.argv[4]))] for _ in range(int(sys.argv[3]))]

        with open("./papers/DSDWA/results/"+sys.argv[5]+".txt", "r") as f:
            for line in f:
                data = line.split()
                if(len(data)): ## To check for empy lines
                    if data[0] == "Time" and (data[9] in list(weight_to_num.keys())) and (int(data[7]) in list(int_to_alg.keys())) : #and data[5]!='0' and data[5]!='1' and data[5]!='6' and data[5]!='8':# and data[9]!='0':
                        table[int(data[7])][weight_to_num[data[9]]] += float(data[13])
                        count_table[int(data[7])][weight_to_num[data[9]]] += 1
                        TheDataSet[int(data[7])][weight_to_num[data[9]]].append(float(data[13]))

        result = np.divide(table, count_table)

        Weights = [1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        # Weights = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]
        # Weights = [1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0]
        
        xpoints = np.array(Weights)
        works = []
        for policy in int_to_alg:
            work_i = []
            for w in list(weight_to_num.keys()):
                work_i.append(result[policy][weight_to_num[w]])
            works.append(work_i)
        for policy in int_to_alg:
            ypoints = []
            for i in works[policy]:
                ypoints.append(i)
            ypoints = np.array(ypoints)

            if int_to_alg[policy]!='DWP' and int_to_alg[policy]!='MAP':
                plt.plot(xpoints, ypoints, linestyle=linestyles[list(linestyles.keys())[5]], linewidth = '2.5', marker=markers[policy], ms='10', markeredgecolor="k", label=int_to_alg[policy], color=colours[policy], alpha=0.5)
            else:
                plt.plot(xpoints, ypoints, linestyle=linestyles[list(linestyles.keys())[5]], linewidth = '2.5', marker=markers[policy], ms='10', markeredgecolor="k", label=int_to_alg[policy], color=colours[policy])
            
            if showErrorBar:
                for w in range(len(xpoints)):
                    d = np.array(TheDataSet[policy][w])
                    l, u = st.norm.interval(confidence=0.95, loc=np.mean(d), scale=st.sem(d))
                    plt.errorbar(x=xpoints[w], y=ypoints[w], yerr=[[np.mean(d)-l], [u-np.mean(d)]], color=colours[policy], capsize=8)
            
        font = {'family':'serif','color':'darkred','size':12}
        plt.ylabel("RunTime per Node", fontdict=font)
        plt.xlabel("Weights", fontdict=font)
        plt.title(sys.argv[5])
        # plt.legend(fontsize="25")
        plt.xticks(xpoints) 
        plt.yscale('log')
        # plt.xscale('log')
        plt.grid(axis='y', color='0.80', which='major')
        # plt.savefig(sys.argv[5]+"E4.pdf", format="pdf", bbox_inches="tight")
        plt.show()

        ##Experiment 6.2: Creates the TotalRunTime/weight plot
        ##arg[3] is the number of algs and arg[4] is the number of weights
        table = np.zeros((int(sys.argv[3]), int(sys.argv[4])))
        count_table = np.zeros((int(sys.argv[3]), int(sys.argv[4])))
        TheDataSet = [[[] for _ in range(int(sys.argv[4]))] for _ in range(int(sys.argv[3]))]

        with open("./papers/DSDWA/results/"+sys.argv[5]+".txt", "r") as f:
            for line in f:
                data = line.split()
                if(len(data)): ## To check for empy lines
                    if data[0] == "Time" and (data[9] in list(weight_to_num.keys())) and (int(data[7]) in list(int_to_alg.keys())) : #and data[5]!='0' and data[5]!='1' and data[5]!='6' and data[5]!='8':# and data[9]!='0':
                        table[int(data[7])][weight_to_num[data[9]]] += float(data[15])
                        count_table[int(data[7])][weight_to_num[data[9]]] += 1
                        TheDataSet[int(data[7])][weight_to_num[data[9]]].append(float(data[15]))

        result = np.divide(table, count_table)

        Weights = [1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        # Weights = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]
        # Weights = [1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0]
        
        xpoints = np.array(Weights)
        works = []
        for policy in int_to_alg:
            work_i = []
            for w in list(weight_to_num.keys()):
                work_i.append(result[policy][weight_to_num[w]])
            works.append(work_i)
        for policy in int_to_alg:
            ypoints = []
            for i in works[policy]:
                ypoints.append(i)
            ypoints = np.array(ypoints)

            if int_to_alg[policy]!='DWP' and int_to_alg[policy]!='MAP':
                plt.plot(xpoints, ypoints, linestyle=linestyles[list(linestyles.keys())[5]], linewidth = '2.5', marker=markers[policy], ms='10', markeredgecolor="k", label=int_to_alg[policy], color=colours[policy], alpha=0.5)
            else:
                plt.plot(xpoints, ypoints, linestyle=linestyles[list(linestyles.keys())[5]], linewidth = '2.5', marker=markers[policy], ms='10', markeredgecolor="k", label=int_to_alg[policy], color=colours[policy])
            
            if showErrorBar:
                for w in range(len(xpoints)):
                    d = np.array(TheDataSet[policy][w])
                    l, u = st.norm.interval(confidence=0.95, loc=np.mean(d), scale=st.sem(d))
                    plt.errorbar(x=xpoints[w], y=ypoints[w], yerr=[[np.mean(d)-l], [u-np.mean(d)]], color=colours[policy], capsize=8)
            
        font = {'family':'serif','color':'darkred','size':12}
        plt.ylabel("TotalRunTime", fontdict=font)
        plt.xlabel("Weights", fontdict=font)
        plt.title(sys.argv[5])
        # plt.legend(fontsize="25")
        plt.xticks(xpoints) 
        plt.yscale('log')
        # plt.xscale('log')
        plt.grid(axis='y', color='0.80', which='major')
        # plt.savefig(sys.argv[5]+"E4.pdf", format="pdf", bbox_inches="tight")
        plt.show()