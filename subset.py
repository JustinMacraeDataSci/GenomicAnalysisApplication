# -*- coding: utf-8 -*-
"""
Created on Thu May 31 12:50:54 2018

@author: Justin Macrae
"""
from tkinter.filedialog import askopenfilename
import numpy as np
import matplotlib.pyplot as plt


def takefourth(elem):
    return elem[3]
    
def pta1v2():
    """
    This function takes in two text files. One being a textfile of bedgraph data.
    The other being a textfile of strands of MRNA.
    The output is a textfile a1.txt containing the MRNA strand info, the number
    of data points it intersected with in the bedgraph data and finally the average
    score across all intersected points.
    """
   # f = open("C:/Users/Justin Macrae/Downloads/rop1/norm.bedgraph",'r')
   #  fmrna = open("C:/Users/Justin Macrae/Downloads/rop1/mrnasub.txt",'r')
   
    filename1 = askopenfilename(title="Select a bedgraph file to parse!")
    filename2 = askopenfilename(title="Select a mrna file to parse!")
    
    cnt=0 #This is used to keep track of the orginal position in the text file.
    fmrna=open(filename2,'r')
    f=open(filename1,'r')
    valuel=[] #An array for holding score values.
    rangel=[] #An array for holding region boundaries.
    for line in f:
        if "scaffold" not in line:
            linesplit= line.split()
            grange= float(linesplit[1])
            totalvalue=float(linesplit[3])
            valuel.append((grange,totalvalue))
    print("done step 1: Aquire values from Data")
    sortedvalue=sorted(valuel,key=lambda x:x[0])
    print("done step 2: Sort values in Data")
        
    for line in fmrna:
        linesplit2= line.split()
        mrangeh= float(linesplit2[4])
        mrangel= float(linesplit2[3])
        minfo= linesplit2[8]
        rangel.append((mrangel,mrangeh,minfo,cnt))
        cnt+=1
    print("done step 3:Aquire values from mRNA")
    sortedrange=sorted(rangel,key=lambda x:x[0])
    print("done step 4: Sort values for mRNA")
    f.close()
    fmrna.close()
    #############################################
    #this section of code has now sorted each lists by range
    #we can access lower bounds by rangel[n][0] and upper rangel[n][1]
    #we can use the sorted ranges to do comparisons
    wf = open("a1.txt", 'w')
    
    cnt=0 #cnt, counts the number of hits on bedgraph data
    
    totalval=0 #variable to hold total score across region
    
    g1=[] #array to hold information that we will graph on a scatter plot
    
    gnum=0 #a way to count the index of our data, I chose to plot one every 50 points to speed up graphing
    
    save=0 #this variable saves the last cut off point for each region, look in presentation under "Time Complexity"
    
    firstflag=True #This is a boolean flag to write into our file the legend of our text file
    
    for n in sortedrange:
        minvalue=n[0]
        maxvalue=n[1]
        info=n[2]
        #we have here the min and max values of the ranges
        for n2 in range(save,len(sortedvalue)):
            if sortedvalue[n2][0]<= maxvalue and sortedvalue[n2][0]>= minvalue:
                totalval+=sortedvalue[n2][1]
                cnt+=1
            elif sortedvalue[n2][0]>maxvalue:
                save=n2
                break
        if cnt!=0:
            if firstflag:
                wf.write("MRNA INFO"+"    "+"# OF DATA POINTS HIT"+"    "+"AVERAGE SCORE"+"    "+"ID"+"    "+"DIVISION NUMBER"+"    "+"DIRECTION"+"\n")
                firstflag=False
            wf.write(info+" "+str(cnt)+" "+str(int(totalval/cnt))+" "+ str(n[3]) +"\n")
            if gnum==50:
                g1.append((int(n[3]),int(totalval/cnt)))
                gnum=0
            gnum+=1
            
        cnt=0
        totalval=0
        
    print("All done! Check your directory for a1.txt!")
    print("Creating Graph")
    x,y = zip(*g1)
    plt.scatter(x,y)
    plt.xlabel("mRNA position")
    plt.ylabel("Score of mRNA")
    plt.show()
        
def threeandfive():
    """
    This function is used to create 2 seperate text files for a given input file.
    The 2 text files filter hits x<3 and x<5 respectively.
    """
    filename = askopenfilename(title="Select a file to parse!")
    wf = open(filename,'r')
    name3= input("Please type the name for the file that parses by greater than 3:    ")
    name5= input("Please type the name for the file that parses by greater than 5:    ")
    f3 = open(name3+".txt",'w')
    f5= open(name5+".txt",'w')
    for line in wf:
        splitline= line.split()
        if int(splitline[1])>5:
            f5.write(line)
        if int(splitline[1])>3:
            f3.write(line)
    print("Subset done!")
    
    
        

def ptamm():
    """
    Gets the mean median and mode, for part A.
    """
    wf = open("C:/Users/Justin Macrae/Downloads/rop1/a1.txt", 'r')
    wf2 = open("C:/Users/Justin Macrae/Downloads/rop1/a2five.txt", 'r')
    wf3 = open("C:/Users/Justin Macrae/Downloads/rop1/a2three.txt", 'r')
    total=0
    total1=0
    total2=0
    l1=[]
    l3=[]
    l5=[]
    for line in wf:
        linesplit=line.split()
        total+=float(linesplit[2])
        l1.append(float(linesplit[2]))
    print("Mean: "+ str(total/11132))
    sortedl1= sorted(l1)
    print("Median: " + str(sortedl1[11133//2]))
    for line in wf2:
        linesplit=line.split()
        total1+=float(linesplit[2])
        l5.append(float(linesplit[2]))
    print("Mean for fives: "+ str(total1/10451))
    sortedl5= sorted(l5)
    print("Median: " + str(sortedl5[10452//2]))

    for line in wf3:
        linesplit=line.split()
        total2+=float(linesplit[2])
        l3.append(float(linesplit[2]))
    print("Mean for threes: "+ str(total2/10681))
    sortedl3=sorted(l3)
    print("Median: " + str(sortedl3[10682//2]))
    
def ptb():
    """
    This function takes two text files as input.
    One being our bedgraph data file, our other being our subset data file for MRNA and exon info.
    The output is a textfile b1weightedaverage.txt, which contains the average across all exons in the MRNA.
    A text file is also created for use, titled b1.txt which contains the average scores for each exon in the MRNA.
    After taking each average score it adds them and divides them by the n scores.
    """
    #f = open("C:/Users/Justin Macrae/Downloads/rop1/norm.bedgraph",'r')
    #fmrna = open("C:/Users/Justin Macrae/Downloads/rop1/subsetdata",'r')
    
    filename1 = askopenfilename(title="Select a bedgraph file to parse!")
    filename2 = askopenfilename(title="Select a file containing both mrna and exons to parse!")
    
    f=open(filename1,'r')
    fmrna= open(filename2,'r')
    
    valuel=[] # An array for holding score values.
    rangel=[] # An array for holding region boundaries.
    
    placement=0   # A variable to hold the section placement (ie: 3 sections could have placement=1 for first section)
    for line in f:
        linesplit= line.split()
        grange= int(linesplit[1])
        totalvalue=float(linesplit[3])
        valuel.append((grange,totalvalue))
    print("done step 1: Aquire values from Data")
    sortedvalue=sorted(valuel,key=lambda x:x[0])
    print("done step 2: Sort values in Data")
    ##########################
    #Here we are sorting the bedgraph data
    for line in fmrna:
        linesplit2= line.split()
        if "longest=1" in line and "scaffold" not in line:
            minfo= linesplit2[8]
            mplace=placement
            placement+=1
        elif "exon" in line and "scaffold" not in line:
            mrangeh= float(linesplit2[4])
            mrangel= float(linesplit2[3])
            rangel.append((mrangel,mrangeh,minfo,mplace))
            #Here I have it as exon lower,upper,mrna info, mrna placement so I can re-sort after.
    print("done step 3:Aquire values from Exons")
    sortedrange=sorted(rangel,key=lambda x:x[1])
    print("done step 4: Sort values for Exons")
    f.close()
    fmrna.close()
    ########################################
    wf = open("b1.txt", 'w')
    cnt=0 #Keep track of # of hits
    totalval=0 #Total score value per section
    save=0 #Look at presentation "Time Complexity" or function A
    listofzero=[] #Array to keep track of exons which had zero hits.
    
    for n in sortedrange:
        minvalue=n[0]
        maxvalue=n[1]
        info=n[2]
        #we have here the min and max values of the ranges
        for n2 in range(save,len(sortedvalue)):
            if sortedvalue[n2][0]<= maxvalue and sortedvalue[n2][0]>= minvalue:
                totalval+=sortedvalue[n2][1]
                cnt+=1
            elif sortedvalue[n2][0]>maxvalue:
                save=n2
                break
        if cnt==0:
            listofzero.append(n[3])
        if cnt!=0:
            wf.write(info+" "+str(cnt)+" "+str(totalval/cnt)+" "+str(n[3])+"\n")
        #Please remember we may lose exons due to 0 hits which is     
        cnt=0
        totalval=0
    #####################################################
    #Now we must sort by mrna-placement or line[3]
    wf.close()
    wf = open("b1.txt", 'r')
    finallist=[] #Our final list of information, now to be sorted by original position to keep true to the original file.
    
    for line in wf:
        linesplit=line.split()
        info=linesplit[0]
        count=int(linesplit[1])
        average=float(linesplit[2])
        position=int(linesplit[3])
        finallist.append((info,count,average,position))
    sortedfinal=sorted(finallist,key=takefourth)
    wf.close()
    wf = open("b1weightedaverage.txt", 'w')
    
    flag=True # This flag is triggered when we need to set the position of the exons.
    
    avg=0 #Keep track of total score.
    
    cnt=0 # Keep track of total # of exons in the strand. Remember we are outputting the average across all exons in the mRNA.
    
    graphcount=0 #Graphcount is  a variable used to keep track of the index, we wants one plot per 50 values.
    
    firstflag=True #A boolean flag to write down our textfile legend at the start of the textfile.
    
    g1=[] #Array to hold graphing information
    
    for n in sortedfinal:
        if flag==True:
            temp=n[3]
            flag=False #iniitally set temp
            
        if temp==n[3]: #check if still the same position
            cnt+=1
            avg+=float(n[2])
            info=n[0]
            
        elif temp!=n[3]:
            if firstflag:
                wf.write("MRNA INFO"+"    "+"# OF DATA POINTS HIT"+"    "+"AVERAGE SCORE"+"    "+"ID"+"    "+"DIVISION NUMBER"+"    "+"DIRECTION"+"\n")
                firstflag=False
                
            wf.write(info+"    "+str(cnt)+"    "+str(int(avg/cnt))+" "+str(temp)+"\n")
            if graphcount==50:
                g1.append((int(n[3]),int(n[2])))
                graphcount=0
            graphcount+=1
                
            temp=n[3]
            cnt=0
            avg=0
            cnt+=1
            avg+=float(n[2])
            info=n[0]
    wf.close()
    
    print("All done, check your directory for text file b1weightedaverage.txt")
    print("Creating graph for averages of Exons per mRNA")
    x,y = zip(*g1)
    plt.scatter(x,y)
    plt.xlabel("mRNA position")
    plt.ylabel("Score of mRNA")
    plt.show()
    




def ptbmm():
    '''
    Here we are checking the mean median and mode of part B.
    '''
    wf = open("C:/Users/Justin Macrae/Downloads/rop1/cv2.txt", 'r')
    wf2 = open("C:/Users/Justin Macrae/Downloads/rop1/bprimefivev2.txt", 'r')
    wf3 = open("C:/Users/Justin Macrae/Downloads/rop1/bprimethreev2.txt", 'r')
    wf4 = open("C:/Users/Justin Macrae/Downloads/rop1/bprimev2.txt", 'r')
    total=0
    total1=0
    total2=0
    total3=0
    l1=[]
    l22=[]
    l33=[]
    l3=[]
    l5=[]
    lp=[]
    tot1=0
    tot2=0
    tot3=0
    for line in wf:
        linesplit=line.split()
        if int(linesplit[4])==1:
            tot1+=float(linesplit[2])
            l1.append(float(linesplit[2]))
        if int(linesplit[4])==2:
            tot2+=float(linesplit[2])
            l22.append(float(linesplit[2]))
        if int(linesplit[4])==3:
            tot3+=float(linesplit[2])
            l33.append(float(linesplit[2]))
        
    print("Mean1: "+ str(tot1/len(l1)))
    sortedl1= sorted(l1)
    print("Median: " + str(sortedl1[len(l1)//2]))
    
    print("Mean2: "+ str(tot2/len(l22)))
    sortedl2= sorted(l22)
    print("Median: " + str(sortedl2[len(l22)//2]))
    
    print("Mean3: "+ str(tot3/len(l33)))
    sortedl3= sorted(l33)
    print("Median: " + str(sortedl3[len(l33)//2]))
    for line in wf2:
        linesplit=line.split()
        total1+=float(linesplit[2])
        l5.append(float(linesplit[2]))
    print("Mean for fives: "+ str(total1/10882))
    sortedl5= sorted(l5)
    print("Median: " + str(sortedl5[10883//2]))

    for line in wf3:
        linesplit=line.split()
        total2+=float(linesplit[2])
        l3.append(float(linesplit[2]))
    print("Mean for threes: "+ str(total2/17301))
    sortedl3=sorted(l3)
    print("Median: " + str(sortedl3[17302//2]))
    
    for line in wf4:
        linesplit=line.split()
        total3+=float(linesplit[2])
        lp.append(float(linesplit[2]))
    print("Mean for prime: "+ str(total3/38197))
    sortedlp=sorted(lp)
    print("Median: " + str(sortedlp[38198//2]))
    


   
def parsetextgene():
    '''
    Takes a gene exons text file, and parse it for exons and mrna.
    '''
    #f = open("C:/Users/Justin Macrae/Downloads/rop1/geneexons.gff3",'r')
    filename = askopenfilename(title="Select a file to parse!")
    f=open(filename,'r')
    wf = open("subsetdata.txt", 'w')
    mrna=False
    gene=True
    for line in f:
        
        if gene: #logic to retrieve the MRNA with the longest strand
            if "longest=1" in line:
                wf.write(line)
                mrna=True
                gene=False
        elif mrna:
            if "exon" in line:
                wf.write(line)
                
            if "gene" in line or "mRNA" in line:#logic to obtain all exons
                gene=True
                mrna=False
    f.close()
    wf.close()            
    print("done!")
            
def parsetextmrna():
    '''
    Creates a text file containing only MRNA.
    '''
    #f = open("C:/Users/Justin Macrae/Downloads/rop1/geneexons.gff3",'r')
    filename = askopenfilename(title="Select a file to parse!")
    f=open(filename,'r')
    wf = open("mrnasub.txt", 'w')
    for line in f:
        if "longest=1" in line:
            wf.write(line)
                
    f.close()
    wf.close()            
    print("done!")      
    
def partc():
    '''
    NOTE: Documentation was originally written to suffice a n=3 base case. Code is updated for just arbitrary n, just a note.
    Takes in bedgraph and subset files as input, outputs a text file containing each mrna strand (Exons Only) made up of exons added together
    and divided into 3 sections. Each section has its average score dispayed.
    '''
    #f = open("C:/Users/Justin Macrae/Downloads/rop1/norm.bedgraph",'r')
    #fmrna = open("C:/Users/Justin Macrae/Downloads/rop1/subsetdata",'r')
    
    n=input("Please enter the desired divisions: ")
    if n=="" or n=="0":
        n=1 
    n= int(n)
    #Code above is for command line taking in divisions as input.
    
    filename1 = askopenfilename(title="Select a bedgraph file to parse!")
    filename2 = askopenfilename(title="Select a file as input!")
    
    f = open(filename1,'r')
    fmrna = open(filename2,'r')
    
    valuel=[] #array for score values
    rangel=[]#array for range boundaries
    
    placement=0 #Keep track of placement in original file.
    
    direction="" #initialize direction
    
    for line in f:
        linesplit= line.split()
        grange= int(linesplit[1])
        totalvalue=float(linesplit[3])
        valuel.append((grange,totalvalue))
    print("done step 1: Aquire values from Data")
    sortedvalue=sorted(valuel,key=lambda x:x[0])
    print("done step 2: Sort values in Data")
    
    cnt=0 #Used to index the exons in our list tempinfo
    
    flagdiv=False #To help us ignore the first iteration
    
    tempinfo=[] #To store our temporary exon values
    
    add=0 #our temporary value used to keep track of left over genes
    
    dividernum=1#starts at 1, max is 3
    
    addflag=False#used to signal a need to use add
    
    temprange=0#placeholder for temp values in range
    
    check=False
    
    temptotal=0
    for line in fmrna:
        linesplit2= line.split()
        if "longest=1" in line and "scaffold" not in line:
            
            if flagdiv:
                divider=temptotal//n  # n This gives us the exact amount of dividers for each section
                sortedtemp=sorted(tempinfo,key=lambda x:x[0])
                
                if direction=="+":
                    dividernum=1
                    while dividernum != n+1: #n+1
                        exon=sortedtemp[cnt]
                        lbound=exon[0]
                        ubound=exon[1]
                        
                        if check:
                            lbound=temprange
                            ubound=savedstate
                            check=False
                            
                        length=ubound-lbound
                        if addflag:#if we noticed that add was changed we apply the math to get our second range
                            temprange=lbound+add
                            
                            if temprange<=exon[1]:#Iddeally what we want
                                rangel.append((lbound,temprange,exon[2],exon[3],dividernum,direction))
                                dividernum+=1
                                addflag=False
                                check=True
                                savedstate=exon[1]
                                if cnt+1!=len(sortedtemp):
                                    cnt+=1
                        
                            else:#If not repeat the process
                                rangel.append((lbound,ubound,exon[2],exon[3],dividernum,direction))
                                add=temprange-exon[1]
                                if cnt+1!=len(sortedtemp):
                                    cnt+=1
                        
                        elif length == divider: #case (1), if they exactly equal copy exon to rangel
                            rangel.append((lbound,ubound,exon[2],exon[3],dividernum,direction))
                            dividernum+=1
                            if cnt+1!=len(sortedtemp):
                                    cnt+=1
                        
                        elif length < divider:
                            add=divider - length
                            rangel.append((lbound,ubound,exon[2],exon[3],dividernum,direction))
                            #we now have an additional add value which we will use to compute the next range
                            addflag=True
                        elif length > divider:
                            temprange=lbound+divider
                            savedstate=exon[1]
                            rangel.append((lbound,temprange,exon[2],exon[3],dividernum,direction))
                            #we have to set up another if nest to let them know where to start
                            check=True
                            dividernum+=1
                            if cnt+1!=len(sortedtemp):
                                    cnt+=1
                elif direction=="-":
                    dividernum=n
                    while dividernum != 0: #n+1
                        exon=sortedtemp[cnt]
                        lbound=exon[0]
                        ubound=exon[1]
                        
                        if check:
                            lbound=temprange
                            ubound=savedstate
                            check=False
                            
                        length=ubound-lbound
                        if addflag:#if we noticed that add was changed we apply the math to get our second range
                            temprange=lbound+add
                            
                            if temprange<=exon[1]:#Iddeally what we want
                                rangel.append((lbound,temprange,exon[2],exon[3],dividernum,direction))
                                dividernum-=1
                                addflag=False
                                check=True
                                savedstate=exon[1]
                                if cnt+1!=len(sortedtemp):
                                    cnt+=1
                            
                            else:#If not repeat the process
                                rangel.append((lbound,ubound,exon[2],exon[3],dividernum,direction))
                                add=temprange-exon[1]
                                if cnt+1!=len(sortedtemp):
                                    cnt+=1
                        
                        elif length == divider: #case (1), if they exactly equal copy exon to rangel
                            rangel.append((lbound,ubound,exon[2],exon[3],dividernum,direction))
                            dividernum-=1
                            if cnt+1!=len(sortedtemp):
                                    cnt+=1
                        
                        elif length < divider:
                            add=divider - length
                            rangel.append((lbound,ubound,exon[2],exon[3],dividernum,direction))
                            #we now have an additional add value which we will use to compute the next range
                            addflag=True
                        elif length > divider:
                            temprange=lbound+divider
                            savedstate=exon[1]
                            rangel.append((lbound,temprange,exon[2],exon[3],dividernum,direction))
                            #we have to set up another if nest to let them know where to start
                            check=True
                            dividernum-=1
                            if cnt+1!=len(sortedtemp):
                                    cnt+=1
                    
                                
            tempinfo.clear()
            dividernum=1
            minfo= linesplit2[8]
            direction=linesplit2[6]
            mplace=placement
            placement+=1
            temptotal=0
            check=False
            cnt=0
                        
                    
                
            
        elif "exon" in line and "scaffold" not in line:
            mrangeh= float(linesplit2[4])
            mrangel= float(linesplit2[3])
            temptotal+=(mrangeh-mrangel)
            tempinfo.append((mrangel,mrangeh,minfo,mplace))
            flagdiv=True
            #rangel.append((mrangel,mrangeh,minfo,mplace))
            #Here I have it as exon lower,upper,mrna info, mrna placement so I can re-sort after.
    print("done step 3:Aquire values from Exons")
    sortedrange=sorted(rangel,key=lambda x:x[1])
    print("done step 4: Sort values for Exons")
    f.close()
    fmrna.close()
    #######################################################
    #ok now we have a list of data, structured as (lowbound,upbound,info,placement,divider(Either 1,2,3))
    #The goal is to have our averages for each mrna, and its divisions average of 1,2,3 
    
    wf = open("c.txt", 'w')
    
    cnt=0 #Keep track of # of hits.
    
    totalval=0 # Total score value.
    
    save=0 #Keep track of where last section left off at, to improve time complexity.
    
    listofzero=[] #List of sections hitting 0 places.
    
    for n in sortedrange:
        minvalue=n[0]
        maxvalue=n[1]
        info=n[2]
        #we have here the min and max values of the ranges
        for n2 in range(save,len(sortedvalue)):
            if sortedvalue[n2][0]<= maxvalue and sortedvalue[n2][0]>= minvalue:
                totalval+=sortedvalue[n2][1]
                cnt+=1
            elif sortedvalue[n2][0]>maxvalue:
                save=n2
                break
        if cnt==0:
            listofzero.append(n[3])
        if cnt!=0:
            wf.write(info+" "+str(cnt)+" "+str(totalval)+" "+str(n[3])+" "+str(n[4])+ " "+n[5]+"\n")
        #Please remember we may lose exons due to 0 hits.    
        cnt=0
        totalval=0
    #####################################################
    #Now we must sort by mrna-placement or line[3]
    wf.close()
    wf = open("c.txt", 'r')
    
    finallist=[] #This is to hold our information to eventually resort into originl file order.
    
    for line in wf:
        linesplit=line.split()
        info=linesplit[0]
        count=int(linesplit[1])
        average=float(linesplit[2])
        position=int(linesplit[3])
        finallist.append((info,count,average,position,linesplit[4],linesplit[5]))
    sortedfinal=sorted(finallist,key=takefourth)
    #We have it structured right
    wf.close()
    wf = open("cv2.txt", 'w')
    #Remember divide total value/ # of data points
    
    flag=True #Used to initially set the mRNA transcript we are working with.
    
    firstflag=True #Used to write the textfile legend at top of textfile.
    
    avg=0 # total score
    
    cnt=0 # # of hits/sections
    
    #The code initially, provides a way to average by number, i need to do that but now alsp average by division
    for n in sortedfinal:
        if firstflag:
            wf.write("MRNA INFO"+"    "+"# OF DATA POINTS HIT"+"    "+"AVERAGE SCORE"+"    "+"ID"+"    "+"DIVISION NUMBER"+"    "+"DIRECTION"+"\n")
            firstflag=False
        if flag==True:
            temp=n[3]
            divis=n[4]
            direct=n[5]
            flag=False #iniitally set temp
            
        if temp==n[3] and divis==n[4]: #check if still the same position, now also check for same number
            cnt+=float(n[1])
            avg+=float(n[2])
            info=n[0]
            
        elif temp!=n[3] or divis!=n[4]:
            wf.write(info+"    "+str(cnt)+"    "+str(avg/cnt)+"    "+str(temp)+"    "+ str(divis)+"    "+direct+"\n")
            temp=n[3]
            divis=n[4]
            direct=n[5]
            cnt=0
            avg=0
            cnt+=float(n[1])
            avg+=float(n[2])
            info=n[0]
        
              
    wf.close()
    print("All done! Check your directory for text file cv2.txt")
    
def promoter():
    '''
    Program designed to calculate the averages of promoter and downstream regions
    both as a whole and in mini batches of n.
    
    Takes in bedgraph and mrna subset as input files as well as desired # of divisions.
    Outputs 4 files including promoter and downstream regions, as well as them divided in n sections.
    
    Also included commented out code for an example when n=3
    
    '''
    n=input("Please enter the desired divisions: ")
    if n=="" or n=="0":
        n=1 
    #Above is a safety check for bad input.
    
    n= int(n)
    filename1 = askopenfilename(title="Select a bedgraph file to parse!")
    filename2 = askopenfilename(title="Select a file to parse!")
    f = open(filename1,'r')
    fmrna = open(filename2,'r')
    
    divn=2000/n #This is to create our division range for our promoter/downstream region.
    
    valuel=[] #An array to hold bedgraph info.
    
    rangepromote=[] #An array to hold information on our promoter reigon
    
    rangedown=[] #An array to hold info on our downstream region
    
    rangebatchd=[] #An array to hold info on our downstream region divided in n sections
    
    rangebatchp=[] #An array to hold info on our promoter region divided in n sections.
    
    rpr2=[] #Our second array to hold sorted values of each.
    rds2=[] #Our second array to hold sorted values of each.
    rbp2=[] #Our second array to hold sorted values of each.
    rbd2=[] #Our second array to hold sorted values of each.
    
    cnt=0 #Keep track of original position in file.
    
    for line in f:
        if "scaffold" not in line:
            linesplit= line.split()
            grange= float(linesplit[1])
            totalvalue=float(linesplit[3])
            valuel.append((grange,totalvalue))
    print("done step 1: Aquire values from Data")
    sortedvalue=sorted(valuel,key=lambda x:x[0])
    print("done step 2: Sort values in Data")
        
    for line in fmrna:
        if "scaffold" not in line:
            linesplit2= line.split()
            mrangeh= float(linesplit2[4])
            mrangel= float(linesplit2[3])
            minfo= linesplit2[8]
            direction= linesplit2[6]
            if direction=="+":
                rangepromote.append((mrangel-2000,mrangel,minfo,direction,cnt))
                
                for a in range(1,n+1):
                    rangebatchp.append(((mrangel-2000)+(divn*(a-1)),(mrangel-2000)+(divn*a),minfo,direction,a,cnt))
                
                #rangebatchp.append((mrangel-2000,(mrangel-2000)+666,minfo,direction,1,cnt))
                #rangebatchp.append(((mrangel-2000)+667,(mrangel-2000)+1333,minfo,direction,2,cnt))
                #rangebatchp.append(((mrangel-2000)+1334,mrangel,minfo,direction,3,cnt))
                
                rangedown.append((mrangeh,mrangeh+2000,minfo,direction,cnt))
                
                for a in range(1,n+1):
                    rangebatchd.append(((mrangeh)+(divn*(a-1)),(mrangeh)+(divn*a),minfo,direction,a,cnt))
                    
                #rangebatchd.append((mrangeh,mrangeh+2000-1334,minfo,direction,1,cnt))
                #rangebatchd.append((mrangeh+2000-1333,mrangeh+2000-667,minfo,direction,2,cnt))
                #rangebatchd.append((mrangeh+2000-666,mrangeh+2000,minfo,direction,3,cnt))
                
            elif direction == "-":
                    rangedown.append((mrangel-2000,mrangel,minfo,direction,cnt))
                    
                    for a in range(1,n+1):
                        rangebatchd.append(((mrangel-2000)+(divn*(a-1)),(mrangel-2000)+(divn*a),minfo,direction,a,cnt))
                    
                    
                    #rangebatchp.append((mrangeh,mrangeh+2000-1334,minfo,direction,1,cnt))
                    #rangebatchp.append((mrangeh+2000-1333,mrangeh+2000-667,minfo,direction,2,cnt))
                    #rangebatchp.append((mrangeh+2000-666,mrangeh+2000,minfo,direction,3,cnt))
                    
                    rangepromote.append((mrangeh,mrangeh+2000,minfo,direction,cnt))
                    
                    for a in range(1,n+1):
                         rangebatchp.append(((mrangeh)+(divn*(a-1)),(mrangeh)+(divn*a),minfo,direction,a,cnt))
                    
                    #rangebatchd.append((mrangel-2000,(mrangel-2000)+666,minfo,direction,1,cnt))
                    #rangebatchd.append(((mrangel-2000)+667,(mrangel-2000)+1333,minfo,direction,2,cnt))
                    #rangebatchd.append(((mrangel-2000)+1334,mrangel,minfo,direction,3,cnt))
        cnt+=1
            
            #I think this is good, just clarify

    print("done step 3:Aquire values from mRNA")
    sortedrange=sorted(rangepromote,key=lambda x:x[0])
    sortedrangedown=sorted(rangedown,key=lambda x:x[0])
    sortedrangebp=sorted(rangebatchp,key=lambda x:x[0])
    sortedrangebd=sorted(rangebatchd,key=lambda x:x[0])
    print("done step 4: Sort values for mRNA")
    f.close()
    fmrna.close()
    #############################################
    #this section of code has now sorted each lists by range
    #we can access lower bounds by rangel[n][0] and upper rangel[n][1]
    #we can use the sorted ranges to do comparisons
    wf = open("promoter.txt", 'w')
    
    cnt=0 #Keep track of number of hits
    
    countallzero=0 #keep track of regions that hit 0 values in bedgraph.
    
    totalval=0 #Total score.
    
    save=0 #Keep track of last place left off, see "Time Complexity" in presentation.
    
    for n in sortedrange:
        minvalue=n[0]
        maxvalue=n[1]
        info=n[2]
        position=n[4]
        #we have here the min and max values of the ranges
        for n2 in range(save,len(sortedvalue)):
            if sortedvalue[n2][0]<= maxvalue and sortedvalue[n2][0]>= minvalue:
                totalval+=sortedvalue[n2][1]
                cnt+=1
            elif sortedvalue[n2][0]>maxvalue:
                save=n2
                break
        if cnt==0:
            countallzero+=1
        if cnt!=0:
            rpr2.append((info,str(cnt),str(int(totalval/cnt)),position))
            #wf.write(info+" "+str(cnt)+" "+str(totalval/cnt)+" "+str(position)+"\n")
            
        cnt=0
        totalval=0
    wf.close()
    wf = open("downstream.txt", 'w')
    cnt=0
    countallzero=0
    totalval=0
    save=0
    for n in sortedrangedown:
        minvalue=n[0]
        maxvalue=n[1]
        info=n[2]
        position=n[4]
        #we have here the min and max values of the ranges
        for n2 in range(save,len(sortedvalue)):
            if sortedvalue[n2][0]<= maxvalue and sortedvalue[n2][0]>= minvalue:
                totalval+=sortedvalue[n2][1]
                cnt+=1
            elif sortedvalue[n2][0]>maxvalue:
                save=n2
                break
        if cnt==0:
            countallzero+=1
        if cnt!=0:
            rds2.append((info,str(cnt),str(int(totalval/cnt)),position))
            #wf.write(info+" "+str(cnt)+" "+str(totalval/cnt)+" "+ str(position)+"\n")
            
        cnt=0
        totalval=0
    wf.close()
    
    wf = open("promoter3batch.txt", 'w')
    cnt=0
    countallzero=0
    totalval=0
    save=0
    for n in sortedrangebp:
        minvalue=n[0]
        maxvalue=n[1]
        info=n[2]
        place=n[4]
        position=n[5]
        #we have here the min and max values of the ranges
        for n2 in range(save,len(sortedvalue)):
            if sortedvalue[n2][0]<= maxvalue and sortedvalue[n2][0]>= minvalue:
                totalval+=sortedvalue[n2][1]
                cnt+=1
            elif sortedvalue[n2][0]>maxvalue:
                save=n2
                break
        if cnt==0:
            
            countallzero+=1
        if cnt!=0:
            rbp2.append((info,str(cnt),str(int(totalval/cnt)),position,str(place)))
            #wf.write(info+" "+str(cnt)+" "+str(totalval/cnt)+" "+str(position)+ " " + str(place)+"\n")
            
        cnt=0
        totalval=0
    wf.close()
    
    wf = open("downstream3batch.txt", 'w')
   
    cnt=0 #Keep track of # of hits.
    
    countallzero=0 #Keep track of # of 0 hits.
    
    totalval=0 #Total score value.
    
    save=0 
    
    for n in sortedrangebd:
        minvalue=n[0]
        maxvalue=n[1]
        info=n[2]
        place=n[4]
        position=n[5]
        #we have here the min and max values of the ranges
        for n2 in range(save,len(sortedvalue)):
            if sortedvalue[n2][0]<= maxvalue and sortedvalue[n2][0]>= minvalue:
                totalval+=sortedvalue[n2][1]
                cnt+=1
            elif sortedvalue[n2][0]>maxvalue:
                save=n2
                break
        if cnt==0:
            countallzero+=1
        if cnt!=0:
            rbd2.append((info,str(cnt),str(int(totalval/cnt)),position,str(place)))
            #wf.write(info+" "+str(cnt)+" "+str(totalval/cnt)+" "+str(position)+" " + str(place) +"\n")
            
        cnt=0
        totalval=0
    wf.close()
    sort1=sorted(rpr2,key=lambda x:x[3])
    sort2=sorted(rds2,key=lambda x:x[3])
    sort3=sorted(rbp2,key=lambda x:x[3])
    sort4=sorted(rbd2,key=lambda x:x[3])
    
    wf1 = open("promoter.txt", 'w')
    wf2 = open("downstream.txt", 'w')
    wf3 = open("promoterNbatch.txt", 'w')
    wf4 = open("downstreamNbatch.txt", 'w')
    
    g1=[] #Array to hold graph values for promoter
    
    g2=[]# Array to hold graph values for downstream
    
    g3=[]# Array to hold graph values for promoter with n divisions
    
    g4=[]# Array to hold graph values for downstream with n divisions
    
    cnt=0 #Keep track of index of our values, when cnt=50 we add to our graph data.
    
    firstflag=True #Used to write down legend as first line.
    
    
    for n in sort1:
        if firstflag:
            wf1.write("MRNA INFO"+"    "+"# OF DATA POINTS HIT"+"    "+"AVERAGE SCORE"+"    "+"ID"+"    "+"DIVISION NUMBER"+"    "+"DIRECTION"+"\n")
            firstflag=False
            
        wf1.write(n[0]+" "+n[1]+" "+n[2]+" "+str(n[3])+"\n")
        if cnt==50:
            g1.append((int(n[3]),int(n[2])))
            cnt=0
        cnt+=1
        
    firstflag=True
    
    for n in sort2:
        if firstflag:
            wf2.write("MRNA INFO"+"    "+"# OF DATA POINTS HIT"+"    "+"AVERAGE SCORE"+"    "+"ID"+"    "+"DIVISION NUMBER"+"    "+"DIRECTION"+"\n")
            firstflag=False
            
        wf2.write(n[0]+" "+n[1]+" "+n[2]+" "+str(n[3])+"\n")
        if cnt==100:
            g2.append((n[3],n[2]))
            cnt=0
        cnt+=1

    for n in sort3:
        if firstflag:
            wf3.write("MRNA INFO"+"    "+"# OF DATA POINTS HIT"+"    "+"AVERAGE SCORE"+"    "+"ID"+"    "+"DIVISION NUMBER"+"    "+"DIRECTION"+"\n")
            firstflag=False
            
        wf3.write(n[0]+" "+n[1]+" "+n[2]+" "+str(n[3])+" "+n[4]+"\n")
        if cnt==100:
            g3.append((n[3],n[2]))
            cnt=0
        cnt+=1
        
    firstflag=True
    
    for n in sort4:
        if firstflag:
            wf4.write("MRNA INFO"+"    "+"# OF DATA POINTS HIT"+"    "+"AVERAGE SCORE"+"    "+"ID"+"    "+"DIVISION NUMBER"+"    "+"DIRECTION"+"\n")
            firstflag=False
        wf4.write(n[0]+" "+n[1]+" "+n[2]+" "+str(n[3])+" "+n[4]+"\n")
        if cnt==100:
            g4.append((n[3],n[2]))
            cnt=0
        cnt+=1
    print("All done check your directory for the textfiles promoter.txt, downstream.txt,promoterNbatch.txt,downstreamNbatch.txt")
    print("Displaying Graphs")
    print("Standard Promoter")
    x,y = zip(*g1)
    plt.scatter(x,y)
    plt.xlabel("mRNA position")
    plt.ylabel("Score of mRNA")
    plt.show()
    
    
