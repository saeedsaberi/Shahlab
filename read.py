#!python2.7.1
#coding: utf8 
import sys
sys.path.append('/Users/ssaberim/epigenomics/code')
import os , gzip   
import commands
import read
import analyse
import numpy as np

def tabless(st):
       i=str(st)
       try:
         while '\n' in i:
             i=i[0:i.index('\n')]+i[i.index('\n')+1:]
         while '\t' in i:
             i=i[0:i.index('\t')]+i[i.index('\t')+1:]
         while '\r' in i:
             i=i[0:i.index('\r')]+i[i.index('\r')+1:]       
             #print i
         return i
       except:
         return st 

def read_dat(file,*spl):
    if len(spl)!=1:
       spl=" "
    f=open(file,'r')
    a=[]
    if len(spl)>0:
        spl=str(spl[0])
    ind=0
    for  ln  in f:
       s = ln.split(spl)
       ind+=1
       if ind>100000000:
              break
       for i in xrange(len(s)):
            s[i]=s[i].strip()
            try :
               try:
                 s[i]=int(s[i])
               except:
                 s[i]=float(s[i])
            except:
               pass
       a.append(s)
    f.close()
    return  a

def readall_coverage(dir,mark,libs):

    ext='.coverage'
    try:
        direct=os.popen('ls '+dir +'/*'+mark+'*'+ext)
    except:
         direct=os.popen('ls '+dir +'/'+mark+'*'+ext)
    data=[]
    files=[]
    for fl in direct.readlines():
			#print fl
			tmpdata=read.read_dat(fl[:-1],'\t')
			data.append(tmpdata)
			files.append(fl.split('/')[-1][0:-len('.coverage')-1])
    genes=[]
    for i in data[0]:
        genes.append(i[3])

    mat=[]
    #\genes=np.array(genes)
    for i in data:
        mat.append(np.zeros((len(genes)),np.float))
        for j in i:
            try:
                mat[-1][genes.index(j[3])]+=np.float(j[4])
                #print j[2]
            except:
                pass
    mat=np.array(mat)
    print 'ending'
    libs=read.read_dat(libs,'\t')
    for i in libs:
        for j in range(len(files)):
          if i[-1] in files[j]:
              files[j]=i[0]
          if 'hg19v69_genes.TSS_2000.pc.' in files[j]:
              ind=files[j].index('hg19v69_genes.TSS_2000.pc.')
              files[j]=files[j][ind]+files[j][ind+len('hg19v69_genes.TSS_2000.pc.'):]


              #break

    return  files,mat,genes

def read_dir(dir,mark,loc,ext,*spl):
	
	try:
              direct=os.popen('ls '+dir +'/*'+mark+'*'+ext)
	except:
              direct=os.popen('ls '+dir +'/'+mark+'*'+ext)
	data=[]
	files=[]
	for fl in direct.readlines():
		if loc in fl:
			#print fl
			if len(spl)>0:
				tmpdata=read.read_dat(fl[:-1],str(spl[0]))
			else:
				tmpdata=read.read_dat(fl[:-1])
			data.append(tmpdata)
			#location=fl[::-1].index('/')
			libst=['A','HS','E']
			ind=0
			run=True
			while run and ind<len(libst):
				try:
					Ai=fl.index(libst[ind])
					run=False
					Aend=Ai+fl[Ai:].index('.')
				except:
					ind+=1
			fl=read.tabless(fl)
			files.append(fl[Ai:Aend])

	return  files, data




def read_alldir(dir,ext,*include):
#h3k27files, allh3k27=read.read_alldir('rhabdoid/coverage/TSS_2000_all','coverage','H3K27me3')
       data=[]
       files=[]
       spl='\t'
       if len(include)>0:
              include='*'+str(include[0])+'*'
       else:
              include='*'
       direct=os.popen('ls '+dir +'/'+include+ext)
       for fl in direct.readlines():
            tmpdata=read.read_dat(fl[:-1],str(spl))
            data.append(tmpdata)
            files.append(fl[-fl[::-1].index('/'):][:-1])

       
       return files, data

def readall_bed(dir,mark,libs,*flagstat):
       
       direct=os.popen('ls '+dir +'/*'+mark+'*.coverage')
       data=[]
       files=[]
       for fl in direct.readlines():
       	fl=fl.strip()
       	for alib in libs:
       		if alib in fl and '#' not in alib:
       			tmpdata=read.read_bed(fl)
       			data.append(tmpdata)
       		#location=fl[::-1].index('/')
       			libst=['A','HS','E']
       			ind=0
       			run=True
       			while run and ind<len(libst):
       				try:
       					Ai=fl.index(libst[ind])
       					run=False
       					Aend=Ai+fl[Ai:].index('.')
       				except:
       					ind+=1
       					Ai=fl.index('/')
       					Aend=len(fl)
       					#print libst[ind-1] , 'tryed';
       			fl = fl[Ai:Aend]
       			#print fl
       			if mark in fl:
       				fl=fl[fl.index(mark):]
       				fl=fl[:len(mark)+1+fl[len(mark)+1:].index('.')]
       				#print fl
       			files.append(fl.strip())
       
       if len(flagstat)>0:
                flagstat=str(flagstat[0])
                normalization=read.read_dat(flagstat,'\t')
                for fl in range(len(files)):
                     bool=False
                     for i in normalization:
                            if files[fl] in i[0]:
                                   data[fl]=norm_bed(data[fl],i[-1])
                                   #print i[-1]
                                   bool=True
                                   break
                     if not bool:
                            print files[fl], 'not found'
                            
       return  files, data
def norm_bed(bed,norm):
       
       for chrom in bed.keys():
                     bed[chrom][:,2]*=50000000./norm


       return bed
def norm_wig(wig,norm):
       
       for chrom in wig.keys():
                     wig[chrom][:,2]*=50000000./norm


       return wig
	
def readall_genescore(dr,mark,lib,genes):
	
	direct=os.popen('ls '+dr +'*'+mark+'*.coverage')
	data=[]
	mylib=[]
	vs=[]
	for i in xrange(len(lib)):
          if '#' not in lib[i]:
		if '-' in lib[i]:
			mylib.append(lib[i].split('-')[0])
		else:
			mylib.append(lib[i])
	data={}
	ind=-1
	for fl in direct.readlines():
		ind+=1
		for alib in mylib:
			if alib in fl:
				found=True
				tmpdata=read.read_dat(fl[:-1],'\t')
				data[alib]=tmpdata
				break

	mydata=[]
	for alib in mylib:
          #if '#' not in lib:
		mydata.append(data[alib])

	
	genenum=len(genes)
	libnum=len(mylib)
	mat=np.zeros((genenum,libnum),np.float)
	for genei in xrange(genenum):
		for libi in xrange(libnum):
			found=False
			for line in mydata[libi]:
				if genes[genei] in line[3]:
					mat[genei,libi]=np.float(line[-2])
					found=True
					break
			#if found:
			#	print 'found'
			#else:
			#	print genes[genei],mylib[libi] 
	return   mat




def readnames_dir(dir,loc):
	
	
	direct=os.popen('ls '+dir +'/*'+loc)

	files=[]
	for fl in direct.readlines():
		if loc in fl:


			files.append(fl.strip())

	return  files
	
	
def read_out(file):
    f=open(file,'r')
    s = f.readline()
    a=[]
    while s != '':
       s = s.split()
       tmp=map(float,s)
       a.append(tmp)
       s = f.readline()
    f.close()
    return  a

def read_plot(file):
    f=open(file,'r')
    s = f.readline()
    a=[]
    while s != '':
       #s = s.split()
       tmp=float(s)
       a.append(tmp)
       s = f.readline()
    f.close()
    return  a

def read_dist(dr):

	lib=[]
	dist=[]
	direct=os.popen('ls '+dr +'/*.dist')
	#print direct,direct[0], len(direct)
	for fl in direct.readlines():
		#print fl
		data=read.read_dat(fl[:-1])
		tmp=[]
		for line in data:
			tmp.append(float(line[-2]))
		dist.append(np.array(tmp))
		lib.append(fl[1+len(dr):len(dr)+7])
		
	return lib,dist
	

		

def read_bed(file):

    elts = {}
    if 'narr' in file:
       num=6
    elif 'coverage' in file:
       num=-2
    else:
       num=3

    f = open(file,'r')
    s = f.readline()

    addition=''
    if 'chr' not in str(s):
       addition='chr'
    f.close()
    bool=False
    if len(s.split('\t'))>3:
       bool=True
    f = open(file,'r')
    for  s in f:
        s = s.split('\t')
        try:
          name = str(s[3])
        except:
              pass
        chrom = str(s[0])
        begin = int(eval(s[1]))
        end = int(eval(s[2]))

        if bool:
              
              try:
                score = float(eval(s[num]))
              except:
                score = float(s[num].strip().split(','))
        else:
              score=0
        if not elts.has_key(chrom):
            elts[chrom] = []
        if ((len(elts[chrom])==0) or (begin!=elts[chrom][-1][0])):
           elts[chrom].append([begin,end,score])



    #print '#num of read bed peaks:', ind

    for key in elts.keys():
       elts[key]=np.array(elts[key])
       if 'chr' not in key:
              elts['chr'+key] = elts.pop(key)

    f.close()
    return elts

def read_gene_pos(file):

    elts = {}
    f = open(file,'r')
    s = f.readline()
    ind=0

    while s != '':
        s = s.split('\t')

        
        
        chrom = str(s[0])
        if 'chr' not in chrom:
            chrom = 'chr'+ chrom
        begin = int(eval(s[1]))
        end = int(eval(s[2]))

        if not elts.has_key(chrom):
            elts[chrom] = []
        try: 
          strand=int(eval(s[3]))
          name = str(s[-1]).strip()          
          elts[chrom].append([int(begin),int(end),int(strand),name])
        except:
          name = str(s[3]).strip()              
          elts[chrom].append([int(begin),int(end),name])

        
        #print name
        
        s = f.readline()

    f.close()
    #print '#num of read bed peaks:', ind
    for key in elts.keys():
       elts[key]=np.array(elts[key])
    return elts
def read_multibed(file):

    elts = {}
    f = open(file,'r')
    s = f.readline()

    while s != '':
        try:
          s = s.split('\t')
          chrom = str(s[0])
          begin = float(s[1])
          end = int(s[2])
          overlap = int(s[3])
          ls=map(int,s[5:])
          #print ls[0:10], sum(ls[0:])
          number = len(ls)
          if not elts.has_key(chrom):
              elts[chrom] = []

          elts[chrom].append([begin,end,float(sum(ls[0:]))/float(number),ls])
        except:
          s=s

        s = f.readline()
    f.close()

    return elts


def read_wig(file,*flagstat):

    elts = {}
    
    f = open(file,'r')
    if '.gz' in file:
       f = gzip.open(file)
    else:
       f = open(file,'r')
    s = f.readline()
    print s
    s = f.readline()
    ind=0
    cover1=int(s[:-1])
    f.close()
    if '.gz' in file:
       f = gzip.open(file)
    else:
       f = open(file,'r')
    count=0
    data=f.read()
    data=data.splitlines()
    print 'reading done'
    for s in data: 
        #count+=1
        #if count >10000000: 
        #      break
        try:
           cover=int(s)
           if cover==cover1:
             reg+=step
           else:
             begin=start+ind*step
             end=begin+reg
             start=end
             ind=0
             elts[chrom].append([begin,end,cover1])
             cover1=cover
             reg=step
             #ind+=step
        except:
           ind=0
           s = s.split()
           s=map(str,s)
           tmp=s[1]
           chrom=tmp[6:]
           tmp=s[2]
           start=int(tmp[6:])
           tmp=s[3]
           step=int(tmp[5:])
           reg=step
           if not elts.has_key(chrom):
                  elts[chrom] = []
        #s = f.readline()
    f.close()
    for key in elts.keys():
       elts[key]=np.array(elts[key])
    try:
      if len(flagstat)>0:
              flagstat=str(flagstat[0])
              normalization=read.read_dat(flagstat,'\t')
              fl=file[0:6]
              print fl
              bool=False
              for i in normalization:
                            if fl in i[0]:
                                   elts=norm_wig(elts,i[-1])
                                   #print i[-1]
                                   bool=True
                                   break
              if not bool:
                            print fl, 'not found'
    except:
       pass
    return elts
            
            



def read_encode(file):

    elts = {}


    f = open(file,'r')
    for  s in f:
        s = s.split('\t')

        name = str(s[3])
        chrom = str(s[0])
        begin = int(eval(s[1]))
        end = int(eval(s[2]))
        score= int(eval(s[5]))
        if not elts.has_key(name):
            elts[name] = {}
        if not elts[name].has_key(chrom):
               elts[name][chrom] = []
        elts[name][chrom].append([begin,end,score])
    for name in elts.keys():
       for chrom in elts[name].keys():
              elts[name][chrom]=np.array(elts[name][chrom])

    f.close()
    #print '#num of read bed peaks:', ind
    return elts








def read_encode_intersect(file):

    elts = {}


    f = open(file,'r')
    for  s in f:
        s = s.split('\t')

        name = str(s[3])
        chrom = str(s[0])
        begin = int(eval(s[1]))
        end = int(eval(s[2]))
        score= int(eval(s[-1]))
        xi=int(eval(s[-3]))
        xf=int(eval(s[-2]))
        if not elts.has_key(name):
            elts[name] = {}
        if not elts[name].has_key(chrom):
               elts[name][chrom] = []
        if (len(elts[name][chrom])>0) and (begin==elts[name][chrom][-1][0]):
              elts[name][chrom][-1][-1]=max(elts[name][chrom][-1][-1],score)
              elts[name][chrom][-1][-2]=max(elts[name][chrom][-1][-2],xf)
              elts[name][chrom][-1][-3]=min(elts[name][chrom][-1][-3],xi)
        else:
              elts[name][chrom].append([begin,end,xi,xf,score])
    for name in elts.keys():
       for chrom in elts[name].keys():
              elts[name][chrom]=np.array(elts[name][chrom])

    f.close()
    #print '#num of read bed peaks:', ind
    return elts



def read_bed_intersect(file):

    elts = {}


    f = open(file,'r')
    for  s in f:
        s = s.split('\t')
        chrom = str(s[0])
        begin = int(eval(s[1]))
        end = int(eval(s[2]))
        coverage=np.float(eval(s[3]))
        score= np.float(eval(s[-2]))
        cpgnum= int(eval(s[-1]))


        if not elts.has_key(chrom):
               elts[chrom] = []
        elts[chrom].append([begin,end,coverage,score,cpgnum])
    for chrom in elts.keys():
              elts[chrom]=np.array(elts[chrom])
    f.close()
    #print '#num of read bed peaks:', ind
    return elts


def readall_bedintersect(dir,mark,*flagstat):
       
       direct=os.popen('ls '+dir +'/*'+mark+'*.bed')
       data=[]
       files=[]
       for fl in direct.readlines():
          fl=fl.strip()
          tmpdata=read.read_bed_intersect(fl)
          data.append(tmpdata)
          ind=0
          run=True
          fl=fl[-fl[::-1].index('/'):fl.index('.')]
          print fl
          files.append(fl)
       
       if len(flagstat)>0:
                flagstat=str(flagstat[0])
                normalization=read.read_dat(flagstat,'\t')
                for fl in range(len(files)):
                     bool=False
                     for i in normalization:
                            if files[fl] in i[0]:
                                   data[fl]=norm_bed(data[fl],i[-1])
                                   bool=True
                                   break
                     if not bool:
                            print files[fl], 'not found'
                            
       return  files, data




def read_mark_encode_intersect(file):

    elts = {}


    f = open(file,'r')
    for  s in f:
        s = s.split('\t')
        chrom = str(s[0])
        begin = int((s[5]))
        end = int((s[6]))
        score=np.float(eval(s[3]))
        TF= str(s[-5])
        if not elts.has_key(TF):
               elts[TF] = {}
        if not elts[TF].has_key(chrom):
               elts[TF][chrom] = []
        elts[TF][chrom].append([begin,end,score])
    for TF in elts.keys():
       for chrom in elts[TF].keys():
           elts[TF][chrom]=np.array(elts[TF][chrom])
    f.close()
    #print '#num of read bed peaks:', ind
    return elts




def read_nuc_bed(file):

    elts = {}

    num=3
    f = open(file,'r')
    s = f.readline()

    addition=''
    if 'chr' not in str(s):
       addition='chr'
    f.close()
    bool=False
    if len(s.split('\t'))>3:
       bool=True
    f = open(file,'r')
    for  s in f:
        s = s.split('\t')
        try:
          name = str(s[3])
        except:
              pass
        chrom = str(s[0])
        begin = int(eval(s[1]))
        end = int(eval(s[2]))

        if bool:
                score = map(float,s[num:])
        tmp=[begin,end]
        for i in score:
              tmp.append(i)
        if not elts.has_key(chrom):
            elts[chrom] = []
        if ((len(elts[chrom])==0) or (begin!=elts[chrom][-1][0])):
           elts[chrom].append(tmp)



    for key in elts.keys():
       elts[key]=np.array(elts[key])
       if 'chr' not in key:
              elts['chr'+key] = elts.pop(key)

    f.close()
    return elts









