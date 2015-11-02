#!python2.7.1
import sys, os
sys.path.append('/Users/ssaberim/epigenomics/code')
import os 
import commands
import read
import write
import numpy as np





def write_QC1(fl,ordered,targetid,ref):
	
	#fread=open(fl,'r')
	lable=['H3K4me1','H3K4me3','H3K9me3','H3K27me3','H3K36me3','H3K27ac','Input DNA']
	data=read.read_dat(fl)
	num=len(ordered)/len(lable)
	ordered=list(ordered)
	ln=data[0]
	i=0
	for tmpi in range(1,len(data)):
		tmpdata=data[tmpi]
		libid=str(tmpdata[2])
		if libid in ordered:
			ind=ordered.index(libid)
			tmpdata.append(lable[int(ind/float(num))])
			idi=0
			for idtmp in range(1,len(targetid)):
				for libtmp in targetid[idtmp]:
					if str(libid[-4:]) in str(libtmp):
						idi=idtmp
			tmp=targetid[idi][-1] 
			for itm in ref:
				if tmp in itm[0]:
					tmp=itm[-1]
			if idi==0:
				print 'aha'
			tmpdata.append(str(tmp))
		if float(tmpdata[5])<20000000:
			tmpdata.append('failed')
		else:
			tmpdata.append('passed')
				
				
		#data[tmpi]=tmpdata
	fwrite=open('table.txt','w')
	for i in ln:
		print >> fwrite, str(i),
	print '\n'
	for tmpdata in data:
		for i in tmpdata:
			
			try:
				print >> fwrite, float(i),'\t',
			except:
				print >> fwrite, str(i),'\t',
		print >> fwrite, '\n',

	fwrite.close()
	return data
	
def write_hub(table,bwfiles):
	
	lable=['H3K4me1','H3K4me3','H3K9me3','H3K27me3','H3K36me3','H3K27ac','V5','input for V5','input for histon']
	num=(len(table)-1)/len(lable)
	data=[]
	c=[]
	proi=[]

# H3K4me1
	c.append( [ 255,102,51])
# H3K4me3
	c.append( [255,0,0])
# H3K9me3
	c.append( [ 0,0,102])
# H3K27me3
	c.append( [ 102,51,0])
# H3K36me3
	c.append( [ 153,0,153])
# H3K9ac
	c.append( [153,0,153])
# H3K27ac 
	c.append( [51,102,255])
# V5
	c.append( [151,10,255])
#input for hiton 
	c.append([0,0,0])
#input for V5
	c.append([100,100,100])
	ind=1
	for tmp in table:
		libid=str(tmp[3])
		libindex=str(tmp[2])
		libE=str(tmp[1])
		h=str(tmp[-2])
		solution=str(tmp[-1])[:-1]
		vector=str(tmp[-3])
		txt='\n\n'
		rep=int(tmp[-4])
		true=0
		for ibw in bwfiles:
			bw=str(ibw[0])
			for histon in lable:
				bwind=bw.split('_')[-1]
				
				bwind=bwind[:bwind.index('.')]
				if libid in bw and libindex == bwind  and histon in h:
					true+=1
					#print libindex
					p_f=''
					#if 'fail' in tmp[-1]:
					#	p_f='(<20M aligned reads)'

					txt+="track "+libid+'_'+libindex+'_'+vector+"\n"
					txt+="shortLabel "+libid+'_'+libindex+"\n"
					
					txt+="longLabel "+h+', '+vector+', '+ libE+', rep'+str(rep)+"\n"
					txt+="type bigWig\nvisibility full\nmaxHeightPixels 70:70:32"
					txt+="\nconfigurable on\nautoScale on\nalwaysZero on\n"
				
					txt+="priority "+ str(ind)+"\n"
					ind+=1
					txt+="bigDataUrl "+bw+"\n"
					colortmp=c[lable.index(histon)]
					#print lable.index(histon), histon
					ctxt=''
					txt+="color "+str(colortmp)+"\n"
					txt+="\n"
					txt+="\n"
					if  true>1:
						print str(true) +'_'+libid+'_'+libindex 
					else:
						data.append(txt)

	print (ind)
	out='trackDb'
	f=open(out+'.txt','w')
	for tmp in data:
		print >> f, tmp
	f.close()
	
def write_data(fl,data):
	fwrite=open(fl,'w')

	for i in data:
		tmpstr=''
		for tmpdata in i:
			tmpstr+= str(tmpdata)+'\t'
	#	print tmpstr
		print >> fwrite, tmpstr[:-1]

	fwrite.close()


def write_bed(intersect,file,*strand):
	
	if len(strand)>0:
		strand=strand[0]
	f=open(file,'w')
	ind=0
	for key in intersect.keys():
		for elt in intersect[key]:
			stout=key+'\t'
			for j in range(2):
				stout+=str(int(elt[j]))+'\t'
			if strand:
				stout+='+1'+'\t'
				
			try:
				stout=stout+str(float(elt[-1]))
			except:
				stout+='peak_'+str(ind)+'\t'
			ind+=1
			print >> f, stout
		
	print stout
	print ind
	f.close()


def write_bed_loc(projlist,tresh,coverfile):
	
	
	data=read.read_dat(coverfile,'\t')
	
	f=open('out.bed','w')
	ind=0
	print len(data), len(data[0])
	genes=[]
	for i in xrange(len(projlist)):
		
		if abs(projlist[i])>tresh:
			ind+=1
			stout=''
			genes.append(str(data[i][3]))
			for j in range(3):
				stout+=str(data[i][j])+'\t'
			print >> f, stout[:-1]
	print 'num of peaks: ',ind
	f.close()
	
	return genes



	
def write_mat(xheader,yheader,data,fl):
	fwrite=open(fl,'w')
	
	tmpstr='\t'
	if len(xheader)>0:
		for i in xheader:
			tmpstr=tmpstr+str(i)+'\t'
	print >> fwrite, tmpstr[:-1]

	ind=0
	for i in list(data):
		tmpstr=yheader[ind]+'\t'
		ind += 1
		for tmpdata in list(i):
				tmpstr = tmpstr+str(tmpdata)+'\t'
		print >> fwrite, tmpstr[:-1]
		#print  tmpstr[:-1]
		

	fwrite.close()


def write_encode(encode,file,):
	
	f=open(file,'w')
	ind=0
	for tf in encode.keys():
		
		elt=encode[tf]
		for chrom in elt.keys():
		    for i in range(len(elt[chrom])):
			stout=chrom+'\t'
			for j in xrange(2):
				
				try:
					stout+=str(int(elt[chrom][i,j]))+'\t'
				except:
					print elt[chrom][i,j]
			stout+=str(tf)
			print >> f, stout
		
	print stout
	f.close()
