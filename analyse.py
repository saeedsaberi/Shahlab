#!/usr/local/bin/env python
import sys, os
sys.path.append('/Users/ssaberim/epigenomics/code')
import  bisect, itertools
import scipy
import scipy.stats as ss
import numpy as np
import random
import read, analyse, write
import matplotlib.pyplot as plt
import os
import numpy as np
import fastcluster
import numpy.linalg as LA
import scipy.cluster.hierarchy  as sch
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import numpy
#color=iter(cm.rainbow(np.linspace(0,1,n)))
#c=next(color)
#plt.plot(x,y,c=c)

import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as colors


# NUM_COLORS = 20
#
# cm = plt.get_cmap('gist_rainbow')
# cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
# scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
# fig = plt.figure()
# ax = fig.add_subplot(111)
# # old way:
# #ax.set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
# # new way:
# ax.set_color_cycle([scalarMap.to_rgba(i) for i in range(NUM_COLORS)])
# for i in range(NUM_COLORS):
#     color = cm(1.*i/NUM_COLORS)  # color will now be an RGBA tuple
#     ax.plot(np.arange(10)*(i+1))
#
# fig.savefig('moreColors.png')
# plt.show()
def plt_data(data):
   lib=[]
   total=[]
   for i in range(len(data)):
        try:
            tmp=data[i]
            
            total.append(float(tmp[-1][6:]))
            tmp=data[i-1]
            lib.append(str(tmp[1][:][:6]))
            #lib.append(str(tmp[0][len('/projects/analysis/analysis12/IX1983/C2NR9ACXX_5/'):][:6]))
        except:
            j=0
   return lib,total
         
def saveplt():
	savefig('coverage.pdf', bbox_inches='tight')
	

def order(total,lib,marks):
    
    
    lable=['H3K4me1','H3K4me3','H3K9me3','H3K27me3','H3K36me3','H3K27ac','Input DNA']
    num=len(lable)
    ordered=[]
    marks_tot=len(marks)
    
    for i in range(num):
       for j in range(marks_tot/num):
         ordered.append(marks[i+num*j][0])
    tot=[]
    for tmp in ordered:
       tot.append(total[lib.index(tmp)])
    return np.array(tot),np.array(ordered)



	
def pl_coverage(tot1,tot2,tot3,ordered):
	
	lable=['H3K4me1','H3K4me3','H3K9me3','H3K27me3','H3K36me3','H3K27ac','Input DNA']
	
	num=len(tot1)
	num=num/len(lable)
	for i in range(len(lable)):
		fig, ax = plt.subplots()
		plt.plot(tot1[i*num:(i+1)*num]/1000000,'gs-')
		plt.plot(tot2[i*num:(i+1)*num]/1000000,'bs-')
		plt.plot(tot3[i*num:(i+1)*num]/1000000,'rs-')
		plt.title(lable[i])
		plt.xticks(range(len(ordered[i*num:(i+1)*num])),ordered[i*num:(i+1)*num])
		plt.xticks(rotation=90)
		plt.xticks(fontsize=14)
		plt.yticks(range(0,161,20))
		plt.ylabel('Milion reads')
		plt.legend(('total reads after filtering','total reads','total aligned'),loc='upper right',prop={'size':14})
		fig.savefig('figs/'+lable[i]+'.pdf', bbox_inches='tight')
		
		
		
					
	return libs
def pl_PET(dr,table):
	
	def pl(tmpdata):
		ls=[]
		for i in tmpdata:
			ls.append(float(i[1]))
		return np.array(ls)
	
	lable=['H3K4me1','H3K4me3','H3K9me3','H3K27me3','H3K36me3','H3K27ac']
	tab=read.read_dat(table)
	libs=[]
	for i in range(len(lable)):
		for j in range(len(tab[0])):
			if lable[i] in tab[0][j]:
				libs.append([])
				for k in range(1,17):
					libs[-1].append(tab[k][j])
    
	num=len(libs[0])
	data=[]
	ind=0
	dr=os.popen('ls '+ dr +'*'+'A'+'*.dist')
	for fl in dr.readlines():
		for i in range(len(lable)):
			fig=plt.figure(i)
			plt.title(lable[i])
			for lib in libs[i]:
					#print lib, lable[i]
					if lib in fl[:-1]:
						#print lib,fl[:-1],lable[i]
						tmpdata=pl(read.read_dat(fl[:-1]))
						data.append(tmpdata)
						plt.plot(range(50,550),tmpdata[50:550]/np.max(tmpdata[50:550]),'s-',label=lib)
						#libs.append(lib)
			plt.legend(prop={'size':14})
	print i
	for i in range(len(lable)):
		fig=plt.figure(i)
		fig.savefig(lable[i]+'.pdf', bbox_inches='tight')
		

	return data

	
	
	
		
def pl_qc1(qc,columns,ordered):
	lable=['H3K4me1','H3K4me3','H3K9me3','H3K27me3','H3K36me3','H3K27ac','Input DNA']

	qc0=[]
	lib=[]
	for col in columns:
		col_i=qc[0].index(col)
		qc0.append([])
		lib.append([])

		for el in ordered:
			for tmpqc in qc:
				if el in tmpqc[1]:
					try:
						qc0[-1].append(float(tmpqc[col_i]))
					except:
						qc0[-1].append(float(tmpqc[col_i+2]))
					lib[-1].append(str(el))
				
		num=len(lib[-1])
		num=num/len(lable)
		qc0[-1]=np.array(qc0[-1])
		for i in range(len(lable)):
			fig=plt.figure(i)
			plt.plot(qc0[-1][i*num:(i+1)*(num)],'s-',label=col)
			plt.xticks(range(num),lib[-1][i*num:(i+1)*(num)])
			plt.ylabel('percentage%')
			plt.xlabel('Library')
		#plt.xticks(range(len(lib)),lib)
			plt.xticks(rotation=90)
			plt.xticks(fontsize=14)
			plt.title(lable[i])
			if i==6:
				plt.legend(prop={'size':14},loc=0)
			fig.savefig('QC1/QC1-'+lable[i]+'map_percent.pdf', bbox_inches='tight')

	return (qc0),lib
	
	
	
				
				
				
				
	
def mk_cluster(dr,mark,loc,*perc):
	#writes the tanle for clustering
	
	header=[]
	if len(perc)>0:
		f=open('cluster-'+mark+'-'+loc+'-'+str(perc[0])+'.txt','w')
	else:
		f=open('cluster-'+mark+'-'+loc+'.txt','w')
	direct=os.popen('ls '+dr +'/*'+mark+'*.coverage')
	for fl in direct.readlines():
		#print fl
		libst=['A','HS','E']
		ind=0
		run=True
		while run and ind<=len(libst):
			try:
				Ai=fl.index(libst[ind])
				run=False
			except:
				ind+=1
		Aend=Ai+fl[Ai:].index('.')
		print >> f, '\t',fl[Ai:Aend],
		header.append(fl[Ai:Aend])
		
	print >> f, '\n',
	direct=os.popen('ls '+dr +'/*'+mark+'*.coverage')
	data=[]
	for fl in direct.readlines():
		if loc in fl:
			tmpdata=read.read_dat(fl[:-1],'\t')
			if  len(tmpdata)<1:
				print len(tmpdata),fl[:-1]
			else:
				data.append(tmpdata)
	print len(data),len(data[0])
	enrich=[]
	lable=[]
	table=[]
	
	for line in range(len(data[0])):

		table.append([])
		#lable.append([])
		for tmpdata in data:
			#print  (tmpdata[0])
			table[-1].append(float(tmpdata[line][-2]))
		table[-1]=np.array(table[-1])
		tmp2=''
		for i in range(4):
				tmp2=tmp2+str(tmpdata[line][i])+'_'
		tmp2=tmp2[:-1]+'_1'
		lable.append(tmp2)
	table=np.array(table)
	print len(table[0,:]), len(table[:,0])
	if len(perc)>0:
		for row in range(len(table[0,:])):
			tmp=np.percentile(table[row,:],float(perc[0]))
			print tmp,
			table[row,:][table[row,:]<tmp]=0.
	enrich=[]
	for i in range(len(table[:,0])):
		enrich.append('')
		for j in range(len(table[0,:])):
			enrich[-1]+=str(table[i,j])+'\t'
	for line in range(len(data[0])):
			print >> f, lable[line],enrich[line]
	#print lable[:3],enrich[:3]
	f.close()
	
	return lable,header,table
	
def dist_coverage(dr,mark,loc,*perc):
	# generates the data for clustering

	files,data=read.read_dir(dr,mark,loc,'\t')


	#print len(data),len(data[0])
	enrich=[]
	lable=[]
	for line in range(len(data)):
		enrich.append([])
		tmp=[]
		for tmpdata in data[line]:
			tmp.append(float(tmpdata[-2]))
		enrich[-1]=np.array(tmp)
	enrich = np.array(enrich)
	print len(enrich[0,:]) , len(enrich[:,0])
	if len(perc)>0:
		for row in range(len(enrich)):
			#print len(enrich[0,:])
			tmp= np.percentile(enrich[row,:],float(perc[0]))
			#print tmp,
			enrich[row,:][enrich[row,:]<tmp]=0.
	return enrich, files
	
	
def PCA(data):

    Ncol = len(data[0])
    Nrow = len(data)

    if Nrow > Ncol:
        cov_data = np.cov(np.transpose(data - np.mean(data,axis=0)) )
    else:
        cov_data = np.cov(data-np.mean(data,axis=1))

    w, V = LA.eigh(cov_data)

    return w, V
    
def PCA2(cij):
	
    cij=np.array(cij)
    cavg = np.mean(cij, axis=0)
    cij_avg = cij - cavg

    w, V = analyse.PCA(cij_avg)

    #pc1=(V[:,-1])
    proj1 = np.dot( cij_avg, V)
	
    #cij_cor = cij_avg - np.outer(proj1, V[:,-1])
    #pc=np.array(np.outer(proj1*pc1))
    #print(len(V[:,-1]))
    return w, V, proj1, cij_avg, cavg
    
    

def stats(dists):
	ln=len(dists)
	heat=np.zeros((ln,ln),np.float)
	for i in range(ln):
		for j in range(i,ln):
			tmp=np.log(ks_2samp(dists[i],dists[j])[1])
			heat[i,j]=tmp
			heat[j,i]=tmp
			
			
	return heat
	
	
def mean_err(dists):
	
	mean=[]
	st=[]
	for i in dists:

		mean.append(np.average(i))
		st.append(np.std(i))
	print len(i)
	
	
	return np.array(mean),np.array(st)
	

def all_corr(vec):

   ln=len(vec[0])
   print ln, ln
   print len(vec), len(vec[0])
   vec=np.array(vec)


   norm=[]
   for ii in xrange(ln):
      norm.append(ss.zscore(vec[:,ii]))
      #norm.append((vec[:,ii]))
   norm=analyse.data2arr(norm)
   print 'Spearman ...'
   cor_mat=ss.spearmanr(norm)[0]
   #cor_mat=[]
   #for ii in range(len(vec)):
   #   cor_mat.append([])
   #   for jj in xrange(len(vec)):
   #      tmp=ss.pearsonr(vec[ii],vec[jj])
   #      cor_mat[-1].append(tmp[0])
   #print 'pearson ...'
   cor_mat=np.array(cor_mat)
   ####
   dist=1-cor_mat
   distance=[]
   for i in xrange(len(dist)):
      distance=distance+list(dist[i][i+1:])
   distance=np.array(distance)
   return norm.T,cor_mat,distance
 
 
def heatmap_cor( x, vec, minval, maxval,*track ):

  preprocessCore = importr('preprocessCore')

# Compute and plots heatmap & dendrogram. it shows the correlation in the right panel
  #norm,corr,dist=analyse.all_corr(vec)

  print 'statrting to cluster...'
  fig = plt.figure(figsize=(8,8))
  ax1 = fig.add_axes([0.09,0.1,0.2,0.8])
  #z=fastcluster.linkage(dist, method='ward')
  v = robjects.FloatVector([ element for col in vec.T for element in col ])
  m = robjects.r['matrix'](v, ncol = len(vec.T), byrow=False)
  Rnormalized_matrix = preprocessCore.normalize_quantiles(m)
  normalized_matrix = numpy.array( Rnormalized_matrix)
  norm,corr,dist=analyse.all_corr(vec)
  z=fastcluster.linkage(norm,metric='correlation', method="ward")
  #z=fastcluster.linkage(dist, method='complete')
  print 'clustering done, drawing the dendogram'
  Z1 = sch.dendrogram(z, labels=x,orientation='right')


  plt.yticks(fontsize=8)
  #ax1.set_yticks([])
  ticks = ax1.get_xticks() #/ max(ax1.get_xticks())
  ticks=map(float,ticks)
  ticks = ['%.2f' % (a/2.) for a in ticks]
  ax1.set_xticklabels(ticks)
  
# Plot distance matrix.
  axmatrix = fig.add_axes([0.4,0.1,0.5,0.8])
  axmatrix.set_xticks([])
  axmatrix.set_yticks([])

  axmatrix.xaxis.tick_top()
  #axmatrix.set_frame_on(False)
  idx1 = Z1['leaves']
  idx2 = Z1['leaves']
  xx=[]
  for i in idx1:
      xx.append(x[int(i)]) 
  D = corr[idx1,:]
  D = D[:,idx2]

  im = axmatrix.pcolor(D,  cmap=plt.cm.RdYlBu)
  #im = axmatrix.pcolor(D,  cmap=plt.cm.RdYlBu,edgecolor='w',)
  plt.xticks(fontsize=5)
  plt.yticks([])
  del norm
  del dist
  del corr
  del D
  xx=[]
  for i in idx1:
      xx.append(x[int(i)])   
   
  
  #plt.yticks(np.arange(len(x)),xx,fontsize = 12)
  plt.xticks(np.arange(len(x)),xx)
  plt.xticks(rotation=90)
  plt.xticks(fontsize=8)

  axcolor = fig.add_axes([0.91,0.1,0.02,0.1])
  plt.colorbar(im, cax=axcolor)
  fig.show()  
  if len(track)>0:
		track=track[0]
		labels=track[0]
		track=track[1]
		for i in xrange(len(labels)):
			labels[i]=labels[i].replace('_','')
		category=[]
		print len(xx), len(labels)
		#print labels
		for i in xx:
		   try:
			   ind=labels.index(i)
		   except:
			   print i, 'not in labels'
		   category.append(track[ind])
		category=np.array(category)

		category=np.array([category,category])
		axmatrix = fig.add_axes([0.35,0.1,0.02,.87])
		axmatrix.set_frame_on(False)
		im=axmatrix.pcolor(category.T)
		plt.xticks([])
		plt.yticks([])
		return xx , category[0]

  return xx
  
def rep_libid(ids,libid,col):
   xx=[]
   for i in ids:
      try:
         xx.append(i[:i.index('.G.A')])
         #print xx
      except:
         xx.append(i)
   for tmp in libid[1:]:
     if len(tmp[0])>0 :
       for i in xrange(len(ids)):
            if tmp[0] in ids[i]: 
              xx[i]=tmp[col]
              #print tmp[0],xx[i],tmp[col]
              xx[i]=xx[i].strip().strip('"')
              break
   return xx
   
   
def heatmap_vec( x, y , vec ):


# Compute and plot first dendrogram.

  vec=np.array(vec)
  mn=np.mean(vec,axis=1)

  mat=[]
  ylb=[]
  for i in xrange(len(mn)):
    if mn[i]>0.:
      mat.append(vec[i])
      ylb.append(y[i])

  mat=np.array(mat)
  norm,corr,dist=analyse.all_corr(mat.T)
  
  del corr
  fig = plt.figure(figsize=(8,8))
  ax1 = fig.add_axes([0.09,0.1,0.25,0.6])
  print 'fastcluster...'
  z=fastcluster.linkage(dist , method='complete')
  del dist
  #print 'dendogramming...'
  Z1 = sch.dendrogram(z,orientation='right')
  ticks = ax1.get_xticks()
  ticks=np.array(ticks)
  ticks/=2.
  ticks = ['%.1f' % a for a in ticks]
  ax1.set_xticklabels(ticks)
  idx1 = Z1['leaves']
  yy=[]
  for i in idx1:
      yy.append(ylb[int(i)]) 
  ax1.set_yticks(range(len(yy)),yy)
  
  if len(yy)<20:
    ax1.set_yticklabels(yy,fontsize=12)
  elif len(yy)<50:
   ax1.set_yticklabels(yy,fontsize=6)
  elif len(yy)<150:
    ax1.set_yticklabels(yy,fontsize=4)
  elif len(yy)<250:
    ax1.set_yticklabels(yy,fontsize=3)
  elif len(yy)<500:
    ax1.set_yticklabels(yy,fontsize=2)
  elif len(yy)<1500:
    ax1.set_yticklabels(yy,fontsize=1)
  else:
   ax1.set_yticklabels(yy,fontsize=.2)
# Plot distance matrix.
  axmatrix = fig.add_axes([0.4,0.1,0.5,0.6])
  D = norm[idx1,:]
  D=D[::-1,:]
  
  im = axmatrix.matshow(D, aspect='auto', origin='lower',cmap='RdYlBu', alpha=0.8,vmin=0)

  plt.xticks(np.arange(len(x)),x)
  plt.xticks(rotation=90)
  mytemplate(D)
  plt.xticks(fontsize=6)
  axmatrix.set_yticks([])

  #print x
# Plot colorbar.
  axcolor = fig.add_axes([0.91,0.3,0.01,0.4])
  plt.colorbar(im, cax=axcolor)

  return   yy



def heat_rna(fl,genelist,field,libfile):
    
    data=read.read_dat(fl)
    data=data[:-1]
    if [''] in data:
       data.delete([['']])
    lib=[]
    lib2=[]
    gene=[]
    for i in data:
      try:
        if i[2] not in gene:
          gene.append(i[2])
      except:
          print i

    libstable=read.read_dat(libfile,'\t')
    for i in data:
      try: 
        if i[0]+'-'+i[1] not in lib:
          lib.append(i[0]+'-'+i[1])
          lib2.append(i[0]+'-'+i[1])

        if len(libfile)>0:
         for ID in libstable:
              if i[0] in ID:
                lib2[-1]=ID[0]+'-'+ID[field]
                break 
      except:
          print i    

    pheno=[]
    pnum=-1.5
    hist=''
    for i in lib:
       j=i.split('-')
       if j[1] != hist:
               pnum+=1
       pheno.append(pnum)
       hist=j[1]
       #print hist,pnum

    mat=np.zeros((len(gene),len(lib)),np.float)
    for i in data: 
       mat[gene.index(i[2]),lib.index(i[0]+'-'+i[1])]=np.float(i[3])
    proj=np.dot(mat,pheno)
    for i in xrange(len(proj)):
       norm=np.sqrt(np.dot(mat[i,:],mat[i,:]))
       if norm!=0.:
         proj[i]/=norm
    #proj=np.array(proj)
    inds=sortedinds=proj.argsort()
    mat=mat[inds]
    gene=np.array(gene)
    gene=gene[inds]
    inds=np.any(mat != 0, axis=1)
    mat=mat[inds]
    gene=np.array(gene)
    gene=gene[inds]
    #for i in mat:
    #  if np.mean(abs(i))==0:
    #   print i
    nm=norm_max(mat)
    #for i in nm:
    #  if np.mean(abs(i))==0:
     #  print i
    genex=gene
    if len(genelist)>0:
       lb=ens_genes(gene,genelist)
       genex=lb
    else:
      return nm,lib2,gene,genex,lib
    
    b=plt.matshow(nm,aspect='auto',cmap='RdYlBu')
    if len(libfile)>0:
          plt.xticks(range(len(lib2)),lib2)
    else:
         plt.xticks(range(len(lib)),lib)
    plt.colorbar()
    plt.xticks(rotation=90)
    plt.yticks(range(len(genex)),genex)

    plt.yticks(fontsize=8)
    mytemplate(nm)
    plt.xlabel(lb)
    lb=fl.split('/')[-1][:-5]
 

    return nm,lib2,gene,genex,lib


def heat_density(dr,mark,libs,ensemble,*genelist):
    gene=ensemble
    mat=read.readall_genescore(dr,mark,libs,gene)

    nm=norm_max(mat)

    b=plt.matshow(nm,aspect='auto',cmap='RdYlBu', alpha=0.8)
    

    plt.xticks(range(len(libs)),libs)
    plt.colorbar()
    plt.xticks(rotation=90)    
    genex=gene
    #print len(genelist)
    if len(genelist)>0:
       genelist=genelist[0]
       lb=ens_genes(gene,genelist)
       genex=lb

    plt.yticks(range(len(genex)),genex)
    lb=mark+' Promoter Density' 
    plt.xlabel(lb)       


    mytemplate(nm)

    return nm


def ens_genes(gene,genelist):
     
     hugo=[]
     for g in gene:
          hugo.append(g)
     for i in xrange(len(gene)):
         for j in xrange(len(genelist)): 
           try:
            if gene[i] in genelist[j][1]+genelist[j][0]:
              hugo[i]=genelist[j][0]
              #print hugo[i]
              break
           except:
               pass
     
     return hugo


def mytemplate(nm):
    ax = plt.gca()
    for t in ax.xaxis.get_major_ticks():
      t.tick1On = False
      t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
     t.tick1On = False
     t.tick2On = False
    ax.set_frame_on(False)
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    if len(nm[:,0])>150:
       plt.yticks(fontsize=2)
       plt.xticks(fontsize=2)
    elif len(nm[:,0])>100:
       plt.yticks(fontsize=4)
       plt.xticks(fontsize=4)
    elif len(nm[:,0])>50:
       plt.yticks(fontsize=5)
       plt.xticks(fontsize=5)
    elif len(nm[:,0])>40:
       plt.xticks(fontsize=6)              
       plt.yticks(fontsize=6)
    elif len(nm[:,0])>15:
       plt.xticks(fontsize=12)              
       plt.yticks(fontsize=12)       
       
       

     

def TSS_profile_genes(tssdata,genes):
   
   list=[]
   for i in tssdata:
      for gene in genes:
        if gene[0] in i[0] and len(i)==len(tssdata[0]):
          list.append(np.array(i[1:]))
          break
   list=np.array(list)
   err=np.std(list,axis=0)
   list=np.mean(list,axis=0)
   a=list
   tmp=np.mean(a[0:10])
   a=a/tmp
   err=err/(tmp)
   
   return a,err
   
def TSS_profile_all(tssdata):
   
   list=[]
   for i in tssdata:
      if len(i)==len(tssdata[0]):
          list.append(np.array(i[1:]))

   list=np.array(list)
   err=np.std(list,axis=0)
   list=np.mean(list,axis=0)
   a=list
   tmp=np.mean(a[0:10])
   a=a/tmp
   err=err/(tmp)
   
   return a,err

def rpkm_pd_corr(rpkmmat,pdmat,librpkm,libpd):
  

  corr=[]

  for ind in range(len(rpkmmat[:,0])):
    list1=[]
    list2=[]
    for j in libpd:
      if j in librpkm:
        indlibrpkm=librpkm.index(j)
        indlibpd=libpd.index(j)
        list1.append(rpkmmat[ind,indlibrpkm])
        list2.append(pdmat[ind,indlibpd]) 
        #print len(list1),len(list2)
        try:
         if len(list1)>1:
           print list1
           tmp=ss.pearsonr(list1,list2)
           corr.append(tmp[1])
        except:
           print list1
           print list2
   
  if len(rpkmmat)==0:
     return 0.,0.
  #print 'aha'
  return np.mean(corr), np.std(corr)


def norm_max(mat):

    if len(mat)==0.:
      return mat
    nm=0.*mat
    
    for i in range(len(mat[:,0])):
       a=1.*mat[i,:]
       if len(a)>0 and len(a[a==0.]) != len(a):
          zscore=ss.zscore(a)
          nm[i,:]=zscore
 
    return nm

def data2arr(data):

   
   for row in data:
      for i in row:
		i=float(i)

      row=np.array(row)
   data=np.array(data)
   
   return data



def rpkm_mat(rpkmdata):
   
   x=len(rpkmdata)
   y=len(rpkmdata[0])
   
   mat=np.zeros((x,y),np.float)
   
   rpkmmat=[]
   genes=[]
   for gene in rpkmdata[0]:
      genes.append(gene[0])
   genes=np.array(genes)

   i=-1
   for rpkmlib in rpkmdata:
      i+=1
      for rpkmtmp in rpkmlib:
           #print rpkmtmp[0]
           try:
              j=bisect.bisect_left(genes,rpkmtmp[0])
              mat[i,j]=rpkmtmp[2]
              #print genes[j], rpkmtmp[0]
           except:
            pass
        

   return mat, genes

def gene_rpkm_lineage(hugo,rpkmdata,rpkmfiles,allgenes,libs,*plotbool):
   
   genelist=[]
   for i in range(len(hugo)):
      genelist.append(hugo[i])
   genelist=ens_genes(genelist,allgenes)
   rpkmls=[]
   files=[]
   for gene in genelist:
     rpkmls.append([])
     for ind in xrange(len(rpkmdata)):
       for i in rpkmdata[ind]:
          if gene in i[0]:
            rpkmls[-1].append(i[2])
            if len(rpkmls)==1:
              files.append(rpkmfiles[ind][:rpkmfiles[ind].index('.')])
            break

   files=np.array(files)
   rpkmls=data2arr(rpkmls)
   sortedinds=rpkmls[0].argsort()
   files=files[sortedinds]
   files=files[::-1]
   for i in range(len(genelist)):
      rpkmls[i]=rpkmls[i][sortedinds]
      rpkmls[i]=rpkmls[i][::-1]
   normalfiles=[]
   normalrpkms=[]
   RTfiles=[]
   RTrpkms=[] 
   libs=read.read_dat(libs,'\t')
   for ind in range(len(files)): 
     for lib in libs:
        if files[ind] in lib[0]:
            files[ind]=lib[1]+'-'+lib[-1]
            break
   for i in range(len(genelist)):
     rpkmls[i]=np.log(rpkmls[i]+0.001)/np.log(10)
     
   for i in range(len(genelist)):
      normalrpkms.append([])
      RTrpkms.append([])
   for ind in range(len(files)):
      if 'Cancer' not in files[ind]  and 'ES' not in files[ind]:
         print files[ind]
         normalfiles.append(files[ind])
         for i in range(len(genelist)):
            normalrpkms[i].append(rpkmls[i][ind])
      elif 'RT'  in files[ind]:
         RTfiles.append(files[ind])
         for i in range(len(genelist)):
            RTrpkms[i].append(rpkmls[i][ind])      
   
   plt.figure()
   for i in range(len(genelist)):
     if i<7:
       plt.plot(xrange(len(list(files))),rpkmls[i],'s-',label=hugo[i])
     else:
       plt.plot(xrange(len(list(files))),rpkmls[i],'^-',label=hugo[i])
   plt.legend(loc=1)
   plt.xticks(range(len(files)),files)
   plt.xticks(rotation=90)
   plt. xticks(fontsize=8)
   plt.ylabel('$log_{10}$ $RPKM$')
   plt.grid()
   plt.show()

   def heatmap(rpkmls,hugo,txt):
     rpkmls=data2arr(rpkmls)
     norm,matrix,dist=analyse.all_corr(rpkmls.T)
     plt.matshow((matrix),cmap='RdYlBu',vmax=1,vmin=-1)
     plt.yticks(range(len(hugo)),hugo)
     plt.xticks(range(len(hugo)),hugo)
     plt.xticks(rotation=90)
     plt.xlabel('Correlation-'+txt)
     plt.colorbar()
     #plt.matshow(np.log(pval),cmap='PuBu',vmax=-5,vmin=-25)
     #plt.yticks(range(len(hugo)),hugo)
     #plt.xticks(range(len(hugo)),hugo)
     #plt.xlabel('P-Value $\log$ -'+txt)
     #plt.xticks(rotation=90)
     #plt.colorbar()
   if len(plotbool)>0:
     if plotbool[0]:
       heatmap(rpkmls,hugo,'All Samples')
       heatmap(normalrpkms,hugo,'Normal Samples')
       heatmap(RTrpkms,hugo,'RT Samples')
       for i in range(len(hugo)):
        print hugo[i],ss.ks_2samp(normalrpkms[i],RTrpkms[i])[1]*np.float(len(allgenes)), 'KS'
        print hugo[i],ss.ttest_ind(normalrpkms[i],RTrpkms[i])[1]*np.float(len(allgenes)), 'T-test'
   
   
   #fig=plt.figure()
   #plt.plot()
   return rpkmls,files,data2arr(normalrpkms),np.array(normalfiles),genelist
   
   
   
def gene_rpkm_corr_compare(thegene,thersh,rpkmdata,rpkmfiles,libs,*enslist):
   #usage:
   #genescompared,genesfc,pvals=\
   #analyse.gene_rpkm_corr_compare('ENSG00000108799',0.000001,rpkm,files,'rhabdoid/RNA-ALL-key.txt','resources/list.genes2.txt')
   genelist=[]
   rpkmgene=[]
   geneimp=[]
   pvals=[]
   files=[]
   for ind in xrange(len(rpkmdata)):
       for i in rpkmdata[ind]:
          if thegene in i[0]:
            rpkmgene.append(i[2])
            files.append(rpkmfiles[ind][:rpkmfiles[ind].index('.')])
            break



   libs=read.read_dat(libs,'\t')
   for ind in range(len(files)):
     for lib in libs:
        if files[ind] in lib[0] :
            files[ind]=lib[1]+'-'+lib[-1]
            break

   rpkmgeneRT=[]
   rpkmgenenormal=[]
   for ind in xrange(len(rpkmgene)):
          #if 'Blood' not in files[ind]:
            #print files[ind]
            if 'RT' in files[ind]:
               rpkmgeneRT.append(rpkmgene[ind])
            else:
               rpkmgenenormal.append(rpkmgene[ind]) 
   rpkmgeneRT=np.array(rpkmgeneRT)
   rpkmgenenormal=np.array(rpkmgenenormal)
   print len(rpkmgeneRT), len(rpkmgenenormal)

   for igene in range(len(rpkmdata[0])):
         rpkmRT=[]
         rpkmnormal=[]
         for ind in xrange(len(rpkmdata)):
            if 'Blood' not in files[ind]:
               if 'RT' in files[ind] :
                 rpkmRT.append(rpkmdata[ind][igene][2])
               else:
                 rpkmnormal.append(rpkmdata[ind][igene][2])
         
         #print len(rpkmRT)
         #print len(rpkmnormal), rpkmnormal[0]
         rpkmRT=np.array(rpkmRT)
         rpkmnormal=np.array(rpkmnormal)
         a=ss.pearsonr(rpkmRT,rpkmgeneRT)
         
         b=ss.pearsonr(rpkmnormal,rpkmgenenormal)
         mn=np.mean(rpkmRT)/np.mean(rpkmnormal)

         if b[1]<thersh  and (a[1]/b[1])>thersh:
            genelist.append(rpkmdata[0][igene][0])
            if mn>1. or mn <1.:
               #print rpkmdata[0][igene][0], 'correlation for normal=%.3f, pvalue normal=%.5f, pvalue RT=%.5f, foldchange=%.3f'%(b[0],b[1],a[1],mn)
               geneimp.append(mn)
               pvals.append(np.log(a[1]/b[1]))
   
        
   return genelist , geneimp, pvals
   
   
def gene_rpkm_compare(thersh,rpkmdata,rpkmfiles,libs,*enslist):

   genelist=[]
   files=[]
   for ind in xrange(len(rpkmdata)):
            files.append(rpkmfiles[ind][:rpkmfiles[ind].index('.')])

   libs=read.read_dat(libs,'\t')
   for ind in range(len(files)):
     for lib in libs:
        if files[ind] in lib[0]:
            files[ind]=lib[1]+'-'+lib[-1]            
            break
   pval=[]
   foldchange=[]
   for igene in range(len(rpkmdata[0])):
         rpkmRT=[]
         rpkmnormal=[]
         for ind in xrange(len(rpkmdata)):
            if 'RT' in files[ind]:
               rpkmRT.append(rpkmdata[ind][igene][2])
            elif 'Cancer' not in files[ind]:
               rpkmnormal.append(rpkmdata[ind][igene][2])
         rpkmRT=np.array(rpkmRT)
         #rpkmRT=np.log(rpkmRT)
         rpkmnormal=np.array(rpkmnormal)
         #rpkmnormal=np.log(rpkmnormal)
         a=ss.ks_2samp(rpkmRT,rpkmnormal)

         if a[1]<thersh/np.float(20000) and np.mean(rpkmnormal)+np.mean(rpkmRT)>1.:
            genelist.append(rpkmdata[0][igene][0])
            #print rpkmdata[0][igene][0], 'pvalue RT=',(a[1])
            pval.append(-np.log(a[1]*np.float(20000)))
            #if np.mean(rpkmnormal)<0:
               #print rpkmnormal
            #rpkmRT=np.array(rpkmRT)
            #rpkmnormal=np.log(rpkmnormal)
            foldchange.append(np.log(np.median(rpkmnormal)/np.median(rpkmRT)))
   print 'RT sample:',len(rpkmRT), ', Normal sample:',len(rpkmnormal)
   if len(enslist)>0:
         enslist=read.read_dat(enslist[0])
         genelist2=ens_genes(genelist,enslist)
         #print genelist2
   
   return genelist, foldchange, pval



def beds_corr(beds,motif,flank):
   corr=[]
   pval=[]
   motifarr={}
   for chrom in motif:
      motifarr[chrom]=[]
      for i in motif[chrom]:
         motifarr[chrom].append(i[0])
      motifarr[chrom]=np.array(motifarr[chrom])
         
         
   for samp in beds:
     x=[]
     for chrom in samp:
       for i in samp[chrom]: 
          #print chrom
          #print motif.keys()
          ind = bisect.bisect_left(motifarr[chrom],int(i[0]))
          j=motif[chrom][ind]
          dist=abs(i[0]-j[0]) +abs(i[1]-j[1])
          if dist<flank:
                x.append([i[-1],j[-2]])

     #print i
     #print j
     x=np.array(x)
     #print len(x)
     #print len(x[0])
     corr.append(ss.pearsonr(x[:,0],x[:,1])[0])
     pval.append(ss.pearsonr(x[:,0],x[:,1])[1])
   
   return corr,pval
   
def methyl_bed(bed):
   
   methylbed=[]
   i=0
   while i <(len(bed)):
        try:
         gene=bed[i][3]
        except:
         print bed[i]
        cpg=0
        methyl=0.
        cpg=0
        methyl=0
        tmpgene=gene
        while tmpgene==gene and i<len(bed)-1:
             i+=1
             try:
               tmpgene=bed[i][3]
               cpg+=1
               methyl+=np.float(bed[i][-3])
             except:
               print bed[i]
        try:
          methylbed.append(bed[i-1][:4])
          methylbed[-1].append(cpg)
          methylbed[-1].append(methyl/np.float(cpg))
        except:
          pass
        i+=1
   return methylbed


def marks_corr(methyl,mark):
   
   
   
   cpg=[]
   fraction=[]
   markvec=[]
   for i in methyl:
      for j in mark:
         if len(i)>3 and len(j)>3 and i[3]==j[3]:
            cpg.append(i[-2])
            fraction.append(i[-1])
            markvec.append(j[-2])
            break

   cpg=np.array(cpg)
   fraction=np.array(fraction)
   markvec=np.array(markvec)
   return ss.pearsonr(cpg,markvec)[0],ss.pearsonr(fraction,markvec)[0],ss.pearsonr(cpg,fraction)[0]
   

def marksbed_corr(methyl,markbed):
   
   
   
   cpg=[]
   fraction=[]
   markvec=[]
   for i in methyl:
      chrom=str(i[0])
      if 'chr' not in chrom:
         chrom='chr'+chrom

      if chrom in markbed.keys():
         ind=bisect.bisect_left(markbed[chrom][:,0],i[1])
         xi=markbed[chrom][ind][0]
         xf=markbed[chrom][ind][1]

         xi=max(xi,i[1])
         xf=min(xf,i[2])
         if xi<=xf:
             cpg.append(i[-2])
             fraction.append(i[-1])
             markvec.append(markbed[chrom][ind][-1])


   cpg=np.array(cpg)
   fraction=np.array(fraction)
   markvec=np.array(markvec)
   print len(cpg)
   return ss.pearsonr(cpg,markvec)[0],ss.pearsonr(fraction,markvec)[0],ss.pearsonr(cpg,fraction)[0]

def allsamps_marks_corr(allmethyl,allmark,methylfiles,markfiles,mark,libsdata):
   
   kfiles=[]
   mfiles=[]
   col=libsdata[0].index(mark)
   for i in range(len(markfiles)):
      #print i, markfiles[i]
      kfiles.append(markfiles[i])
      for j in libsdata[1:]:
        if  j[col]!='' and j[col] in  markfiles[i]:
            kfiles[i]=j[0]
            break


   col=libsdata[0].index('methylation')
   for i in range(len(methylfiles)):
      mfiles.append(methylfiles[i])
      for j in libsdata[1:]:
        if j[col]!='' and j[col] in methylfiles[i] :
            mfiles[i]=j[0]
            break
   ans=[]
   fls=[]
   methylbed=allmethyl
   for i in range(len(mfiles)):
    for j in range(len(kfiles)):
        if kfiles[j] in mfiles[i]:
            #methylbed.append(analyse.methyl_bed(allmethyl[i]))
            #try:
              print kfiles[j]
              ans.append(analyse.marksbed_corr(methylbed[i],allmark[j]))
              fls.append(kfiles[j])
            #except:
            #   pass
              break
   
   
   
   plt.plot(ans,'.-')
   plt.xticks(range(len(fls)),fls)
   plt.xticks(rotation=90)
   plt.ylabel('correlation')
   plt.legend(('CpG, '+mark,'fractional methylation, '+mark,'CpG, methylation'),loc=0)

   return ans, fls

    
def plot_bed(bed,c,*chrom):
    
    i=0.
    y=[]
    x=[]
    if len(chrom)>0:
       chrom=chrom[0]
       bed[chrom]=bed[chrom][bed[chrom][:,0]<10000000]
       x=list(np.ones((len(bed[chrom][:,1])),np.float))
       y=list((bed[chrom][:,1]+bed[chrom][:,0])/2)
       plt.plot(y,x,c+'.')
       plt.yticks(range(3),['',chrom,''])
    else:
      for chrom in bed.keys():
        i+=1.
        x+=list(i*np.ones((len(bed[chrom][:,1])),np.float))
        y+=list((bed[chrom][:,1]+bed[chrom][:,0])/2)
      plt.plot(y,x,c+'.')
         
      plt.yticks(range(len(bed.keys())+2),['']+bed.keys()+[''])
    

def stat_bed(bed):
   ave=0
   vec2=0
   sum=0
   for chrom in bed.keys():
      a=bed[chrom][:,1]-bed[chrom][:,0]
      ave+=np.dot(a,bed[chrom][:,2])
      vec2+=np.dot(a,bed[chrom][:,2]**2)
      sum+=np.sum(a)
   ave/=sum
   vec2/=sum
   vec2-=ave**2

   
   return ave,np.sqrt(vec2)
   

def bedintersect_corr(bed):

   cpg=[]
   fraction=[]

   for chrom in bed.keys():
      cpg+=list(bed[chrom][:,-1])
      fraction+=list(bed[chrom][:,-2])

   cpg=np.array(cpg)
   fraction=np.array(fraction)

   print len(cpg)
   return np.mean(fraction)/100.,ss.pearsonr(cpg,fraction)[0]


def corr_vec_rpkm(rpkmmat,ensmble,vec,col):
   
   data=[]
   idvec=[]
   
   for i in vec:
     try:
       data.append(np.int(i[1]))
     except:
       print i
       #data.append(np.float(i[1]))
     idvec.append(i[0])
   data=np.array(data)
   idvec=np.array(idvec)
   sorting=data.argsort()
   idvec=idvec[sorting][::-1]
   data=data[sorting][::-1]
   genenum=len(rpkmmat[0])
   genes=[]
   for igene in xrange(genenum):
         genes.append(rpkmmat[0][igene])
   samples=len(rpkmmat)
   ids=[]
   id_i=[]

   for tmp in range(len(idvec)):
      for i in range(1,samples):
         if idvec[tmp][:16] in rpkmmat[i][0]:
           ids.append(rpkmmat[i][0])
           id_i.append(i)
           idvec[tmp]=idvec[tmp][:16]+'_'+str(data[tmp])
           break
   rpkms=[]

   cors=[]
   for itmp in xrange(len(ensmble)):
      igene=genes.index(ensmble[itmp])
      rpkm=[]
      for i in id_i:
        rpkm.append(rpkmmat[i][igene])
   #   print rpkm
      #print data
      cor=ss.pearsonr(rpkm,data)
      rpkms.append(rpkm)
      cors.append(cor)
         

   rpkms=np.array(rpkms)
   
   
   return rpkms, idvec, data, cors
   
   
   
   

def CPG_RPKM(RPKM, CPG, lim):


  CGI=[]
  cpg=read.read_dat(CPG,'\t')
  for i in cpg:
    CGI.append(i[0])
  genes=RPKM[0][1:]
  libs=[]
  mat=[]
  for i in RPKM[1:]:
      mat.append(i[1:])
      libs.append(i[0])
  mat=analyse.data2arr(mat)
  genes=np.array(genes)
  allave=[]
  cgiave=[]
  print len(mat)
  for i in range(len(mat)):
      gntmp=genes[mat[i,:]>lim]

      tmp=mat[i,:][mat[i,:]>lim]
      allave.append(np.mean(tmp))
      cgirpkm=[]
      for j in xrange(len(gntmp)):
          if gntmp[j] in CGI:
              cgirpkm.append(tmp[j])
      cgiave.append(np.mean(cgirpkm))
      print len(gntmp), len(cgirpkm)
  allave=np.array(allave)
  cgiave=np.array(cgiave)
  return allave,cgiave, libs





def data2matrix(T,NA):
	x=T[0][1:]
	y=[]
	mat=[]
	for i in xrange(len(T)):
		for j in xrange(len(T[0])):
			if T[i][j]=='NA' or T[i][j]=='N/A':
				T[i][j]=NA
	for i in T[1:]:
		y.append(i[0])
		mat.append(map(float,i[1:]))
	#print mat[0]
	mat=analyse.data2arr(mat)


	return x,y,mat










