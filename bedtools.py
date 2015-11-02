#!python2.7.1
import  bisect, heapq, itertools,random
import scipy.stats as ss
import numpy as np
#from joblib import Parallel, delayed
#from cloud.serialization.cloudpickle import dumps
from operator import itemgetter

def intersect_beds(bed1, bed2, flank):
    
    intrsct = {}
    nonint= {}
    for chrom in bed1.keys():
        if bed2.has_key(chrom):
            intrsct[chrom] = []
            nonint[chrom] = []
            for loc1 in list(bed1[chrom]):
                begin = loc1[0]
                end = loc1[1]-1
                sc1 = loc1[2]
                
                tmp=len(bed2[chrom])
                indloc2=bisect.bisect_left(bed2[chrom][:,0],begin-flank)
                if indloc2<tmp-1:
                  loc2=bed2[chrom][indloc2,:]
                  xi = loc2[0]-flank
                  xf = loc2[1]+flank-2
                  sc2 = loc2[2]
                  xL = max(begin, xi)
                  xR = min(end, xf)
                  #print xL,xR
                  if xR >= xL:
                       intrsct[chrom].append([xL,xR,max(sc1,sc2)])
                  else:
                       nonint[chrom].append([begin,end,0])
            #print chrom, len(intrsct[chrom]), ch,true
            intrsct[chrom]=np.array(intrsct[chrom])
            nonint[chrom]=np.array(nonint[chrom])
            print chrom

    return  intrsct, nonint


def bin_wig(wig,res):
   
   out = {}
   for chrom in wig.keys():
       tmpwig=wig[chrom]
       minn=tmpwig[0,0]
       maxx=tmpwig[-1,0]
       ln=(maxx)/res
       tmp=np.zeros((ln,3),np.float)
       if not out.has_key(chrom):
          out[chrom] = tmp
       for i in xrange(maxx/res-1):
          begin=i*res
          end=i*res+res
          pos1=bisect.bisect_left(tmpwig[:,0],begin)
          pos2=bisect.bisect_left(tmpwig[:,1],end)
          correction=(begin-tmpwig[pos1+1,0])*tmpwig[pos1+1,2]+(end-tmpwig[pos1,0])*tmpwig[pos1,2]
          last=max(1,(np.sum(tmpwig[pos1:pos2,2])+correction)/np.float(res))
          tmp[i,:]=np.array([begin,end,last])
       if tmp[-1,-1]==0.:
         tmp[-1,-1]=1.
       out[chrom]=tmp
   return out
   
def len_wig_thresh(wig,thresh):



    ln=0
    bed={}
    wig=filter_bed_chrom(wig)
    for chrom in wig.keys():
       bed[chrom]=[]
       inds=wig[chrom][:,2]>thresh
       tmp=wig[chrom][inds]
       sm=tmp[:,1]-tmp[:,0]
       ln+=np.sum(sm)
       tmp=wig[chrom][0,2]
       
       itmp=-1
       print len(inds), chrom
       while itmp<len(inds)-1 :
          itmp+=1
          bool=False
          bin=0
          score=0
          if inds[itmp]:
              bed[chrom].append([])
              begin=wig[chrom][itmp,0]
              bool=True
          while inds[itmp] and itmp<len(inds)-1:
              end=wig[chrom][itmp,1]
              itmp+=1
              xi=wig[chrom][itmp,0]
              score+=wig[chrom][itmp,2]*(end-xi)
          if bool:
              bed[chrom][-1]=[begin,end,score/np.float(end-begin)]
       bed[chrom]=np.array(bed[chrom])
          
          
    return ln, bed


def filter_bed_chrom(bed):

   ls=range(23)
   chrom=map(str,ls)
   for i in ls:
      chrom[i]='chr'+chrom[i]
   chrom.append('chrX')
   chrom.append('chrx')
   chrom.append('chry')
   chrom.append('chrY')
   for key in bed.keys():
     if key not in chrom:
         del bed[key]
   
   return bed


def filter_bed(bed,score):

    out = {}
    rank=[]
    for chrom in bed.keys():
          if not out.has_key(chrom):
              out[chrom] = bed[chrom][bed[chrom][:,2]>score]
              rank+=list(bed[chrom][:,2])
    rank=np.array(rank)
    inds=rank.argsort()
    rank=rank[inds][::-1]

    return out, rank

    
def filter_bed_num(bed,num):

    out={}
    a=[]
    tup=[]
    for key in bed.keys():
       out[key]=[]
       a=bed[key][:,2]
       tmp=(heapq.nlargest(num, zip(a, itertools.count())))
       for i in xrange(num):
           try:
             tmp[i]=tmp[i]+(key,)
             tup.append(tmp[i])
           except:
              a=a


    tup.sort(key=lambda tup:tup[0])
    tup=tup[::-1][:num]
    indices=(zip(*tup)[1])
    allkeys=(zip(*tup)[2])

    for i in xrange(num): 
       try:
         ind=indices[i]
         key=allkeys[i]
         out[key].append(bed[key][ind,:])
       except:
         pass
    for key in out.keys():
       out[key]=np.array(out[key])
    return out
    
def filter_bed_pi(bed,pi):

    out={}
    a=[]
    tup=[]
    for key in bed.keys():
       out[key]=[]
       a=bed[key][:,2]
       tmp=(heapq.nlargest(num, zip(a, itertools.count())))
       for i in xrange(num):
           try:
             tmp[i]=tmp[i]+(key,)
             tup.append(tmp[i])
           except:
              pass


    tup.sort(key=lambda tup:tup[0])
    tup=tup[::-1][:num]
    indices=(zip(*tup)[1])
    allkeys=(zip(*tup)[2])

    for i in xrange(num): 
       try:
         ind=indices[i]
         key=allkeys[i]
         out[key].append(bed[key][ind,:])
       except:
         a=a
    for key in out.keys():
       out[key]=np.array(out[key])
    return out





    
def compare_bed(bed1,*bed2):
    
    
    out={}

    if len(bed2)>0:
      bed2=bed2[0]
      for key in bed2.keys():
       tmp=bed1[key]
       tmp2=bed2[key]
       div=0.*tmp[:,2]
       #div+=1.
       a=tmp[:,2]
       b=tmp2[:,2]
       div=a*0.
       div[(a!=0.) * (b==0.)]=1.
       div[(a==0.) * (b!=0.)]=1.
       div[(a!=0.) * (b!=0.)]=a[(a!=0.) * (b!=0.)]/b[(a!=0.) * (b!=0.)]
       out[key]=np.array([tmp[:,0],tmp[:,1],div]).T
    else:
      for key in bed1.keys():
         tmp=bed1[key]     
         div= tmp[:,2]*0.
         a=(tmp[:,2]>0.)
         div[a]=np.log(tmp[:,2][a])
         out[key]=np.array([tmp[:,0],tmp[:,1],div]).T
    #print len(out[key][:,1])
    return out

def Compare_all_bed(beds1,*beds2):
    
   if len(beds2)>0:
     beds2=beds2[0]
     compared=[]
     for j in beds2:
      compared.append([])
      #print len(beds1)
      for i in beds1:
        compared[-1].append(compare_bed(i,j))
        
       
   else:
     compared=[]
     for j in beds1:
        compared.append([])
        compared[-1].append(compare_bed(j))      
   
   
   return compared
   
   

def Compare_all_bed_pi(beds1,beds2,pival):
    

     compared={}
     for key in beds1[0].keys():
        compared[key]=[]
        for i in range(len(beds1[0][key])):
          test2=[]
          test1=[]
          #print beds1[0][key][i][-1]
          for j in range(len(beds1)):
            test1.append(beds1[j][key][i][-1])
          for j in range(len(beds2)):
            test2.append(beds2[j][key][i][-1])
          #print test1 , test2
          test1=np.array(test1)
          test2=np.array(test2)
          
          if len(test1[test1!=0])>0 or len(test2[test2!=0])>0:
             p=ss.ttest_ind(test1,test2)[1]
             #print len(test1[test1==0]), len(test2[test2==0])
          #except:
           #  print test1,test2
          #print p

          if p<pival:
             #print test1 , test2
             compared[key].append(np.zeros(3,(np.float)))
             compared[key][-1][:-1]= beds1[0][key][i][:-1]
             compared[key][-1][-1]=p
             compared[key][-1]=np.array(compared[key][-1])
        compared[key]=np.array(compared[key])
     return compared



def bedsall_prod(beds1,beds2):
    
   prod=[]
   for j in xrange(len(beds1)):
       #print len(beds1), len(beds2), 'prod'
       prod.append(bed_product(beds1[j],beds2[j]))
   #print len(prod), 'prod'
   return prod
   
   
def bed_product(bed1,bed2):
   
   prod={}
   for key in bed1.keys():
       prod[key]=1.*bed1[key]
       prod[key][:,2]*=bed2[key][:,2]
   return prod
   
   
def merge_bed(bed,gap,domain_size):
   
   sz=domain_size
   domains={}
   for key in bed.keys():
       tmp=bed[key]
       res=tmp[0,1]-tmp[0,0]
       ln=len(tmp[:,1])
       tmpd=[[0,0,0.]]
       for ind in xrange(ln):
           st  = tmp[ind,0]
           end = tmp[ind,1]
           sc = tmp[ind,2]
           
           add = True
           for j in xrange(len(tmpd)):
                stdomain  = tmpd[j][0]
                enddomain = tmpd[j][1]
                
                if stdomain<=end+gap and stdomain>=st:
                  tmpd[j][0]=st
                  tmpd[j][2]+=sc
                  add=False
                elif enddomain+gap>=st and enddomain<=end:
                  tmpd[j][1]=end
                  tmpd[j][2]+=sc
                  add=False

           if add:
                 tmpd.append([st,end,sc])
       domains[key]=tmpd          
   num=0
   for key in domains.keys():
       tmp=domains[key]
       for elt in tmp:
         diff=elt[1]-elt[0]
         if diff<=sz-gap:
            
            domains[key].remove(elt)
         else:
           elt[2]/=np.float(diff/res)
           num+=1
   #print 'total number of domains is:', num, ' resolution=', res
   return domains, num
       


def FDR_bed(bed,num,gap,region):
   random.seed()
   def foo(tmpbed,gap,region):
        shuffle_bed(tmpbed)
        dom_rand,peaks=merge_bed(tmpbed,gap,region)
        return peaks
   
   tmpbed=bed.copy()
   domain,ns=merge_bed(bed,gap,region)
   nr=[]
   for i in range(10):
        shuffle_bed(tmpbed)
        dom_rand,peaks=merge_bed(tmpbed,gap,region)
        nr.append(peaks)
   
   nr=np.mean(nr)
   FDR=(nr+1)/(ns+1)
   print 'FDR is ',  FDR 
   
   return domain, FDR
   
def shuffle_bed(bed):   
   
   
   for key in bed.keys():
      maxx=int(max(bed[key][:,1]))
      minn=int(min(bed[key][:,0]))
      ln=len(bed[key][:,0])
      res=int(np.mean(bed[key][:,1]-bed[key][:,0]))
      x=np.array(random.sample(xrange(minn,maxx,res),ln))
      bed[key][:,0]=x
      bed[key][:,1]=x+int(res)
      #print res, key
   return bed
   
def print_bed(bed,*boolean):
	
	if len(boolean)<1:
		boolean=False
	s=0
	sm=0
	for chrom in bed.keys():
		try:
			if boolean:
				print chrom, len(bed[chrom][:,0])
			s+=len(bed[chrom][:,0])
			sm+=np.sum(bed[chrom][:,1]-bed[chrom][:,0])
		except:
			pass
	
	print 'total: ' ,s
	
	
	
	return s ,sm   


def intersect_bed_gene(bed,gene):
	
	genes={}
	for chrom in gene.keys():
		if len(bed[chrom])>0:
			for j in bed[chrom]:
				xi=np.int(j[0])
				xf=np.int(j[1])
				sc=np.int(j[2])
				a=gene[chrom][:,0]
				#b=gene[chrom][:,1]
				a=map(np.int,a)
				#b-map(np.int,b)
				i=bisect.bisect_left(a,xi)
				i=gene[chrom][i,:]
				st=np.int(i[0])
				end=np.int(i[1])
				try:
				    geneid=(i[2][:i[2].index('_')])
				except:
				    geneid=(i[2])
				xL=max(xi,st)
				xR=min(xf,end)
				if xL<=xR:
				    genes[geneid]=[st,end,sc]
				    #print i
					#	break
			#print xi,st, i
	return genes



def intersect_allgenes(ranked):
	elts={}
	for i in range(len(ranked)): 
		for geneid in ranked[i].keys():
			if not elts.has_key(geneid):
				elts[geneid]=[]
			elts[geneid].append(i)
	return elts
	
	
	
def find_common_genes(genes,num):
	
	tup=[]
	for i in genes.keys():
		tup.append((i,len(genes[i])))
	
	tmp=heapq.nlargest(num, tup, key=(itemgetter(1)))
	a=np.array(zip(*tmp)[1])
	#print min(a)
	ranked={}
	for i in tmp:
		ranked[i[0]]=genes[i[0]]
	
	

	return ranked



def filter_common_genes(genes,thresh):
	
	filtered={}
	for key in genes.keys():
		 if len(genes[key])>= thresh:
			filtered[key]=genes[key]
	
	

	return filtered


def intersect_gene_names(genes1,genes2):
	
	intersect={}
	for key1 in genes1.keys():
		for key2 in genes2.keys():
			#print key1, key2
			if key2 == key1:
				intersect[key1]=genes2[key2]+genes1[key1]
				#print key1
				break
	
	return intersect
	
	
	
def makebed_genepos(ensIDs,genespos,dist,chromsizes):
   sizes={}
   for i in chromsizes:
       if 'chr' not in i[0]:
          i[0]='chr'+i[0]
       sizes[i[0]]=np.int(i[1])
       
   ens=ensIDs
   def widen(a,b,dist,chrom):
     b=int(b)
     a=int(a)
     nowdist=b-a
     plus=(dist-nowdist)/2
     a=max(0,a-plus)
     b=min(b+sizes[chrom],plus)

   genesbed=[]
   for ens_i in ensIDs:
     found=False
     for chrom in genespos.keys():
      for j in genespos[chrom]:
           if found:
                break
           if ens_i in j[-1]:
              widen(j[0],j[1],dist,chrom)
              genesbed.append([chrom,np.int(j[0]),np.int(j[1]),np.int(j[2]),ens_i])
              found=True
              break
   return genesbed

def dist_bed(bed1,bed2):
    dist=[]
    for chrom in bed1.keys():
        if chrom in bed2.keys():
            xvec1=bed1[chrom][:,0]
            yvec1=bed1[chrom][:,1]
            xvec1+=yvec1
            xvec1/=2.
            xvec2=bed2[chrom][:,0]
            yvec2=bed2[chrom][:,1]
            xvec2+=yvec2
            xvec2/=2.
            for x1 in xvec1:
                i2=bisect.bisect_left(xvec2,x1)
                if i2>=len(xvec2):
                  x2=max(xvec2)
                else:
                  x2=xvec2[i2]

                dist.append(abs(x1-x2))
    
    
    return dist, xvec1
    
    
    
    
def randomize_bed(bed):
    
    
    bed1={}
    for chrom in bed.keys():
	xvec=bed[chrom][:,0]
	yvec=bed[chrom][:,1]
	
	res= xvec-yvec
	maxx=max(xvec)
	minn=min(xvec)
	ln=len(xvec)
	tmp1=np.random.randint(maxx,size=ln)
	tmp2=tmp1+res
	tmp3=tmp1*0.
	bed1[chrom]=np.array([tmp1,tmp2,tmp3]).T
	
	
    return bed1
	
	

def genedensity(bed,genes):
	
	intersect={}
	for gene in genes:
		chrom=gene[0]
		if chrom not in intersect.keys():
			intersect[chrom]=[]
		if chrom in bed.keys():
				xst=gene[1]
				xend=gene[2]
				name=gene[3]
				indst=bisect.bisect_left(bed[chrom][:,0],xst)
				indend=bisect.bisect_right(bed[chrom][:,1],xend)
				sc=indend-indst
				sc=np.float(sc)/np.float(xend-xst)
				intersect[chrom].append([xst,xend,sc,name])

	return intersect


def filter_multibed(multibed,frac):
    
    
    outbed={}

    for chrom in multibed.keys():
        print chrom
        for i in multibed[chrom]:
           if i[2]>frac:
              if chrom not in outbed.keys():
                outbed[chrom]=[]
              if len(outbed[chrom])>0 and np.int(i[0])<=outbed[chrom][-1][1]:
                 outbed[chrom][-1][1]=np.int(i[1])
              else:
                  outbed[chrom].append([np.int(i[0]),np.int(i[1]),np.float(i[2])])
        if chrom in outbed.keys():
          outbed[chrom]=np.array(outbed[chrom])
    
    
    return outbed
    
    
    
def chr_bed(bed):

   #print bed.keys()
   for chrom in bed.keys():
     #print chrom
     if 'chr' not in chrom:
        #print chrom
        bed['chr'+chrom]=bed.pop(chrom)
    
   return bed
    
    
    

    
    
def plot_bed(bed):
    
    i=-1.
    for chrom in bed.keys():
       i+=1.
       plt.plot(i*np.ones(len(bed[chrom][:,1]),np.float()),(bed[chrom][:,1]+bed[chrom][:,0])/2,'b.')
    plt.xticks(bed.chrom())
    
    
    
def encode_bed_tf(encode,bed):
    
    
    ch=0
    true=0
    intrsct = {}
    for chrom in encode.keys():
        if bed.has_key(chrom):
            for loc1 in encode[chrom]:
                begin = loc1[0]
                try:
		    end = loc1[1]-1
		except:
		    print loc1[1], chrom
                sc1 = loc1[2]
                tf=str(loc1[3])
                intrsct[tf]={}
                intrsct[tf][chrom]={}
                try:
                  tmp=bed[chrom][:,0]
                  i=bisect.bisect_left(tmp,begin)
                  true+=1
                  loc2=bed[chrom][i,:]
                  ch+=1
                  #print loc2
                  xi = loc2[0]-flank
                  xf = loc2[1]+flank-2
                  sc2 = loc2[2]
                  xL = max(begin, xi)
                  xR = min(end, xf)
                  #print xL,xR
                  if xR >= xL:
                       intrsct[tf][chrom].append([xL,xR])

                  #else:
                       #print xR, xL
                except:
                     pass
            #print chrom, len(intrsct[chrom]), ch,true
            intrsct[tf][chrom]=np.array(intrsct[chrom])

    return  intrsct


def print_tf_intersect(intrsct):
    
    for tf in intrsct.keys():
        print tf
        print_bed(intrsct[tf])

    
def eq_bed(bed,bp):
    
    out={}
    ipeak=0
    for chrom in bed.keys():
       out[chrom]=[]
       for i in bed[chrom]:
          ipeak+=1
          out[chrom].append([((i[0]+i[1])/2)-bp,((i[0]+i[1])/2)+bp,ipeak])
       out[chrom]=np.array(out[chrom])
    return out

       
def add_bed_chr(bed):
    
    out={}

    for chrom in bed.keys():
       if 'chr' not in chrom:
          bed['chr'+chrom]=bed.pop(chrom)

    
def print_encode(encode):

    s={}
    for tf in encode.keys():
            s[tf]=0
            for chrom in encode[tf].keys():
                    xi=encode[tf][chrom][:,0]
                    xf=encode[tf][chrom][:,1]
                    a=xf-xi
                    s[tf] += int(np.sum(a))
    return s   


def makebed_intersect(data):
    #to make a bed out of bedtools intersect average over all cpg ... 
    out=[]
    ind=-1
    while ind<len(data)-1:
        ind+=1
        try:
         line=data[ind]
         num=0
         methyl=0
         tmp=line[1]
        
         while line[1]==tmp and ind<len(data)-2:
            num+=1
            ind+=1
            methyl+=line[-3]
            ln=data[ind]
            tmp=ln[1]
         if num!=0:
          line=line[:3]
          line.append(methyl/np.float(num))
          line.append(num)
          if 'chr' not in str(line[0]):
             line[0]='chr'+str(line[0])
          out.append(line)
        except:
         pass
    
    return out

def bin_bed(bed,size,filter):
    
    nums=[]
    bins=[]
    out={}
    for chrom in bed.keys():
       out[chrom]=[]
       for i in range(0,int(bed[chrom][-1,1])-size,size):
         indi=bisect.bisect_left(bed[chrom][:,0],i)
         indj=bisect.bisect_left(bed[chrom][:,0],i+size)
         num=0
         for j in range(indi,indj):
              num+=bed[chrom][j,1]-bed[chrom][j,0]
         nums.append(num)
         if num>filter:
            out[chrom].append([i,i+size,num])
         #bins.append([chrom,i,i+size,num])
       out[chrom]=np.array(out[chrom])
    nums=np.array(nums)
    sortind=nums.argsort()
    #plot(nums[sortind][::-1])
    
    
    return nums[sortind][::-1],out
       
    





def scaleup_beds(bed,flank,*key):
	
	out={}
	if len(key)>0:
	    tmp={}
	    tmp[key[0]]=bed[key[0]]*1.
	else:
	    tmp=bed
	for chrom in tmp.keys():
		out[chrom]=[]
		if len(tmp[chrom])>0:
			j=0
			a=tmp[chrom][:,0]
			b=tmp[chrom][:,1]
			i1=j
			while j <(len(a))-1:
				xf,xi=a[j+1], b[j]
				i2=j
				if xf-xi<flank:
					contin=True
				else:
					contin=False
				if i1!=i2:
				    sc=np.mean(tmp[chrom][i1:i2,2])
				else:
					sc=tmp[chrom][i1,2]
				j=i2+1
				end=b[i2]
				st=a[i1]
				if not contin:
					out[chrom].append([st,end,sc])
					i1=j
		st=a[i1]
		end=b[j]
		out[chrom].append([st,end,sc])
		if len(out[chrom])>0:
		    out[chrom]=np.array(out[chrom])
	return out

def sort_coverage(bed):
    
    cov=[]
    loci=[]
    for chrom in bed.keys():
       cov+=list(bed[chrom][:,-1])
       ls=[]
       for i in bed[chrom][:,:2]:
            ls.append([i,chrom])
       loci+=list(ls)
    cov=np.array(cov)
    inds=cov.argsort()
    loci=np.array(loci)
    print len(loci),len(cov)
    return cov[inds][::-1],loci[inds][::-1]


def bed_loc(loci):
    
    bed={}
    for i in loci:
       if i[-1] not in bed.keys():
          bed[i[-1]]=[]
       bed[i[-1]].append(i[0])
    for chrom in bed.keys():
       bed[chrom]=np.array(bed[chrom])
    return bed
    
def filter_nuc_bed(bed,col,pval):

    out = {}
    rank=[]
    for chrom in bed.keys():
          if not out.has_key(chrom):
              out[chrom] = bed[chrom][bed[chrom][:,col]<pval]
              rank+=list(bed[chrom][:,col])
    rank=np.array(rank)
    inds=rank.argsort()
    rank=rank[inds][::-1]

    return out, rank


def intersect_genes(bed,gene_pos):
    
    genes=[]
    out_pos={}
    for chrom in gene_pos:
      if chrom in bed.keys():
        for  i in gene_pos[chrom]:
	    xigene=np.int(i[0])
	    xfgene=np.int(i[1])
	    ind=bisect.bisect_left(bed[chrom][:,0],xigene)
	    ind=min(ind,len(bed[chrom])-1)
	    xi=bed[chrom][ind][0]
	    xf=bed[chrom][ind][1]
	    xmin=max(xi,xigene)
	    xmax=min(xf,xfgene)
	    if xmin<xmax:
		if chrom not in out_pos.keys():
		    out_pos[chrom]=[]
		out_pos[chrom].append(i)
		genes.append(i[-1])
	    
	
	
    for chrom in out_pos:
	out_pos[chrom]=np.array(out_pos[chrom])
	

    return genes,out_pos






