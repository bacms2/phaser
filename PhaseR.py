#!/usr/bin/env python
# coding=utf-8
"""
beta v1.1 29 April 2014
PhaseR: predicting small RNA components of regulatory networks
@author: bacms2

PhaseR: Predicting Small RNA Components of Regulatory Networks.
Currently, PhaseR is a set of Python scripts.

PhaseR can identify phased small RNA loci using high throughput genomic
sequencing datasets. Its input is the alignment of a small RNA dataset to a
reference sequence. The alignment can be done using one of many freely available
alignment software programs e.g.  PatMaN, bowtie. The alignment file can be 
passed to the PhaseR algorithm to identify potentially phased small RNA loci.
The PhaseR algorithm uses an extension of the probability calculation proposed
 by Chen at al. (2007) to distinguish likely occurrences of phasing from random
events. In Chen’s original algorithm each position along an sRNA locus was 
treated as a binary variable. In other words it could only have two states,
either it was occupied by the 5’end of an sRNA or it was not. The new algorithm 
innovates by using the counts of the number of sRNAs with a 5’ end in each 
position when calculating the probability for the locus. The algorithm also 
considers every possible sRNA in the dataset to be the start of the sRNA locus
for which the end is unknown. So every possible length for that locus is tested 
as long as it is not bigger than a certain number of nucleotides for which there
are no matching sRNAs.

A matrix of probabilities is then built with one dimension corresponding to all
possible phased segments ("locus") and the other to all possible registers. 
Each position in the matrix is filled with the minimum probability from the 
hypergeometric tests performed for the different abundance thresholds. The last
step consists of determining the segment, register, abundance and probability 
for the element in the matrix with the lowest probability. These probabilities
are calculated using the hypergeometric distribution as proposed by Chen at al.
(2007).
Phaser is implemented in the Python programming language and it is released 
under the GPL v3, See http://www.gnu.org/copyleft/gpl.html for more details. 
Please contact Bruno Santos on bacms2@cam.ac.uk for further support.

Version 1.1

"""
from matplotlib import use
use('Agg')

from alignment_library import alignment_file
from numpy import zeros,hstack,arange,log,unique
from time import time,ctime
from sys import argv,exit
from os import path
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch
from collections import defaultdict
import pylab
from math import log as mathlog
try:from mpmath import binomial,hyper
except ImportError:pass
try: import psyco
except ImportError:pass
try: psyco.full()
except NameError: pass
try:
    import rpy2.rinterface as rinterface
    bolean_rpy2 = True
except ImportError:bolean_rpy2 = False

if bolean_rpy2: #Initialize the R interface if present
    rinterface.initr()
    phyper = rinterface.globalenv.get("phyper")
    SexpVector = rinterface.SexpVector
    myparams = {'lower.tail':SexpVector([False,],rinterface.LGLSXP),'log.p': SexpVector([True,],rinterface.LGLSXP)}
    dtPhyper = {}

class alignment_node:
    def __init__(self,name,counts,reads,coordinate):
        """
        self.__init__(self,zName,nCounts,zTag,coordinate)
        Creates a new instance of srna_alignment.
        A srna_alignment correspond to a position on the genome with a 5' end of a sRNA matching to it
        zName <- Unique name to identify the smallRNA, its normally obtained by joining the Chromosome
        name with its the coordinate on the genome (ie chr1_11210)
        nCounts <- Tuple of integers corresponding to the count numbers of both strands
        zTags <- Tuple of strings corresponding to the tag sequences of both strands
        coordinate <- Coordinate of the srna_alignment
        """
        self.name = name 
        self.ls_counts = counts
        self.ls_reads = reads
        self.coordinate = int(coordinate)
    
    def get_parameters(self):
        """
        self.getPar(self)
        Return the name and coordinate of the signature used to create the srna_alignment
        """
        lsAux = self.zname.split(',')
        return (lsAux[0],self.coordinate)
    
    def __str__(self):
        """
        self.__str__(self)
        Changes the default representation of the object to a human readable form.
        """
        return ('Name:%s Coord:%d Counts:%d|%d Tags:%s|%s\n' % (self.zname,self.coordinate,self.count[0],self.counts[1],self.lssequences[0],self.lssequences[1]))
        
class alignment_graph():
    def __init__(self,dt_reads):
        """
        self.__init__(self,rna_signatures)
        Creates a new instance of Graph
        A Graph corresponds to a graph representation of all the 
        rna_signatures <- Dictionary containing all the sRNA signatures in the library
        """
        self.dt_alignment_nodes = {} #store all the nodes in the graph
        self.dt_rna_signatures = dt_reads #Store all sRNA matches 
        self.unique_alignments = 0 #Size of the graph
    
    def add_alignment_node(self,zname):
        """
        self.add_alignment_node(self,zName)
        add a Node to the graph
        zName <- Unique name to identify the smallRNA, its normally obtained by joining the Chromosome
        name with its the coordinate on the genome
        eg. chr1_11210 
        """
        try:
            #If node already exists return
            self.dt_alignment_nodes[zname]
            return
        except KeyError:pass
        #Get the parameters of the node
        coordinate = zname.split(',')[1]
        #If the coordinates are outside the chromosome return
        if coordinate <= 0:return
        
        #First get the sense signature
        counts_forward,tag_forward = (0.0,0.0)
        try:
            counts_forward = self.dt_rna_signatures[zname+',1'].count
            tag_forward = self.dt_rna_signatures[zname+',1'].read
        except KeyError:pass
        #Repeat the same for the opposite strand
        counts_reverse,tag_reverse = (0.0,0.0)
        try:
            counts_reverse = self.dt_rna_signatures[zname+',-1'].count
            tag_reverse = self.dt_rna_signatures[zname+',-1'].read
        except KeyError:pass
        self.dt_alignment_nodes[zname]=alignment_node(zname,(counts_forward,counts_reverse),(tag_forward,tag_reverse),coordinate)
        #Increase the number of srna_alignments in the graph
        self.unique_alignments+=1


    def cardinal(self):
        """
        self.card(self)
        return the number of nodes in the graph
        """
        return self.unique_alignments
       
class locus:
    def __init__(self,start_coord,end_coord,locus_length,read_size,chromosome = 'None'):
        """
        locus(self,start_coord,end_coord,locus_length,read_size)
        Create a new instance of the locus class
        start_coord <- The start coordinate of the locus on the genome
        end_coord <- The end coordinate of the locus on the genome
        locus_length <- The length of the locus
        read_size <- The size class of sRNAs being tested
        chromosome <- Optional argument to define the chromosome name
        Note that the length of the locus will always be higher than the distance between
        their start and end coordinates because two rows of zeros are added to the rigth
        and left extremes of the locus for computation purposes
        """      
        self.start_coord = start_coord
        self.end_coord = end_coord
        self.chromosome = chromosome
        self.forward_reads = locus_length*[0]
        self.forward_counts = zeros(locus_length,dtype='Float64')
        self.reverse_reads = locus_length*[0]
        self.reverse_counts = zeros(locus_length,dtype='Float64')
        self.locus_length = locus_length
        self.read_size = read_size
        self.p_value = 0.0
        self.phased_region_indexes = (0,0)
        self.phased_counts = 0
        
    def mask_repeats(self):
        """
        Function to mask the tags that appear more than once on segment.
        """
        #get repeated position on forward strand
        items = defaultdict(list)
        for i,item in enumerate(self.forward_reads):items[item].append(i)
        indexes = []
        for item, locus in items.iteritems():
            if len(locus) > 1 and item!=0:indexes.extend(locus)
        self.forward_counts[indexes]=0
        
        #get repeated position on anti sense strand  
        for key in items.iterkeys():items[key] = []
        for i,item in enumerate(self.reverse_reads):items[item].append(i)
        indexes =[]
        for item, locus in items.iteritems():
            if len(locus) > 1 and item!=0:indexes.extend(locus)
        self.reverse_counts[indexes]=0
        return
        
    def phaser(self,left_index,right_index):
        """
        self.hypergeometric(self,left_index,rigth_index)
        Performs the hypergeometric test for the locus between [left_index,rigth_index[ for
        all the count number in phase positions and return the minimum p-Value
        left_index<-Integer indicating the left position on the array
        rigth_index <-Integer indicating the right position on the array
        """
        #Create the subarray to test
        reverse_counts = self.reverse_counts[left_index-2:right_index-2]
        #Test to make sure no mistake is happening with the indexes
        if left_index<0 or right_index<0: raise Exception('Indexes cannot be lower than zero')
        if len(reverse_counts)!= len(self.forward_counts[left_index:right_index]):raise Exception('Sense and Anti sense have different lengths')
        #Stack both strands
        locus_strand = hstack([self.forward_counts[left_index:right_index],reverse_counts[::-1]])
        
        #Get the length of the new locus_strand
        locus_length = len(locus_strand)
        
        #Get the indexes of Phased positions
        phased_index=arange(0,locus_length,self.read_size)
        #Get and array with only the Phased Positions
        phased_strand =locus_strand[phased_index]
        #Get the unique values in the array of Phased counts
        phased_counts = unique(phased_strand[phased_strand>0])

        #Caculate the number of phase and non phase positions in the locus
        phased_positions_total = locus_length/self.read_size
        non_phased_positions_total = (locus_length-1)-(locus_length/self.read_size-1) 
        
        #Do a standard call of lphyper to increase speed
        if bolean_rpy2:r_phyper = self.r_phyper
        else: mhyper = self.mhyper
        
        #Remove the zero values from locus_strand to reduce its size
        locus_strand = locus_strand[locus_strand>0]
        locus_total_counts = locus_strand.sum()
        
        #Initialize the p-Value to zero
        p_value = 0
        phased_positions_occupied,total_positions_occupied=0,0
        for count in phased_counts:
            #Calculate matches in Phased positions
            phased_positions_occupied = (phased_strand>=count).sum()
            #Calculate matches in non Phased positions and reduce the array to values higher than r
            locus_strand = locus_strand[locus_strand>=count]
            total_positions_occupied = len(locus_strand)
            #Get the p-value from the hypergeometric distribution
            if bolean_rpy2:temp_p_value = r_phyper(phased_positions_occupied-1, phased_positions_total, non_phased_positions_total,total_positions_occupied)
            else: temp_p_value = mhyper(phased_positions_occupied-1, phased_positions_total, non_phased_positions_total,total_positions_occupied) 
            p_value = min(p_value,temp_p_value)

        #Return p-Value
        return p_value,locus_total_counts
    
    
    def mhyper(self,q,m,n,k):
        u = min(self.read_size,n)-1
        A = binomial(m,q)*binomial(n,k-q)
        B = hyper([1,q-k,q-m],[1+q,1-k+n+q], 1)
        C = -binomial(m,1+u)*binomial(n,k-1-u)
        D = hyper([1,1-k+u,1-m+u],[2+u,2-k+n+u], 1)
        return mathlog(float((A*B + C*D) / binomial(m+n, k)))
    
    def r_phyper(self,q,m,n,k):
        """
        self.phyper(self,q,m,n,k)
        Calculate p-value using R function phyper from rpy2 low-level
        interface. 
        "R Documentation
        phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
        q: vector of quantiles representing the number of white balls
            drawn without replacement from an urn which contains both
            black and white balls.
        m: the number of white balls in the urn.
        n: the number of black balls in the urn.
        k: the number of balls drawn from the urn.
        log.p: logical; if TRUE, probabilities p are given as log(p).
        lower.tail: logical; if TRUE (default), probabilities are P[X <= x],
            otherwise, P[X > x].
        "
        """
        phyper_q = SexpVector([q,], rinterface.INTSXP)
        phyper_m = SexpVector([m,], rinterface.INTSXP)
        phyper_n = SexpVector([n,], rinterface.INTSXP)
        phyper_k = SexpVector([k,], rinterface.INTSXP)
        return phyper(phyper_q,phyper_m,phyper_n,phyper_k,**myparams)[0]
        
    def plot(self,screen_display=False,whole_locus=True):
        """
        self.plot(self,display=False,whole=True)
        whole_locus <- Flag can be used to plot the whole locus or just the most significant part
        Plot the number of signatures in the sense and anti sense 
        strand for the current locus.
        """
        forward = True
        
        #Create a figure object
        fig = pylab.figure(figsize=(16,9),dpi=160)
        
        
        for locus_strand in [self.forward_counts,self.reverse_counts]:        
            
            #Define Patch objects for the legend
            non_phased_positions = Patch(edgecolor='b', facecolor='b')        
            phased_positions = Patch(edgecolor='r', facecolor='r')
            
            
            #Plot strand
            if whole_locus:#Plot the whole strand
                #Define the title to be shown on the plot
                pylab.suptitle('%s:%d..%d with log p-value = %.2f' % (self.chromosome,self.start_coord,self.end_coord,self.p_value),fontsize=18,x=0.39)
                #Remove the extra elements from the sense strand
                locus_strand = locus_strand[self.read_size:-self.read_size]
                #get the first position in phase
                first_position = self.phased_region_indexes[0]-self.phased_region_indexes[0]/self.read_size*self.read_size
                if not forward: first_position=first_position-3
                #Get the positions that should contain phased sRNAs
                phased_positions_index = arange(first_position,len(locus_strand),self.read_size)
                
            else: #Plot only the most significant region
                #Define the title to be shown on the plot
                pylab.suptitle('Most significant region on phased sRNA loci in %s:%d..%d with log p-value = %.2f' % (self.chromosome,self.start_coord+self.phased_region_indexes[0]-self.read_size,self.start_coord+self.phased_region_indexes[1]-self.read_size,self.p_value),fontsize=12,x=0.45)
                #Remove the extra elements from the sense strand
                locus_strand = locus_strand[self.phased_region_indexes[0]:self.phased_region_indexes[1]]
                #Get the positions that should contain phased sRNAs
                if forward: phased_positions_index = arange(0,len(locus_strand),self.read_size)
                else: phased_positions_index = arange(self.read_size-3,len(locus_strand),self.read_size)
            
            #Plot title and legend for the figure
            pylab.figlegend([phased_positions,non_phased_positions],['Phased positions','Out of phase positions'],'upper right',prop={'size':8})
            #Produce the actual plot
            ##Divide the screen in two
            if forward:pylab.subplot(211)
            else:pylab.subplot(212)
            #Adjust the margins
            pylab.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.92, hspace=0.16)
            #Display the labels on the axis 
            pylab.xlabel('Relative position',fontsize= 14)
            ylabel = 'Number of sRNAs'
                        
            #If the values are large plot in log scale
            try:
                phased_positions_occupied = phased_positions_index[locus_strand[phased_positions_index]>0]
                if max(locus_strand)>=1000: 
                    ylabel = 'log(Number of sRNAs)'
                    locus_strand = log(locus_strand)
            except ValueError: pass #May not contain any significant phased region
            
            #Display the ylabel
            pylab.ylabel(ylabel,fontsize= 14)
            
            #Plot the bars corrresponding to whole the locus
            pylab.bar(range(len(locus_strand)),locus_strand,color='b',width=1.0,edgecolor='b')
            #Plot the phased positions bars
            pylab.bar(phased_positions_index,locus_strand[phased_positions_index],color='r',width=2,edgecolor='r',alpha=1,align='center')
            #Compute the position where to place the markers of phased positions
            try:yaxis_pos = max(locus_strand)/3.0*2 
            except ValueError:yaxis_pos = 0
            #Plot the expected phased positions
            try:pylab.scatter(phased_positions_index, [yaxis_pos]*len(phased_positions_index), marker='d', color='white',edgecolor='black')
            except ValueError:pass
            #Plot the occupied phased positions
            pylab.scatter(phased_positions_occupied, [yaxis_pos]*len(phased_positions_occupied),s=40,marker='d',color='r',edgecolor='black')
            #Try to set the axis 
            try:pylab.axis([0,len(locus_strand),0,max(locus_strand)])
            except ValueError:pass
            
            forward=False
            
            
        if screen_display:
            pylab.show()
            pylab.close()
            return
        return fig

class alignment_file_manager:
    def __init__(self,rna_signatures,read_size):
        """
        libraryManager(self,rna_signatures,read_size)
        rna_signatures <- Dictionary containing signature objects
        read_size <- Size class of sRNAs
        Create a new instance of the libraryManager
        """
        self.graph = alignment_graph(rna_signatures) #Instantiate the graph
        self.scc = [] #Will be used to store the Strong Connected Components
        self.read_size = read_size
        self.alignment_results = rna_signatures#Dictionary containing the srna alignments 
        self.dt_locus = {} #store the results to be used in iPython 
        self.p_values_list = [] #store the p-values to be used for debug
        
    def create_graph(self,maximum_gap):
        """
        self.createNodes(self,maximum_gap)
        maximum_gap <- Maximum gap allowed between two signatures
        Creates all the nodes based on the signatures on the rna_signatures.
        A node is defined by its chromosomal location and contains information
        about both strands. For the extension step only the presence or absence
        of the node is taken in consideration.
        """
        #Create a temporary dictionary to store the nodes per chromosome 
        dt_nodes = {}
        
        #Create the nodes for all the signatures
        for signature in self.alignment_results.values():
            signature_name = '%s,%d' % (signature.chromosome,signature.coordinate)
            #Add the node to the graph
            self.graph.add_alignment_node(signature_name)
            #Add the node to the dictionary of its chromosome
            try: dt_nodes[signature.chromosome].append(signature.coordinate)
            except KeyError:dt_nodes[signature.chromosome]=[signature.coordinate]
        
        #Get all the loci
        for chromosome in dt_nodes:
            scc = []#Define a list to store a locus
            #Sort nodes by their chromosome and then by their coordinate
            dt_nodes[chromosome] = dict.fromkeys(dt_nodes[chromosome]).keys()
            dt_nodes[chromosome].sort()
            
            #Create the loci
            for i in xrange(0,len(dt_nodes[chromosome])-1):
                signature_name = '%s,%d' % (chromosome,dt_nodes[chromosome][i])
                gap_between_signatures = dt_nodes[chromosome][i+1]-int(dt_nodes[chromosome][i])+1
                scc.append(signature_name)
                if gap_between_signatures>maximum_gap:#If the gap between two nodes is bigger than the maximum_gap store the current locus and start a new one
                    if len(scc)>1:self.scc.append(scc[:])
                    scc=[]
            
            #Append the last scc since there is no node after to test the distance
            signature_name = '%s,%d' % (chromosome,dt_nodes[chromosome][len(dt_nodes[chromosome])-1])
            scc.append(signature_name)
            if len(scc)>1:self.scc.append(scc[:])
        del dt_nodes#Call garbage collection for dt_nodes to keep a low memory footprint
        return
    
    def build_locus(self,scc):
        """
        self.buildLocus(self,scc)
        scc <- List containing the name of the nodes to create a Locus
        Create and return a Locus object, based on the list of Nodes.
        """
        #Calculate the length of the region with sRNAs span on it
        locus_length = int(scc[-1].split(',')[1])-int(scc[0].split(',')[1])+1
        
        #Create the locus object 
        temp_locus = locus(int(scc[0].split(',')[1]),int(scc[-1].split(',')[1]),locus_length+2*self.read_size,self.read_size)
        temp_locus.chromosome,coordinate = scc[0].split(',')
        coordinate = int(coordinate)
        for node_name in scc:
            node = self.graph.dt_alignment_nodes[node_name]
            temp_locus.forward_counts[node.coordinate-coordinate+self.read_size],temp_locus.reverse_counts[node.coordinate-coordinate+self.read_size] = node.ls_counts
            temp_locus.forward_reads[node.coordinate-coordinate+self.read_size],temp_locus.reverse_reads[node.coordinate-coordinate+self.read_size] = node.ls_reads
        temp_locus.mask_repeats()
        return temp_locus
        
        
    def compute_significance(self,fThreshold,locus_minimum_length,output_file=None,write_pdf=False):
        """
        M.computeSignificance(self,fThreshold,nMilocus_length,zfOutput,bWrite=False)
        Function to compute hypergeometric tests for all the possible paths 
        connecting the nodes.
        """
        #If an output location is specified create the files for them
        if output_file:handle = open(output_file+'.txt','w')
        if write_pdf:pdf = PdfPages(output_file+'.pdf')
        
        #Initialize the values to write for the locus
        locus_values = ['0.0','None','None','None','0']
        
        if len(self.alignment_results)==0:return locus_values
        
        for scc in self.scc:
            #If the segment has less than 2 nodes break
            if len(scc)<=2:continue 
            
            #Build the locus object by calling the buildLocus method
            loc = self.build_locus(scc)
            
            #Control test to ensure that the scc does not has more sRNAs than 
            if len(scc)>loc.locus_length:exit('Error length of scc cannot be higher than loc.length')
            
            #If locus is smaller than the minimum length ignore it, mainly for speed purposes 
            if loc.locus_length<locus_minimum_length:continue
            
            #Store the locus for posterior analysis under an ipython session
            self.dt_locus[scc[0]]=loc
            
            #Initiate minimum p-value at 0.0 and phased index at -1,-1
            phased_p_value = 0.0
            phased_index = (-1,-1)
            
            #Creates a boolean variable to flag if the locus is significant or not
            locus_is_phased = False 
            #Start iterating over each position on the locus as a potential start site for phasing
            for left_index in xrange(2,loc.locus_length-self.read_size+2):
                #If the first node has no reads on phased positions ignore this walk
                if loc.forward_counts[left_index]==0 and left_index+self.read_size-3<loc.locus_length and loc.reverse_counts[left_index+self.read_size-3]==0:continue
                #Start iterate over potential end positions for phasing
                for right_index in xrange(left_index+self.read_size,loc.locus_length,self.read_size):
                    #if the locus length is smaller than the minimum size for the locus ignore it
                    if right_index-left_index+1<locus_minimum_length:continue
                    #If the new locus as no sRNAs in the end positions ignore this walk 
                    if loc.forward_counts[right_index-self.read_size]==0 and loc.reverse_counts[right_index-3]==0:continue
                    
                    #Get the p_value for this segment of the locus
                    p_value = loc.phaser(left_index,right_index)[0]
                    if p_value<phased_p_value:
                        phased_p_value=p_value
                        phased_index = (left_index,right_index)
            loc.phased_region_indexes = phased_index
            loc.p_value = phased_p_value
            
            #If the p-value is higher than the threshold save it on the list
            if loc.p_value<=float(locus_values[0]):
                
                coordinate= scc[0].split(',')[1]    
                counts = loc.phased_counts
                locus_values = [str(loc.p_value),loc.chromosome,str(loc.phased_region_indexes[0]+int(coordinate)-int(self.read_size)),str(loc.phased_region_indexes[1]+int(coordinate)-int(self.read_size)),counts]
            elif loc.p_value>=float(locus_values[0]):locus_values=[str(loc.p_value),str(loc.chromosome),str(loc.start_coord),str(loc.end_coord)]
            
            if loc.p_value<=fThreshold:locus_is_phased=True
            
            if locus_is_phased and output_file:
                locus_is_phased = False
                self.writeOutput(handle,scc,loc)
                if write_pdf:
                    #Plot the all locus
                    fig = loc.plot()
                    pdf.savefig(fig)
                    pylab.close(fig)
                    del(fig)
                    #Plot the phased bit
                    fig = loc.plot(whole_locus=False)
                    pdf.savefig(fig)
                    pylab.close(fig)
                    del(fig)
        if output_file:handle.close()
        if write_pdf:pdf.close()
        return locus_values
                
    def writeOutput(self,file_handle,scc,loc):
        """ 
        writeOutput(self,handle,scc,loc)
        Write a matrix with the p-values for the locus and a file 
        with a summary of the results
        """
        chromosome,coordinate= scc[0].split(',')
        file_handle.write('Identified sRNA locus in %s from position %s to %i\n' % (chromosome,loc.start_coord,loc.end_coord))
        file_handle.write('\tPhased detected from position %s to %i\n' % (loc.phased_region_indexes[0]+int(coordinate)-int(self.read_size),loc.phased_region_indexes[1]+int(coordinate)-int(self.read_size)))
        file_handle.write('\tLog p-Value = %f\n' % loc.p_value)
        file_handle.write('\tTo plot this locus type: manager.dtlocus[\'%s\'].plot()\n' % scc[0])
        file_handle.write('\tThe coordinates in the graph will be:%i to %i\n' % (loc.phased_region_indexes[0]-self.read_size,loc.phased_region_indexes[1]-int(self.read_size)))
        file_handle.write('%s\n' % ('#'*74))
        return
        
        
if __name__=='__main__':
    print ctime()
    print 'Running Phased RNA identification method (PhaseR) '
    start = time()
    #check whther the user input enough parameters 
    if len(argv)<7:exit('Not enough parameters for the run\nPlease refer to the documentation for more information')
    
    #Input file
    if path.isfile(argv[1]):input_file_path = argv[1]
    else: exit('Input file not found %s' % argv[1])
    
    #Format of the alignment file
    if argv[2] in ('default','patman','gff'):file_format = argv[2]
    else: exit('File format not recognized. File needs to be \'default\', \'gff\' or \'patman\'')
    
    #P-value to be used
    try:p_value_treshold = float(argv[3])
    except ValueError: exit("###ERROR### log p-value needs to a negative float")
    if p_value_treshold>0:exit("###ERROR### log p-value needs to a negative float")
    
    #The individual size of each siRNA
    try:srna_size = int(argv[4])
    except ValueError: exit("###ERROR### sRNA size needs to be a positive integer")
    if srna_size<0:exit("###ERROR### sRNA size needs to be a positive integer")
    
    #Minimum length of locus
    try:locus_minimum_length = int(argv[5])
    except ValueError: exit("###ERROR### Minimum length of the locus needs to be apositive integer")
    if locus_minimum_length<0: exit("###ERROR### Minimum length of the locus needs to be a positive integer")
 
    #Maximum gap 
    try:maximum_gap = int(argv[6])
    except ValueError: exit("###ERROR### Maximum srna gap needs to be a positive integer")
    if maximum_gap<0: exit("###ERROR### Maximum srna gap needs to be a positive integer")
    
    
    #Getting the output information
    #Reading the input files
    file_manager = alignment_file(input_file_path,file_format,srna_size)
    print 'Number of rna signatures',len(file_manager.rna_signatures)
    print 'Reading the input file took:%f' % (time()-start)
    sample_name = input_file_path.split('/')[-1]
    
    #Open file to store the phasing loci
    output_file_path = '../results/%s_%d_%s_phasing' % (sample_name,srna_size,p_value_treshold)
    if path.isfile(output_file_path+'.txt'):print '###Warning###:Output file exists and will be overwritten'
    if path.isfile(output_file_path+'.pdf'):print '###Warning###:pdf file exists and will be overwritten'
    
    #Print all the information to the user
    print
    print 'Parameters to be used will be:'
    print 'Input File->%s' % input_file_path
    print 'Output File -> %s' % output_file_path
    print 'Length of Phased siRNAs:%d' % srna_size
    print 'P-value:%f'% p_value_treshold
    print 'Number of distinct sRNAs used:%d' % len(file_manager.rna_signatures)
    print 'Minimum length for locus: %d' % locus_minimum_length
    print 'Maximum gap length between sRNAs: %d' % maximum_gap
    print 
    #Creating the nodes
    manager = alignment_file_manager(file_manager.rna_signatures,srna_size)
    astart = time()
    manager.create_graph(maximum_gap)
    print 'Nodes Created\nCreating nodes took:%fs' % (time()-astart)
    print 
    astart = time()
    manager.compute_significance(p_value_treshold,locus_minimum_length,output_file_path,write_pdf=True)
    #manager.computeSignificance_fixed_lenght(nThreshold,locus_minimum_length,zfOutput,bPdf=True,method='bphase')
    print 'Computing significance took:%fs' % (time()-astart)
    print 'All finished'
    print 'This run took:%f' % (time()-start)
    
    