'''
Modified on April 29, 2014
Implementation of parser for default, gff and patman alignment files
and the structure to store each sRNA aligment tag
@author: bacms2
'''
from HTSeq import BAM_Reader

class rna_signature():
    def __init__(self,chromosome,coordinate,strand,tag,count,tag_size):
        """Class for store the information about each sRNA alignment tag
        in a alignment file
        chromosome -> Target of the alignment
        coordinate -> Start coordinate of the tag
        strand -> Direct or inverted tag
        sequence -> Nucleotide Sequence
        count -> Number of tags from the sequencing
        """
        self.chromosome = chromosome
        self.coordinate = int(coordinate)
        self.strand = int(strand)
        self.tag = tag
        self.count = int(count)
        self.tag_size = tag_size
    
    def output(self):
        """Prints the information stored in each signature in
        a human readable form"""
        return ('%s\t%d\t%d\t%s\t%d')%(self.chromosome,self.coordinate,self.strand,self.sequence,self.count)
    
class alignment_file():
    def __init__(self,input_file,tag_size,file_format='bam'):
        """Read the alignment files and return a python dictionary
        with all the tags. Each value of the map is an  rna_signature 
        object.
        zfInput -> Path to the input file
        nSize -> Size os sRNAs to search for 
        fFormat -> The format of the alignment file, if not given uses
        the same one described in Ho-Ming Chen Paper PNAS2007"""
        self.rna_signatures = {}
        print 
        if file_format=='bam':
            self.__read_bam__(input_file,tag_size,file_format)
        elif file_format=='sam':
            self.__read_sam__(input_file,tag_size,file_format)
        elif file_format=='patman':
            self.__read_patman__(input_file,tag_size,file_format)
        elif file_format=='gff':
            self.__read_gff__(input_file,tag_size,file_format)
        elif file_format=='chen':
            self.__read_chen__(input_file,tag_size,file_format)
        else: exit('Unrecognized file format')
    
    def __read_bam__(self,input_file,file_format,tag_size):
        """
        Use pysam to read a bam file
        """
        convert = {'+':'1','-':'-1'}
        for alignment in BAM_Reader(input_file):
            if alignment.aligned:
                zchr =alignment.iv.chrom
                zstart = alignment.iv.start
                zstrand = convert[alignment.iv.strand]
                ztag = str(alignment.read.seq)
                if len(ztag)==tag_size:
                    try:self.rna_signatures[zchr+','+str(zstart)+','+zstrand].count+=1
                    except KeyError:self.rna_signatures[zchr+','+str(zstart)+','+zstrand]=(zchr,zstart,zstrand,ztag,1,tag_size)
        return 
    
    def __read_sam__(self,input_file,tag_size,file_format):
        convert = {'+':'1','-':'-1'}
        for alignment in BAM_Reader(input_file):
            if alignment.aligned:
                zchr =alignment.iv.chrom
                zstart = alignment.iv.start
                zstrand = convert[alignment.iv.strand]
                ztag = str(alignment.read.seq)
                if len(ztag)==tag_size:
                    try:self.rna_signatures[zchr+','+str(zstart)+','+zstrand].count+=1
                    except KeyError:self.rna_signatures[zchr+','+str(zstart)+','+zstrand]=(zchr,zstart,zstrand,ztag,1,tag_size)
        return
    
    def __read_patman__(self,input_file,tag_size,file_format):
        convert = {'+':'1','-':'-1'}
        handle = open(input_file)
        for zline in handle:
            zline = zline.rstrip('\n')
            ls_parameters = zline.split('\t')
            if (int(ls_parameters[3])-int(ls_parameters[2])+1)==tag_size:#Check if the tag has the correct size
                #Test the start position of the tag
                if ls_parameters[4]=='+':coordinate = ls_parameters[2]
                elif ls_parameters[4]=='-':coordinate = ls_parameters[3]
                #Check whether the file is non redundant
                try:count = int(ls_parameters[1].split()[-1].split(':')[-1])
                except ValueError:count = 1
                sequence = ls_parameters[1].split()[0]
                try:
                    #If tag is redundant check if it is already in the Hash map and update its number 
                    self.rna_signatures[ls_parameters[0]+','+coordinate+','+convert[ls_parameters[4]]].count+=int(count)
                except KeyError:
                    #If not present add tag to the map
                    self.rna_signatures[ls_parameters[0]+','+coordinate+','+convert[ls_parameters[4]]]=rna_signature(ls_parameters[0],coordinate,convert[ls_parameters[4]],sequence,int(count),tag_size)
        handle.close()
        return
    
    def __read_gff__(self,input_file,tag_size,file_format):
        convert = {'+':'1','-':'-1'}
        handle = open(input_file)
        for zline in handle:
            zline = zline.rstrip('\n')
            ls_parameters = zline.split('\t')
            if int(ls_parameters[4])-int(ls_parameters[3])+1==tag_size:
                if ls_parameters[6]=='+':coordinate = ls_parameters[3]
                elif ls_parameters[6]=='-':coordinate = ls_parameters[4]
                try:self.rna_signatures[ls_parameters[0]+','+coordinate+','+convert[ls_parameters[6]]].count+=int(ls_parameters(5))
                except KeyError:self.rna_signatures[ls_parameters[0]+','+coordinate+','+convert[ls_parameters[6]]]=rna_signature(ls_parameters[0],coordinate,convert[ls_parameters[6]],ls_parameters[8].split()[-1],ls_parameters[5],tag_size)
        handle.close()
        return
    
    def __read_chen__(self,input_file,tag_size,file_format):
        handle = open(input_file)
        for zline in handle:
            zline = zline.rstrip('\n')
            ls_parameters = zline.split('\t')
            if len(ls_parameters[3])==tag_size:
                try:self.rna_signatures[ls_parameters[0]+','+ls_parameters[1]+','+ls_parameters[2]].count+=1
                except KeyError:self.rna_signatures[ls_parameters[0]+','+ls_parameters[1]+','+ls_parameters[2]]=rna_signature(ls_parameters[0],ls_parameters[1],ls_parameters[2],ls_parameters[3],1,tag_size)
        handle.close()
        return
    

            
