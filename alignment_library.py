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
        self.__read_file__(input_file,tag_size,file_format)
        
    def __read_file__(self,input_file,file_format,tag_size,sep='\t'):
        """The read file function itself, it is called in the instantiation
        of the class and should not be called directly"""
        #Creates map to convert - and + to -1 and 1
        convert = {'+':'1','-':'-1'}
        if file_format in ('chen','patman','gff'):
            for zline in input_file:
                ls_parameters = zline.split(sep)
                #Procedure to deal with the default format used by Chen et al 2007
                if file_format=='default':
                    try:self.rna_signatures[ls_parameters[0]+','+ls_parameters[1]+','+ls_parameters[2]].count+=int(ls_parameters[4])
                    except KeyError:self.rna_signatures[ls_parameters[0]+','+ls_parameters[1]+','+ls_parameters[2]]=rna_signature(ls_parameters[0],ls_parameters[1],ls_parameters[2],ls_parameters[3],ls_parameters[4],tag_size)
                elif file_format=='patman':
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
                        
                elif file_format=='gff':
                    if int(ls_parameters[4])-int(ls_parameters[3])+1==tag_size:
                        if ls_parameters[6]=='+':coordinate = ls_parameters[3]
                        elif ls_parameters[6]=='-':coordinate = ls_parameters[4]
                        try:self.rna_signatures[ls_parameters[0]+','+coordinate+','+convert[ls_parameters[6]]].count+=int(ls_parameters(5))
                        except KeyError:self.rna_signatures[ls_parameters[0]+','+coordinate+','+convert[ls_parameters[6]]]=rna_signature(ls_parameters[0],coordinate,convert[ls_parameters[6]],ls_parameters[8].split()[-1],ls_parameters[5],tag_size)
        elif file_format=='bam':
            handle = BAM_Reader(input_file)
            for alignment in handle:
                if alignment.aligned:
                    zchr =alignment.iv.chrom
                    zstart = alignment.iv.start
                    zend = alignment.iv.end
                    zstrand = convert[alignment.iv.strand]
                    ztag = str(alignment.read.seq)
                    if len(ztag)==tag_size:
                        try:self.rna_signatures[zchr+','+str(zstart)+','+zstrand].count+=1
                        except KeyError:self.rna_signatures[zchr+','+str(zstart)+','+zstrand]=(zchr,zstart,zstrand,ztag,1,tag_size)
        handle.close()
        return
            
