'''
Modified on April 29, 2014
Implementation of parser for default, gff and patman alignment files
and the structure to store each sRNA aligment read
@author: bacms2
'''
from sys import exit

class rna_signature():
    def __init__(self,chromosome,coordinate,strand,read,count,read_size):
        """Class for store the information about each sRNA alignment read
        in a alignment file
        chromosome -> Target of the alignment
        coordinate -> Start coordinate of the read
        strand -> Direct or inverted read
        sequence -> Nucleotide Sequence
        count -> Number of tags from the sequencing
        """
        self.chromosome = chromosome
        self.coordinate = int(coordinate)
        self.strand = int(strand)
        self.read = read
        self.count = int(count)
        self.read_size = read_size
    
    def output(self):
        """Prints the information stored in each signature in
        a human readable form"""
        return ('%s\t%d\t%d\t%s\t%d')%(self.chromosome,self.coordinate,self.strand,self.sequence,self.count)
    
class alignment_file():
    def __init__(self,input_file,read_size,file_format='default'):
        """Read the alignment files and return a python dictionary
        with all the reads. Each value of the map is an  rna_signature 
        object.
        zfInput -> Path to the input file
        nSize -> Size os sRNAs to search for
        fFormat -> The format of the alignment file, if not given uses
        the same one described in Ho-Ming Chen Paper PNAS2007"""
        self.rna_signatures = {}
        self.read_file(input_file,read_size,file_format)
        
    def read_file(self,input_file,file_format,read_size,sep='\t'):
        """The read file function itself, it is called in the instantiation
        of the class and should not be called directly"""
        handle = open(input_file)
        if file_format not in ('default','gff','patman'):exit('###ERROR###\nFile format not recognised please refer to the documentation')
        try:int(read_size)
        except ValueError:exit('###ERROR###\nRead size needs to be a positive integer')
        #Creates map to convert - and + to -1 and 1
        convert = {'+':'1','-':'-1'}
        for zline in handle:
            ls_parameters = zline.split(sep)
            #Procedure to deal with the default format used by Chen et al 2007
            if file_format=='default':
                try:self.rna_signatures[ls_parameters[0]+','+ls_parameters[1]+','+ls_parameters[2]].count+=int(ls_parameters[4])
                except KeyError:self.rna_signatures[ls_parameters[0]+','+ls_parameters[1]+','+ls_parameters[2]]=rna_signature(ls_parameters[0],ls_parameters[1],ls_parameters[2],ls_parameters[3],ls_parameters[4],read_size)
            elif file_format=='patman':
                if (int(ls_parameters[3])-int(ls_parameters[2])+1)==read_size:#Check if the read has the correct size
                    #Test the start position of the read
                    if ls_parameters[4]=='+':coordinate = ls_parameters[2]
                    elif ls_parameters[4]=='-':coordinate = ls_parameters[3]
                    #Check whether the file is non redundant
                    try:count = int(ls_parameters[1].split()[-1].split(':')[-1])
                    except ValueError:count = 1
                    sequence = ls_parameters[1].split()[0]
                    try:
                        #If read is redundant check if it is already in the Hash map and update its number 
                        self.rna_signatures[ls_parameters[0]+','+coordinate+','+convert[ls_parameters[4]]].count+=int(count)
                    except KeyError:
                        #If not present add read to the map
                        self.rna_signatures[ls_parameters[0]+','+coordinate+','+convert[ls_parameters[4]]]=rna_signature(ls_parameters[0],coordinate,convert[ls_parameters[4]],sequence,int(count),read_size)
                    
            elif file_format=='gff':
                if int(ls_parameters[4])-int(ls_parameters[3])+1==read_size:
                    if ls_parameters[6]=='+':coordinate = ls_parameters[3]
                    elif ls_parameters[6]=='-':coordinate = ls_parameters[4]
                    try:self.rna_signatures[ls_parameters[0]+','+coordinate+','+convert[ls_parameters[6]]].count+=int(ls_parameters(5))
                    except KeyError:self.rna_signatures[ls_parameters[0]+','+coordinate+','+convert[ls_parameters[6]]]=rna_signature(ls_parameters[0],coordinate,convert[ls_parameters[6]],ls_parameters[8].split()[-1],ls_parameters[5],read_size)
        handle.close()
        return
            
