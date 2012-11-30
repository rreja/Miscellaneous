import sys, os, operator
from optparse import OptionParser , IndentedHelpFormatter
import numpy as np


def process_file(infile,options):
    dict_fimo = {}
    results = []
    variants = []
    outfile = os.path.join(os.path.dirname(infile),"motif_variation_information.txt")
    out = open(outfile,"w")
    out.write("FIMO_motif_chr\tFIMO_motif_start\tFIMO_motif_strand\tMotif_from_fimo\tpvalue\tMotif_from_BAM\tCIGAR\tQuality\n")
    input = open(infile,'rt')
    for line in input:
        if line.startswith("#"):
            continue
        cols = line.rstrip().split("\t")
        fimo_motif,pval = extract_fimo_motif(cols[8])
        # Appending the chr, start,strand of fimo motif to be used later to find out which motifs were pulled up.
        motif_extra_info = cols[0]+"\t"+cols[3]+"\t"+cols[6]
        motif_start = int(cols[7].split("-")[0])
        motif_end = int(cols[7].split("-")[1])
        motif_length = motif_end - motif_start  #int(cols[4]) - int(cols[3]) # Storing the actual motif co-ordinates in the 8th col.
        if motif_length < 2:
            print "You are using motif mid-point as input. Your motif length should be greater than 1."
            sys.exit(1)
        if int(cols[12]) == -1 or abs(int(cols[18])) > 30: # If the tag does not lie near ref point or the d
            continue
        else:
            # dict key is composed of tag coordinates that matches the motif vicinity. value is composed of fimo_motif start, length,(with if the motif is ++,+-,-+,--),ref_motif,pval and extra info
            dict_fimo[cols[9]+"\t"+cols[12]+"\t"+cols[15]] = str(motif_start)+"\t"+str(motif_length)+"\t"+cols[6]+cols[15]+"\t"+fimo_motif+"\t"+pval+"\t"+motif_extra_info 
            
    results = retrieve_seq_from_SAM(dict_fimo,options,out,outfile)
    variants = find_variation(results)
    for t in variants:
        out.write(t)
    print "Sucessfully written the output. Your output is in"+os.path.dirname(outfile)

def find_variation(list1):
    list2 = []
    tmp = []
    for i in list1:
        tmp = i.split("\t")
        variant_info = report_mismatch(tmp[3],tmp[5])
        list2.append(i+"\t"+variant_info)
    return(list2)
        
            
def report_mismatch(refMotif,newMotif):
    counter = 1
    for i, j in zip(refMotif, newMotif):
        if i == j:
            print i,j
            continue
        else:
            print i,j,counter
        counter = counter +1
    # All I ned to do here is return the position and the allele and then delete sys.exit from below.
    sys.exit(1)
    
def extract_fimo_motif(line):
    cols = line.rstrip().split(";")
    pval = ""
    motif = ""
    for attr in cols:
        tmp = attr.split("=")
        if tmp[0] == "pvalue":
            pval = tmp[1]
        if tmp[0] == "sequence":
            motif = tmp[1]
    return(motif,pval)
            
def retrieve_seq_from_SAM(dict_fimo,options,out,outfile):
    input = open(options.SAMfile,'rt')
    listx = []
    for line in input:
        if line.startswith("#"):
            continue
        cols = line.rstrip().split("\t")
        if int(cols[1]) == 16:
            strand = "+"   
        elif int(cols[1]) == 0:
            strand = "-"    
        else:
            continue
        key = cols[2]+"\t"+cols[3]+"\t"+strand
        tag_start = int(cols[3])
        tag_end  = tag_start + len(cols[9])
        seq      = cols[9]
        if key in dict_fimo:
            strandedness = (dict_fimo[key].split("\t"))[2] 
            motif_start        = int((dict_fimo[key].split("\t"))[0])
            motif_end          = motif_start + abs(int((dict_fimo[key].split("\t"))[1]))
            fimo_motif  = (dict_fimo[key].split("\t"))[3]
            pval = (dict_fimo[key].split("\t"))[4]
            extra_info = dict_fimo[key].split("\t")
            #print motif_start, motif_end, tag_start, tag_end, strandedness
            motif = check_strand_and_perform_operation(strandedness,motif_start,motif_end,tag_start,tag_end,seq)
            if motif != 0:
                listx.append(extra_info[5]+"\t"+extra_info[6]+"\t"+extra_info[7]+"\t"+fimo_motif+"\t"+pval+"\t"+motif+"\t"+cols[5]+"\t"+cols[4])
                #out.write(extra_info[5]+"\t"+extra_info[6]+"\t"+extra_info[7]+"\t"+fimo_motif+"\t"+pval+"\t"+motif+"\t"+cols[5]+"\t"+cols[4]+"\n")
    return(listx)
              
def check_strand_and_perform_operation(strandedness,motif_start,motif_end,tag_start,tag_end,seq):
    motif_length = abs(motif_end - motif_start)+1
    if strandedness == "++":
        if tag_start < motif_start:
            s = motif_start - tag_start
            e = s + abs(motif_end-motif_start)
            if len(seq[s:e+1]) != motif_length: # Checking if the read covered the whole motif or not.
                return(0)
            else:
                return(seq[s:e+1]) # 's' will be the distance to start from. Motif starts from s+1. So motif starts from s+1 to (s+1)+ e
        else:
            return(0)
    if strandedness == "+-":
        if tag_start > motif_start:
            return(0)
        else:
            s = motif_start - tag_start
            e = s + abs(motif_end-motif_start)
            if len(seq[s:e+1]) != motif_length:
                return(0)
            else:
                return(seq[s:e+1])
    # for both cases below motif start would be motif end and motif end would be motif start.
    if strandedness == "--":
        if tag_start > motif_start:
            return(0)
        else:
            s = abs(motif_start - tag_start)
            e = s + abs(motif_end - motif_start)
            if len(seq[s:e+1:]) != motif_length:
                return(0)
            else:
                # Need to return the reverse compliment of this motif here.
                #print seq[s:e+1:]
                return(reverseComplement(seq[s:e+1:]))
        
    if strandedness == "-+":
        if tag_start > motif_start:
            return(0)
        else:
            s = motif_start - tag_start
            e = s + abs(motif_end-motif_start)
            if len(seq[s:e+1]) != motif_length:
                return(0)
            else:
                # Need to return the reverse compliment of this motif.
                return(reverseComplement(seq[s:e+1:]))


def reverseComplement(sequence):
    complement = {'a':'t','c':'g','g':'c','t':'a','n':'n','A':'T','C':'G','G':'C','T':'A'}
    return "".join([complement.get(nt.upper(), '') for nt in sequence[::-1]])


usage = '''
input_paths may be:
- a single file.

example usages:
python 2_retrieve_seq_from_SAM.py  -s SAM <intersect_output>
'''.lstrip()

# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description

  

def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    #parser.add_option('-p', action='store', type='string', dest='PPdir',
    #                  help='Dir cotaining all the peak-pair file in gff format.')
    ##parser.add_option('-r', action='store',  dest='reffile',
    ##                  help='Reference genome in fasta format. Chromosme should match with input file')
    #parser.add_option('-u', action='store', type='int', dest='up_dist',default=20,
    #                  help='No of bases to subtract from the start coordinate, Default=20')
    #parser.add_option('-d', action='store',  type='int', dest='down_dist',default=20,
    #                  help='No of bases to add to the end coordinate, Default=20')
    parser.add_option('-s', action='store',  type='string', dest = 'SAMfile',
                      help='input SAM file')
    
    (options, args) = parser.parse_args()
    
    # Check if all the required arguments are provided, else exit     
    if not args:
        parser.print_help()
        sys.exit(1)
        
    #try:
    #    from scipy.stats import stats
    #except ImportError:
    #    print "You need to install Scipy(http://www.scipy.org/Download) before you can run this script."
    #    sys.exit(1)
    # 
                 
    if not os.path.isdir(args[0]):
        process_file(args[0],options)
        
 
if __name__ == "__main__":
    run() 
