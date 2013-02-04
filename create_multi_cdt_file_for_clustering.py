import sys, os, operator, pybedtools
from optparse import OptionParser , IndentedHelpFormatter
from operator import itemgetter
from itertools import izip, cycle, tee

def process_files(idxdir,options,outfile):
    # First creating an interval file with up and down distance
    encData = {}
    filehash = {}
    count = 1
    idxData = {}
    alldata = {}
    sortdataby = {}
    
    gff_file = pybedtools.BedTool(options.enrcfile).slop(g=options.gfile,l=options.up,r=options.down)
    out = open(outfile,"w")
   

    # Reading all idx files in the directory one by one
    for fname in os.listdir(idxdir):
        if fname.endswith("idx") or fname.endswith("tab"):
            filehash[count] = fname
            input = open(os.path.join(idxdir,fname),"r")
            for line in input:
                if line.startswith("chrom") or line.startswith("#"):
                    continue;
                tmp = line.rstrip().split("\t")
                ttag = int(tmp[2]) + int(tmp[3])
                #ttag = float(tmp[4])
                chrom = tmp[0]
                start = tmp[1]
                idxData[str(count)+":"+chrom+":"+start] = ttag
            count = count + 1
    #print_header(out,options,len(filehash.keys()))
    
    for line in gff_file:
        start = int(line[3])
        end = int(line[4])
        if (abs(end - start)-1) != (options.up+options.down):
            
            print "This region has a small interval "+chrom+":"+line[3]+"-"+line[4]+"\t."
            continue
        chrom = line[0]
        pline = chrom+":"+line[3]+"-"+line[4]+"\t."
        # big list including binned values for all the factors 
        binnedlist = []
        for key,val in filehash.items():
            # list for one factor at a time.
            tmpbinlist = []
            if key == options.sortby:
                totaltag = 0
            for j in range(start,end+1):
                keyname = str(key)+":"+chrom+":"+str(j)
                if keyname in idxData:
                    pline = pline +"\t"+ str(idxData[keyname])
                    tmpbinlist.append(idxData[keyname])
                    if key == options.sortby:
                        totaltag = totaltag + idxData[keyname]
                else:
                    pline = pline + "\t0"
                    tmpbinlist.append(0)
            binnedlist = binnedlist + bindata(tmpbinlist,options)     
        ##out.write(pline+"\n")
        alldata[chrom+":"+line[3]+"-"+line[4]] = binnedlist
        sortdataby[chrom+":"+line[3]+"-"+line[4]] = totaltag
        ##out.flush()
    ##out.close()
    sort_and_print(alldata,sortdataby,out)

    for k,v in filehash.items():
        print k,v

def sort_and_print(alldata,sortdataby,out):
    sorted_x = sorted(sortdataby.iteritems(), key=itemgetter(1),reverse=True)
    for j in sorted_x:
        line = j[0]
        for m,n in enumerate(alldata[j[0]]):
            line = line +"\t"+str(n)
        out.write(line+"\n")
        
def print_header(out,options,l):
    line = "ID\tWT"
    for i in range((-1)*options.up*l,(options.down)*l+1):
        line = line+"\t"+str(i)
    out.write(line+"\n")
    
def bindata(taglist,options):
    bintaglist = []
    summation = 0
    tmplist = range(0,len(taglist)+1,options.bins)
    for elem , next_elem in pairwise(tmplist):
        for j in taglist[elem:next_elem:1]:
            summation = summation+j
        bintaglist.append(summation)
        summation = 0
    return(bintaglist)
    
    
def pairwise(seq):
    a, b = tee(seq)
    next(b)
    return izip(a, b)
    
usage = '''
input_paths may be:
- a single file.

example usages:
python create_multi_cdt_file_for_clustering [options] /dir_to_idx_files/
'''.lstrip()

# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description

  
def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-e', action='store',  dest='enrcfile',
                      help='File containing the enriched regions.')
    parser.add_option('-u', action='store', type='int', dest='up',default=30,
                      help='Upstream distance to go .Default=30')
    parser.add_option('-d', action='store', type='int', dest='down',default=30,
                      help='Downstream distance to go.Default=30')
    parser.add_option('-g', action='store',  type='string', dest = 'gfile',
                      help='text file containing the size of each chromosome in the format: chr start end, ex: chr1 1 123456.Required Parameter.')
    parser.add_option('-s', action='store',  type='int', dest = 'sortby',default=1,
                      help='Sort by occupancy of this factor.Default = 0.')
    parser.add_option('-b', action='store',  type='int', dest = 'bins',default=5,
                      help='Bin size, Default= 5.')
    (options, args) = parser.parse_args()
    
    # Check if all the required arguments are provided, else exit     
    if not args:
        parser.print_help()
        sys.exit(1)
    outfile = os.path.join(os.path.dirname(options.enrcfile),"all_factors_merged.cdt")
    if not os.path.isdir(args[0]):
        print "Input data direcotry containing all idx files."
    else:
        process_files(args[0],options,outfile)
    
if __name__ == "__main__":
    run() 