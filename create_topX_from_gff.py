import sys, os, re, difflib
from optparse import OptionParser , IndentedHelpFormatter
from operator import itemgetter
from subprocess import call
from collections import defaultdict


def process_files(infile,options,outdir):
    input = open(infile,"r")
    data = {}
    outfile = os.path.join(outdir,os.path.basename(infile))
    out = open(outfile,"w")
    for line in input:
        if line.startswith("chrom") or line.startswith("#"):
            out.write(line)
            continue
        chrom,junk,junk,start,end,score,strand,junk,junk = line.rstrip().split("\t")
        data[line.rstrip()] = float(score)
    

    sorted_x = sorted(data.iteritems(), key=itemgetter(1),reverse=True)
    count = 0
    for j in sorted_x:
        if not count > options.N:
            out.write(j[0]+"\n")
        count = count + 1
        
    out.close()
        
    print "INFO: Complete executing file: "+os.path.basename(infile)





usage = '''
input_paths may be:
- a directory containing gff files.


example usages:
python create_topX_from_gff.py -n <number>
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-n', action='store', type='int', dest='N',default=1000,
                      help='Give top "N" lines from a file.')
    
    
    (options, args) = parser.parse_args()
        
    if not args:
        parser.print_help()
        sys.exit(1)
        
    output_folder = os.path.join(args[0],'top_'+str(options.N)) 
    if not os.path.exists(output_folder): os.makedirs(output_folder)
    
    
    if not os.path.exists(args[0]):
        parser.error('Path %s does not exist.' % args[0])
    if os.path.isdir(args[0]):
        for fname in os.listdir(args[0]):
            if fname.endswith("gff"):
                # intellignetly join paths without worrying about '/'
                fpath = os.path.join(args[0], fname)
                process_files(fpath,options,output_folder)   

         
    
if __name__ == "__main__":
    run() 
    