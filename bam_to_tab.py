import sys, os, re
from optparse import OptionParser , IndentedHelpFormatter


def process_file(fname,samfile,outfile,tmpfile):
    fwd = {}
    rev = {}
    list_all = []
    out = open(outfile,"w")
    os.system(" /usr/local/bin/bamToBed -i "+fname+" >"+samfile)
    input = open(samfile,"rt")
    for line in input:
        tmp = line.rstrip().split("\t")
        if int(tmp[4]) == 0:
            continue
        if tmp[5] == "+":
            # BED file hence add to the start only.
            start = int(tmp[1]) + 1
            string = tmp[0]+"\t"+str(start)
            if string in fwd:
                fwd[string] = fwd[string] + 1
                
            else:
                fwd[string] = 1
                list_all.append(string)
        elif tmp[5] == "-":
            string = tmp[0]+"\t"+tmp[2]
            if string in rev:
                rev[string] = rev[string] + 1
                
            else:
                rev[string] = 1
                list_all.append(string)            
    for val in list_all:
        if val in fwd and val in rev:
            total = fwd[val]+rev[val]
            out.write(val+"\t"+str(fwd[val])+"\t"+str(rev[val])+"\t"+str(total)+"\n")
        elif val in fwd:
            total = fwd[val]
            out.write(val+"\t"+str(fwd[val])+"\t0\t"+str(total)+"\n")
        elif val in rev:
            total = rev[val]
            out.write(val+"\t0\t"+str(rev[val])+"\t"+str(total)+"\n")
        out.flush()
    os.system("sort -k 1,1 -k 2,2n "+outfile+" >"+tmpfile)
    os.system("mv "+tmpfile+" "+outfile)
    print "Converted "+fname+" to "+outfile
    #out.write("chrom\tindex\tforward\treverse\tvalue\n")


usage = '''
input_paths may be:
- a directory only.


example usages:
python convert_to_tab.py <dir_with_data>
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    
    (options, args) = parser.parse_args()
    
    if not args:
        parser.print_help()
        sys.exit(1)
    
    outdir = os.path.join(args[0], "_tab")
    if not os.path.exists(outdir): os.makedirs(outdir)
    
    if not os.path.isdir(args[0]):
        print "The input to this script is a directory."
    else:
                
        if not os.path.exists(args[0]):
            parser.error('Path %s does not exist.' % args[0])
        for name in os.listdir(args[0]):
            if not name.endswith(".bam"):
                continue
            fname = os.path.join(args[0], name)
            samfile = os.path.join(args[0], os.path.splitext(name)[0]+".bed")
            outfile = os.path.join(outdir,os.path.splitext(name)[0]+".tab")
            tmp = os.path.join(outdir,os.path.splitext(name)[0]+".tmp")
            process_file(fname,samfile,outfile,tmp)
    
    
if __name__ == "__main__":
    run() 