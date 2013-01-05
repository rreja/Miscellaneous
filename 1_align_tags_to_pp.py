import sys, os, re, difflib
from optparse import OptionParser , IndentedHelpFormatter
from operator import itemgetter
from subprocess import call
from collections import defaultdict

def process_files(idx,pp,options,output_folder,tmp_folder):
    print "INFO: Input idx file = "+idx+" and input REF file = "+pp
    # Creatin temp idx file.
    tmpidx = os.path.join(tmp_folder,"tmpidx.gff")
    outidx = open(tmpidx,"w")
    # Reding tag file and converting to gff.
    idxRead = open(os.path.join(options.idxDir,idx),"r")
    for line in idxRead:
        if line.startswith("chrom") or line.startswith("#"):
            continue
        chrom,start,fwtag,rvtag,ttag = line.rstrip().split("\t")
        if int(start) == 0:
            print "WARNING: You should re-check your idx files. BAM files were not properly converted to tab!"
            print line
            continue
        if int(fwtag) > 0:
            end = int(start)+1
            outidx.write(chrom+"\t.\t.\t"+start+"\t"+str(end)+"\t"+fwtag+"\t+\t.\t.\n")
        if int(rvtag) > 0:
            end = int(start)+1
            outidx.write(chrom+"\t.\t.\t"+start+"\t"+str(end)+"\t"+rvtag+"\t-\t.\t.\n")
    outidx.close()
    # Creating a temporary peak-pair file.
    tmpPP = os.path.join(tmp_folder,"tmpPP.gff")
    outpp = open(tmpPP,"w")
    ppRead = open(os.path.join(options.ppDir,pp),"r")
    for line in ppRead:
        if line.startswith("#"):
            continue
        chrom,junk,junk,start,end,score,junk,junk,dist = line.rstrip().split("\t")
        nstart = int(start) - options.up
        nend = int(end) + options.down #-1 # since end is alsways midpoint+1
        if nstart > 0:
            outpp.write(chrom+"\t.\t.\t"+str(nstart)+"\t"+str(nend)+"\t"+score+"\t.\t"+start+"-"+end+"\t"+dist+"\n")
        else:
             outpp.write(chrom+"\t.\t.\t1\t"+str(nend)+"\t"+score+"\t.\t"+start+"-"+end+"\t"+dist+"\n")
    outpp.close()
    print "INFO: Running intersectBed on these files."
    tmpInter = os.path.join(tmp_folder,"tmpintersect.txt")
    intersect = os.path.join(options.bedDir,"intersectBed")
    call(intersect+" -wo -a "+tmpPP+" -b "+tmpidx+" >"+tmpInter,shell=True)
    
    print "INFO: Creating CDT file"
    # Creating CDT file from this tempintersect file.
    outfile = os.path.join(output_folder,os.path.basename(idx).split(".")[0]+".cdt")
    cdtout = open(outfile,"w")
    print_header(cdtout,options)
    
    input = open(tmpInter,"r")
    cdt_dict = defaultdict(list)
    for line in input:
        tmp_dict = []
        data = line.rstrip().split("\t")
        key = data[0]+":"+data[7]+"\t"+data[8]
        if int(data[12]) == -1:
            cdt_dict[key].append("NA")
            continue
        refStart = int(data[7].split("-")[0])
        tagStart = int(data[12])
        distance = tagStart - refStart
        cdt_dict[key].append(str(distance)+":"+data[14])
    
    dict_with_ttag = {}
    # In order to print total tags uncomment the below lines at three places:
    for k,v in cdt_dict.items():
        tmpdict = {}
        # Comment this line for total tag
        line = k
        #line = ""
        totaltag = 0
        if v == "NA":
            for i in range((-1)*options.up,options.down+1):
                line = line+"\t0"
                totaltag = totaltag + 0
                
        for val in v:
            pos = int(val.split(":")[0])
            tag = (val.split(":"))[1]
            tmpdict[pos] = tag
        for i in range((-1)*options.up,options.down+1):
            if i in tmpdict:
                line = line+"\t"+tmpdict[i]
                totaltag = totaltag + int(tmpdict[i])
                
            else:
                line = line+"\t0"
        # Comment this line for total tags.
        cdtout.write(line+"\n")
        #dict_with_ttag[k.split("\t")[0]+"\t"+str(totaltag)+line] = totaltag
        
    #sort_and_print(dict_with_ttag,cdtout)
    

def process_files_tagStrand(idx,pp,options,output_folder,tmp_folder):
    
    print "INFO: Input idx file = "+idx+" and input REF file = "+pp
    # Creatin temp idx file.
    sense_tmpidx = os.path.join(tmp_folder,"Stmpidx.gff")
    anti_tmpidx = os.path.join(tmp_folder,"Atmpidx.gff")
    Soutidx = open(sense_tmpidx,"w")
    Aoutidx = open(anti_tmpidx,"w")
    # Reding tag file and converting to gff.
    idxRead = open(os.path.join(options.idxDir,idx),"r")
    for line in idxRead:
        if line.startswith("chrom") or line.startswith("#"):
            continue
        chrom,start,fwtag,rvtag,ttag = line.rstrip().split("\t")
        if int(start) == 0:
            print "WARNING: You should re-check your idx files. BAM files were not properly converted to tab!"
            print line
            continue
        if int(fwtag) > 0:
            end = int(start)+1
            Soutidx.write(chrom+"\t.\t.\t"+start+"\t"+str(end)+"\t"+fwtag+"\t+\t.\t.\n")
        if int(rvtag) > 0:
            end = int(start)+1
            Aoutidx.write(chrom+"\t.\t.\t"+start+"\t"+str(end)+"\t"+rvtag+"\t-\t.\t.\n")
    Soutidx.close()
    Aoutidx.close()
    # Reading peak-pair and aligning tags .
    tmpPP = os.path.join(tmp_folder,"tmpPP.gff")
    outpp = open(tmpPP,"w")
    ppRead = open(os.path.join(options.ppDir,pp),"r")
    for line in ppRead:
        if line.startswith("#"):
            continue
        chrom,junk,junk,start,end,score,junk,junk,dist = line.rstrip().split("\t")
        # new start and new end
        nstart = int(start) - options.up
        nend = int(end) + options.down #-1 # since end is alsways midpoint+1
        if nstart > 0:
            outpp.write(chrom+"\t.\t.\t"+str(nstart)+"\t"+str(nend)+"\t"+score+"\t.\t"+start+"-"+end+"\t"+dist+"\n")
        else:
             outpp.write(chrom+"\t.\t.\t1\t"+str(nend)+"\t"+score+"\t.\t"+start+"-"+end+"\t"+dist+"\n")
    outpp.close()
    print "INFO: Running intersectBed on these files."
    StmpInter = os.path.join(tmp_folder,"Stmpintersect.txt")
    AtmpInter = os.path.join(tmp_folder,"Atmpintersect.txt")
    intersect = os.path.join(options.bedDir,"intersectBed")
    call(intersect+" -wo -a "+tmpPP+" -b "+sense_tmpidx+" >"+StmpInter,shell=True)
    call(intersect+" -wo -a "+tmpPP+" -b "+anti_tmpidx+" >"+AtmpInter,shell=True)
    
    print "INFO: Creating CDT file"
    # Creating CDT file from this tempintersect file.
    Soutfile = os.path.join(output_folder,"sense_"+os.path.basename(idx).split(".")[0]+".cdt")
    Aoutfile = os.path.join(output_folder,"anti_"+os.path.basename(idx).split(".")[0]+".cdt")
    Scdtout = open(Soutfile,"w")
    Acdtout = open(Aoutfile,"w")
    print_header(Scdtout,options)
    print_header(Acdtout,options)
    
    input = open(StmpInter,"r")
    Scdt_dict = defaultdict(list)
    for line in input:
        tmp_dict = []
        data = line.rstrip().split("\t")
        key = data[0]+":"+data[7]+"\t"+data[8]
        if int(data[12]) == -1:
            Scdt_dict[key].append("NA")
            continue
        refStart = int(data[7].split("-")[0])
        tagStart = int(data[12])
        distance = tagStart - refStart
        Scdt_dict[key].append(str(distance)+":"+data[14])
    
    input = open(AtmpInter,"r")
    Acdt_dict = defaultdict(list)
    for line in input:
        tmp_dict = []
        data = line.rstrip().split("\t")
        key = data[0]+":"+data[7]+"\t"+data[8]
        if int(data[12]) == -1:
            Acdt_dict[key].append("NA")
            continue
        refStart = int(data[7].split("-")[0])
        tagStart = int(data[12])
        distance = tagStart - refStart
        Acdt_dict[key].append(str(distance)+":"+data[14])
    
    Sdict_with_ttag = {}
    Adict_with_ttag = {}
    dict_with_ttag = {}
    # In order to print total tags look at the comments below. Only positive or negative numeber comments should be present at a time.
    for k,v in Scdt_dict.items():
        Stmpdict = {}
        Atmpdict = {}
        # Comment this line for total tag
        Sline = k   #1
        Aline = k    #2
        #Sline = ""     #-1
        #Aline = ""    #-2
        Stotaltag = 0
        Atotaltag = 0
        ## Check for "NA" values.
        if v == "NA":
            for i in range((-1)*options.up,options.down+1):
                Sline = Sline+"\t0"
                Stotaltag = Stotaltag + 0
        if Acdt_dict[k] == "NA":
            for i in range((-1)*options.up,options.down+1):
                Aline = Aline+"\t0"
                Atotaltag = Atotaltag + 0
        #################################################        
        for val in v:
            pos = int(val.split(":")[0])
            tag = (val.split(":"))[1]
            Stmpdict[pos] = tag
        for val in Acdt_dict[k]:
            pos = int(val.split(":")[0])
            tag = (val.split(":"))[1]
            Atmpdict[pos] = tag
        for i in range((-1)*options.up,options.down+1):
            if i in Stmpdict:
                Sline = Sline+"\t"+Stmpdict[i]
                Stotaltag = Stotaltag + int(Stmpdict[i])
                
            else:
                Sline = Sline+"\t0"
        for i in range((-1)*options.up,options.down+1):
            if i in Atmpdict:
                Aline = Aline+"\t"+Atmpdict[i]
                Atotaltag = Atotaltag + int(Atmpdict[i])
                
            else:
                Aline = Aline+"\t0"
        # Comment this line for total tags.
        Scdtout.write(Sline+"\n") #3
        Acdtout.write(Aline+"\n") #4
        #Sdict_with_ttag[k.split("\t")[0]] = str(Stotaltag)+Sline #-3
        #Adict_with_ttag[k.split("\t")[0]] = str(Atotaltag)+Aline  #-4
        
    #sort_and_print(Sdict_with_ttag,Adict_with_ttag,Scdtout,Acdtout) #-5


 
def sort_and_print(Sdicti,Adicti,Sout,Aout):
    tmpdict = {}
    for k,v in Sdicti.items():
        tmpdict[k] = int(v.split("\t")[0])
    sorted_x = sorted(tmpdict.iteritems(), key=itemgetter(1),reverse=True)
    for j in sorted_x:
        Sout.write(j[0]+"\t"+Sdicti[j[0]]+"\n")
        Aout.write(j[0]+"\t"+Adicti[j[0]]+"\n")
                
def print_header(out,options):
    line = "ID\tcw-dist"
    for i in range((-1)*options.up,options.down+1):
        line = line+"\t"+str(i)
    out.write(line+"\n")
    
            

usage = '''
input_paths may be:
- a directory to run on all files in them


example usages:
python 1_align_tags_to_pp.py -i <index_dir> -d <peakPairs_dir>
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-i', action='store', type='string', dest='idxDir',
                      help='Path to index directory.')
    parser.add_option('-r', action='store', type='string', dest='ppDir',
                      help='Path to reference feature directory.')
    parser.add_option('-u', action='store', type='int', dest='up',default=100,
                      help='Distance to subtract from the start coordinate.default=100')
    parser.add_option('-d', action='store', type='int', dest='down',default=100,
                      help='Distance to add to the end coordinate. default=100')
    parser.add_option('-b', action='store', type='string', dest='bedDir',default="/usr/local/bin",
                      help='Directory containing BedTools.Default="/usr/local/bin"')
    parser.add_option('-s', action='store', type='int', dest='tagStrand',default=0,
                      help='1=> Consider tag strandedness only. 0=> Do not consider any strandedness.Default=0')
    #parser.add_option('-S', action='store', type='int', dest='allStrand',default=0,
    #                  help='1=> Consider strandedness of both reference and tags.0=> Do not consider any strandedness.Default=0')
    
    (options, args) = parser.parse_args()
        
    idx_files = []
    pp_files = []
    
    
    output_folder = os.path.join(options.idxDir,'CDT/') 
    if not os.path.exists(output_folder): os.makedirs(output_folder)
    
    tmp_folder = os.path.join(options.idxDir,'tmp/') 
    if not os.path.exists(tmp_folder): os.makedirs(tmp_folder)
    
    if not os.path.exists(options.idxDir) or not os.path.exists(options.ppDir):
        parser.error('Path %s does not exist.')
        
    if os.path.isdir(options.idxDir):
        for fname in os.listdir(options.idxDir):
            if fname.endswith("tab") or fname.endswith("idx"):
                idx_files.append(fname)
        
        for fname in os.listdir(options.ppDir):
            if not fname.endswith('gff'):
                continue
            Oratio = 0
            for file_idx in idx_files:
                Nratio = difflib.SequenceMatcher(None,fname,file_idx).ratio()
                if(Nratio >= Oratio):
                    idx = file_idx
                    pp = fname
                    Oratio = Nratio
            if options.tagStrand == 1:
                process_files_tagStrand(idx,pp,options,output_folder,tmp_folder)
            #elif options.allStrand == 1:
            #    process_files_allStrand(idx,pp,options,output_folder,tmp_folder)
            else:
                process_files(idx,pp,options,output_folder,tmp_folder)
                
        print "Completed creating all CDT.Result in "+output_folder
        os.system("rm -fr "+tmp_folder)
    
if __name__ == "__main__":
    run() 
    