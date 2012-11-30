import sys, os, re, difflib
from optparse import OptionParser , IndentedHelpFormatter
from operator import itemgetter

def process_files(sense, antisense, options,output_folder):
    #print antisense
    #sys.exit(1)
    input1 = open(sense,'rt')
    sense_dict,up,down,bins = parse_file(input1)
    input2 = open(antisense,'rt')
    antisense_dict,up,down,bins = parse_file(input2)
    outfile = os.path.join(output_folder, os.path.basename(sense).split('.')[0] + "_normTag.txt")
    out = open(outfile,"w")
    out.write("Identifier\tNumerator_sum\tDenominator_sum\tRatio\n")
    ranges = open(options.rangeFile, "rt")
    for line in ranges:
        if line.startswith("#"):
            continue
        elif line.startswith("NS"):
            NS = line.rstrip().split(";")
            NS.remove('')
        elif line.startswith("DS"):
            DS = line.rstrip().split(";")
            DS.remove('')
        elif line.startswith("NA"):
            NA = line.rstrip().split(";")
            NA.remove('')
        elif line.startswith("DA"):
            DA = line.rstrip().split(";")
            DA.remove('')
    sum_NS = sum_by_range(sense_dict,NS)
    sum_DS = sum_by_range(sense_dict,DS)
    sum_NA = sum_by_range(antisense_dict,NA)
    sum_DA = sum_by_range(antisense_dict,DA)
    for k,v in sum_NS.items():
        numerator = v + sum_NA[k]
        denominator = sum_DS[k] + sum_DA[k]
        if denominator != 0:
            ratio = float(numerator)/denominator
        else:
            ratio = "NA"
        out.write(k+"\t"+str(numerator)+"\t"+str(denominator)+"\t"+str(ratio)+"\n")
    print "Output of "+os.path.basename(sense)+" and "+os.path.basename(antisense)+" written to the output folder "+output_folder

# Execute this function when only one strand is present. like tags aligned to peak-pair midpoint    
def process_files_nostrand(f,options,output_folder):
    input = open(f,'rt')
    fdict,up,down,bins = parse_file(input)
    outfile = os.path.join(output_folder, os.path.basename(f).split('.')[0] + "_normTag.txt")
    out = open(outfile,"w")
    out.write("Identifier\tNumerator_sum_per_base\tDenominator_sum_per_base\tRatio\n")
    ranges = open(options.rangeFile, "rt")
    for line in ranges:
        if line.startswith("#"):
            continue
        elif line.startswith("NS"):
            NS = line.rstrip().split(";")
            NS.remove('')
        elif line.startswith("DS"):
            DS = line.rstrip().split(";")
            DS.remove('')
    sum_NS = sum_by_range(fdict,NS)
    sum_DS = sum_by_range(fdict,DS)
    for k,v in sum_NS.items():
        numerator = v
        denominator = sum_DS[k]
        if denominator != 0:
            ratio = float(numerator)/denominator
        else:
            ratio = "NA"
        out.write(k+"\t"+str(numerator)+"\t"+str(denominator)+"\t"+str(ratio)+"\n")
    print "Output of "+os.path.basename(f)+" written to the output folder "+output_folder
         
    
        
def sum_by_range(dicti,RG):
    tmpdict = {}
    # Initialization of keys in tmpdict.
    for k,v in dicti.items():
        tmpdict[k] = 0
    # Starting loop over all given ranges.
    for i in RG:
        tmp = i.split("=")
        rang = tmp[1].split(":")
        length = abs(int(rang[1]) - int(rang[0]))
        for key,val in dicti.items():
            s = 0
            # Look at each line of the cdt file where val is a dictionary of the form coord=>tag  -40=>3, 30=>5.
            for k,v in val.items():
                if k > int(rang[0]) and k < int(rang[1]):
                    s = s+ int(v)
                else:
                    s = s+0
            tmpdict[key] = tmpdict[key]+ (float(s)/length)
    return(tmpdict)
                    
         
def parse_file(input):
    file_hash = {}
    first_line = input.readline().rstrip().split("\t")
    # Somce CDT files, the first two columns are space seperated and not tab seperated. This will take care of that.
    try:
        up_distance = int(first_line[2])
    except ValueError:
        up_distance = int(first_line[3])
    # Reading the last column to get the up stream distance
    down_distance = int(first_line[-1])
    # calculating bin
    bins = abs(int(first_line[5]) - int(first_line[4]))
    # Skipping the second line which has no information.
    input.next()
    # Iterating through the rest of the file
    for line in input:
        data = line.rstrip().split("\t")
        header = data.pop(0)
        data.pop(0)
        file_hash[header] = convert_list_to_dict(data,up_distance,down_distance,bins)
        
    return(file_hash,up_distance,down_distance,bins)


def convert_list_to_dict(list1,up,down,bins):
    new_dict = {}
    list2 = map(int, list1)    
    sum = 0
    c=0
    for i in range(up,down+1,bins):
        if list2[c] > 0:
            new_dict[i] = list2[c]
        c+=1    
    return(new_dict)

def print_dict(dict1):
    for k,v in dict1.items():
        if v:
            print k,v

usage = '''
input_paths may be:
- a directory to run on all files in them


example usages:
python local_normalization_in_cdt.py -r <rangeInfo file> <path_to_cdt_folder>
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-r', action='store', type='string', dest='rangeFile',
                      help='A file specifying the range to consider for summing up tags')
    #parser.add_option('-d', action='store', type='int', dest='down_distance',default=0,
    #                  help='Downstream distance to use in plotting. By default, it will take the downstream distance from the CDT file.')


    (options, args) = parser.parse_args()
    
    if not args:
        parser.print_help()
        sys.exit(1)
        
    antisense_files = []
    sense_files = []
    
    output_folder = os.path.join(args[0],'Normalized_CDT/') 
    if not os.path.exists(output_folder): os.makedirs(output_folder)
    
    #outfile =  os.path.join(output_folder, "all_normalized_cdt.txt")
    #out = open(outfile,"w")
    
    if not os.path.exists(args[0]):
        parser.error('Path %s does not exist.' % args[0])
    if os.path.isdir(args[0]):
        for fname in os.listdir(args[0]):
            if fname.endswith("cdt"):
                # intellignetly join paths without worrying about '/'
                fpath = os.path.join(args[0], fname)
                matchobj = re.match(r'(.*)anti(.*)',fname)
                if matchobj:
                    antisense_files.append(fpath)
                else:
                    sense_files.append(fpath)
        if not len(antisense_files) == 0:
            for f1 in antisense_files:
                Oratio = 0
                file1_anti = ""
                file2_sense = ""
                for f2 in sense_files:
                    # Calculating the distance ration between two file names and grouping together based upon the highest ratio
                    # Requirement: Antisense file should have anti word in it and both sense and antisense files should have similar names
                    Nratio = difflib.SequenceMatcher(None,f1,f2).ratio()
                    if(Nratio >= Oratio):
                        file1_anti = f1
                        file2_sense = f2
                        Oratio = Nratio        
                process_files(file2_sense,file1_anti,options,output_folder)
                
        else:
            for f in sense_files:
                process_files_nostrand(f,options,output_folder)    
    
if __name__ == "__main__":
    run() 
    
    