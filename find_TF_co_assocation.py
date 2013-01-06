import sys, os, pybedtools, operator
from optparse import OptionParser , IndentedHelpFormatter
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from optparse import OptionParser , IndentedHelpFormatter
import matplotlib.mlab as mlab
import matplotlib.colors as cl
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from operator import itemgetter
from timeit import Timer
import matplotlib.colorbar as clb

def process_file(infile,PPFiles,options):
    input = open(infile,'rt')
    fullfile = {}
    y_cord = 0
    bedfile = ""
    dict_main = {}
    count = 0
    x = []
    y = []
    vals = []
    outfile = os.path.join(options.PPdir,os.path.splitext(os.path.basename(infile))[0]+"_cp-association.png")
    # Start reading the peak-pair file and sort it.
    for line in input:
        if line.startswith("#"):
            continue
        if(len(line.split("\t")) == 9):
            #countLines = countLines+1
            chrom,junk,junk,start,end,tag,strand,junk,attr = line.split("\t")
            fullfile[line] = float(tag)
    sortedfile = sorted(fullfile.iteritems(),key=operator.itemgetter(1),reverse=True)
    # Create a string of the sorted gff file to create bedtool object.
    for k in sortedfile:
        bedfile = bedfile+k[0]
        row = k[0].split("\t")
        dict_main[row[0]+"\t"+row[3]+"\t"+row[4]] = count
        count +=1 
    # Create peak-pair bedtool object from the string
    a = pybedtools.BedTool(bedfile, from_string=True)
    # Run through all the fimo output in the FIMO directory conatining fimo output for various motifs corresponding to the same peak-pair file
    for pp in PPFiles:
        intersect = {}
        linecount = 0
        matchcount = 0
        try:
            
            # Calculate intersection.
            b = a.intersect(pybedtools.BedTool(pp).slop(g=options.gfile,l=options.up_dist,r=options.down_dist),u=True)
        except pybedtools.helpers.BEDToolsError:
            print "No intersection found with peak-pairs from "+fimo
            continue
            
        # Create a dictionary of intersection results.
        for i in b:
            key = i[0]+"\t"+i[3]+"\t"+i[4]
            if key in dict_main:
                x.append(dict_main[key])
                y.append(y_cord)
                vals.append(float(i[5]))
                #vals.append(1.0)
        
        y_cord +=1
        
    draw_cluster_plot(x,y,vals,outfile)
    
def draw_cluster_plot(x,y,vals,outfile):
    
    fig = Figure(figsize=(9,6))
    canvas = FigureCanvas(fig)
    ax = fig.add_subplot(111,frame_on=False)    
    ax.set_aspect('equal',adjustable='box',anchor='C')
    ax.set_xlim(0,max(x))
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    cmap = cm.get_cmap('Blues')
    cmap.set_over(color='b')
    cmap.set_under(color='w')
    bounds = [100,150,200,250,300,350,400,450]
    norm = cl.BoundaryNorm(bounds, ncolors=256, clip = False)
   
    
    
    img = ax.scatter(x,y,c=vals,marker='o',cmap=cmap,norm=norm,s=.25,edgecolors='none')
    #clb.ColorbarBase(ax,cmap=cmap, norm=norm, boundaries=bounds, ticks=bounds)    
    canvas.print_figure(outfile,dpi=200,bbox_inches='tight')



usage = '''
input_paths may be:
- a directory to run on all files in them
- a single file.

example usages:
python find_TF_co_association.py -d <dir_with_all_TF> -g sg07.txt -u 50 -d 50  -g <genome_length_file>
'''.lstrip()

# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description

  

def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-p', action='store', type='string', dest='PPdir',
                      help='Dir cotaining all the peak-pair file in gff format.')
    #parser.add_option('-r', action='store',  dest='reffile',
    #                  help='Reference genome in fasta format. Chromosme should match with input file')
    parser.add_option('-u', action='store', type='int', dest='up_dist',default=20,
                      help='No of bases to subtract from the start coordinate, Default=20')
    parser.add_option('-d', action='store',  type='int', dest='down_dist',default=20,
                      help='No of bases to add to the end coordinate, Default=20')
    parser.add_option('-g', action='store',  type='string', dest = 'gfile',
                      help='text file containing the size of each chromosome in the format: chr start end, ex: chr1 1 123456')
    
    (options, args) = parser.parse_args()
    
    # Check if all the required arguments are provided, else exit     
    if not args:
        parser.print_help()
        sys.exit(1)
        
    try:
        from scipy.stats import stats
        import pybedtools
    except ImportError:
        print "You need to install Scipy(http://www.scipy.org/Download) before you can run this script."
        print "If you have Scipy installed, then you should check if you have pybedtools installed"
        sys.exit(1)
     
    # Create array to store all the fimo output files corresponding to different motifs.
    PPFiles = []
    for name in os.listdir(options.PPdir):
            if name.endswith('.gff'):
                PPFiles.append(os.path.join(options.PPdir, name))
                     
    if not os.path.isdir(args[0]):
        process_file(args[0],PPFiles,options)
        
 
if __name__ == "__main__":
    run() 
