import sys, os, re, difflib, csv
from optparse import OptionParser , IndentedHelpFormatter
from PIL import Image


def process_files(sense,anti,output_folder):
    
    outfile = os.path.join(output_folder,os.path.basename(sense).split(".")[0]+".png")
    background = Image.open(sense)
    overlay = Image.open(anti)
    
    background = background.convert("RGBA")
    overlay = overlay.convert("RGBA")
    
    new_img = Image.blend(background, overlay, 0.5)
    new_img.save(outfile,"PNG")








usage = '''
input_paths may be:
- a directory to run on all files in them
- "." to run in the current directory

example usages:
python overlay_two_images.py  <DIR containing both sense and antisense cluster plots>
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
        
    antisense_files = []
    sense_files = []
    
    output_folder = os.path.join(args[0],'Combined_images/') 
    if not os.path.exists(output_folder): os.makedirs(output_folder) 
    
    
    if not os.path.exists(args[0]):
        parser.error('Path %s does not exist.' % args[0])
    if os.path.isdir(args[0]):
        for fname in os.listdir(args[0]):
            if fname.endswith("png") or fname.endswith("jpg"):
                # intellignetly join paths without worrying about '/'
                fpath = os.path.join(args[0], fname)
                matchobj = re.match(r'(.*)anti(.*)',fname)
                if matchobj:
                    antisense_files.append(fpath)
                else:
                    sense_files.append(fpath)
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
            process_files(file2_sense,file1_anti,output_folder)
            print "Completed. Output in  "+output_folder
            
    
if __name__ == "__main__":
    run() 
    