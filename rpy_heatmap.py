print "Loading molecular biology results,",
#Load the molecular biology results (second column)
#for each patient into a Python dictionary, using the
#patient ID as the key (first column)
input = open("ALL_mol_biol.tsv")
lines = input.readlines()
input.close()
mol_biol = dict()
for line in lines :
    #Strip any trailing newline
    if line[-1] == "\n" : line = line[:-1]
    #Break up the line using the tabs
    parts = line.split("\t")
    #Check there are exactly two fields...
    assert len(parts)==2
    #Store the patients mol biol results:
    mol_biol[parts[0]] = parts[1]
print "Done"

print "Loading expression data,",
input = open("ALL_exprs.tsv")
lines = input.readlines()
input.close()
#The genes are the rows, patients are columns.
#Deal with the header row...
line = lines[0]
#Strip any trailing newline
if line[-1] == "\n" : line = line[:-1]
#Break up the line using the tabs
parts = line.split("\t")
#Check first entry is blank
assert parts[0]==""
#Store the patient (column) names
col_names = parts[1:]
row_names = []
exprs = []

for line in lines[1:] :
    #Strip any trailing newline
    if line[-1] == "\n" : line = line[:-1]
    #Break up the line using the tabs
    parts = line.split("\t")
    #Check there are exactly n+1 fields...
    assert len(parts)==len(col_names)+1
    #Store the expression levels for this gene
    row_names.append(parts[0])
    exprs.append(map(float, parts[1:]))
    
assert len(exprs)==len(row_names)
print "Done"
print "%i patients, %i genes" % (len(col_names), len(row_names))

print "Converting expression data into Numeric array,",
import Numeric
exprs_as_array = Numeric.array(exprs, Numeric.Float)
print "Done"

def patient_colour(patient_id) :
    assert patient_id in mol_biol, \
           "Patient ID of '%s' not in list!" % patient_id
    if mol_biol[patient_id] == "ALL1/AF4" :
        return "#FF0000" # Red
    else :
        return "#0000FF" # Blue
patient_colours = map(patient_colour, col_names)

from rpy import *
print "Heatmap as PNG,",
r.png("heatmap_from_python.png", width=600, height=589)
r.heatmap(exprs_as_array,
          cexRow=0.5,
          labRow=row_names, labCol=col_names,
          ColSideColors = patient_colours,
          col = r.topo_colors(50))
r.dev_off()
print "Done"

print "Heatmap as PDF,",
r.pdf("heatmap_from_python.pdf")
r.heatmap(exprs_as_array,
          cexRow=0.3,
          labRow=row_names, labCol=col_names,
          ColSideColors = patient_colours,
          col = r.topo_colors(50))
r.dev_off()
print "Done"

#Note
#In R, the exprs object would include the row names
#and column names.  This isn't possible with an array
#in Python, so we must explicitly provide them.
