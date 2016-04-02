import sys
import csv
import math
from igraph import *
from scipy import spatial

edgeList = []  
attrList = []
g = Graph()
memArray = []
sim = []


if len(sys.argv) != 2:
	print "python sac1.py <alpha as 0 | 0.5 | 1 >"
	exit(1)

alpha = float(sys.argv[1])

def main():
	load_data()
	g.add_vertices(range(324))
	g.add_edges(edgeList)
	global sim
	sim = createSimilarityMatrix(attrList)
	d = sac1(g, sim)
	writeTofile(d)


# Writes the final communities to the communities.txt file. Here d is the dictionary of communities.
def writeTofile(d):
	filename = "communities.txt"
	fo = open(filename, "w+")

	for key in d:
		if len(d[key]) != 0:
			v = ",".join(map(str, d[key]))
			fo.write(v)
			fo.write("\n")
	
	fo.close()		



# Creates similarity matrix based on attr list
def createSimilarityMatrix(attrList):
	tempSim = [[0 for x in range(len(attrList))] for x in range(len(attrList))]
	for i in range(len(tempSim)):
		for j in range(len(tempSim)):
			tempSim[i][j] = round(1 - spatial.distance.cosine(attrList[i], attrList[j]), 2)
	return tempSim 




# SAC-1 Algorithm
def sac1(g, sim):
	iterations = 1
	convergence = True
	d = {} # this dictionary stores a dictionary of communities at every stage. Key has the community id and value is a list of vertices that belong to that community

	# initialize the dictionary where each vertex is its own community
	for i in range(g.vcount()):
		d[i] = [i]

	# until convergence is happening and max iterations less than 15
	while iterations <= 15 and convergence:
		numberOfCommunities = g.vcount()
		vcount = g.vcount()
		memArray = range(g.vcount()) # Membership array. The value at each index indicates the community it belongs to
		print "Number of communities: ", numberOfCommunities

		#phase 1
		improvement = True
		phase1Iterations = 0
		while improvement and phase1Iterations < 15:
			for i in range(vcount):
				maxDelta = 0
				tmp = memArray[i]
				modularity = g.modularity(memArray)
				captureJ = i # could be = i to begin with
				for j in range(vcount):
					memArray[i] = memArray[j]
					newModularity = g.modularity(memArray)
					memArray[i] = tmp
					verticesInJ = [m for m, x in enumerate(memArray) if x == memArray[j]] # gives the vertices that are present in j's community
					deltaQAttr = sum([sim[i][v] for v in verticesInJ])
					l = len(verticesInJ)
					uniqueCommunities = set(memArray)
					if l!=0:
						deltaQAttr = float(deltaQAttr)/((l**2)*len(uniqueCommunities))
					else:
						deltaQAttr = float(deltaQAttr)/len(uniqueCommunities)

					delta = (newModularity - modularity) * float(alpha) + (1-alpha) * deltaQAttr # Composite modularity calculated
					
					
					if maxDelta < delta:
						maxDelta = delta
						captureJ = j

				if maxDelta > 0:
					memArray[i] = memArray[captureJ]
				else:
					improvement = False

			setCommunities = set(memArray)
			phase1Iterations = phase1Iterations + 1

		#phase 2
		uniq = list(set(memArray))

		# Simplify the membership array such that communities follow zero based indexing because igraph needs simplified zero based indexing for contracting vertices later. 
		cnt = 0
		for i, elem in enumerate(uniq):
			memArray = [cnt if (elem == x) else x for x in memArray]
			cnt = cnt+1

		
		# Since membership array will change, we need to store the actual vertices and communities mapping
		temp = {}
		for i in range(len(memArray)):
			temp[i] = []
		for i,m in enumerate(memArray):
			temp[m] = temp[m] + d[i]
		d = temp


		# Check for convergence i.e. whether we need to end phase 2 here
		if numberOfCommunities == len(uniq):
			convergence = False
			break	

		g.contract_vertices(memArray)

		# Now that graph has been compressed by combining vertices belonging to the same community, we need to combine attribute lists as well by taking 
		# sum of the attributes. This is like comibined attribute value for a community and then recalculate similiarity matrix
		newAttrList = []	
		for i, x in enumerate(memArray):
			t = [i if (elem == x) else elem for elem in memArray]
			selectedAttributesOfCombinedVertices = []
			for r in t:
				selectedAttributesOfCombinedVertices.append(attrList[r])
			newAttrList.append(map(sum, zip(*selectedAttributesOfCombinedVertices)))
			#print "********",len(newAttrList)

		sim = createSimilarityMatrix(newAttrList)	
		
		iterations = iterations + 1

	return d
	



# Loads the data from the given files. Creates edge list and attribute list. These lists are global

def load_data():
	with open("data/fb_caltech_small_edgelist.txt", "r") as f:
		for e in enumerate(f):
			edge = e[1].strip().split()
			edgeList.append((int(edge[0]),int(edge[1])))
	f.close()
	
	with open('data/fb_caltech_small_attrlist.csv', 'r') as f:
		next(f)
		for line in f:
			temp = line.split(',')
			temp = [int(v.strip()) for v in temp]
			attrList.append(temp)    
	f.close()
	   
				  
					
if __name__ == "__main__":
	main()