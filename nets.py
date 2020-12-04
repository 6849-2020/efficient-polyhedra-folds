"""

h = 3**0.5/2
triangle = Polygon([(0,0), (0.5,h), (1,0)])
Tetrahedron = Polyhedron([(triangle,"ABD"),(triangle,"BAC"),
                          (triangle,"CDB"),(triangle,"ADC")],
                         {(v,w):1 for v,w in permutations("ABCD",2)})

# the second argument is vertex distances along the polyhedron's surface
# (for the purpose of determining faithfulness)

polys = GeneralNet(Tetrahedron, radius=3.5,
                   checkFaithful=True, debug=False)
# this is a list of (Polygon, labeled vertex string like "ABD")

showPolys(polys) # pretty picture

"""

import shapely
import matplotlib.pyplot as plt
from shapely.geometry import Point, LinearRing, Polygon
from shapely.ops import cascaded_union
from useful import window
import numpy as np
from random import choice
from itertools import permutations
import math

def showPolys(polygons, colorDict=None):
    """ requires type shapely.Polygon """
    for poly, vtxs in polygons:
        if colorDict and vtxs in colorDict:
            color = "#"+colorDict[vtxs]
        else:
            color = ((hash(vtxs)%10)/11,
                     (hash(vtxs)%11)/12,
                     (hash(vtxs)%13)/14)
        x, y = poly.exterior.xy
        if len(vtxs)>1:
            plt.fill(x,y,color=color,lw=0) # plt.plot for not filling
            # shrink slightly around the centroid
            #centroidShrink = 0.99
            #x = [centroidShrink*r+(1-centroidShrink)*poly.centroid.x for r in x]
            #y = [centroidShrink*r+(1-centroidShrink)*poly.centroid.y for r in y]
        else:
            plt.plot(x,y,color=color,lw=0.5)
        # add text!
        '''char = "A"
        dx, dy = x[-2]-x[0], y[-2]-y[0]
        rotation = math.atan2(dy, dx)*360/math.tau
        midx, midy = (x[0]+x[-2])/2, (y[0]+y[-2])/2
        normalx, normaly = -dy, dx
        plt.text(midx+0.3*normalx, midy+0.3*normaly,
                 char, size=10, rotation=rotation,
                 ha="center", va="center", color=color)'''
    plt.show()

def GeneralNet(polyhedron, radius, checkFaithful=False, debug=False):
    """ Returns a list of (Polygon, vertex sequence). Sticks in the given
    radius from the origin. TODO: Some kind of optimization for using
    different faces more in close regions?? """
    basePoly, baseVtxs = polyhedron.faces[0]
    polys = [(basePoly, baseVtxs)]
    existingArea = basePoly
    unusedEdges = [] # edges of polygon
    for (vtx1, x_y1, vtx2, x_y2) in edges(basePoly, baseVtxs):
        unusedEdges.append((vtx1, x_y1, vtx2, x_y2))
    while unusedEdges:
        # choose a side from the unused edges of the figure
        # check for appropriate distance from origin
        index, (vtx1, x_y1, vtx2, x_y2) = choice(list(enumerate(unusedEdges)))
        if norm(x_y1) > radius:
            unusedEdges = unusedEdges[:index]+unusedEdges[index+1:]
            continue
        # assign the polygon on the other side of it
        facePoly, faceVtxs = polyhedron.faceWith(vtx2, vtx1)
        v1,p1,v2,p2 = [(v1, p1, v2, p2) for (v1, p1, v2, p2)
                       in edges(facePoly, faceVtxs)
                       if v1==vtx2 and v2==vtx1][0]
        T = SegmentTransform((p1,p2), (x_y2,x_y1))
        facePolyT = T(facePoly)
        # check for intersection
        if facePolyT.intersection(existingArea).area > 1e-10:
            unusedEdges = unusedEdges[:index]+unusedEdges[index+1:]
            continue
        # ((check for faithfulness))
        if checkFaithful:
            if not faithful((facePolyT, faceVtxs), polys,
                            polyhedron.vertexDistances):
                unusedEdges = unusedEdges[:index]+unusedEdges[index+1:]
                continue
        # either (1) discard polygon and mark edge as used, or
        # (2) add polygon to list of polys and update edges ((and vtx map))
        polys.append((facePolyT, faceVtxs))
        for (vtx1, x_y1, vtx2, x_y2) in edges(facePolyT, faceVtxs):
            unusedEdges.append((vtx1, x_y1, vtx2, x_y2))
        existingArea = existingArea.union(facePolyT)
        if debug:
            showPolys([poly for poly,vtxs in polys])
            showPolys([existingArea])
    return polys

def faithful(newPoly, polys, vertexDistances, tol=1e-6):
    """ Can the new polygon be faithfully added? Or would this create a
    pair of vertices closer on the polyhedron than on the paper? """
    for label, vertex in vertices(*newPoly):
        for oldPoly in polys:
            for oldLabel, oldVertex in vertices(*oldPoly):
                if oldLabel != label and (oldLabel, label) in vertexDistances:
                    d = norm((vertex[0]-oldVertex[0], vertex[1]-oldVertex[1]))
                    if d + tol < vertexDistances[oldLabel, label]:
                        return False
    return True
    

class Polyhedron:
    # polyhedron.faces --> list of faces
    # polyhedron.faceWith(v1, v2) --> face
    def __init__(self, faces, vertexDistances=None):
        self.faces = faces
        self.vertexDistances = vertexDistances
        self.colorDict = {vtxs:color for ((poly, vtxs), color) in
                          zip(faces, kelly_colors[1:])}

    def faceWith(self, v1, v2):
        """ e.g. P.faceWith("A","B") could return (polygon, "BDA") """
        for face in self.faces:
            poly, vtxs = face
            for a,b in window(vtxs, 2):
                if a==v1 and b==v2:
                    return face
            if vtxs[-1]==v1 and vtxs[0]==v2: return face
        print(v1,v2)
        crash


kelly_colors = ['F2F3F4', '222222', 'F3C300', '875692', 'F38400', 'A1CAF1',
                'BE0032', 'C2B280', '848482', '008856', 'E68FAC', '0067A5',
                'F99379', '604E97', 'F6A600', 'B3446C', 'DCD300', '882D17',
                '8DB600', '654522', 'E25822', '2B3D26']

### MAKE PLATONIC SOLIDS ###

def edgeDistance(P,v,w):
    """ How many edges to walk from v to w on P? """
    edges = {v:[] for poly,vtxs in P.faces for v in vtxs}
    for poly, vtxs in P.faces:
        for i,j in window(list(vtxs)+[vtxs[0]], 2):
            edges[i].append(j)
    horizon, seen = {v}, {v}
    steps = 0
    while True:
        if w in seen: return steps
        horizon = {j for i in horizon for j in edges[i] if j not in seen}
        for i in horizon: seen.add(i)
        steps += 1
    

h = 3**0.5/2
triangle = Polygon([(0,0), (0.5,h), (1,0)])
Tetrahedron = Polyhedron([(triangle,"ABD"),(triangle,"BAC"),
                          (triangle,"CDB"),(triangle,"ADC")],
                         {(v,w):1 for v,w in permutations("ABCD",2)})
Octahedron = Polyhedron(
    [(triangle,vtxs) for vtxs in
     ("ACD","ADE","EDF","DCF","ABC","AEB","CBF","BEF")])
Octahedron.vertexDistances = {(v,w):[0,1,2*h][edgeDistance(Octahedron,v,w)]
                              for v,w in permutations("ABCDEF",2)}

square = Polygon([(0,0),(0,1),(1,1),(1,0)])
Cube = Polyhedron(
    [(square,vtxs) for vtxs in
     ("ABDC","CDHG","DBFH","BAEF","ACGE","FEGH")])
Cube.vertexDistances = {(v,w):[0,1,2**0.5,5**0.5][edgeDistance(Cube,v,w)]
                        for v,w in permutations("ABCDEFGH",2)}

def cis(angle):
    return np.array([math.cos(angle),math.sin(angle)])
pentagon = Polygon([(0,0),cis(math.tau*3/10),
                    cis(math.tau*3/10)+cis(math.tau*1/10),
                    cis(math.tau*1/5)+(1,0),
                    (1,0)])
Dodecahedron = Polyhedron(
    [(pentagon,vtxs) for vtxs in
     ("ABGHF","CBAED","AFPQE","FHLRP","GJKLH","LKMSR","PRSTQ",
      "GBCIJ","JINMK","ICDON","NOTSM","ODEQT")])
dodecahedronDistances = [0, 1, phi := np.linalg.norm(cis(0)+cis(math.tau/5)),
                         np.linalg.norm(cis(0)+2*cis(math.tau/5)),
                         np.linalg.norm(cis(math.tau/5)+(phi+1)*cis(0)),
                         np.linalg.norm((2+phi)*cis(0)+cis(math.tau/5))]
Dodecahedron.vertexDistances = {
    (v,w):dodecahedronDistances[edgeDistance(Dodecahedron,v,w)]
    for v,w in permutations("ABCDEFGHIJKLMNOPQRST", 2)
}

Icosahedron = Polyhedron(
    [(triangle,vtxs) for vtxs in
     ("ADF","AFE","AEB","FDK","IFK","EFI","GEI","BEG","GIL","IKL",
      "ACD","ABC","CJD","JKD","HJC","BHC","BGH","GLH","HLJ","JLK")])
icosahedronDistances = [0,1,2*h,7**0.5]
Icosahedron.vertexDistances = {
    (v,w):icosahedronDistances[edgeDistance(Icosahedron,v,w)]
    for v,w in permutations("ABCDEFGHIJKL", 2)
}

### HELPER FUNCTIONS ###

def edges(basePoly, baseVtxs):
    # list of (vtx1, x_y1, vtx2, x_y2)
    x, y = basePoly.exterior.xy
    points = zip(x,y)
    labels = list(baseVtxs)+[baseVtxs[0]]
    for ((p1,l1),(p2,l2)) in window(zip(points, labels), 2):
        yield (l1,p1,l2,p2)

def vertices(basePoly, baseVtxs):
    # list of (vtx1, x_y1)
    x, y = basePoly.exterior.xy
    points = zip(x,y)
    labels = list(baseVtxs)
    for point, label in zip(points, labels): yield (label, point)

def SegmentTransform(segment1, segment2):
    def transformFromBase(P1, P2):
        """ Give a transform that maps (0,0) to P1 and (1,0) to P2. """
        Pd = (P2[0]-P1[0], P2[1]-P1[1])
        return np.array([[Pd[0],-Pd[1],P1[0]],
                         [Pd[1],Pd[0],P1[1]],
                         [0,0,1]])
    transform1 = np.linalg.inv(transformFromBase(*segment1))
    transform2 = transformFromBase(*segment2)
    transform = transform2 @ transform1
    return lambda shape: shapely.affinity.affine_transform(
        shape,
        (transform[0,0], transform[0,1],
         transform[1,0], transform[1,1],
         transform[0,2], transform[1,2])
        )

def norm(x_y):
    x,y=x_y
    return (x**2+y**2)**0.5


