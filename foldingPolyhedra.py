
import shapely
import numpy
from shapely.geometry import Point, LinearRing, Polygon
from shapely.ops import cascaded_union
from scipy.optimize import minimize
import random

from nets import (GeneralNet, Tetrahedron, Octahedron, Cube, Dodecahedron,
                  Icosahedron, showPolys)

def save(obj, filename):
    with open("C:/Users/linyk/Desktop/"+filename+".txt", "w") as f:
        f.write(repr(obj))

def PolygonRepr(P):
    return f"Polygon({list(zip(*P.exterior.coords.xy))[:-1]})"

def oneNet(polyhedron, paper, radius=8):
    # Step 1: move `paper` to be centered at the origin
    paper = PutPolygon(paper, (-paper.centroid.x, -paper.centroid.y), (1,0))
    def fuzz(r): return random.random()*2*r-r
    polygons = GeneralNet(polyhedron, radius=radius, checkFaithful=True)
    if not len(set(b for a,b in polygons)) == len(polyhedron.faces):
        print("Not all faces appeared in the unfolded net. Use a larger "\
              "radius or (shudder) assemble a net by hand.")
    # First, try like 25 random positions for the paper to see which is best
    testPapers = [((fuzz(radius/2),fuzz(radius/2)),(fuzz(1),fuzz(1)))
                  for i in range(25)] # 25 chosen empirically for speed
    o,x = min(testPapers,
              key=lambda s:ScaleFactor(paper,s[0],s[1],polygons,
                                       precision=1e-2))
    #showPolys(polygons+[(PutPolygon(paper,o,x),"")], Tetrahedron.colorDict)##
    #global Y; Y = o,x,polygons
    x = ScaleToLength(x,ScaleFactor(paper,o,x,polygons,precision=1e-2))
    #showPolys(polygons+[(PutPolygon(paper,o,x),"")], Tetrahedron.colorDict)##
    # Now plug it into scipy for gradient descent
    o,x = TryMinimizingPolygon(paper, o, x, polygons)
    #showPolys(polygons+[(PutPolygon(paper,o,x),"")], Tetrahedron.colorDict)##
    score = ScaleFactor(paper, o, x, polygons)
    x = ScaleToLength(x, score)
    #showPolys(polygons+[(PutPolygon(paper,o,x),"")], Tetrahedron.colorDict)##
    return (polygons, PutPolygon(paper,o,x), score)

def TryMinimizingPolygon(paper, origin, xaxis, polys):
    """ polys are a list of (Polygon, labeled vertex string like "ABD").
    paper is a polygon.
    origin and xaxis are tuples like (1,2) and (4,5) defining a square.
    Gradient descents to find a smaller square """
    epsilon = 1e-10
    baseDict = {l:p for p,l in polys} # label string: polygon
    areaToCapture = sum(p.area for p in baseDict.values())
    def Nonneg(pt):
        return Score(PutPolygon(paper,pt[:2],pt[2:]),polys)\
               - areaToCapture + epsilon
    result = minimize(lambda pt:(pt[2]**2+pt[3]**2)**0.5,
                      origin+xaxis,
                      constraints=[{"type":"ineq",
                                    "fun":Nonneg}])
    origin, xaxis = result.x[:2], result.x[2:]
    return origin, xaxis

def Score(placedPaper, polygons):
    """ Total area of polygons covered by the paper (moved to (origin, xaxis))?
    (optionally, scale the square first.) """
    score = 0
    baseDict = {l:p for p,l in polygons}
    D = Coverings(placedPaper, polygons, baseDict)
    for label, polys in D.items():
        score += cascaded_union(polys).area
    return score

def Coverings(placedPaper, polygons, baseDict):
    D = {label:[] for polygon,label in polygons}
    for polygon, label in polygons:
        I = polygon.intersection(placedPaper)
        nuqu = transformToBase(*vertices(polygon)[:2])(I)
        if not nuqu.is_empty: D[label].append(nuqu)
    return D

def ScaleFactor(paper, origin, xaxis, polys, epsilon=1e-10, precision=1e-8):
    """  """
    baseDict = {l:p for p,l in polys} # label string: polygon
    areaToCapture = sum(p.area for p in baseDict.values())
    def okFn(origin, xaxis):
        #print(xaxis, Score(origin, xaxis, polys))
        return Score(PutPolygon(paper, origin, xaxis), polys)\
               + epsilon >= areaToCapture
    return _ScaleFactor(origin, xaxis, okFn, epsilon=precision)

def _ScaleFactor(origin, xaxis, okFn, epsilon=1e-8):
    def ok(scale): return okFn(origin, ScaleToLength(xaxis,scale))
    if not ok(50): return 1e20
    lo, hi = 0.0, 50.0
    while hi - lo > epsilon:
        mid = (hi+lo)/2
        if ok(mid):
            hi = mid
        else:
            lo = mid
    return hi

def ScaleToLength(vector, length):
    currentLength = (vector[0]**2 + vector[1]**2)**0.5
    return (vector[0]/currentLength*length, vector[1]/currentLength*length)

'''
def Scale(polygon, scaleFactor):
    """ scales the polygon up around its centroid. """
    cx, cy = polygon.centroid.x, polygon.centroid.y
    def scaled(point): return ((point[0]-cx)*scaleFactor + cx,
                               (point[1]-cy)*scaleFactor + cy)
    return PointMap(scaled, polygon)'''

def PutPolygon(polygon,o,x):
    """ move the polygon mapping (0,0) -> o and (1,0) -> x """
    return PointMap(
        lambda p: (o[0]+p[0]*x[0]-p[1]*x[1], o[1]+p[0]*x[1]+p[1]*x[0]),
        polygon
        )

def PointMap(f, P):
    return Polygon([f(v) for v in vertices(P)])

def vertices(P):
    return list(zip(*P.boundary.coords.xy))

def transformToBase(P1, P2):
    """ Give a transform that maps P1 to (0,0) and P2 to (1,0). """
    Pd = (P2[0]-P1[0], P2[1]-P1[1])
    Minv = numpy.array([[Pd[0],-Pd[1],P1[0]],
                        [Pd[1],Pd[0],P1[1]],
                        [0,0,1]])
    transform = numpy.linalg.inv(Minv)
    return lambda shape: shapely.affinity.affine_transform(
        shape,
        (transform[0,0], transform[0,1],
         transform[1,0], transform[1,1],
         transform[0,2], transform[1,2])
        )

def add(*pts):
    return [*map(sum, zip(*pts))]
