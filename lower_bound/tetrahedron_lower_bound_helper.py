
import numpy as np
from random import random, choice, shuffle
from itertools import product, permutations
from math import cos, sin, tau

import shapely
from shapely.geometry import Point, LinearRing, Polygon
from shapely.ops import cascaded_union

from scipy.optimize import minimize

h = 3**0.5/2

def baseSquare(sideLength):
    return Polygon([(0,0), (0,sideLength),
                    (sideLength,sideLength), (sideLength,0)])

def aimEpsilon(stage, pts, target):
    for P in permutations(range(len(pts))):
        if _bestEpsilon(stage, pts, P, target) >= target:
            return True
    return False

def _bestEpsilon(stage, pts, permutation, target):
    """ We're just optimizing in 2D x {24} now ok
    and the `permutation` arg takes care of the {24} part """
    # first of all, is a point just outside the stage?
    if any(stage.distance(Point(pt)) >= target for pt in pts):
        return max(stage.distance(Point(pt)) for pt in pts)
    def badness(xy, tryHard=True):
        distances = tetrahedronDistances(xy)
        distances = [distances[i] for i in permutation]
        rArea = realizable(stage, pts, distances, resolution=100)
        if rArea > 0: return rArea
        epsilon = badDistance(stage, pts, distances, resolution=100,
                              iterations=40 if tryHard else 10)
        return -epsilon
    firstTry = min(((x/10,y/10) for x in range(11) for y in range(11)
                   if x+y <= 10), key=lambda xy:badness(xy, tryHard=False))
    if badness(firstTry) >= 0: return 0
    return -badness(firstTry)
    # result = minimize(badness, firstTry, bounds=((0,1),(0,1)))
    # TODO: ensure x+y <= 1 as well (maybe by just mapping (x,y) to (1-x,1-y))
    # return result

def Diamond(xytheta):
    """ Returns locations of vertices A,B,C,D locked in the diamond configuration
    (two edge-adjacent triangles) """
    def rotateThenTranslate(pt, theta, dx, dy):
        x,y = pt
        return (x*cos(theta)-y*sin(theta)+dx, x*sin(theta)+y*cos(theta)+dy)
    dx, dy, theta = xytheta
    return [rotateThenTranslate(pt, theta, dx, dy)
            for pt in [(0,0.5), (0,-0.5), (h,0), (-h,0)]]

def realizable(stage, points, distanceLowerBounds,
               resolution = 100):
    """ Is there a point in `stage` with distances from `points`
    at least `distanceLowerBounds` ? """
    badZones = [circle(center, radius, resolution)
                for center, radius in zip(points, distanceLowerBounds)]
    badZone = cascaded_union(badZones)
    return (stage - badZone).area

def badDistance(stage, points, distanceLowerBounds,
                resolution = 100, iterations=40):
    """ How much must you adjust `points` so that there's a point in `stage`
    satisfying `distanceLowerBounds` ? """
    # If you can adjust all points by epsilon and win,
    # then shrinking buffers by epsilon makes stuff realizable.
    minEpsilon, maxEpsilon = 0, 1.0
    for i in range(iterations):
        midEpsilon = (minEpsilon + maxEpsilon)/2
        if realizable(stage, points,
                      [d-midEpsilon for d in distanceLowerBounds],
                      resolution):
            maxEpsilon = midEpsilon
        else:
            minEpsilon = midEpsilon
    return minEpsilon

def tetrahedronDistances(tetPt):
    """ tetPt is nonnegative (x,y) where x+y<=1, which we interpret
    as a point in a certain triangle occupying one-sixth of the area
    of one face of the tetrahedron.
    Returns distances along the surface to the four vertices. """
    vertices = np.array(((0.5, h), (0.75, h/2), (0.5, h/3)))
    xyz = np.array([tetPt[0], tetPt[1], 1-tetPt[0]-tetPt[1]])
    pt = np.sum(vertices * xyz.T[:,None], axis=0)
    tetVertices = np.array(((0, 0), (0.5, h), (1, 0), (1.5, h)))
    distances = np.linalg.norm(pt[None, :] - tetVertices, axis=1)
    return distances

def circle(center, radius, resolution):
    """ actually makes a polygon that's strictly inside the circle """
    return Point(center).buffer(radius, resolution=resolution)

def randomTetPt():
    """ returns nonnegative (x,y) where x+y<=1. """
    x, y = random(), random()
    if x + y > 1: x, y = 1-x, 1-y
    return x,y

