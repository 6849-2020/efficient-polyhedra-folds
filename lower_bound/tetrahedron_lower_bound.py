

from tetrahedron_lower_bound_helper import (baseSquare, aimEpsilon, Diamond,
                                    realizable, badDistance,
                                    tetrahedronDistances, randomTetPt)

from math import tau
from itertools import product
from collections import namedtuple

# with side length 1.3, subdividing the whole thing with factor 16 works!

def ProveOnePointFiveFive():
    """ This takes a while. It prints partial progress accounts. Be patient. """
    x = 1.55
    stage = baseSquare(x)
    # we'll search configuration space box from (0,0,0) to (tau/2, tau/2, tau/2)
    # this is fine since tau/2 > x
    configurationSpace = Cubelet((tau/4,tau/4,tau/4),tau/4)
    unsolvedCubelets = subdivide(configurationSpace, 16)
    i = 0
    while unsolvedCubelets:
        print(len(unsolvedCubelets), "cubelets still remain")
        stillUnsolved = []
        for c in unsolvedCubelets:
            i += 1
            if i%200==0: print(i,end=" cubelets examined in total")
            try:
                assert aimEpsilon(stage, Diamond(c.center),
                                  c.boxRadius*diamondScaleFactor)
            except:
                for subcubelet in subdivide(c, 2):
                    stillUnsolved.append(subcubelet)
        unsolvedCubelets = stillUnsolved
    return stillUnsolved
        
    

h = (3**0.5) / 2
diamondScaleFactor = h + 2**0.5

Proof = namedtuple("Proof", "x y permutation")
Cubelet = namedtuple("Cubelet", "center boxRadius")

def DiamondProves(sidelength, cubelet, proof):
    points = Diamond(cubelet.center)
    epsilon = cubelet.boxRadius * diamondScaleFactor
    distanceBounds = tetrahedronDistances((proof.x, proof.y))
    distanceBounds = [distanceBounds[i]-epsilon for i in permutation]
    stage = baseSquare(sidelength)
    return not realizable(stage, points, distanceBounds, resolution=100)
    

def DiamondProof(sidelength):
    # this means that if we get an epsilon of 0.1, then params can each
    # change by at most 0.1/scaleFactor
    paramDimension = max(tau/4, sidelength)
    cubelet = Cubelet((paramDimension/2,paramDimension/2,paramDimension/2),
                      paramDimension/2)

def subdivide(cubelet, factor=2):
    coords, boxRadius = cubelet.center, cubelet.boxRadius
    newRadius = boxRadius/factor
    subCordOptions = [[x-boxRadius+i*newRadius for i in range(1,2*factor,2)]
                      for x in coords]
    return [Cubelet(c,newRadius) for c in product(*subCordOptions)]

