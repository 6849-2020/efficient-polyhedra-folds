requirements: scipy, numpy, matplotlib, shapely

Let's document by example:

# What is the smallest 1 x 2 rectangle that perfectly folds around the surface of a tetrahedron?

```
from foldingPolyhedra import Tetrahedron, Polygon, oneNet, display
rectangle = Polygon([(0,0), (0,1), (2,1), (2,0)]) # give points in clockwise order
```

The following function call either crashes (only occasionally, and I can't control this -- shapely gives errors sometimes) or outputs a valid folding. It will take several seconds.

`netPolygons, paperPlacement, scaleFactor = oneNet(Tetrahedron, rectangle, radius=5)`

Now we can visualize:

`display(netPolygons, paperPlacement, polyhedron=Tetrahedron)`

By repeatedly running the function and keeping track of the result with the smallest scaleFactor, we can answer our question.

# Defining new polyhedra

```
import itertools
from nets import Polyhedron
triangle = Polygon([(0,0), (0.5,3\*\*0.5/2), (1,0)])
Tetrahedron = Polyhedron([(triangle,"ABD"),(triangle,"BAC"),    # give points in clockwise order!!!
                          (triangle,"CDB"),(triangle,"ADC")],
                          {(v,w):1 for v,w in itertools.permutations("ABCD",2)})
```

The final dict argument is a map from vertex pairs (v,w) to the distance between v and w on the surface of the polyhedron. This is the trickiest part of making your own polyhedron. The code uses this argument to make sure the polyhedron's nets actually fold properly onto it, without having to stretch or tear the paper in between. If you don't care about this, you can just use {(v,w):0 for v,w in itertools.permutations(<your vertex list>,2)}.
  
Email linyks@gmail.com with questions/bugs/etc.
