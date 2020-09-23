"""
# =============================================================================
# Implements pathfinding behavior in order to generate basic topologies for
# stress testing. Algorithm used is A*, but adapted for the grid style used in
# the program (square or cube elements, nodes in corners, boundary conditions
# placed on nodes).
#
# Author: Tarcísio L. de Oliveira
# Copyright (C) 2020, Tarcísio L. de Oliveira
# =============================================================================
"""

from queue import PriorityQueue
import heapq
import itertools

import numpy as np

__all__ = ['split', 'get_path', 'active_points', 'passive_points']

def split(array, nelx, nely, nelz):
    """
    Splits an array of boundary conditions into an array of collections of
    elements. Boundary conditions that are more than one node in size are
    grouped together. From the nodes, the function returns the neighboring
    elements inside the array.
    """
    if len(array) == 0:
        return []

    array.sort()
    
    nlist = []
    tmp = _get_elem(array[0], nelx, nely, nelz)
    for i in range(1, len(array)):
        if _neighbors_node(array[i-1], array[i], nelx, nely, nelz):
            tmp = tmp.union(_get_elem(array[i], nelx, nely, nelz))
        else:
            nlist.append(list(tmp))
            tmp = _get_elem(array[i], nelx, nely, nelz)

    nlist.append(list(tmp))
    return nlist

def active_points(array, nelx, nely, nelz):
    """
    Groups an array of elements into an array of collections of element points.
    Neighboring nodes are grouped together. Used for passive and active
    elements.
    """
    if len(array) == 0:
        return []

    array.sort()

    elist = []
    tmp = [elem2point(array[0], nelx, nely, nelz)]
    for i in range(1, len(array)):
        if _neighbors_elem(array[i-1], array[i], nelx, nely, nelz):
            tmp.append(elem2point(array[i], nelx, nely, nelz))
        else:
            elist.append(list(tmp))
            tmp = [elem2point(array[i], nelx, nely, nelz)]

    elist.append(list(tmp))
    return elist

def passive_points(array, nelx, nely, nelz):
    """
    Converts element numbers into points.
    """
    if len(array) == 0:
        return []

    elist = []
    for i in range(len(array)):
        elist.append(elem2point(array[i], nelx, nely, nelz))

    return elist

def dist_elem(e1, e2, nelx, nely, nelz):
    """
    Returns the distance between two elements.
    """
    p1 = elem2point(e1, nelx, nely, nelz)
    p2 = elem2point(e2, nelx, nely, nelz)

    return np.linalg.norm(np.asarray(p1) - np.asarray(p2))

def dist_node(n1, n2, nelx, nely, nelz):
    """
    Returns the distance between two nodes. Used to define which active
    elements are associated to which loads.
    """
    p1 = node2point(n1, nelx, nely, nelz)
    p2 = node2point(n2, nelx, nely, nelz)

    return np.linalg.norm(np.asarray(p1) - np.asarray(p2))

def dist(p1, p2):
    """
    Returns the distance between two points.
    """

    return np.linalg.norm(np.asarray(p1) - np.asarray(p2))

def shortest_dist(a1, p1):
    """
    Returns the smallest distance between a point (p1) and an array of points
    (p1).
    """

    d = dist(p1, a1[0])
    for i in range(1, len(a1)):
        d2 = dist(p1, a1[i])
        if d2 < d:
            d = d2

    return d

def elem2point(e1, nelx, nely, nelz):
    """
    Converts an element number to a point.
    """
    if nelz == 0:
        x1 = int(np.floor(e1 / nely))
        y1 = int(e1 % nely)

        p1 = (y1, x1)
    else:
        plane = nelx*nely

        z1 = int(np.floor(e1 / plane))
        rest1 = int(e1 % plane)

        x1 = int(np.floor(rest1 / nely))
        y1 = int(rest1 % nely)

        p1 = (z1, y1, x1)

    return p1

def node2point(n, nelx, nely, nelz):
    """
    Converts node's DOF data to an element point (equivalent to the element to
    its lower right).
    """
    if nelz == 0:
        e = int(np.floor(n/2))
        return elem2point(e, nelx + 1, nely + 1, 0)
    else:
        e = int(np.floor(n/3))
        return elem2point(e, nelx + 1, nely + 1, nelz + 1)


def point2elem(x, y, z, nelx, nely):
    """
    Converts a point to an element number.
    """
    return y + nely*x + nelx*nely*z

def get_path(mesh, supports, active, passive, load):
    """
    Gets the optimum path from a load to every support, passing through the
    active elements and avoiding passive elements.
    """
    class PathNode:
        def __init__(self, elem, prev):
            self.elem = elem
            self.prev = prev

        def __lt__(self, n):
            return self.elem < n.elem

        def __str__(self):
            return "Current: "+str(self.elem)+"; Previous: "+str(self.prev)
        def __repr__(self):
            return "Current: "+str(self.elem)+"; Previous: "+str(self.prev)

    if len(mesh.shape) == 2:
        nely, nelx = mesh.shape
        nelz = 0
    else:
        nelz, nely, nelx = mesh.shape

    queue_start = []

    if len(active) > 0:
        queue_a = []

        for l in load:
            if l not in passive:
                heapq.heappush(queue_a, (shortest_dist(active[0], l), PathNode(l, None)))

        for i in range(len(active)):
            a = active[i]
            for p in a:
                mesh[p] = 1
            e = heapq.heappop(queue_a)
            while e[1].elem not in a:
                neighbors = _get_elem_neighbors(e[1].elem, nelx, nely, nelz)
                for n in neighbors:
                    if n not in passive:
                        heapq.heappush(queue_a, (shortest_dist(a, n), PathNode(n, e)))
                e = heapq.heappop(queue_a)

            queue_a = []
            if i < len(active)-1:
                heapq.heappush(queue_a, (shortest_dist(active[i+1], e[1].elem), e[1]))

        queue_start.append((shortest_dist(supports[0], n), e[1]))

    else:
        for l in load:
            if l not in passive:
                queue_start.append((shortest_dist(supports[0], l), PathNode(l, None)))

    for s in supports:
        queue_s = queue_start.copy()
        heapq.heapify(queue_s)

        e = heapq.heappop(queue_s)
        while e[1].elem not in s:
            neighbors = _get_elem_neighbors(e[1].elem, nelx, nely, nelz)
            for n in neighbors:
                if n not in passive:
                    heapq.heappush(queue_s, (shortest_dist(s, n), PathNode(n, e)))
            e = heapq.heappop(queue_s)

        while e != None:
            mesh[e[1].elem] = 1
            e = e[1].prev

    return mesh

def _neighbors_node(n1, n2, nelx, nely, nelz):
    """
    Checks if two nodes are neighboring.
    """
    if nelz == 0:
        n1 = int(np.floor(n1/2)) + 1
        n2 = int(np.floor(n2/2)) + 1
    else:
        n1 = int(np.floor(n1/3)) + 1
        n2 = int(np.floor(n2/3)) + 1
    diff = abs(n2 - n1)
    plane = (nelx + 1)*(nely + 1)
    if nelz == 0:
        return diff == 0 or \
               diff == 1 or \
               diff == nely or \
               diff == nely + 1 or \
               diff == nely - 1
    else:
        return diff == 0 or \
               diff == 1 or \
               diff == nely or \
               diff == nely + 1 or \
               diff == nely - 1 or \
               diff == plane or \
               diff == plane + 1 or \
               diff == plane - 1

def _neighbors_elem(e1, e2, nelx, nely, nelz):
    """
    Checks if two elements are neighboring.
    """
    diff = abs(e1 - e2)
    if nelz == 0:
        return diff == 0 or \
               diff == 1 or \
               diff == nely
    else:
        plane = nely*nelx
        return diff == 0 or \
               diff == 1 or \
               diff == nely or \
               diff == plane


def _get_elem(n, nelx, nely, nelz):
    """
    Transforms a node into a set of element points.
    """
    eset = {}
    if nelz == 0:
        e = [0]*4
        e[0] = node2point(n, nelx, nely, nelz) 
        e[1] = (min(e[0][0], nely-1), max(e[0][1]-1,    0))
        e[2] = (max(e[0][0]-1,    0), min(e[0][1], nelx-1))
        e[3] = (max(e[0][0]-1,    0), max(e[0][1]-1,    0))
        e[0] = (min(e[0][0], nely-1), min(e[0][1], nelx-1))

        eset = set(e)
    else:
        e = [0]*8
        e[0] = node2point(n, nelx, nely, nelz) 
        e[1] = (min(e[0][0], nelz-1), min(e[0][1], nely-1), max(e[0][2]-1,    0))
        e[2] = (min(e[0][0], nelz-1), max(e[0][1]-1,    0), min(e[0][2], nelx-1))
        e[3] = (min(e[0][0], nelz-1), max(e[0][1]-1,    0), max(e[0][2]-1,    0))
        e[4] = (max(e[0][0]-1,    0), min(e[0][1], nely-1), min(e[0][2], nelx-1))
        e[5] = (max(e[0][0]-1,    0), min(e[0][1], nely-1), max(e[0][2]-1,    0))
        e[6] = (max(e[0][0]-1,    0), max(e[0][1]-1,    0), min(e[0][2], nelx-1))
        e[7] = (max(e[0][0]-1,    0), max(e[0][1]-1,    0), max(e[0][2]-1,    0))
        e[0] = (min(e[0][0], nelz-1), min(e[0][1], nely-1), min(e[0][2], nelx-1))

        eset = set(e)

    return eset

def _get_elem_neighbors(p, nelx, nely, nelz):
    """
    Gets the elements which surround the element in point p1 (no diagonals).
    """

    eset = {}
    if nelz == 0:
        e1 = (max(p[0]-1, 0), p[1])
        e2 = (min(p[0]+1, nely-1), p[1])
        e3 = (p[0], max(p[1]-1, 0))
        e4 = (p[0], min(p[1]+1, nelx-1))

        eset = {e1, e2, e3, e4}
    else:
        e1 = (max(p[0]-1, 0), p[1], p[2])
        e2 = (min(p[0]+1, nelz-1), p[1], p[2])
        e3 = (p[0], max(p[1]-1, 0), p[2])
        e4 = (p[0], min(p[1]+1, nely-1), p[2])
        e5 = (p[0], p[1], max(p[2]-1, 0))
        e6 = (p[0], p[1], min(p[2]+1, nelx-1))

        eset = {e1, e2, e3, e4, e5, e6}

    return list(eset)
