import random
from sortedcontainers import SortedList
import itertools

class PointInfo:
    """
    """
    def __init__(self, index: int):
        self.idx = index
        self.adj_edges = set()
        
    def __repr__(self):
        return f"p{self.idx}, adj_edges={set((e.p.idx, e.q.idx) for e in self.adj_edges)}"
    
    def __hash__(self):
        return hash(self.idx)
    
   
class PathFragment:
    
    def __init__(self, index: int):
        self.sid = index
        self.seq = [index]
        
    def fst(self):
        return self.seq[0]
    
    def lst(self):
        return self.seq[-1]
        
    def connectOverEdge(self, e, other):
        p, q = e.p.idx, e.q.idx
        if p != self.fst() and p != self.lst():
            p, q = q, p
        
#        print("connect", self, other, e)
        m = PathFragment(min(self.sid, other.sid))

        if p == self.fst():
            if q == other.fst():
                m.seq = self.seq[::-1] + other.seq
            elif q == other.lst():
                m.seq = other.seq + self.seq
            else:
                raise ValueError()
        elif p == self.lst():
            if q == other.fst():
                m.seq = self.seq + other.seq
            elif q == other.lst():
                m.seq = self.seq + other.seq[::-1]
            else:
                raise ValueError()
        else:
            raise ValueError()
        
#        print("->", m)
        
        return m

    def __repr__(self):
        return f"{self.seq}"

class Edge:
    
    def __init__(self, p1, p2, w):
      self.p = p1
      self.q = p2
      self.w = w
      self.score = None
      
    def compute_score(self):
        degp, degq = len(self.p.adj_edges), len(self.q.adj_edges)
        self.score = (degp, degq, self.w) if degp < degq else (degq, degp, self.w)
        
    def __repr__(self):
        return f"({self.p.idx}--{self.q.idx} {self.score})"
    
    def __hash__(self):
        return hash(self.p) + hash(self.q)
    
    def __eq__(self, other):
        return self.p.idx == other.p.idx and self.q.idx == other.q.idx

def set_super(sets, idx, to):
    s = sets[idx]
    if s.sid == idx:
        sets[idx] = to
    else:
        ss = set_super(sets, s.sid, to)
        s.sid = to.sid

def get_set(sets, idx):
    s = sets[idx]
    if s.sid == idx:
        return s
    else:
        ss = get_set(sets, s.sid)
        s.sid = ss.sid
        return ss

class GreedyEdgeConst:
    """ Choose edges that have the least neighbors greedily
        There might still be some bugs in here...
        Also slightly random even if alpha=0 for some reason?
    """
    
    def __init__(self, n, edges):
        self.n = n
        self.orig_edges = edges


    def construct(self, alpha):
        
        ps = [PointInfo(i) for i in range(self.n)]
        edges = [Edge(ps[pi],ps[qi],w) for (pi,qi,w) in self.orig_edges ]
        
        for e in edges:
            pi, qi = e.p.idx, e.q.idx
            e.p.adj_edges.add(e)
            e.q.adj_edges.add(e)
        
        for e in edges:
            e.compute_score()
        
        n = self.n
        edges = SortedList(edges, key=lambda e: e.score)
        setof = [PathFragment(i) for i in range(n)]
        seledges = [[] for i in range(n)]
        
        k = 0
        cur_obj = 0
        
        while edges :
            sel_len = int(len(edges)*alpha)
            rand_sel = random.randrange(max(1, sel_len))
            e = edges[rand_sel]
            
            #update chosen edges, sets
            p,q = e.p, e.q
            seledges[p.idx].append(q.idx)
            seledges[q.idx].append(p.idx)
            sq,sp = get_set(setof, q.idx), get_set(setof, p.idx)
            
#            assert(sq.sid != sp.sid)
            
            spq = sp.connectOverEdge(e, sq)
            set_super(setof, sp.sid, spq)
            set_super(setof, sq.sid, spq)

            #update edges
            #all edges that share a point with current edge plus edges neighboring these again
            
            p.adj_edges.remove(e)
            q.adj_edges.remove(e)
            edges.remove(e)
                      
            #need to look at edges incident to end points of merged fragments
            affected_points = set([sp.fst(),sp.lst(), sq.fst(), sq.lst()])
#            print("affected points", affected_points)
            adjacient_edges = set(itertools.chain.from_iterable(ps[i].adj_edges for i in affected_points))
#            print("affected edges", adjacient_edges)

            #edges that are kept, but must be updated
            affected_edges = set()
            
            for e in adjacient_edges:
                #filter out now illegal edges
                #is k<n-1 correct? want to be able to close cycle on last iteration
                if len(seledges[e.p.idx]) > 1 or len(seledges[e.q.idx]) > 1 or (get_set(setof, e.p.idx).sid == get_set(setof, e.q.idx).sid and k < n-1):
#                    print("remove", e)
                    e.p.adj_edges.remove(e)
                    e.q.adj_edges.remove(e)
                    edges.remove(e)
                    affected_points.add(e.p.idx)
                    affected_points.add(e.q.idx)
                else:
                    affected_edges.add(e)
                    
#            print("aff points", affected_points)
                    
            affected_edges = affected_edges.union(set(itertools.chain.from_iterable(ps[i].adj_edges for i in affected_points)))
            for e in affected_edges:
                edges.remove(e)
                e.compute_score()
                edges.add(e)
#                print("updated",e)


            cur_obj += e.w
            k = k + 1
            
            
        x = []
        for s in set([get_set(setof, i) for i in range(n)]):
            x += s.seq
            
        return x
    
