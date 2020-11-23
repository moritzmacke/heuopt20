import random

class PointInfo:
    """
    """
    def __init__(self, index: int):
        self.idx = index
        self.ns_open = set()
        self.ns_path = set()
        
    def __repr__(self):
        return f"p{self.idx}#{self.visited}, ns_open={set(i.idx for i in self.ns_open)}, ns_path={set(i.idx for i in self.ns_path)}"
    
    def __hash__(self):
        return hash(self.idx)
    
class PathFragment:
    
    def __init__(self, index: int):
        self.fid = index
        self.seq = [index]
        
    def fst(self):
        return self.seq[0]
    
    def lst(self):
        return self.seq[-1]
        
    def connectOverEdge(self, e, other):
        p, q = e[0], e[1]
        if p != self.fst() and p != self.lst():
            p, q = q, p
        
#        print("connect", self, other, e)
        m = PathFragment(min(self.fid, other.fid))

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

def edge_rank(e, ps):
    i1, i2 = e[0], e[1]
    p1, p2 = ps[i1], ps[i2]
    g1, g2 = len(p1.ns_open), len(p2.ns_open)
    return (g1, g2, e[2]) if g1 < g2 else (g2, g1, e[2])

class GreedyEdgeConst:
    
    def __init__(self, n, edges):
        self.n = n
        self.edges = edges
#        self.ps = [PointInfo(i) for i in range(n)]


    def construct(self, alpha):
        
        ps = [PointInfo(i) for i in range(self.n)]
        for (n1, n2, w) in self.edges:
            ps[n1].ns_open.add(ps[n2])
            ps[n2].ns_open.add(ps[n1])
        
        n = self.n
        edges = sorted(self.edges, key=lambda e: edge_rank(e,ps))
        setof = [PathFragment(i) for i in range(n)]
        seledges = [[] for i in range(n)]
        
        k = 0
        cur_obj = 0
        
        while edges:
            cand_edges = []
            for e in edges:
                p, q = e[0], e[1]
                if len(seledges[p]) > 1 or len(seledges[q]) > 1:
                    continue
                if setof[p].fid == setof[q].fid and k < n:
                    continue
                cand_edges.append(e)
            if not cand_edges:
                break
            cand_edges.sort(key=lambda e: edge_rank(e,ps))
#            print("filtered edges", cand_edges)
            sel_len = int(len(cand_edges)*alpha)
            rand_sel = random.randrange(max(1, sel_len))
#            print(rand_sel, len(cand_edges))
            e = cand_edges[rand_sel]
            p,q = e[0], e[1]
            seledges[p].append(q)
            seledges[q].append(p)
            sq,sp = setof[q], setof[p]
#            print(sp,sq)
            spq = sp.connectOverEdge(e, sq)
            for i in range(n):
                if setof[i].fid == sq.fid or setof[i].fid == sp.fid:
                    setof[i] = spq
#            print(["{}:{}".format(i,e) for (i,e) in enumerate(seledges)])
#            print(setof)
            edges = cand_edges
            cur_obj += e[2]
            k = k + 1
            
        x = []
        for s in set(setof):
#            print(s)
            x += s.seq

        return x
    
