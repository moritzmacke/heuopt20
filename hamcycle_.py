import itertools

class PointInfo:
    """
    """
    def __init__(self, index: int):
        self.idx = index
        self.nxt = None
        self.prv = None
        self.ns_open = set()
        self.ns_path = set()
        self.visited = False
        
    def __repr__(self):
        return f"p{self.idx}#{self.visited}, ns_open={set(i.idx for i in self.ns_open)}, ns_path={set(i.idx for i in self.ns_path)}"
    
    def __hash__(self):
        return hash(self.idx)
       
class Path:
    """
    """
    def __init__(self):
        self.fst = None
        self.lst = None
        self.length = 0
        
    def append(self, e : PointInfo):
        if self.lst == None:
            self.fst = self.lst = e
            e.nxt = e.prv = None
        else:
            e.nxt = None
            self.lst.nxt = e
            e.prv = self.lst
            self.lst = e
        self.length += 1
        
    def rev_after(self, p):
        if p == self.lst:
            return
        c = self.lst
        while c != p:
            c.prv, c.nxt = c.nxt, c.prv
            c = c.nxt
        self.lst.prv = p
        c = p.nxt
        p.nxt = self.lst
        c.nxt = None
        self.lst = c
        
    def tolist(self):
        l = []
        c = self.fst
        while c != None:
            l.append(c.idx)
            c = c.nxt
        return l
        
    def pr(self):
        print(self.tolist())
        
    def rpr(self):
        l = []
        c = self.lst
        while c != None:
            l.append(c.idx)
            c = c.prv
        print(l)


class HamCycle:
    
    def __init__(self, n, edges):
        self.n = n
        self.edges = edges
        self.ps = [PointInfo(i) for i in range(n)]
        
    def construct(self):
        self.ps = [PointInfo(i) for i in range(self.n)]
        for (n1, n2, w) in self.edges:
            self.ps[n1].ns_open.add(self.ps[n2])
            self.ps[n2].ns_open.add(self.ps[n1])
            
        es = sorted([p for p in self.ps], key=lambda x: -len(x.ns_open))
        
#        print(es)
                            
        last = es[0]
        path = Path()
        path.append(last)
        tryrotate = True
        iteration = 0
        
        tried_moves = set()
        
        while tryrotate:
            self.greedy_path_extend(path)
            if path.length < self.n:
                print("iteration: {} path length: {}, trying to rotate".format(iteration, path.length))
                tryrotate = self.rotate_path(path, tried_moves)
            else:
                tryrotate = False
            iteration += 1
        
#        path.pr()
        
        if(path.length == self.n):
            self.try_close_cycle(path)
        
        x = path.tolist()
        x += [i.idx for i in self.ps if not i.visited]
        
        return x
        
        
    
    def greedy_path_extend(self, path):
        
        last = path.lst
        while path.length <= self.n:
            last.visited = True
            for i in itertools.chain(last.ns_open, last.ns_path):
                if last in i.ns_open: 
                    i.ns_open.remove(last)
                i.ns_path.add(last)
            
#            print(self.ps)
            nxt = None
            candidates = sorted([p for p in last.ns_open], key=lambda p: (len(p.ns_open),len(p.ns_path))) #neighbors, sort by degree, pick lowest degree

#            print("next candidates", candidates)
            for p in candidates: #does testing neighbors do anything? seems to almost always choose the first point?
#                print("try", p.idx)
                p.visited = True #for sake of neighbor test
                ok = True # what if last?
                for q in p.ns_open: #any neighbor unreachable? does this actually test it?
                    s = sum(1 for r in q.ns_open if not r.visited)
#                    print(q, "->", s)
                    if s == 0 and path.length < self.n-2:  #it's okay if we are at last two points?
                        ok = False
                        break
                p.visited = False #
                if ok:
                    nxt = p
                    break
            if nxt == None:
#                print("no viable neighbors")
                break
            path.append(nxt)
            last = nxt
#            print("choose", last)
            
    def rotate_path(self, path, tried_moves):
        
        e = path.lst
        best_p = None
        best_le = -1
#        print("end", e)
        for p in e.ns_path: #find neighbors in path
 #           print(p)
            new_end = p.nxt
            if new_end == e or (e.idx, new_end.idx) in tried_moves:
                continue
            le = len(new_end.ns_open)
#            print(new_end)
            if le > best_le:
                best_le = le
                best_p = p
        if best_p != None:
            print("cur end: {} -> new end: {}".format(e,best_p.nxt))
            tried_moves.add((e.idx, best_p.nxt.idx))
            print("tried moves", tried_moves)
            path.rev_after(best_p)
            return True
        else:
            return False
            
    
    def try_close_cycle(self, path):
               
        tried = [False for i in range(self.n)]
               
        s, e = path.fst, path.lst
            
        while e != None:
            #with path len == n all neighbors will be in ns_path
            if s in e.ns_path: #is there edge to start?
                print("found cycle!")
                break
                
#            print("cur end", e)
            tried[e.idx] = True
            best_p = None
            best_deg = 0
            
            for p in e.ns_path: #
                new_end = p.nxt
                if tried[new_end.idx]:
                    continue
                if s in new_end.ns_path: #edge to start?
                    best_p = p
                    break
                if len(new_end.ns_path) > best_deg: #otherwise chose point with most neighbors as new end
                    best_p = p
                    best_deg = len(new_end.ns_path)
            if best_p != None:
#                print("new end", best_p.nxt)
                path.rev_after(best_p)
            e = path.lst
                
    
if __name__ == '__main__':
    
    
    ps = [PointInfo(i) for i in range(10)]
    
    path = Path()
    for p in ps:
        path.append(p)
        
    path.rev_after(ps[4])
        
    path.pr()
    path.rpr()
    
    
