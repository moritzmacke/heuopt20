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
        
    def reverse(self):
        c = self.fst
        while c != self.lst:
            c.prv,c.nxt = c.nxt, c.prv
            c = c.prv
        c.prv,c.nxt = c.nxt, c.prv
        self.lst = self.fst
        self.fst = c
        
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
    
    def __init__(self, inst):
        self.inst = inst
        
    def sol_obj(self, x):
        w = 0
        d = self.inst.weights
        for i in range(len(x) - 1):
            w += d[x[i]][x[i + 1]]
        w += d[x[-1]][x[0]]
        return w
          
    def construct(self, find_best_solution = False):
        pts = [PointInfo(i) for i in range(self.inst.n)]
        for (n1, n2, w) in self.inst.edges:
            pts[n1].ns_open.add(pts[n2])
            pts[n2].ns_open.add(pts[n1])
            
        es = sorted([p for p in pts], key=lambda x: -len(x.ns_open))
               
        best_sol = None
        best_obj = None
               
        max_path = None
        max_path_len = 0
        ctr = 0
        for e in es:
#            print("start from", e.idx)
            ps, path = self.try_construct_path_from(e.idx)
            if path.length > max_path_len:
                max_path = path
                max_path_len = path.length
            print("loop: {}, p: {}, len: {}, max: {}".format(ctr, e.idx, path.length, max_path_len))
            if path.length == self.inst.n and self.try_close_cycle(path):
                if find_best_solution:
                    x = path.tolist()
                    x_obj = self.sol_obj(x)
#                    print("obj", x_obj, best_obj)
                    if not best_sol or abs(x_obj) < abs(best_obj):
                        best_sol = x
                        best_obj = x_obj
                        print("best obj", best_obj)
                else:
                    break
            ctr += 1
            if ctr > 100: break; # searching all takes quite long...
        if not best_sol:
            x = max_path.tolist()
            x += [i.idx for i in ps if i.idx not in x]
            best_sol = x

        return best_sol
                
    def try_construct_path_from(self, start: int):
        ps = [PointInfo(i) for i in range(self.inst.n)]
        for (n1, n2, w) in self.inst.edges:
            ps[n1].ns_open.add(ps[n2])
            ps[n2].ns_open.add(ps[n1])
            
        path = Path()
        path.append(ps[start])
        self.greedy_path_extend(path)
        
        tried_moves = set()
        rev = False
        iteration = 0
        while path.length < self.inst.n:
#            print("iteration: {} path length: {}, trying to rotate".format(iteration, path.length))
            if not self.rotate_path(path, tried_moves):
                break
            self.greedy_path_extend(path)
            iteration += 1
        
#        path.pr()      
        
#        print("rem points", [p for p in ps if not p.visited ])
        
        return ps, path
    
    def greedy_path_extend(self, path):
        
        last = path.lst
        while path.length <= self.inst.n:
            last.visited = True
            for i in itertools.chain(last.ns_open, last.ns_path):
                if last in i.ns_open: 
                    i.ns_open.remove(last)
                i.ns_path.add(last)
            
#            print(self.ps)
            
            candidates = sorted([p for p in last.ns_open], key=lambda p: (len(p.ns_open),len(p.ns_path))) #neighbors, sort by degree, pick lowest degree

            if not candidates:
                break
            nxt = candidates[0]
            path.append(nxt)
            last = nxt
#            print("choose", last)
            
    def rotate_path(self, path, tried_moves):
        
#        if len(path.fst.ns_open) > len(path.lst.ns_open):
#            print("reversing...")
#            path.reverse()
        
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
            new_end = best_p.nxt
#            print("cur end: {} -> new end: {}".format(e, new_end))
            tried_moves.add((e.idx, new_end.idx))
#            print("tried moves: {}".format(tried_moves))
            path.rev_after(best_p)
            return True
        else:
            return False
            
    
    def try_close_cycle(self, path):
               
        tried = [False for i in range(self.inst.n)]
               
        s, e = path.fst, path.lst
        cycle = False
        while e != None:
            #with path len == n all neighbors will be in ns_path
            if s in e.ns_path: #is there edge to start?
#                print("found cycle!")
                cycle = True
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
        return cycle
                
