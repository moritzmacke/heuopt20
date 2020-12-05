import itertools
import time

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
    def __init__(self, dists):
        self.fst = None
        self.lst = None
        self.length = 0
        self.weight = 0
        self.d = dists
        
    def append(self, e : PointInfo):
        if self.lst == None:
            self.fst = self.lst = e
            e.nxt = e.prv = None
        else:
            self.weight += self.d[self.lst.idx][e.idx]
            e.nxt = None
            self.lst.nxt = e
            e.prv = self.lst
            self.lst = e
        self.length += 1
        
    def rev_after(self, p):
        if p == self.lst:
            return
        self.weight -= self.d[p.idx][p.nxt.idx]
        c = self.lst       
        while c != p:
            c.prv, c.nxt = c.nxt, c.prv
            c = c.nxt
        self.lst.prv = p
        c = p.nxt
        p.nxt = self.lst
        c.nxt = None
        self.lst = c
        self.weight += self.d[p.idx][p.nxt.idx]
        
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
          
    def construct(self, find_best_solution = False, max_time = 10):
        """
        :param find_best_solution:
        :param max_time: maximum time in second to try for solution
        """
        pts = [PointInfo(i) for i in range(self.inst.n)]
        for (n1, n2, w) in self.inst.edges:
            pts[n1].ns_open.add(pts[n2])
            pts[n2].ns_open.add(pts[n1])
            
        es = sorted([p for p in pts], key=lambda x: -len(x.ns_open))
               
        best_sol = None
        best_obj = None
        max_path = None
        max_path_len = 0
        t_start = time.process_time()
        ctr = 0
        for e in es:
            ps, path = self.try_construct_path_from(e.idx)
            if path.length > max_path_len:
                max_path = path
                max_path_len = path.length
            print("loop: {}, p: {}, len: {}, max: {}".format(ctr+1, e.idx, path.length, max_path_len))
            if path.length == self.inst.n and self.try_close_cycle(path):
                if find_best_solution:
                    obj_val = path.weight + self.inst.weights[path.lst.idx][path.fst.idx]
                    if not best_sol or abs(obj_val) < abs(best_obj):
                        best_sol = path.tolist()
                        best_obj = obj_val
#                        assert(best_obj == self.sol_obj(best_sol))
                        print("best obj", best_obj)
                else:
                    break
            ctr += 1
            elapsed = time.process_time() - t_start
            if elapsed > max_time: break
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
            
        path = Path(self.inst.weights)
        path.append(ps[start])
        self.greedy_path_extend(path)
        
        tried_moves = set()
        while path.length < self.inst.n:
            if not self.rotate_path(path, tried_moves):
                break
            self.greedy_path_extend(path)
            
        return ps, path
    
    def greedy_path_extend(self, path):
        last = path.lst
        while path.length <= self.inst.n:
            last.visited = True
            for i in itertools.chain(last.ns_open, last.ns_path):
                i.ns_open.discard(last)
                i.ns_path.add(last)
            candidates = sorted(last.ns_open, key=lambda p: (len(p.ns_open),len(p.ns_path))) #neighbors, sort by degree, pick lowest degree
            if not candidates:
                break
            last = candidates[0]
            path.append(last)

            
    def rotate_path(self, path, tried_moves):
        e = path.lst
        best_p = None
        best_le = -1
        for p in e.ns_path: #find neighbors in path
            new_end = p.nxt
            if new_end == e or (e.idx, new_end.idx) in tried_moves:
                continue
            le = len(new_end.ns_open)
            if le > best_le:
                best_le = le
                best_p = p
        if best_p != None:
            new_end = best_p.nxt
            tried_moves.add((e.idx, new_end.idx))
            path.rev_after(best_p)
            return True
        else:
            return False
            
    
    def try_close_cycle(self, path):
               
        tried = [False for i in range(self.inst.n)]
               
        s, e = path.fst, path.lst
        cycle = False
        while not tried[e.idx]:
            #with path len == n all neighbors will be in ns_path
            if s in e.ns_path: #is there edge to start?
#                print("found cycle!")
                cycle = True
                break
                
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
                path.rev_after(best_p)
            e = path.lst
        return cycle
                
