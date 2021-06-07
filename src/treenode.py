class Node:
    def __init__(self):
        self.label = ""
        self.length = 0.0
        self.time_length = 0.0
        self.parent = None
        self.children = []
        self.data = {}
        self.istip = False
        self.height = 0
        self.note = ""
        self.number = 0  # this is my mod of Stephen's
    
    def add_child(self,child):
        #make sure that the child is not already in there
        assert child not in self.children
        self.children.append(child)
        child.parent = self
    
    def remove_child(self,child):
        #make sure that the child is in there
        assert child in self.children
        self.children.remove(child)
        child.parent = None
    
    def leaves(self,v=None):
        if v == None:
            v = []
        if len(self.children) == 0:
            v.append(self)
        else:
            for child in self.children:
                child.leaves(v)
        return v

    def leaves_fancy(self):
        return [n for n in self.iternodes() if n.istip ]

    def lvsnms(self):
        return [n.label for n in self.iternodes() if len(n.children) == 0 ]

    def iternodes(self,order="preorder"):
        if order.lower() == "preorder":
            yield self
        for child in self.children:
            for d in child.iternodes(order):
                yield d
        if order.lower() == "postorder":
            yield self

    def is_rooted(self):
        if len(self.children) == 2:
            return True
        elif len(self.children) == 3:
            return False

    def number_tree(self):
        c = 0
        for n in self.iternodes(order="postorder"):
            if n.parent is None:
                continue
            else:
                n.number = c
                c += 1

    def prune(self):
        p = self.parent
        if p != None:
            p.remove_child(self)
        return p
    
    def get_newick_repr_paint(self):
        ret = ""
        painted_children = []
        for i in self.children:
            if "paint" in i.data:
                painted_children.append(i)
        for i in range(len(painted_children)):
            if i == 0:
                ret += "("
            ret += painted_children[i].get_newick_repr_paint()
            if i == len(painted_children)-1:
                ret += ")"
            else:
                ret += ","
        if self.label != None and "paint" in self.data:
            ret += self.label
        return ret

    def get_newick_repr(self, showbl=False, shownum=False):
        ret = ""
        for i in range(len(self.children)):
            if i == 0:
                ret += "("
            ret += self.children[i].get_newick_repr(showbl, shownum)
            if i == len(self.children)-1:
                ret += ")"
            else:
                ret += ","
        if self.label != None:
            ret += self.label
            #ret += str(self.label)
        if shownum == True:
            if self.istip:
                ret += "_" + str(self.number)
            else:
                ret += str(self.number)
        if self.note != "":
            ret += "["+self.note+"]"
        if showbl == True:
            if self.parent != None:
                ret += ":" + format(self.length,'.12f')
        return ret

    def get_newick_repr_data(self,datan,showbl=False):
        ret = ""
        for i in range(len(self.children)):
            if i == 0:
                ret += "("
            ret += self.children[i].get_newick_repr_data(datan,showbl)
            if i == len(self.children)-1:
                ret += ")"
            else:
                ret += ","
        if datan in self.data:
            ret += self.data[datan]
        elif self.label != None:
            ret += self.label
        if self.note != "":
            ret += "["+self.note+"]"
        if showbl == True:
            ret += ":" + str(self.length)
        return ret

    def __repr__(self):
        return self.get_newick_repr(False, False)
