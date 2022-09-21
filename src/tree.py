class TranscriptsTreeNode:
    def __init__(self, kwargs, parent=None, divider_exon=-1):
        self.kwargs = kwargs
        self.parent = parent
        self.divider_exon = divider_exon
        self.left_child = None
        self.right_child = None
        self.node_id = None
        self.df = None
        self.res = {}
        self.tissue_res = {}

    def set_children(self, left=None, right=None):
        self.left_child = left
        self.right_child = right
