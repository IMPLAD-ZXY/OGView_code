class noGeneNameError(Exception):
    def __init__(self, feature):
        self.feature = feature

class noExonNumberError(Exception):
    def __init__(self, feature):
        self.feature = feature

class noAnnotatedGenomeError(Exception):
    def __init__(self):
        return
