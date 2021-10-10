#models.py

class dataData:
    def __init__(self,filepath,func):
        self.filepath = filepath
        self.func = func

class proteinData:
    def __init__(self,id,filepath,func):
        self.id = id
        self.filepath = filepath
        self.func = func

class matchData:
    def __init__(self,id,filepath,res,func):
        self.id = id
        self.filepath = filepath
        self.aa = res[0].upper()
        self.pos = int(res[1:])
        self.func = func

class pmidData:
    def __init__(self,id,filepath,func):
        self.id = id
        self.filepath = filepath
        self.func = func

class digestData:
    def __init__(self,id,protease,filepath,func):
        self.id = id
        self.protease = protease
        self.filepath = filepath
        self.func = func

class genData:
    def __init__(self,filepath,lst,func):
        self.filepath = filepath
        self.func = func
        self.input = lst

class alignData:
    def __init__(self,filepath,infile,func):
        self.filepath = filepath
        self.func = func
        with open(infile,'r') as f:
            self.input = f.read().splitlines()
