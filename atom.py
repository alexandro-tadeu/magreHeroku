#atom.py por Nigel S. Michki; parte de nmwizConvertor.py
#Define o tipo de dados Atom e quaisquer funções necessárias para operar em seus dados.

 

class Atom(object):
    def __init__(self):
        self.resID = ""
        self.atomName = ""
        self.resName = ""
        self.resNumber = ""
        self.A3D = " "
        self.coordinates = ["0.0" ,"0.0" ,"0.0"]
        self.bFactor = ""
        self.tag = ""   #Última coluna do arquivo .pdb, geralmente MAIN ou SOLV
        self.modes = []

    def printCoordinates(self):
        return (str(self.coordinates[0]) + " " + str(self.coordinates[1]) 
                + " " + str(self.coordinates[2])) 

    def printMode(self, modeNumber):
        return (self.modes[modeNumber-1]).printMode()

