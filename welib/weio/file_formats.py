
class FileFormat():
    def __init__(self,fileclass=None):
        self.constructor = fileclass
        if fileclass is None:
            self.extensions = []
            self.isValid = False
            self.name = ''
        else:
            self.extensions  = fileclass.defaultExtensions()
            self.name        = fileclass.formatName()
            self.isValid     = fileclass.isRightFormat

    def __repr__(self):
        return 'FileFormat object: {} ({})'.format(self.name,self.extensions[0])

