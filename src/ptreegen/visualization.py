## @package visualization
# Contains just the ptreegen::visualization::Visualization class.
#
import codecs
import os

from ete2 import TreeStyle

from enums import *


class Visualization:

    def __init__(self, tree, options):
        self._tree = tree
        self._options = options
        self._formatingMethods = []
        self._storageDir = None
        self._filePrefix = None
        self._style = TreeStyle()
        self._style.show_leaf_name = True
        self._style.show_branch_length = False
        self._style.scale = 80
        self._style.branch_vertical_margin = 15
        self._style.rotation = 90
        self._style.arc_start = -180
        self._style.arc_span = 180

        self.parseOptions(options)

    def parseOptions(self, options):
        self._storageDir = os.path.dirname(options["alignment_file"])
        self._filePrefix = os.path.basename(options["alignment_file"])
        if '.' in self._filePrefix:
            self._filePrefix = self._filePrefix[0:self._filePrefix.rfind('.')]
        out_formats = options['out_form'].split(',')
        for form in out_formats:
            if form == OutputForm.NEWICK:
                self._formatingMethods.append(self.saveNewick)
            elif form == OutputForm.IMAGE_PNG:
                self._formatingMethods.append(self.savePNG)
            elif form == OutputForm.IMAGE_SVG:
                self._formatingMethods.append(self.saveSVG)
            elif form == OutputForm.GUI:
                self._formatingMethods.append(self.showGUI)
            elif form == OutputForm.PRINT:
                self._formatingMethods.append(self.printToStdout)
            else:
                raise RuntimeError("Unknown output format: " + form)
        if options['tree_type'] == TreeType.CIRC:
            self._style.mode = 'c'
        elif options['tree_type'] == TreeType.RECT:
            self._style.mode = 'r'
        else:
            raise RuntimeError("Unknown tree type: " + options['tree_type'])

    def saveNewick(self):
        file_path = os.path.join(self._storageDir, self._filePrefix + ".newick")
        with codecs.open(file_path, "w", "utf-8") as outfile:
            outfile.write(self._tree.write(format=5))

    def savePNG(self):
        file_path = os.path.join(self._storageDir, self._filePrefix + ".png")
        self._tree.render(file_path, tree_style=self._style, dpi=300)

    def saveSVG(self):
        file_path = os.path.join(self._storageDir, self._filePrefix + ".svg")
        self._tree.render(file_path, tree_style=self._style, dpi=300)

    def showGUI(self):
        self._tree.show(tree_style=self._style)

    def printToStdout(self):
        print self._tree

    def show(self):
        for func in self._formatingMethods:
            func()
