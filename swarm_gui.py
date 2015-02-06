#!/usr/bin/env python

import sys
from PyQt4.QtCore import *
from PyQt4.QtGui import *


class GasInput(QWidget):
    """A widget to select input gases, their
    fractions and cross section files"""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.nRows = 0
        self.nameL = QLabel("Gas &name")
        self.fracL = QLabel("Gas &frac")
        self.fileL = QLabel("&Input file")
        self.addB = QPushButton("&Add gas")

        self.lo = QGridLayout()
        self.gasList = []
        self.addRow()

        self.lo.addWidget(self.nameL, 0, 0)
        self.lo.addWidget(self.fracL, 0, 1)
        self.lo.addWidget(self.fileL, 0, 2)
        self.lo.addWidget(self.addB, self.nRows+1, 0)
        self.connect(self.addB, SIGNAL("clicked()"), self.addRow)
        self.lo.setAlignment(Qt.AlignTop)
        self.setLayout(self.lo)

    def addRow(self):
        tmp = (QLineEdit(), QLineEdit(), QPushButton("Load"),
               QLabel("No file"), QToolButton(), {})
        tmp[4].setText("x")
        self.connect(tmp[2], SIGNAL("clicked()"), lambda: self.loadFile(tmp))
        self.connect(tmp[4], SIGNAL("clicked()"), lambda: self.rmRow(tmp))
        self.gasList.append(tmp)

        for col, w in enumerate(self.gasList[self.nRows][0:5]):
            self.lo.addWidget(w, self.nRows+1, col)

        self.nameL.setBuddy(self.gasList[self.nRows][0])
        self.fracL.setBuddy(self.gasList[self.nRows][1])
        self.fileL.setBuddy(self.gasList[self.nRows][2])
        self.lo.removeWidget(self.addB)
        self.nRows += 1
        self.lo.addWidget(self.addB, self.nRows+1, 0)

    def rmRow(self, row):
        for w in row[0:5]:
            w.close()
        self.gasList.remove(row)
        self.nRows -= 1

    def loadFile(self, row):
        fname = QFileDialog.getOpenFileName(self, "Load cross sections")
        if fname:
            row[5]['filename'] = fname
            if (len(fname) > 15):
                fname = "..." + fname[-12:]
            row[3].setText(fname)


class SwarmSettings(QWidget):
    """A widget for the swarm settings
    (pressure, temperature, field range etc.)"""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.inputs = {}
        self.fl = QFormLayout()
        self.addSpinBox("&Temperature (K)", 1e-6, 1e99, 300.)
        self.addSpinBox("&Pressure (bar)", 1e-6, 1e99, 1.)
        self.addSpinBox("&Min. field (V/m)", 1e3, 1e12, 1e6)
        self.addSpinBox("Ma&x. field (V/m)", 1e3, 1e12, 1e7)
        self.addSpinBox("N&um. of fields", 1, 1e6, 11)
        self.setLayout(self.fl)

    def addSpinBox(self, label, minval, maxval, defval):
        le = QLineEdit(str(defval))
        self.fl.addRow(label, le)
        self.inputs[label] = (le, minval, maxval, defval)


class SwarmMainWindow(QWidget):
    """The main window for the swarm gui"""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.gasInput = GasInput(self)
        self.gasSettings = SwarmSettings(self)
        self.lo = QGridLayout()

        self.lo.addWidget(self.gasInput, 0, 0)
        self.lo.addWidget(self.gasSettings, 0, 1)
        self.setLayout(self.lo)
        self.setWindowTitle("Swarm GUI")

app = QApplication(sys.argv)
smw = SwarmMainWindow()
smw.show()
app.exec_()
