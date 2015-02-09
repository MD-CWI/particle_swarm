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


class smartLine(QObject):

    def __init__(self, label, defVal, validRange=None, parent=None):
        super(QObject, self).__init__(parent)
        self.t = type(defVal)
        self.w = QLineEdit(str(defVal))
        self.l = QLabel(label)
        self.l.setBuddy(self.w)  # Set widget and label as buddies
        self.v0 = defVal
        if validRange:
            if defVal < validRange[0] or defVal > validRange[1]:
                raise ValueError('Default value not in range')
            self.lbound = validRange[0]
            self.ubound = validRange[1]
        self.w.editingFinished.connect(self.validator)

    def validator(self):
        isValid = False
        try:
            val = self.t(self.w.text())
            if self.lbound <= val <= self.ubound:
                isValid = True
        except valueError:
            pass

        if not isValid:
            box = QMessageBox()
            msg = (self.l.text().replace("&", "") + " should lie between " +
                   str(self.lbound) + " and " + str(self.ubound))
            box.setText(msg)
            box.exec_()
            self.w.setText(str(self.v0))


class SwarmMainWindow(QWidget):
    """The main window for the swarm gui"""

    def __init__(self, parent=None):
        super().__init__(parent)

        self.vdict = {}
        self.vdict['temp'] = smartLine("&Temperature (K)",
                                       300., (1e-6, 1e99))
        self.vdict['pres'] = smartLine("&Pressure (bar)",
                                       1., (1e-6, 1e99))
        self.vdict['minfld'] = smartLine("&Min field (V/m)",
                                         1e6, (1e0, 1e99))
        self.vdict['maxfld'] = smartLine("Ma&x field (V/m)",
                                         1e7, (1e0, 1e99))
        self.vdict['nfld'] = smartLine("&Num fields",
                                       11, (1, 1000*1000))

        self.gasW = GasInput(self)
        self.lo = QGridLayout()
        self.lo.addWidget(self.gasW)

        for i, d in enumerate(self.vdict.values()):
            self.lo.addWidget(d.w, i+1, 1)
            self.lo.addWidget(d.l, i+1, 0)

        self.startB = QPushButton("Start computation")
        self.startB.clicked.connect(self.startSim)
        self.lo.addWidget(self.startB)
        self.setLayout(self.lo)
        self.setWindowTitle("Swarm GUI")

    def startSim(self):
        gasNames, gasFracs, gasFiles = self.gasW.getData()

app = QApplication(sys.argv)
smw = SwarmMainWindow()
smw.show()
app.exec_()
