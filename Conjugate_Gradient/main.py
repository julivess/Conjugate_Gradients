# This is a sample Python script.
import sys
import CG_test_tesk
import tesks
from PyQt5 import uic
import numpy as np
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import QApplication, QMainWindow, QHeaderView, QLabel

from PyQt5.QtWidgets import QTableWidgetItem
# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.



class App(QMainWindow):
    def __init__(self):
        super().__init__()
        uic.loadUi('form.ui', self)
        self.tableWidget.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.button.clicked.connect(self.res) ##описываем обработчик нажатия


    def res(self):  #обработчик нажатия кнопки button
        n = int(self.lineEdit_n.text())
        m = int(self.lineEdit_m.text())
        eps = float(self.lineEdit_eps.text())
        Nmax = int(self.lineEdit_Nmax.text())
        self.tableWidget.setColumnCount(0)
        self.tableWidget.setRowCount(0)
        self.tableWidget.setColumnCount(n + 2)
        self.tableWidget.setRowCount(m + 2)



        o = self.comboBox.currentIndex()
        if o == 0:

            tesks.runTest(self, n, m, Nmax, eps)
        else:
            tesks.runMain(self, n, m, Nmax, eps, Nmax, eps)




if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    ex.show()
    sys.exit(app.exec())
   ## Conjugate_gradient
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
