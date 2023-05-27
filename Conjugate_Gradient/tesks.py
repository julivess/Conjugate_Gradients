import numpy as np
import matplotlib.pyplot as plt
import CG_main_tesk
import CG_test_tesk
from PyQt5.QtWidgets import QApplication, QMainWindow, QHeaderView, QTableWidgetItem, QLabel
from PyQt5.QtGui import QColor
def runMain(self, n, m, N_max, eps, N_max_x2, eps_x2):
    v = CG_main_tesk.Conjugate_gradient(n, m, N_max, eps)
    spravka_main(self, n, m, eps, N_max)
    v_x2 = CG_main_tesk.Conjugate_gradient(2*n, 2*m, N_max_x2, eps_x2)
    eps_tesk = 0
    for i in range(1, m):
        for j in range(1, n):
            if(np.abs(v[i][j] - v_x2[2 * i][2 * j]) > eps_tesk):
                eps_tesk = np.abs(v[i][j] - v_x2[2 * i][2 * j])

    spravka_main2(self, eps, eps_tesk, N_max)
    renderMain(v, v_x2, n, m)
    fillTable(self,v, n, m, 0, 3, 0, 1)

def runTest(self, n, m, N_max, eps):
    self.label_spravka2.setText(' ')
    v = CG_test_tesk.Conjugate_gradient(n, m, N_max, eps)
    renderTest(v, n, m)
    fillTable(self,v, n, m, 0, 3, 0, 1)
    spravka_test(self, n, m, eps, N_max)

def renderMain(v, v_x2, n, m):
    x = np.linspace(0, 3, 2 * n + 1)
    y = np.linspace(0, 1, 2 * m + 1)
    xgrid, ygrid = np.meshgrid(x, y)

    fig = plt.figure()
    V_g = fig.add_subplot(1, 2, 1, projection='3d')

    U_g = fig.add_subplot(1, 2, 2, projection='3d')
    U_g.plot_surface(xgrid, ygrid, v_x2, color="plum")
    U_g.set_title('Численное решение на половинной сетке')
    U_g.set_xlabel('x', fontsize=15)
    U_g.set_ylabel('y', fontsize=15)
    U_g.set_zlabel('z', fontsize=15)

    x = np.linspace(0, 3, n + 1)
    y = np.linspace(0, 1, m + 1)
    xgrid, ygrid = np.meshgrid(x, y)
    V_g.plot_surface(xgrid, ygrid, v, color="coral")
    V_g.set_title('Численное решение')
    V_g.set_xlabel('x', fontsize=15)
    V_g.set_ylabel('y', fontsize=15)
    V_g.set_zlabel('z', fontsize=15)

    mng = plt.get_current_fig_manager()
    mng.window.showMaximized()
    plt.show()
def renderTest(v, n, m):
    x = np.linspace(0, 3, 50)
    y = np.linspace(0, 1, 50)
    xgrid, ygrid = np.meshgrid(x, y)

    z = CG_test_tesk.U_xy(xgrid, ygrid)

    fig = plt.figure()
    V_g = fig.add_subplot(1, 2, 1, projection='3d')

    U_g = fig.add_subplot(1, 2, 2, projection='3d')
    U_g.plot_surface(xgrid, ygrid, z, color="plum")
    U_g.set_title('Точное решение')
    U_g.set_xlabel('x', fontsize=15)
    U_g.set_ylabel('y', fontsize=15)
    U_g.set_zlabel('z', fontsize=15)

    x = np.linspace(0, 3, n + 1)
    y = np.linspace(0, 1, m + 1)
    xgrid, ygrid = np.meshgrid(x, y)
    V_g.plot_surface(xgrid, ygrid, v, color="coral")
    V_g.set_title('Численное решение')
    V_g.set_xlabel('x', fontsize=15)
    V_g.set_ylabel('y', fontsize=15)
    V_g.set_zlabel('z', fontsize=15)

    mng = plt.get_current_fig_manager()
    mng.window.showMaximized()
    plt.show()


def fillTable(self, v, n, m, a, b, c, d):
    h = (b - a) / n
    k = (d - c) / m

    up = ['x_i']
    i = 0
    while i < n + 1:
        up.append(str(a + h * i))
        self.tableWidget.setItem(0, i + 1, QTableWidgetItem(str(i)))
        self.tableWidget.item(0, i + 1).setBackground(QColor(205, 240, 255))
        i += 1
    self.tableWidget.setHorizontalHeaderLabels(up)

    down = ['y_j']
    i = 0
    while i < m + 1:
        down.append(str(c + k * i))
        self.tableWidget.setItem(i + 1, 0, QTableWidgetItem(str(i)))
        self.tableWidget.item(i + 1, 0).setBackground(QColor(205, 240, 255))
        i += 1
    self.tableWidget.setVerticalHeaderLabels(down)
    self.tableWidget.setItem(0, 0, QTableWidgetItem("j\i"))
    self.tableWidget.item(0, 0).setBackground(QColor(205, 240, 255))
    i = 0
    while i < m + 1:
        j = 0
        while j < n + 1:
            self.tableWidget.setItem(i + 1, j + 1, QTableWidgetItem(str(v[i][j])))
            j += 1
        i += 1
def spravka_test(self, n, m, eps, Nmax):
    self.label_spravka.setText('Справка для тестовой задачи' + '\n' +
                               'Для решения тестовой задачи использовалась сетка с числом разбиений по x n= ' + str(n) +
                               ' и числом разбиений по y= ' + str(m) +'\n'+ 'Метод сопряженных градиентов'+'\n' +
                               'критерий оставновки по точности Eps = ' + str(eps) + ' и по числу итераций Nmax=  ' + str(Nmax)+ '\n'
                               'На схемы (СЛАУ) затрачено итераций N=' + str(CG_test_tesk.S)+ ' и достигнута точность итерационного метода ε = ' +str(CG_test_tesk.eps_max)+ '\n'
                               'Схема (СЛАУ) решена с невязкой ||Rn|| = ' + str(CG_test_tesk.max_r)+ ' Для невязки СЛАУ использована норма максимума'+'\n'+
                               'Тестовая задача должна быть решена с погрешностью не более ε = 0.5*10^(-6)'+ '\n'+
                               'Задача решена с погрешностью ε = ' + str(CG_test_tesk.eps_tesk) +'\n'+
                               'Максимальное отклонение точного и численного решения наблюдается в узле \n' +
                               'x= ' + str(CG_test_tesk.xm) +' y= '+ str(CG_test_tesk.ym)+'\n'+
                               'В качвестве начального приближения использовано нулевое приближение\n' +
                               'Невязка СЛАУ на начальном приближении ||Ro||= ' +str(CG_test_tesk.max_ro)+ " Норма максимума ")
def spravka_main(self, n, m, eps, Nmax):
    self.label_spravka.setText('Справка для основной задачи' + '\n' +
                               'Для решения основной задачи использовалась сетка с числом разбиений по x n= ' + str(n) +
                               ' и числом разбиений по y= ' + str(m) +'\n'+ 'Метод сопряженных градиентов'+'\n' +
                               'критерий оставновки по точности ε = ' + str(eps) + ' и по числу итераций Nmax=  ' + str(Nmax)+ '\n'
                               'На схемы (СЛАУ) затрачено итераций N=' + str(CG_main_tesk.S)+ ' и достигнута точность итерационного метода ε = ' +str(CG_main_tesk.eps_max)+ '\n'
                               'Схема (СЛАУ) решена с невязкой ||Rn|| = ' + str(CG_main_tesk.max_r)+ ' Для невязки СЛАУ использована норма максимума'+'\n'+
                               'Максимальное отклонение точного и численного решения наблюдается в узле \n' +
                               'x= ' + str(CG_main_tesk.xm) +' y= '+ str(CG_main_tesk.ym)+'\n'+
                               'В качвестве начального приближения на основной сетке использовано нулевое приближение\n' +
                               'На основной сетке невязка СЛАУ на начальном приближении ||Ro||= ' +str(CG_main_tesk.max_ro)+ " Норма максимума " +'\n\n')
def spravka_main2(self, eps, eps_tesk, Nmax):
    self.label_spravka2.setText('Для контроля точности использована сетка с половинным шагом,'+'\n'+
                               'метод сопряженных градиаентов'+'\n'+
                               'критерий оставновки по точности ε2 = ' + str(eps) + ' и по числу итераций Nmax2=  ' + str(Nmax) + '\n' +
                               'На схемы (СЛАУ) затрачено итераций N2=' + str(CG_main_tesk.S)+ ' и достигнута точность итерационного метода ε2= ' +str(CG_main_tesk.eps_max)+ '\n'+
                               'Схема (СЛАУ) на сетке с половинным шагом решена с невязкой ||Rn2|| = ' + str(CG_main_tesk.max_r)+ ' Для невязки СЛАУ использована норма максимума'+'\n'+
                               'Основная задача должна быть решена с точностью не более ε = 0.5*10^(-6)' + '\n' +
                               'Задача решена с точностью ε = ' + str(eps_tesk) + '\n' +
                               'Максимальное отклонение численного решения на основной сетке и сетке с половинным шагом наблюдается в узле \n' +
                               'x= ' + str(CG_main_tesk.xm) + ' y= ' + str(CG_main_tesk.ym) + '\n' +
                               'В качвестве начального приближения использовано нулевое приближение\n' +
                               'Невязка СЛАУ на начальном приближении ||Ro||= ' + str(CG_main_tesk.max_ro) + " Норма максимума ")