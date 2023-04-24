#pragma once
#include "Header.h"

const int   l = 10;
const int   ip = 2;
const int   N = l * l * l;
const double k0 = 0.5;

const double a = 0.0;
const double b = 1.0;
const double h = (b - a) / ((double)(l));

void test()
{
    // Счётчики
    int i, j, k;

    // 
    double xi, yi, zi;

    // l - размерность куба
    // ip - количесвво точек интегрирования
    // ni - номер куба 1
    // nj - номер куба 2
    int  ni, nj;

    // a - нулевая левая нижняя точка
    // b - конечная точка
    // h - шаг
    // sum - численное решение итеграла
    double sum;

    // col - колличество столбцов
    // row - колличество строк
    // kcol - выбор столбца
    int row, col, kcol;

  
    ni = 1;
    nj = 2;



    sum = 0;

    cout << "Int6 = " << integral6(k0, ni, nj, h, ip, l) << endl;

    cout << "Int3 = " << integral3(k0, ni, h, ip, l) << endl;

}

// Выделение памяти под квадратную матрицу и вектор
void CreateMatrixMemory(int N, double**& A)
{
    // Вход: N - размерноть квадратной матрицы, А - указатель матрицы
    // Выход: А - указатель матрицы

    int i1, i2;

    A = new double* [N];
    for (i1 = 0; i1 < N; i1++) {
        A[i1] = new double[N];
        for (i2 = 0; i2 < N; i2++) {
            A[i1][i2] = 0.0;
        }
    }
}
void CreateMatrixMemory(int N, complex**& A)
{
    // Вход: N - размерноть квадратной матрицы, А - указатель матрицы
    // Выход: А - указатель матрицы

    int i1, i2;

    A = new complex* [N];
    for (i1 = 0; i1 < N; i1++) {
        A[i1] = new complex[N];
        for (i2 = 0; i2 < N; i2++) {
            A[i1][i2] = 0.0;
        }
    }
}
void CreateVectorMemory(int N, double*& V, double fill)
{
    int i1;

    V = new double[N];
    for (i1 = 0; i1 < N; i1++) {
        V[i1] = fill;
    }
}
void CreateVectorMemory(int N, complex*& V, complex fill)
{
    int i1;

    V = new complex[N];
    for (i1 = 0; i1 < N; i1++) {
        V[i1] = fill;
    }
}

// Удаление квадратной матрицы или вектора 
void DeleteMatrixMemory(int N, double** A)
{
    int i1;

    for (i1 = 0; i1 < N; i1++) {
        delete A[i1];
    }
    delete[]A;


}
void DeleteMatrixMemory(int N, complex** A)
{
    int i1;

    for (i1 = 0; i1 < N; i1++) {
        delete A[i1];
    }
    delete[]A;


}
void DeleteVectorMemory(double* V)
{
    delete[]V;
}
void DeleteVectorMemory(complex* V)
{
    delete[]V;
}

// Поэлементное копирование матрицы или вектора
void EqualMatrix(int N, double** A1, double** A2)
{
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            A2[i][j] = A1[i][j];
        }
    }
}
void EqualVector(int N, double* B1, double* B2)
{
    int i;

    for (i = 0; i < N; i++) {
        B2[i] = B1[i];
    }
}
void EqualMatrix(int N, complex** A1, complex** A2)
{
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            A2[i][j] = A1[i][j];
        }
    }
}
void EqualVector(int N, complex* B1, complex* B2)
{
    int i;

    for (i = 0; i < N; i++) {
        B2[i] = B1[i];
    }
}
// Координата X Y Z
double coord_x(int ii, double h)
{
    // ii - координата по X
    // h - длина шага
    return (double)(ii)+h;
}
double coord_y(int jj, double h)
{
    // jj - координата по Y
    // h - длина шага
    return (double)(jj)+h;
}
double coord_z(int kk, double h)
{
    // kk - координата по Z
    // h - длина шага
    return (double)(kk)+h;
}

// Метод Гаусса
double Gauss(int n, double** A, double* B, double* U)
{
    /*
    Вход:
    n - размерность квадратной матрицы
    A - матрица row x col (коэффициенты переменных)
    B - вектор столбец row (свободные члены)
    u - ??

    Выход:
    mult - ??
    */

    int i, j, k, ii, jj, kk;
    double d, s, mult;
    double  bb, cc;
    double eps;
    eps = 0.0000001;

    double** AA;
    CreateMatrixMemory(n, AA);

    double* BB;
    CreateVectorMemory(n, BB, 0.0);

    EqualMatrix(n, A, AA);
    EqualVector(n, B, BB);

    for (k = 0; k < n; k++) {

        if (abs(AA[k][k]) < eps)
        {
            kk = 0;
            for (ii = k + 1; ii < n; ii++)
            {
                if ((abs(AA[ii][k]) > eps) && (kk == 0))
                {
                    kk = 1;
                    cc = BB[k];
                    BB[k] = BB[ii];
                    BB[ii] = cc;

                    for (jj = k; jj < n; jj++)
                    {
                        bb = AA[k][jj];
                        AA[k][jj] = AA[ii][jj];
                        AA[ii][jj] = bb;
                    }
                }
            }
            if (kk == 0) {
                cout << "System error!!" << endl;
                system("pause");
                exit(1);
            }
        }

        for (j = k + 1; j < n; j++)
        {
            d = AA[j][k] / AA[k][k];
            for (i = k; i < n; i++)
            {
                AA[j][i] = AA[j][i] - d * AA[k][i];
            }
            BB[j] = BB[j] - d * BB[k];
        }
    }

    mult = 1.0;
    for (ii = 0; ii < n; ii++) {
        mult *= AA[ii][ii];
    }

    for (k = n - 1; k > -1; k--) {
        d = 0.0;
        for (j = k; j < n; j++) {
            s = AA[k][j] * U[j];
            d = d + s;
        }
        U[k] = (BB[k] - d) / AA[k][k];
    }

    DeleteVectorMemory(BB);
    DeleteMatrixMemory(n, AA);

    return mult;
}

complex Gauss(int n, complex** A, complex* B, complex* U)
{
    /*
    Вход:
    n - размерность квадратной матрицы
    A - матрица row x col (коэффициенты переменных)
    B - вектор столбец row (свободные члены)
    u - ??

    Выход:
    mult - ??
    */

    int i, j, k, ii, jj, kk;
    complex d, s, mult;
    complex  bb, cc;
    double eps;
    eps = 0.0000001;

    complex** AA;
    CreateMatrixMemory(n, AA);

    complex* BB;
    CreateVectorMemory(n, BB, 0.0);

    EqualMatrix(n, A, AA);
    EqualVector(n, B, BB);

    for (k = 0; k < n; k++) {

        if (abs(AA[k][k]) < eps)
        {
            kk = 0;
            for (ii = k + 1; ii < n; ii++)
            {
                if ((abs(AA[ii][k]) > eps) && (kk == 0))
                {
                    kk = 1;
                    cc = BB[k];
                    BB[k] = BB[ii];
                    BB[ii] = cc;

                    for (jj = k; jj < n; jj++)
                    {
                        bb = AA[k][jj];
                        AA[k][jj] = AA[ii][jj];
                        AA[ii][jj] = bb;
                    }
                }
            }
            if (kk == 0) {
                cout << "System error!!" << endl;
                system("pause");
                exit(1);
            }
        }

        for (j = k + 1; j < n; j++)
        {
            d = AA[j][k] / AA[k][k];
            for (i = k; i < n; i++)
            {
                AA[j][i] = AA[j][i] - d * AA[k][i];
            }
            BB[j] = BB[j] - d * BB[k];
        }
    }

    mult = 1.0;
    for (ii = 0; ii < n; ii++) {
        mult *= AA[ii][ii];
    }

    for (k = n - 1; k > -1; k--) {
        d = 0.0;
        for (j = k; j < n; j++) {
            s = AA[k][j] * U[j];
            d = d + s;
        }
        U[k] = (BB[k] - d) / AA[k][k];
    }

    DeleteVectorMemory(BB);
    DeleteMatrixMemory(n, AA);

    return mult;
}
// Ищем координаты левого нижнего кубика по его номеру
void Get_Pos(int ni, int l, double h, double& xi, double& yi, double& zi)
{
    /*
    Вход:
    ni - нумер маленького кубика
    l - размерноть большого куба (состоит из маленькиз кубиков)
    h - шаг или длина грани маленького кубика

    Выход: его координаты xi, yi, zi.
    */

    // Счётчики
    int i, j, k;

    // Счётчик кубиков
    int n = 0;

    // Проходим кубики пока не найдём
    for (k = 0; k < l; k++)
    {
        for (i = 0; i < l; i++)
        {
            for (j = 0; j < l; j++)
            {
                // Номера кубиков совпадают?
                if (ni == n)
                {
                    // Записываем текущие координаты и выходим из j-го цикла 
                    xi = (double)(j)*h;
                    yi = (double)(i)*h;
                    zi = (double)(k)*h;
                }
                n += 1;
            }
        }
    }
}

complex f(double xp, double yp, double zp)
{
    //return exp(xp + yp + zp);
    return  1.0;
}
complex G(double k, double xp, double yp, double zp, double xq, double yq, double zq)
{
    double r = sqrt((xp - xq) * (xp - xq) + (yp - yq) * (yp - yq) + (zp - zq) * (zp - zq));
    return exp(_i*k*r) / r;
}

complex integral6(double k0, int ni, int nj, double h, double ip, int l)
{
    // Счётчики
    int i, j, k;

    // позиции кубиков p_i или p_j
    int i1, i2, i3, j1, j2, j3;
    complex sum;
    double h1 = h / ((double)(ip));
    double h2 = h / ((double)(ip + 1));

    // координаты кубика p_i или p_j и бегующая точка интегррирования
    double xi, yi, zi, xp, yp, zp;
    double xj, yj, zj, xq, yq, zq;

    // xi, yi, zi - Вернуть из функции Get_pos
    Get_Pos(ni, l, h, xi, yi, zi);
    Get_Pos(nj, l, h, xj, yj, zj);
    //cout << xi << yi << zi << endl;
    //cout << xj << yj << zj << endl;
    /*xi = 0.0;    yi = 0.0;    zi = 0.0;
    xj = 2. * h; yj = 2. * h; zj = 2. * h;*/

    sum = 0.;

    for (i1 = 0; i1 < ip; i1++) {
        xp = xi + i1 * h1 + h1 / 2.;
        for (i2 = 0; i2 < ip; i2++) {
            yp = yi + i2 * h1 + h1 / 2.;
            for (i3 = 0; i3 < ip; i3++) {
                zp = zi + i3 * h1 + h1 / 2.;

                for (j1 = 0; j1 < ip + 1; j1++) {
                    xq = xj + j1 * h2 + h2 / 2.;
                    for (j2 = 0; j2 < ip + 1; j2++) {
                        yq = yj + j2 * h2 + h2 / 2.;
                        for (j3 = 0; j3 < ip + 1; j3++) {
                            zq = zj + j3 * h2 + h2 / 2.;

                            sum += G(k0, xp, yp, zp, xq, yq, zq);
                            // cout << xp << yp << zp << endl << xq << yq << zq << endl;;
                            // system("pause");
                        }
                    }
                }
            }
        }
    }

    //cout << "h = "<< h << endl;
    //cout << "h1 = " << h1 << endl;
    return -1.0 * sum * h1 * h1 * h1 * h2 * h2 * h2;
}
complex integral3(double k0, int ni, double h, double ip, int l)
{
    // Счётчики
    int i, j, k;

    // позиции кубиков p_i или p_j
    int i1, i2, i3, j1, j2, j3;
    complex sum;
    double h1 = h / ((double)(ip));

    // координаты кубика p_i или p_j // бегующая точка интегррирования
    double xi, yi, zi, xp, yp, zp;

    // xi, yi, zi - Вернуть из функции Get_pos
    Get_Pos(ni, l, h, xi, yi, zi);

    sum = 0.;
    for (i1 = 0; i1 < ip; i1++) {
        xp = xi + i1 * h1 + h1 / 2.;
        for (i2 = 0; i2 < ip; i2++) {
            yp = yi + i2 * h1 + h1 / 2.;
            for (i3 = 0; i3 < ip; i3++) {
                zp = zi + i3 * h1 + h1 / 2.;

                sum += f(xp, yp, zp);

            }
        }
    }

    return sum * h1 * h1 * h1;
}