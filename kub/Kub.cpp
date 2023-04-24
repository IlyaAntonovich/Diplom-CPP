#include "Header.h"
#include "complex.h"
#include "Main.h"
#include "Out.h"

int main()
{
    // Консоль на кирилице
    setlocale(LC_ALL, "RUS");
    setprecision(4);

    // l3 = l^3
    int i, j, l3;
    double vol;

    cout << "h = " << h << endl;
    cout << "N = " << N << endl;

    // Размерность l3 x l3
    l3 = l * l * l;

    complex** A;
    CreateMatrixMemory(l3, A);
    complex* B;
    CreateVectorMemory(l3, B, 1.0); // integral3 
    complex* U;
    CreateVectorMemory(l3, U, 1.0);
    
    vol = h * h * h;

    cout << "Заполнение матрицы А <" << l3 << " x " << l3 << "> 6-тиричным интегралом" << endl;
    double procent = 100 / (double)(l3); // Процент выполнения
    for (i = 0; i < l3; i++){
        
        if (i % 100 == 0) {cout << i * procent <<" %" << endl; } // Процент выполнения
        
        for (j = 0; j < l3; j++){
            if (i != j) {
                A[i][j] = (-1.0)*integral6(k0, i, j, h, ip, l);
            }
            else {
                A[i][j] = vol - integral6(k0, i, j, h, ip, l);
            }
        }
    }
    cout << 100 << " %" << endl;
    cout << "---------" << endl;

    cout << A[0][0] << endl;
    Gauss(l3, A, B, U);
    cout << "Решение нелинейного уравнения методом Гауса завершено" << endl;

    

    cout <<"---------"<< endl;
    
    // Создаём файл 3d
    Exprint(l, 0, "C:\\Temp\\U.txt", U);

    //Создаём командный файл для Gnuplot
    ComGnuplot("C:\\Temp\\U.txt", "C:\\Temp\\ComU.txt");

    // Создаём график 3d через Gnuplot
    system("C:\\gnuplot\\bin\\wgnuplot.exe -p C:\\Temp\\ComU.txt");

    cout <<"---------"<< endl;

    DeleteVectorMemory(U);
    DeleteVectorMemory(B);
    DeleteMatrixMemory(l3, A);
 
    return 0;
}
