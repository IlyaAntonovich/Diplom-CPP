#pragma once
#include "Header.h"

// ����� � ������� ��� ����

// ����� � �������� ������� ��� �������
void Terminal_Matrix(double** a, int row, int col, int s)
{
    cout << endl;
    cout << "�������:" << endl;
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            cout << setw(s) << a[i][j];
        }
        cout << endl;
    }
}
void Terminal_Vector(double* a, int row)
{
    cout << endl;
    cout << "������:" << endl;
    for (int i = 0; i < row; i++)
    {
        cout << a[i] << endl;
    }
    cout << endl;
}

// ��������� ������� � ����� �������
void num(int l) // ���-�� ��������� �������
{
    // ��������
    int i, j, k;

    // n - ����� �������
    int n = 0;

    for (k = 0; k < l; k++)
    {
        cout << "k = " << k << endl;
        for (i = 0; i < l; i++)
        {
            for (j = 0; j < l; j++)
            {
                cout << setw(l) << n;
                n += 1;
            }
            cout << endl;
        }
        cout << endl;
    }
}

// ������ ������� � ���� .txt
void Exprint(int l, int layer, const char* FileName, complex* U)
{
    int i, j, k, ind;

    ofstream file(FileName);
    ind = 0;
    for (i = 0; i < l; i++) 
    {
        for (j = 0; j < l; j++)
        {
            for (k = 0; k < l; k++)
            {
                if (layer == i)
                {

                    file  << abs(U[ind])<< "\t";
                }
                ind++;
            }
            if (layer == i)
            {
                file << endl;
            }
        }
    }
    file.close();
}

// �������� FileName3dCom 
void ComGnuplot(const char* FileName3d, const char* FileName3dCom) /* C:\\Temp\\ */
{
    const char* com1, * com2;

    com1 = "splot '";
    com2 = "' matrix with pm3d";

    ofstream file(FileName3dCom);
    file << com1 << FileName3d << com2;
    file.close();
}