#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>
#include "complex.h"

using namespace std;

void CreateMatrixMemory(int N, double**& A);
void CreateVectorMemory(int N, double*& V, double fill);
void DeleteMatrixMemory(int N, double** A);
void DeleteVectorMemory(double* V);
void CreateMatrixMemory(int N, complex**& A);
void CreateVectorMemory(int N, complex*& V, complex fill);
void DeleteMatrixMemory(int N, complex** A);
void DeleteVectorMemory(complex* V);
void EqualVector(int N, double* B1, double* B2);
void EqualMatrix(int N, double** A1, double** A2);
void Terminal_Matrix(double** a, int row, int col, int s);
void Terminal_Vector(double* a, int row);
void num(int l);
void test();
void Exprint(int l, int layer, const char *FileName, double* U);
void ComGnuplot(const char* FileName3d, const char* FileName3dCom);

double coord_x(int ii, double h);
double coord_y(int jj, double h);
double coord_z(int kk, double h);
double Gauss(int n, double** A, double* B, double* U);
complex Gauss(int n, complex** A, complex* B, complex* U);
complex integral3(double k0, int ni, double h, double ip, int l);
complex integral6(double k0, int ni, int nj, double h, double ip, int l);
complex f(double xp, double yp, double zp);
complex G(double k, double xp, double yp, double zp, double xq, double yq, double zq);