#include <stdlib.h> // for exit(1)

using namespace std;



void
LU_decomposition(double **L, double **U, int n, int m);
// здесь A(исходная) = L, либо можешь еще 1 аргумент добавить.
// PS функция "возвращает" 2 матрицы L и U
// По хорошему еще должен быть вектор P, задающий перестановки.
// Лидирующий элемент выбирается как max(abs(x)) по всему столбцу.

