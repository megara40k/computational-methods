#include <iostream>
#include <limits>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>

using std::cout;
using std::cin;
using std::endl;

class Matrix
{
private:
    double** M;
    int m, n;
   

public:
    Matrix() //Конструктор по умолчанию 
    {
        n = m = 0;
        M = nullptr;
    }
    Matrix(int n) // для единичной матрицы
    {
        this->m = n;
        this->n = n;

        M = (double**) new double* [m];

        for (size_t i = 0; i < m; i++) //Выделение памяти
        {
            M[i] = (double*) new double[n];
        }

        for (size_t i = 0; i < m; i++) //Заполнение диагонали единицами
        {
            for (size_t j = 0; j < n; j++)
            {
                if (i != j) M[i][j] = 0;
                else M[i][i] = 1;
            }
        }
    }
    Matrix(int m, int n) //Конструктор для прямоугольной матрицы
    {
        this->m = m; //Строки
        this->n = n; //Столбцы

        M = (double**) new double* [m];

        for (size_t i = 0; i < m; i++) //Выделение памяти
        {
            M[i] = (double*) new double[n];
        }

        for (size_t i = 0; i < m; i++) //Заполнение нулями
        {
            for (size_t j = 0; j < n; j++)
            {
                M[i][j] = 0;
            }
        }
    }
    Matrix(const Matrix& Matr) //Коснструктор копирования
    {
        this->m = Matr.m;
        this->n = Matr.n;

        this->M = (double**) new double* [m];

        for (size_t i = 0; i < m; i++)
        {
            this->M[i] = (double*) new double[n];
        }

        for (size_t i = 0; i < m; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                this->M[i][j] = Matr.M[i][j];
            }
        }


    }

    double GetElem(unsigned int i, unsigned int j) //Получение элемента 
    {
        if ((this->m > 0 && this->m > i) && (this->n > 0 && this->n > j))
        {
            return M[i][j];
        }
        else
        {
            return -0.000001;
        }
    }
    bool isSymmetrical()
    {
        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < i; j++)
            {
                if (M[i][j] != M[j][i]) return false;
            }
        }
        return true;
    }

    void SetElem(unsigned int i, unsigned int j, double a) //Изменение элемента
    {
        M[i][j] = a;
    }

    void multiplOnElem(int line, double elem) //Умножение строки на элемент
    {
        determ /= elem;
        for (size_t i = 0; i < n; i++)
        {
            M[line][i] = M[line][i] * elem;
        }

    }

    void sumLine(int second, int first) //Сумма строк
    {
        double result;
        for (size_t i = 0; i < n; i++)
        {
            result = M[second][i] + M[first][i];
            M[second][i] = result;
        }
    }

    void upperTriang()
    {
        for (size_t i = 0; i < m - 1; i++) //Верхнетреугольный вид
        {
            for (size_t j = 0; j < m - 1; j++)
            {
                if (j + 1 != i && M[i][i] != 0 && M[j + 1][i] != 0) //проверка на совпадение столбца и строки
                { //(если её не сделать, то i строка умножится на себя и сложится сама с собой(т.е. умножится на 2))
                    multiplOnElem(i, -(M[j + 1][i] / M[i][i]));
                    sumLine(j + 1, i);
                }
                else if (j + 1 != i && M[i][i] == 0 && M[j + 1][i] != 0)
                {
                    lineInterchange(j + 1, i);
                }

            }
        }
    }
    void lowerTriang()
    {
        for (size_t i = m - 1; i > 0; i--) //Нижнетреугольный вид
        {
            for (size_t j = m - 1; j > 0; j--)
            {
                if (j - 1 != i && M[i][i] != 0 && M[j - 1][i] != 0) //проверка на совпадение столбца и строки
                { //(если её не сделать, то i строка умножится на себя и сложится сама с собой(т.е. умножится на 2))
                    multiplOnElem(i, -(M[j - 1][i] / M[i][i]));
                    sumLine(j - 1, i);
                }
                else if (j - 1 != i && M[i][i] == 0 && M[j - 1][i] != 0) lineInterchange(j - 1, i);
            }
        }
    }

    void diag() //Делает матрицу диагональной
    {
        upperTriang();
        lowerTriang();
    }

    void print() //Вывод матрицы
    {
        for (size_t i = 0; i < m; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                cout << std::fixed << std::setprecision(8) << GetElem(i, j) << " ";
            }
            cout << endl;
        }
    }
    void print(int f, int s) //Вывод матрицы
    {
        for (size_t i = 0; i < f; i++)
        {
            for (size_t j = 0; j < s; j++)
            {
                if (GetElem(i, j) == 0) cout << std::fixed << std::setprecision(8) << GetElem(i, j) << " ";
                else if (GetElem(i, j) == -0) cout << std::fixed << std::setprecision(8) << -GetElem(i, j) << " ";
                else cout << std::fixed << std::setprecision(8) << -GetElem(i, j) << " ";
            }
            cout << endl;
        }
    }

    void lineInterchange(unsigned int first, unsigned int second) //Перестановка строк
    {
        determ *= -1;
        double elem;
        for (size_t i = 0; i < n; i++)
        {
            elem = M[first][i];
            M[first][i] = M[second][i];
            M[second][i] = elem;
        }
    }
    void columnInterchange(unsigned int first, unsigned int second) //Перестановка столбцов
    {
        determ *= -1;
        double elem;
        for (size_t i = 0; i < m; i++)
        {
            elem = M[i][first];
            M[i][first] = M[i][second];
            M[i][second] = elem;
        }
    }

    double determinant()
    {
        diag();
        for (size_t i = 0; i < m; i++)
        {
            determ *= M[i][i];
        }
        return determ;
    }

    double norm() //нормирование матрицы
    {
        double max = -INT16_MAX;
        for (size_t i = 0; i < m; i++)
        {
            double currSum = 0;
            for (size_t j = 0; j < n; j++)
            {
                currSum += abs(M[i][j]);
            }
            if (currSum > max) max = currSum;
        }
        return max;
    }

    Matrix operator=(const Matrix& Matr) //Оператор копирования
    {
        this->m = Matr.m;
        this->n = Matr.n;

        this->M = (double**) new double* [m];
        for (size_t i = 0; i < m; i++)
        {
            this->M[i] = (double*) new double[m];
        }

        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                this->M[i][j] = Matr.M[i][j];

        return *this;
    }
    Matrix operator+(const Matrix& matrix)
    {
        if ((matrix.n == 0 && n == 0) || (matrix.n != n))
        {
            Matrix result = Matrix();
            return result;
        }

        Matrix result = Matrix(m, n);
        for (size_t i = 0; i < m; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                result.M[i][j] = M[i][j] + matrix.M[i][j];
            }
        }

        return result;
    }
    Matrix operator-(const Matrix& matrix)
    {
        if ((matrix.n == 0 && n == 0) || (matrix.n != n))
        {
            Matrix result = Matrix();
            return result;
        }

        Matrix result = Matrix(matrix.m, matrix.n);
        for (size_t i = 0; i < m; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                result.M[i][j] = M[i][j] - matrix.M[i][j];
            }
        }

        return result;
    }
    Matrix operator*(const Matrix& matrix) // перемножение матриц
    {
        if ((matrix.n == 0 && n == 0) || (n != matrix.m)) //Совпадают ли строки левой матрицы со столбцами правой
        {
            Matrix result = Matrix();
            return result;
        }
        Matrix result = Matrix(matrix.m, matrix.n);
        {
            for (size_t i = 0; i < matrix.m; i++)
            {
                for (size_t j = 0; j < matrix.n; j++)
                {
                    result.M[i][j] = 0;
                    for (size_t k = 0; k < matrix.m; k++)
                    {
                        if (abs(M[i][k]) < 0.0000001)
                        {
                            result.M[i][j] += 0;
                        }
                        else if (abs(matrix.M[k][j]) < 0.0000001)
                        {
                            result.M[i][j] += 0;
                        }
                        else
                        {
                            result.M[i][j] += M[i][k] * matrix.M[k][j];
                        }

                    }

                }

            }
        }
        return result;
    }
    Matrix operator*(double a) // умножение вектора на скаляр
    {
        Matrix result = Matrix(m, n);
        double elem;
        for (size_t i = 0; i < m; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                elem = M[i][j] * a;
                result.SetElem(i, j, elem);
            }
        }
        return result;
    }

    double* operator()(int column)
    {
        double* result = new double[m];
        for (size_t i = 0; i < m; i++)
        {
            result[i] = M[i][column];
        }
        return result;
    }

};

double dotProduct(double* A, double* B, int n)
{
    double result = 0;
    for (size_t i = 0; i < n; i++)
    {
        result += A[i] * B[i];
    }
    return result;
}

void main()
{
    std::ifstream in("input.txt");
    setlocale(LC_ALL, "Russian");

    int n;
    double elem, epsilon;

    in >> n;//размерность

    Matrix matrix = Matrix(n, n); //Создал матрицу

    for (size_t i = 0; i < n; i++) //Ввод элементов матрицы (матрица должна быть симметричной!)
    {
        for (size_t j = 0; j < n; j++)
        {
            in >> elem;
            matrix.SetElem(i, j, elem);
        }
    }

    Matrix stolbec = Matrix(n, 1);
    for (size_t i = 0; i < n; i++) //Ввод свободных коэффициентов
    {
        in >> elem;
        stolbec.SetElem(i, 0, elem);
    }

    Matrix x = Matrix(n, 1); //начальное приближение
    Matrix r = Matrix(n, 1); //невязка
    Matrix z = Matrix(n, 1); //ошибка
    Matrix rPrev = Matrix(n, 1);//невязка предыдущая (для рекурсии)
    Matrix xPrev = Matrix(n, 1);//решение предыдущее (для рекурсии)
    Matrix zPrev = Matrix(n, 1);//ошибка предыдущая (для рекурсии)
    for (size_t i = 0; i < n; i++)
    {
        in >> elem;
        x.SetElem(i, 0, elem);
    }

    in >> epsilon;//погрешность


    //Считаем начальную невязку
//подготовка к итерационному процессу
    r = stolbec - matrix * x;
    z = r;

    double stolbNorm = stolbec.norm();
    double a, b;
    int k = 0;
//итерационный процесс
    while (epsilon < (r.norm()) / (stolbNorm)) //пока относительная невязка не станет >= заданному числу
    {
        k++;
        rPrev = r;
        xPrev = x;
        zPrev = z;

        a = (dotProduct(rPrev(0), rPrev(0), n)) / (dotProduct((matrix * zPrev)(0), zPrev(0), n));
        x = xPrev + zPrev * a;
        r = rPrev - (matrix * zPrev) * a;

        b = (dotProduct(r(0), r(0), n)) / (dotProduct(rPrev(0), rPrev(0), n));
        z = r + zPrev * b;
    }

    cout << "Ответ" << endl;
    x.print();
    cout << "Количество итераций: " << endl;
    cout << k;
}
