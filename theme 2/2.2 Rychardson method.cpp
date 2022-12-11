#include <iostream>
#include <limits>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>

using std::cout;
using std::cin;
using std::endl;

std::ifstream fin("input.txt");
std::ofstream fout("output.txt");

class Matrix
{
private:
    double** M;
    int m, n;
    double determ = 1; //Коэффициент, необходимый для отслеживания влияния преобразования матрицы на определитель

public:

    Matrix() //Конструктор по умолчанию 
    {
        n = m = 0;
        M = nullptr;
    }
    Matrix(int n)
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
                fout << std::fixed << std::setprecision(8) << GetElem(i, j) << " ";
            }
            fout << endl;
        }
        fout << endl;
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
    double norm()
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
    Matrix operator*(const Matrix& matrix)
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
    Matrix operator*(double a)
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
    Matrix operator/(double a)
    {
        if (a == 0)
        {
            return 0;
        }
        Matrix result = Matrix(m, n);
        double elem;
        for (size_t i = 0; i < m; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                elem = M[i][j] / a;
                result.SetElem(i, j, elem);
            }
        }
        return result;
    }
};

double dotProduct(double* A, double* B, int n) //скалярное произведение двух векторов
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


    int dimension;
    double elem, epsilon, lambdaMax;

    fin >> dimension;

    Matrix A = Matrix(dimension, dimension); //Создал матрицу
    for (size_t i = 0; i < dimension; i++) //Ввод элементов матрицы
    {
        for (size_t j = 0; j < dimension; j++)
        {
            fin >> elem;
            A.SetElem(i, j, elem);
        }
    }

    Matrix b = Matrix(dimension, 1);
    for (size_t i = 0; i < dimension; i++) //Ввод свободных коэффициентов
    {
        fin >> elem;
        b.SetElem(i, 0, elem);
    }
    Matrix xPrevPrev = Matrix(dimension, 1);
    for (size_t i = 0; i < dimension; i++)// начальное приближение
    {
        fin >> elem;
        xPrevPrev.SetElem(i, 0, elem);//x0
    }
    fin >> epsilon;
    size_t k = 0;

    // степенной метод
    
    Matrix prevEVector = Matrix(dimension, 1);
    for (size_t i = 0; i < dimension; i++)//приближение
    {
        prevEVector.SetElem(i,0,xPrevPrev.GetElem(i,0));
    }
    Matrix firstEVector = prevEVector;
    double prevEValue = prevEVector.norm();
    Matrix currEVector = A * (prevEVector / prevEVector.norm());
    double currEValue = currEVector.norm();
   
    //подсчет максимального собственного значения
    //условие выхода: т.к. последовательность из рассчитанных с.з. стремиться к макс. с.з, из цикла выйдем, когда разность с.з станет < epsilon
    while (abs(currEValue - prevEValue) > epsilon)
    {
        k++;
        prevEVector = currEVector; //сохраняем предыдущий с.в
        currEVector = A * currEVector * (1 / currEVector.norm()); //считаем новый по формуле
        prevEValue = currEValue;  //сохраняем предыдущее с.з
        currEValue = currEVector.norm(); //новое текущее с.з - норма собственного вектора
    }
    lambdaMax = currEValue;
    fout << "Максимальное собственное значение: " << std::fixed << std::setprecision(8) << lambdaMax << endl;
    fout << "Количество итераций:" << k << endl << endl;

    //минимального собственного значения
    k = 0; //обнулил счетчик итераций
    prevEVector = firstEVector;
    Matrix I = Matrix(dimension);
    currEVector = (I * A.norm() - A) * prevEVector * (1 / prevEVector.norm());
    prevEValue = A.norm() - prevEVector.norm();
    currEValue = A.norm() - currEVector.norm();
    //условие выхода: аналогично, только последовательность стремится к минимальному с.з
    while (abs(prevEValue - currEValue) > epsilon)
    {
        //все аналогично
        k++;
        prevEVector = currEVector;
        currEVector = (I * A.norm() - A) * prevEVector * (1 / prevEVector.norm());
        prevEValue = currEValue;
        currEValue = A.norm() - currEVector.norm();
    }
    double lambdaMin = currEValue;
    fout << "Минимальное собственное значение: " << std::fixed << std::setprecision(8) << lambdaMin << endl;
    fout << "Количество итераций:" << k << endl;

    // трехчленный метод ричардсона

    double stolbNorm = b.norm();
    const double tao = 2 / (lambdaMin + lambdaMax);
    double omegaFirst = -((lambdaMax - lambdaMin) / (lambdaMax + lambdaMin)); //w1
    Matrix r = A * xPrevPrev - b;
    Matrix xPrev = xPrevPrev - r * tao;
    r = A * xPrev - b;
    double omegaPrevPrev = omegaFirst;
    double omegaPrev = 1 / (2 * (1 / omegaFirst) - omegaFirst); //w2
    Matrix xCurr = xPrev + (xPrev - xPrevPrev) * (omegaFirst * omegaPrev) - r * (tao * (1 + omegaFirst * omegaPrev)); //x2
    r = A * xCurr - b; //r2
    double omegaCurr;

    while ((r.norm() / stolbNorm) > epsilon) //условие выхода - когда норма невязки / норму столбца станет < epsilon
    {
        k++;
        omegaCurr = 1 / (2 * (1 / omegaFirst) - omegaPrev);
        xPrevPrev = xPrev;
        xPrev = xCurr;
        xCurr = xPrev + (xPrev - xPrevPrev) * (omegaPrev * omegaCurr) - r * (1 + omegaPrev * omegaCurr) * tao; 
            r = A * xCurr - b;
        omegaPrev = omegaCurr;
        double norma = r.norm();
    }
    fout << endl;
    xCurr.print();
    fout << "Количество итераций: " << k << endl;
}
