#include <stdio.h>
#include <math.h>

#define m 12  // размер области по горизонтали
#define n 10  // размер области по вертикали
#define R 2   // радиус круглого выреза
#define T 100 // постоянная температура в правом верхнем углу
#define T3 50 // начальная температура внутри области

#define EPSILON 0.0001  // Допустимая погрешность
#define MAX_ITER 1000   // Максимальное количество итераций

void initialize(double temperature[m][n]);
void calculateTemperature(int iterations, double temperature[m][n]);
void computeTemperature(double temperature[m][n]);

int main() {
    double temperature[m][n];

    initialize(temperature);

    printf("Before:\n");

    // Вывод распределения температур
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (pow(i-round((double)(m/2-1)), 2) + pow(j-round((double)(n/2)), 2) <= R*R) printf("   \t");
            else printf("%.2f\t", temperature[i][j]);
        }
        printf("\n");
    }

    computeTemperature(temperature);

    printf("\nAfter:\n");

    // Вывод распределения температур
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (pow(i-round((double)(m/2-1)), 2) + pow(j-round((double)(n/2)), 2) <= R*R) printf("   \t");
            else printf("%.2f\t", temperature[i][j]);
        }
        printf("\n");
    }
    
    return 0;
}

void initialize(double temperature[m][n]) {
    // Инициализация значений температур в начальный момент времени
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (j == n - 1) {
                temperature[i][j] = T - (double)T/(m-1)*i;  // равномерно убывающая температура на правой границе
            } 
            else if (j == 0 || i == m-1) {
                temperature[i][j] = 0;  // нулевая температура на левой и нижней границах
            } else {
                temperature[i][j] = T3; // начальная температура внутри области
            }
        }
    }
}

double getAverageTemp(double temperature[m][n], int row, int col) {

    if (col == 0 || col == n-1 || row == m-1) return temperature[row][col]; // постоянная температура

    double newTemp = 0.0;
    int count = 0;

    for (int i = row-1; i<=row+1;i++) {
        if (i < 0 || i > m-1) continue;
        for(int j = col-1;j<=col+1;j++) {
            if (j < 0 || j > n-1) continue;
            if (pow(i-round((double)(m/2-1)), 2) + pow(j-round((double)(n/2)), 2) > R*R) {
                if (i != row && j != col) newTemp += temperature[i][j] * 0.6;
                else if (i == row && j == col) newTemp += temperature[i][j] * 0.9;
                else newTemp += temperature[i][j];
                count++;
            }
        }
    }

    if (count != 0) return newTemp/count;
    else return temperature[row][col];

}

void computeTemperature(double temperature[m][n]) {

    // Метод последовательных приближений (метод прогонки)
    double diff = 1.0;  // Разница между текущей и предыдущей итерациями
    int iter = 0;       // Счетчик итераций

    double tempNew[m][n];

    while (diff > EPSILON) {
        diff = 0.0;

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                
                tempNew[i][j] = getAverageTemp(temperature, i,j);
                //printf("[%d, %d]: temp=%lf\n", i, j, tempNew[i][j]);

                
            }
        }

        for (int i = 0; i < m; i++) {
            for (int j = 0;j < n; j++) {
                diff += fabs(tempNew[i][j] - temperature[i][j]);
                temperature[i][j] = tempNew[i][j];
                
            }
        }

        //printf("diff = %lf\n", diff);

        iter++;
    }

    printf("Count: %d", iter);

    
}