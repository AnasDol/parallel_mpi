#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define T 100 // постоянная температура в правом верхнем углу
#define T3 50 // начальная температура внутри области

#define EPSILON 0.0001  // Допустимая погрешность
#define MAX_ITER 1000   // Максимальное количество итераций

#define min(a,b) (((a) < (b)) ? (a) : (b))

void initialize(double** temperature, int m, int n);
double get_new_temp(double** temperature, int rows, int cols, int current_row, int current_col, double R);
void compute_temperature(double** temperature, int m, int n, double R);

int main() {

    printf("Row number (>=2): ");
    int m; // размер области по вертикали
    if (!(scanf("%d", &m)==1 && m>=2)) {
        printf("Input error\n");
        return 0;
    }

    printf("Column number (>=2): ");
    int n; // размер области по горизонтали
    if (!(scanf("%d", &n)==1 && n>=2)) {
        printf("Input error\n");
        return 0;
    }

    printf("Cut radius (1<R<%.2lf): ", (double)min(m,n)/2);
    double R; // радиус круглого выреза
    if (!(scanf("%lf", &R)==1 && R>=1 && R<=(double)min(m,n)/2)) {
        printf("Input error\n");
        return 0;
    }

    double** temperature = (double**)malloc(m * sizeof(double*));
    for (int i = 0; i < m; i++) {
        temperature[i] = (double*)malloc(n * sizeof(double));
    }

    initialize(temperature, m, n);

    compute_temperature(temperature, m, n, R);

    if (m <= 10 && n <= 10) {
        for (int i = 0; i < n; i++) {
            printf("%.2lf\t", temperature[0][i]);
        }
    }

    // Вывод распределения температур
    // for (int i = 0; i < m; i++) {
    //     for (int j = 0; j < n; j++) {
    //         if (pow(i-round((double)(m/2-1)), 2) + pow(j-round((double)(n/2)), 2) <= R*R) printf("   \t");
    //         else printf("%.2f\t", temperature[i][j]);
    //     }
    //     printf("\n");
    // }


    for (int i = 0;i<m;i++) {
        free(temperature[i]);
    }
    free(temperature);
    
    return 0;
}

void initialize(double** temperature, int rows, int cols) {
    // Инициализация значений температур в начальный момент времени
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (j == cols - 1) {
                temperature[i][j] = T - (double)T/(rows-1)*i;  // равномерно убывающая температура на правой границе
            } 
            else if (j == 0 || i == rows-1) {
                temperature[i][j] = 0;  // нулевая температура на левой и нижней границах
            } else {
                temperature[i][j] = T3; // начальная температура внутри области
            }
        }
    }
}

double get_new_temp(double** temperature, int rows, int cols, int current_row, int current_col, double R) {

    if (current_col == 0 || current_col == cols-1 || current_row == rows-1) return temperature[current_row][current_col]; // постоянная температура

    double newTemp = 0.0;
    int count = 0;

    for (int i = current_row-1; i<=current_row+1;i++) {
        if (i < 0 || i > rows-1) continue;
        for(int j = current_col-1;j<=current_col+1;j++) {
            if (j < 0 || j > cols-1) continue;
            if (pow(i-round((double)(rows/2-1)), 2) + pow(j-round((double)(cols/2)), 2) > R*R) {
                newTemp += temperature[i][j];
                count++;
            }
        }
    }

    if (count != 0) return newTemp/count;
    else return temperature[current_row][current_col];

}

void compute_temperature(double** temperature, int rows, int cols, double R) {

    // Метод последовательных приближений
    double diff = 1.0;  // Разница между текущей и предыдущей итерациями
    int iter = 0;       // Счетчик итераций

    double** new_temp = (double**)malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        new_temp[i] = (double*)malloc(cols * sizeof(double));
    }

    while (diff > EPSILON) {
        diff = 0.0;

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                
                new_temp[i][j] = get_new_temp(temperature, rows, cols, i, j, R);
            }
        }

        for (int i = 0; i < rows; i++) {
            for (int j = 0;j < cols; j++) {
                diff += fabs(new_temp[i][j] - temperature[i][j]);
                temperature[i][j] = new_temp[i][j];
                
            }
        }

        iter++;
    }

    printf("Count: %d", iter);

    for (int i = 0;i<rows;i++) {
        free(new_temp[i]);
    }
    free(new_temp);

    
}