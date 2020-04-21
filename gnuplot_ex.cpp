#include <iostream>
#include <cstdio>
#include <stdio.h>
#include <math.h>
#include <bits/stdc++.h>

using namespace std;


#define casualgnuplot "C:\\gnuplot\\bin\\gnuplot -persist"

#define f first
#define s second
#define pb push_back

using namespace std;

const double EPS = 1e-9;

const int N = 1e5 + 10;

int m;

int t[N], bb[N];

int n;

class Matrix {
public:
    double a[100][100];
    int N;
    int M;

    Matrix() {
    }

    Matrix(int N, int M) {
        this->N = N;
        this->M = M;
    }

    friend istream & operator >> (istream &in, Matrix &F) {
        in >> F.N;
        in >> F.M;
        for (int i = 1; i <= F.N; i++) {
            for (int j = 1; j <= F.M; j++) {
                in >> F.a[i][j];
            }
        }
        return in;
    }

    friend ostream & operator << (ostream &out, const Matrix &F) {
        for (int i = 1; i <= F.N; i++) {
            for (int j = 1; j <= F.M; j++) {
                out << fixed << setprecision(2) << F.a[i][j];
                if (j < F.M) out << ' ';
            }
            out << endl;
        }
        return out;
    }


    Matrix operator+(Matrix b) {
        Matrix c = Matrix(N, M);
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= M; j++) {
                c.a[i][j] = a[i][j]+b.a[i][j];
            }
        }
        return c;
    }

    Matrix operator*(Matrix b) {
        Matrix c = Matrix(N, b.M);
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= b.M; j++) {
                c.a[i][j] = a[i][1] * b.a[1][j];
                for (int k = 2; k <= M; k++) {
                    c.a[i][j] = c.a[i][j] + a[i][k] * b.a[k][j];
                }
            }
        }
        return c;
    }

    Matrix transponse() {
        Matrix c = Matrix(M, N);
        for (int j = 1; j <= M; j++) {
            for (int i = 1; i <= N; i++) {
                c.a[j][i] = a[i][j];
            }
        }
        return c;
    }
};

class EliminationMatrix: public Matrix {
public:
    EliminationMatrix(int N, pair <int, int> coor, double val){
        this->N = N;
        this->M = N;
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                this->a[i][j] = 0.0;
                if (i == j) this->a[i][j] = 1.0;
            }
        }
        this->a[coor.f][coor.s] = val;
    }
};

class IdentityMatrix: public Matrix {
public:
    IdentityMatrix(int N) {
        this->N = N;
        this->M = N;
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                this->a[i][j] = 0;
                if (i == j) this->a[i][j] = 1;
            }
        }
    }
};

class PermutationMatrix: public Matrix {
public:
    PermutationMatrix(int N, int row, int torow){
        this->N = N;
        this->M = N;
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                this->a[i][j] = 0;
            }
        }

        vector <int> permutation;
        int cur = 0;
        for (int i = 1; i <= N; i++) permutation.pb(i);
        swap(permutation[torow - 1], permutation[row - 1]);

        for (int i = 1; i <= N; i++) {
            this->a[i][permutation[i - 1]] = 1.0;
        }
    }
};

void EliminationProcesstoFindInverseMatrix(Matrix K, Matrix& I) {
    int curcolumn = 1;
    while (curcolumn < K.N) {
        int mx = curcolumn;
        for (int i = curcolumn; i <= K.N; i++) {
            if (abs(K.a[i][curcolumn]) > abs(K.a[mx][curcolumn])) mx = i;
        }
        if (mx != curcolumn) {
            Matrix P = PermutationMatrix(K.N, mx, curcolumn);
            K = P * K;
            I = P * I;
        }
        for (int i = curcolumn + 1; i <= K.N; i++) {
            Matrix E = EliminationMatrix(K.N, {i, curcolumn}, -(K.a[i][curcolumn] / K.a[curcolumn][curcolumn]));
            K = E * K;
            I = E * I;
        }
        curcolumn++;
    }
    curcolumn = K.N;
    while (curcolumn > 1) {
        for (int i = curcolumn - 1; i >= 1; i--) {
            double diff = -(K.a[i][curcolumn] / K.a[curcolumn][curcolumn]);
            K.a[i][curcolumn] = 0;
            for (int j = 1; j <= K.N; j++) {
                I.a[i][j] += diff * I.a[curcolumn][j];
            }
        }
        curcolumn--;
    }
    for (int i = 1; i <= K.N; i++) {
        for (int j = 1; j <= K.N; j++) {
            I.a[i][j] /= K.a[i][i];
            I.a[i][j] += EPS;
        }
    }
}

Matrix b;
Matrix A;
Matrix A_T;

void input() {
    cin >> m;
    b = Matrix(m, 1);
    for (int i = 1; i <= m; i++) {
        cin >> t[i] >> bb[i];
    }
    cin >> n;
    A = Matrix(m, n + 1);
    for (int i = 1; i <= m; i++) {
        int cur = 1;
        for (int j = 1; j <= n + 1; j++) {
            A.a[i][j] = cur * 1.0;
            cur = cur * t[i];
        }
    }
    for (int i = 1; i <= m; i++) {
        b.a[i][1] = bb[i] * 1.0;
    }
}


int main() {
    input();
    cout << "A:" << '\n';
    cout << A;
    A_T = A.transponse();
    cout << "A_T*A:" << '\n';
    Matrix K = A_T*A;
    cout << K;
    cout << "(A_T*A)^-1:" << '\n';
    Matrix I = IdentityMatrix(K.N);
    EliminationProcesstoFindInverseMatrix(K, I);
    cout << I;
    cout << "A_T*b:" << '\n';
    Matrix G = A_T*b;
    cout << G;
    cout << "x~:" << '\n';
    Matrix X = I*G;
    cout << X;


    FILE *pipe = _popen(casualgnuplot, "w");

    if (pipe != NULL) {
        const int npoints = 200;
        const double step = (t[m] - t[1])*1.0 / npoints;
        fprintf(pipe, "%s\n", "plot '-' using 1:2 title 'appr' with lines, '-' using 1:2 title 'exp' w p ls 1");

        for (int i = 0; i < npoints + 1; i++) {
            double x = t[1] + i * step;
            double curdegreex = 1.0;
            double y = 0;
            for (int j = 0; j <= n; j++) {
                y = y + curdegreex * X.a[j + 1][1];
                curdegreex = curdegreex * x;
            }
            fprintf(pipe, "%f\t%f\n", x, y);
        }
        fprintf(pipe, "%s\n", "e");

        for (int i = 1; i <= m; i++) {
            double x = t[i];
            double y = bb[i];
            fprintf(pipe, "%f\t%f\n", x, y);
        }
        fprintf(pipe, "%s\n", "e");
        fflush(pipe);
        _pclose(pipe);
    } else {
        cout << "Could not open pipe" << endl;
    }
    return 0;
}
