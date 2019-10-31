#include <cstdio>
#include <cmath>
#include <valarray>
#include <vector>
#include <string>
#include "DoubleMass.h"
#include "RungeKutta.h"
#include <omp.h>

const char * gnuplot = "\"/usr/bin/gnuplot\" -persist";

using namespace std;

void saveToFile(const char * fName, valarray<double> & x, valarray<double> & y, 
        valarray<valarray<double> > & T);

void startMenuLoop();

void IntegrateDoubleMass();
void multitreadingOptions();
void clearStream();

int main(int argc, char * argv[]) {
    startMenuLoop();
    
    return 0;
}

void clearStream() {
    while (getchar() != '\n') {}
}

void multitreadingOptions() {
    int nThreads = omp_get_max_threads();
    printf("Number of thread [%d]: ", nThreads);
    scanf("%d", &nThreads);
    omp_set_num_threads(nThreads);
}

void startMenuLoop() {
    vector<string> items;
    items.push_back("1 - Integrate DoubleMass");
    items.push_back("9 - Multi-threading options");
    items.push_back("0 - Exit");

    int selected = -1;
    while (selected != 0) {
        selected = -1;
        system("clear");
        for (int n = 0; n < (int)items.size(); n++) {
            printf("%s\n", items[n].c_str());
        }
        printf("Select: ");
        
        int res = scanf("%d", &selected);
        if (res != 1) {
            clearStream();
            continue;
        }
        clearStream();
        switch (selected) {
        case 1:
            IntegrateDoubleMass();
            break;
        case 5:
            multitreadingOptions();
            break;
        }

    }

}

void saveToFile(const char * fName, valarray<double> & x, valarray<double> & y, 
        valarray<valarray<double> > & T) {
    FILE *f = fopen(fName, "w");
    int N_x = (int)x.size();
    int N_y = (int)y.size();

    for (int i = 0; i < N_x; i++) {
        fprintf(f, "\t%.15lg", x[i]);
    }
    fprintf(f, "\n");

    for (int j = 0; j < N_y; j++) {
        fprintf(f, "%.15lg", y[j]);
        for (int i = 0; i < N_x; i++) {
            fprintf(f, "\t%.15lg", T[i][j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

void inputMesh(const char * par1Name, const char * par2Name, int * par1N, double * par1min, double * par1max, 
                                                             int * par2N, double * par2min, double * par2max) {
    printf("%s, nodes: ", par1Name);
    scanf("%d", par1N);
    printf("%s, min  : ", par1Name);
    scanf("%lf", par1min);
    printf("%s, max  : ", par1Name);
    scanf("%lf", par1max);

    printf("%s, nodes: ", par2Name);
    scanf("%d", par2N);
    printf("%s, min  : ", par2Name);
    scanf("%lf", par2min);
    printf("%s, max  : ", par2Name);
    scanf("%lf", par2max);
}

void createMesh(int N, double left, double right, valarray<double> & mesh, bool logScale = false) {
    double step;
    mesh.resize(N);
    if (logScale) {
        double p1 = log10(left);
        double p2 = log10(right);
        step = (p2 - p1) / (N - 1);
        for (int i = 0; i < N; i++)
            mesh[i] = pow(10, p1 + i * step);
    }
    else {
        step = (right - left) / (N - 1);
        for (int i = 0; i < N; i++)
            mesh[i] = left + i * step;
    }
}

void allocate2d(int N1, int N2, valarray<valarray<double> > & A) {
    A.resize(N1);
    for (int i = 0; i < N1; i++) {
        A[i].resize(N2, 0.0);
    }
}

void IntegrateDoubleMass() {
    DoubleMass odes;
    odes.setParameter(0, 1.0);//A
    odes.setParameter(1, 1.0);//B
    odes.setParameter(2, 5.0);//J
    odes.setParameter(3, 1.0);//mp
    odes.setParameter(4, 1.0);//r
    odes.setParameter(5, 2.0);//R
    odes.setParameter(6, 2.0 * M_PI);//Omega
    odes.setParameter(7, 0.0);//Gamma0
    odes.setParameter(8, 0.1);//epsGamma
    odes.setParameter(9, 1.6);//delta
    odes.setParameter(10, 0.1);//mu1
    odes.setParameter(11, 0.1);//mu2
    odes.setParameter(12, 0.01);//mu3

    RungeKutta method;
    method.init(&odes);

    const int N = 100001;
    double dt = 0.001; 
    double t[N];
    double q[N][6];

    t[0] = 0.0; q[0][0] = 0.0; q[0][1] = 0.0; q[0][2] = 0.0; q[0][3] = 0.0; q[0][4] = 0.0; q[0][5] = 0.0;
    for (int i = 0; i < N - 1; i++) {
        method.makeStep(&odes, t[i], q[i], q[i + 1], dt);
        t[i + 1] = t[i] + dt;
    }
    
    FILE * f = fopen("DoubleMassIntegrated.csv", "w");
    for (int i = (N - 1) / 10 * 9; i < N; i++) {
        fprintf(f, "%.15lg\t%.15lg\t%.15lg\t%.15lg\t%.15lg\t%.15lg\t%.15lg\n", 
                    t[i], q[i][0], q[i][1], q[i][2], q[i][3], q[i][4], q[i][5]);
    }
    fclose(f);

    FILE * pipe = popen(gnuplot, "w");
    fprintf(pipe, 
            "set terminal pngcairo enhanced size 800, 600 font \'Times, 16\'\n"
            "set style line 12 dashtype 2 lc \'black\' \n"
            "set grid xtics,ytics ls 12\n"
            "set output \'p1.png\'\n"
            "plot \'DoubleMassIntegrated.csv\' using 1:2 w l lw 2 lc \'black\'\n"
            "set output \'p2.png\'\n"
            "plot \'DoubleMassIntegrated.csv\' using 1:3 w l lw 2 lc \'black\'\n"
            "set output \'M.png\'\n"
            "plot \'DoubleMassIntegrated.csv\' using 1:4 w l lw 2 lc \'black\'\n"
            "set output \'xy.png\'\n"
            "plot \'DoubleMassIntegrated.csv\' using 5:6 w l lw 2 lc \'black\'\n"
            "set output \'phi.png\'\n"
            "plot \'DoubleMassIntegrated.csv\' using 1:7 w l lw 2 lc \'black\'\n");

    pclose(pipe);
}

