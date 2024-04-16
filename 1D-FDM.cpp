#include <iostream>
#include <cmath>
#include <vector>
#include <sciplot/sciplot.hpp>


double timestep(double dx, double dt, double lambda, std::vector<double> T, int j){
    double factor = lambda * dt / pow(dx, 2);
    return 
        T[j] + factor * (T[j+1] - 2 * T[j] + T[j-1]);
}

std::vector<double> T_setup(int N_1, int N_2, double L, double t_0, double t_1){
    double dx = L / N_1;
    double dt = (t_1 - t_0) / N_2;
    double T_x_0 = 100.0;
    double T_x_L = 0.0;
    double lambda = 1;

    std::vector<std::vector<double>> T(N_2, std::vector<double>(N_1));

    for (int i = 0; i < N_2; i++){
        for (int j = 0; j < N_1; j++){
            if (j == 0){
                T[i][j] = T_x_0;
            }
            else if (j == N_1 - 1){
                T[i][j] = T_x_L;
            }
            else {
                T[i][j] = timestep(dx, dt, lambda, T[i], j);
            }
        }
    }
    return T[N_2 - 1];
}


std::vector<double> linspace(double start, double end, int num){
    std::vector<double> x(num);
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; i++){
        x[i] = start + i * step;
    }
    return x;
}


void plot_T(std::vector<double> T, double L){
    std::vector<double> x_range = linspace(0, L, T.size());
    sciplot::Plot2D plot;
    plot.drawCurve(x_range, T).label("T");
    plot.legend().hide();
    plot.xlabel("x");
    plot.ylabel("T");
    sciplot::Figure figure = {{plot}};
    sciplot::Canvas canvas = {{figure}};
    canvas.show();
}

int main(){
    int N_1 = 50;
    int N_2 = 10000;
    double L = 1.0;
    double t_0 = 0.0;
    double t_1 = 1.0;
    std::vector<double> T = T_setup(N_1, N_2, L, t_0, t_1);
    std::cout << "T: " << T[1] << std::endl;
    plot_T(T, L);
    return 0;
}