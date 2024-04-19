#include <iostream>
#include <cmath>
#include <vector>
#include <sciplot/sciplot.hpp>

//Forward methods - deprecated!
double timestep_forward(double dx, double dt, double lambda, std::vector<double> T, int j){
    double factor = lambda * dt / pow(dx, 2);
    return 
        T[j] + factor * (T[j+1] - 2 * T[j] + T[j-1]);
}

std::pair<std::vector<std::vector<double>>, std::vector<double>> forward_elimination(std::vector<std::vector<double>> tridiag, std::vector<double> T){
    for (int i = 1; i < tridiag.size(); i++) {
        double m = tridiag[i][i-1] / tridiag[i-1][i-1];
        tridiag[i][i] -= m * tridiag[i-1][i];
        T[i] -= m * T[i-1];
    }
    return {tridiag, T};
}

std::pair<std::vector<std::vector<double>>, std::vector<double>> back_substitution(std::vector<std::vector<double>> tridiag, std::vector<double> T, double T_x_0, double T_x_L){
    std::vector<double> T_next(T.size());
    int n = T.size() - 1;
    T_next[n] = T[n] / tridiag[n][n]; // Start by solving the last row
    for (int i = n-1; i >= 0; i--) {
        T_next[i] = (T[i] - tridiag[i][i+1] * T_next[i+1]) / tridiag[i][i];
    }
    return {tridiag, T_next};
}
std::vector<std::vector<double>> tridiagonal_constructor(const int N, double factor){
    std::vector<std::vector<double>> tridiag(N, std::vector<double>(N, 0));
    for (int i = 0; i < N; i++){
        tridiag[i][i] = 1 + 2 * factor;
        if (i > 0){
            tridiag[i][i-1] = -factor;
        }
        if (i < N - 1){
            tridiag[i][i+1] = -factor;
        }
    }
    return tridiag;
}

//constructs a tridiagonal matrix to be used in the backward time step via the Thomas algorithm
std::vector<double> timestep_backward(double factor, std::vector<double> T, double T_x_0, double T_x_L){
    int N = T.size();
    std::vector<std::vector<double>> tridiag = tridiagonal_constructor(N, factor);
    std::pair<std::vector<std::vector<double>>, std::vector<double>> new_system = forward_elimination(tridiag, T);
    std::vector<double> T_next = back_substitution(new_system.first, new_system.second, T_x_0, T_x_L).second;
    T_next[0] = T_x_0;
    T_next[N - 1] = T_x_L;
    return T_next;
}



std::vector<double> T_setup(int N_1, int N_2, double L, double t_0, double t_1, double lambda, std::vector<double> boundary_conditions){
    double dx = L / N_1;
    double dt = (t_1 - t_0) / N_2;
    double T_x_0 = boundary_conditions[0];
    double T_x_L = boundary_conditions[1];

    std::vector<std::vector<double>> T(N_2, std::vector<double>(N_1, 0));
    T[0][0] = T_x_0;
    T[0][N_1 - 1] = T_x_L;
    for (int i = 1; i < N_2; i++){

        double factor = lambda * dt / pow(dx, 2);
        T[i] = timestep_backward(factor, T[i - 1], T_x_0, T_x_L);
        T[i][0] = T_x_0;
        T[i][N_1 - 1] = T_x_L;
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
    int N_1 = 400; //lengthwise discretization
    int N_2 = 300; //temporal discretization
    double L = 1.0; //length of rod
    double t_0 = 0.0; //initial time
    double t_1 = 1.0; //final time
    double lambda = 0.01; //thermal diffusivity

    std::vector<double> boundary_conditions = {100.0, 0.0}; //Dirichlet boundary conditions
    std::vector<double> T = T_setup(N_1, N_2, L, t_0, t_1, lambda, boundary_conditions);
    plot_T(T, L);
    return 0;
}