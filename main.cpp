#include <iostream>
#include <cmath>
#include <list>

using namespace std;

class Resultado {

public:
    Resultado();

    bool solucao{};
    double d{};
    double a3{};
    double a2{};
    double lambda{};
    string metodo;
    int num_iter{};

    void imprimir_resultado();
};

void Resultado::imprimir_resultado() {
    cout << "f(d) = (" << a3 << ")d^3 - 9(" << a2 << ")d + 3" << endl;
    if (metodo == "FL") cout << "lambda = " << lambda << endl;
    cout << "d = " << d << endl;
    cout << "metodo = " << metodo << endl;
    cout << "numero de iteracoes = " << num_iter << endl;
}

Resultado::Resultado() {
}


class Quadro {
public:
    list<Resultado> resultados;

    void imprimir_quadro();
};

void Quadro::imprimir_quadro() {
    int tam = resultados.size();
    for (int i = 0; i < tam; ++i) {
        cout << "Resultado " << i + 1 << ":\n\n";
        resultados.front().imprimir_resultado();
        resultados.pop_front();
        cout << "~~~~~~~~~~~~~~~~~~~~~~~\n\n";
    }
}


double a3 = 1;
double a2 = 1;
double precisao = 0.001;
double d0 = 0.5;
double lambda = 0.05;

double fpendulo(double d) {
    return a3 * pow(d, 3) - 9 * a2 * d + 3;
}

double derivada_fpendulo(double d) {
    return 3 * a3 * pow(d, 2) - 9 * a2;
}

double nDeriv(double x, double (*function)(double)) {
    double h = 0.001;
    return (function(x + h) - function(x - h)) / (2 * h);
}

Resultado newton_original() {
    double x = d0;
    int maxIteracoes = 100;
    int numIteracoes;
    bool solucao = false;

    for (int i = 0; i < maxIteracoes; ++i) {
        double y = fpendulo(x);
        double dy = derivada_fpendulo(x);

        if (abs(y) < precisao) {
            solucao = true;
            numIteracoes = i + 1;
            break;
        }
        x = x - y / dy;
    }
    Resultado resultado;
    resultado.solucao = solucao;
    if (solucao) {
        resultado.d = x;
        resultado.num_iter = numIteracoes;
    }
    resultado.a2 = a2;
    resultado.a3 = a3;
    resultado.metodo = "original";
    return resultado;
}

Resultado newton_FL() {
    double x = d0; // palpite inicial (derivada do ponto inicial deve ser >= lambda)
    int maxIteracoes = 20;
    bool solucao = false;
    double xw;
    int numIteracoes;

    for (int i = 0; i < maxIteracoes; ++i) {
        double y = fpendulo(x);
        double dy = derivada_fpendulo(x);

        if (abs(y) < precisao) {
            solucao = true;
            numIteracoes = i + 1;
            break;
        }

        double FL;
        if (abs(dy) > lambda) {
            xw = x;
            FL = dy;
        } else {
            FL = derivada_fpendulo(xw);
        }

        x = x - y / FL;
    }
    Resultado resultado;
    resultado.solucao = solucao;
    if (solucao) {
        resultado.d = x;
        resultado.num_iter = numIteracoes;
    }
    resultado.lambda = lambda;
    resultado.a2 = a2;
    resultado.a3 = a3;
    resultado.metodo = "FL";
    return resultado;
}

Resultado newton_derivada_numerica() {
    double x = d0; // palpite inicial
    int maxIteracoes = 20;
    bool solucao = false;
    int numIteracoes;

    for (int i = 0; i < maxIteracoes; ++i) {
        double y = fpendulo(x);
        double dy = nDeriv(x, fpendulo);

        if (abs(y) < precisao) {
            solucao = true;
            numIteracoes = i + 1;
            break;
        }

        x = x - y / dy;
    }
    Resultado resultado;
    resultado.solucao = solucao;
    if (solucao) {
        resultado.d = x;
        resultado.num_iter = numIteracoes;
    }
    resultado.a2 = a2;
    resultado.a3 = a3;
    resultado.metodo = "derivada num.";
    return resultado;
}

void encontrarRaizes(Quadro &quadro) {
    quadro.resultados.push_back(newton_original());
    quadro.resultados.push_back(newton_derivada_numerica());
    quadro.resultados.push_back(newton_FL());
}

int main() {
    Quadro quadro;
    int n_lambda;

    cout << "Insira o numero de opcoes para lambda: ";
    cin >> n_lambda;

    cout << "Insira a precisao desejada para os calculos: ";
    cin >> precisao;

    for (int i = 0; i < n_lambda; ++i) {

        cout << "Lambda numero " << i + 1 << endl;
        cout << "Insira os seguintes valores: \n";

        cout << "Valor de lambda: ";
        cin >> lambda;

        cout << "Valor de a3: ";
        cin >> a3;

        cout << "Valor de a2: ";
        cin >> a2;

        encontrarRaizes(quadro);
    }

    quadro.imprimir_quadro();
    return 0;
}


