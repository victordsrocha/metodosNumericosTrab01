#include <iostream>
#include <cmath>
#include <list>
#include <fstream>

using namespace std;

double a3 = 1;
double a2 = 1;
double precisao = 0.001;
double d0 = 0.5;
double lambda = 0.05;
int n_lambda = 1;
int n_parametros = 1;
list<double> lambda_list;
list<double> a3_list;
list<double> a2_list;

class Resultado {

public:
    Resultado();

    //bool solucao{};
    double raiz_1;
    double raiz_2;
    double raiz_3;
    double f_d_1;
    double f_d_2;
    double f_d_3;
    double a3;
    double a2;
    double lambda;
    string metodo;
    int num_iter;

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

    void gerar_quadro_csv();
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

void Quadro::gerar_quadro_csv() {
    ofstream myFile;
    myFile.open("resultados.csv");

    myFile << "id,a3,a2,metodo,lambda,d,numero de iteracoes" << endl;

    int tam = resultados.size();
    for (int i = 0; i < tam; ++i) {

        Resultado r = resultados.front();

        myFile << i + 1 << ","
               << r.a3 << ","
               << r.a2 << ","
               << r.metodo << ",";

        if (r.solucao) {
            if (r.metodo == "FL") {
                myFile << r.lambda << ","
                       << r.d << ","
                       << r.num_iter << endl;
            } else {
                myFile << "" << ","
                       << r.d << ","
                       << r.num_iter << endl;
            }
        } else {
            if (r.metodo == "FL") {
                myFile << r.lambda << ","
                       << "" << ","
                       << "" << endl;
            } else {
                myFile << "" << ","
                       << "" << ","
                       << "" << endl;
            }
        }

        resultados.pop_front();
    }
}


double fpendulo(double d) {
    return a3 * pow(d, 3) - 9 * a2 * d + 3;
}

bool bolzano_fpendulo(double a, double b) {
    return fpendulo(a) * fpendulo(b) < 0;
}

// Preenche vetor isolamentos com os isolamentos das 3 raizes procuradas no turno atual
// utilizando o teorema de bolzano no intervalo [-100,100]
double *gerar_vetor_isolamentos(double *isolamentos) {
    double min = -100; // valor minimo para a busca de isolamentos
    double max = 100; // valor maximo para a busca de isolamentos
    double resolucao = 0.1; // resolucao da busca

    // a e b são os valores testados a cada iteracao
    double a = min;
    double b = min + resolucao;

    int i = 0; // posicao do vetor isolamentos
    while (i < 6 || b > max) {
        // cria um par de isolamentos sempre que o teorema de bolzano retorna verdadeiro
        if (bolzano_fpendulo(a, b)) {
            isolamentos[i] = a;
            isolamentos[i + 1] = b;
            i += 2;
        }
        a = a + resolucao;
        b = b + resolucao;
    }
    return isolamentos;
}

// gera os pontos inicias com base na média dos valores a e b de isolamento da raiz
double *gerar_vetor_pontos_iniciais(double *iniciais, const double *isolamentos) {
    for (int i = 0; i < 3; ++i) {
        iniciais[i] = (isolamentos[2 * i] + isolamentos[2 * i + 1]) / 2;
    }
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

void encontrar_raizes(){
    for (int i = 0; i < n_parametros; ++i) {
        a3 = a3_list.front();
        a2 = a2_list.front();



    }

}

void print_cabecalho() {
    cout << "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-++-+-+-+-+-+-+-+-+-+-+" << endl;
    cout << "+    Calculadora de raizes da funçao f(d) = (a3)d^3 - 9(a2)d + 3       +" << endl;
    cout << "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-++-+-+-+-+-+-+-+-+-+-+" << endl;
    cout << endl;
}

void input_precisao() {
    system("clear");
    print_cabecalho();
    cout << "Digite um valor para a precisao: ";
    cin >> precisao;

    while ((!cin) || precisao <= 0 || precisao > 0.1) {
        cin.clear();
        cin.ignore(numeric_limits<int>::max(), '\n');
        cout
                << "Precisao não está no formato correto. A Precisao deve ser um valor real no intervalo ]0,0.1].\nDigite novamente: ";
        cin >> precisao;
    }
}

void input_parametros() {
    system("clear");
    print_cabecalho();

    cout << "Digite a quantidade de pares para os parametros a3 e a2: ";
    cin >> n_parametros;

    while ((!cin) || n_parametros <= 0 || n_parametros > 10) {
        cin.clear();
        cin.ignore(numeric_limits<int>::max(), '\n');
        cout << "Quantidade não está no formato correto. "
                "A quantidade deve ser um número inteiro maior ou igual a 1 e menor do que 10.";
        cout << "\nDigite novamente: ";
        cin >> n_lambda;
    }

    for (int i = 0; i < n_parametros; ++i) {
        double a3_, a2_;
        cout << "\nPar (" << i + 1 << "). Digite um valor para o parametro a3: ";
        cin >> a3_;

        while ((!cin)) {
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            cout << "Formato incorreto. O parametro a3 deve ser um numero real.\nDigite novamente: ";
            cin >> a3_;
        }

        cout << "Par (" << i + 1 << "). Digite um valor para o parametro a2: ";
        cin >> a2_;

        while ((!cin)) {
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            cout << "Formato incorreto. O parametro a2 deve ser um numero real.\nDigite novamente: ";
            cin >> a2_;
        }

        a3_list.push_back(a3_);
        a2_list.push_back(a2_);
    }


}

void input_lambda() {
    system("clear");
    print_cabecalho();
    cout << "Digite a quantidade de opcoes para lambda: ";
    cin >> n_lambda;

    while ((!cin) || n_lambda <= 0 || n_lambda > 10) {
        cin.clear();
        cin.ignore(numeric_limits<int>::max(), '\n');
        cout
                << "Quantidade não está no formato correto. A quantidade deve ser um número inteiro maior ou igual a 1 e menor do que 10.\nDigite novamente: ";
        cin >> n_lambda;
    }

    for (int i = 0; i < n_lambda; ++i) {
        double lambda_;
        cout << "Digite o valor de lambda da opcao " << i + 1 << ": ";
        cin >> lambda_;
        while ((!cin || lambda_ <= 0 || lambda_ > 0.1)) {
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            cout << "Formato incorreto. Lambda deve ser um valor real no intervalo ]0,0.1].\nDigite novamente: ";
            cin >> lambda_;
        }
        lambda_list.push_back(lambda_);
    }
}

void input() {
    input_parametros();
    input_lambda();
    input_precisao();
    system("clear");
}

int main() {
    Quadro quadro;

    input();

    /*
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

    //quadro.imprimir_quadro();
    quadro.gerar_quadro_csv();

    double isolamentos[6];
    gerar_vetor_isolamentos(isolamentos);
    double iniciais[3];
    gerar_vetor_pontos_iniciais(iniciais, isolamentos);
    */

    return 0;
}


