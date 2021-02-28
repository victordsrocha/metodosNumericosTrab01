#include <iostream>
#include <cmath>
#include <list>
#include <fstream>

using namespace std;

int registro_id = 1;
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
double isolamentos[6];
double iniciais[3];
int qtd_raizes;

struct resultado_geral {
    bool solucao_1{};
    bool solucao_2{};
    bool solucao_3{};
    double isolamento_1[2]{};
    double isolamento_2[2]{};
    double isolamento_3[2]{};
    double raiz_1{};
    double raiz_2{};
    double raiz_3{};
    double f_d_1{};
    double f_d_2{};
    double f_d_3{};
    int num_iter_1{};
    int num_iter_2{};
    int num_iter_3{};
    double a3{};
    double a2{};
    double lambda{};
    string metodo;
};

struct resultado_individual {
    bool solucao{};
    double isolamento[2]{};
    double raiz{};
    double f_d{};
    int num_iter{};
};

double fpendulo(double d) {
    return a3 * pow(d, 3) - 9 * a2 * d + 3;
}

bool bolzano_fpendulo(double a, double b) {
    return fpendulo(a) * fpendulo(b) < 0;
}

// Preenche vetor isolamentos com os isolamentos das 3 raizes procuradas no turno atual
// utilizando o teorema de bolzano no intervalo [-100,100]
double *gerar_vetor_isolamentos(double *isolamentos_) {
    double min = -100; // valor minimo para a busca de isolamentos
    double max = 100; // valor maximo para a busca de isolamentos
    double resolucao = 0.1; // resolucao da busca

    // a e b são os valores testados a cada iteracao
    double a = min;
    double b = min + resolucao;

    int i = 0; // posicao do vetor isolamentos
    while (i < 6 && b <= max) {
        // cria um par de isolamentos sempre que o teorema de bolzano retorna verdadeiro
        if (bolzano_fpendulo(a, b)) {
            isolamentos_[i] = a;
            isolamentos_[i + 1] = b;
            i += 2;
        }
        a = a + resolucao;
        b = b + resolucao;
    }

    // Evita que o isolamento seja mostrado na tabela utilizando notacao exponencial
    for (int j = 0; j < 6; ++j) {
        if (isolamentos_[j] < 0.01 && isolamentos_[j] > -0.01){
            isolamentos_[j] = 0;
        }
    }

    qtd_raizes = i / 2;
    return isolamentos_;
}

// gera os pontos inicias com base na média dos valores a e b de isolamento da raiz
double *gerar_vetor_pontos_iniciais(double *iniciais_, const double *isolamentos_) {
    for (int i = 0; i < qtd_raizes; ++i) {
        iniciais_[i] = (isolamentos_[2 * i] + isolamentos_[2 * i + 1]) / 2;
    }
}

double derivada_fpendulo(double d) {
    return 3 * a3 * pow(d, 2) - 9 * a2;
}

double nDeriv(double x, double (*function)(double)) {
    double h = 0.001;
    return (function(x + h) - function(x - h)) / (2 * h);
}

resultado_individual newton_original(int i_x, int i_a) {

    double x = iniciais[i_x];

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
    resultado_individual res;
    res.solucao = solucao;
    if (solucao) {
        res.raiz = x;
        res.f_d = fpendulo(x);
        res.num_iter = numIteracoes;
    }

    res.isolamento[0] = isolamentos[i_a];
    res.isolamento[1] = isolamentos[i_a + 1];
    return res;
}

resultado_geral newton_original_geral() {

    resultado_individual r1 = newton_original(0, 0);

    resultado_geral rg;

    rg.a3 = a3;
    rg.a2 = a2;
    rg.metodo = "original";
    rg.isolamento_1[0] = r1.isolamento[0];
    rg.isolamento_1[1] = r1.isolamento[1];

    rg.solucao_1 = r1.solucao;

    if (rg.solucao_1) {
        rg.raiz_1 = r1.raiz;
        rg.num_iter_1 = r1.num_iter;
        rg.f_d_1 = r1.f_d;
    }

    if (qtd_raizes == 3) {
        resultado_individual r2 = newton_original(1, 2);
        resultado_individual r3 = newton_original(2, 4);

        rg.isolamento_2[0] = r2.isolamento[0];
        rg.isolamento_2[1] = r2.isolamento[1];
        rg.isolamento_3[0] = r3.isolamento[0];
        rg.isolamento_3[1] = r3.isolamento[1];

        rg.solucao_2 = r2.solucao;
        rg.solucao_3 = r3.solucao;

        if (rg.solucao_2) {
            rg.raiz_2 = r2.raiz;
            rg.num_iter_2 = r2.num_iter;
            rg.f_d_2 = r2.f_d;
        }

        if (rg.solucao_3) {
            rg.raiz_3 = r3.raiz;
            rg.num_iter_3 = r3.num_iter;
            rg.f_d_3 = r3.f_d;
        }
    }

    return rg;
}

resultado_individual newton_derivada(int i_x, int i_a) {

    double x = iniciais[i_x];

    int maxIteracoes = 100;
    int numIteracoes;

    bool solucao = false;

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
    resultado_individual res;
    res.solucao = solucao;
    if (solucao) {
        res.raiz = x;
        res.f_d = fpendulo(x);
        res.num_iter = numIteracoes;
    }

    res.isolamento[0] = isolamentos[i_a];
    res.isolamento[1] = isolamentos[i_a + 1];
    return res;
}

resultado_geral newton_derivada_geral() {

    resultado_individual r1 = newton_derivada(0, 0);

    resultado_geral rg;

    rg.a3 = a3;
    rg.a2 = a2;
    rg.metodo = "derivada num.";
    rg.isolamento_1[0] = r1.isolamento[0];
    rg.isolamento_1[1] = r1.isolamento[1];

    rg.solucao_1 = r1.solucao;

    if (rg.solucao_1) {
        rg.raiz_1 = r1.raiz;
        rg.num_iter_1 = r1.num_iter;
        rg.f_d_1 = r1.f_d;
    }

    if (qtd_raizes == 3) {
        resultado_individual r2 = newton_derivada(1, 2);
        resultado_individual r3 = newton_derivada(2, 4);

        rg.isolamento_2[0] = r2.isolamento[0];
        rg.isolamento_2[1] = r2.isolamento[1];
        rg.isolamento_3[0] = r3.isolamento[0];
        rg.isolamento_3[1] = r3.isolamento[1];

        rg.solucao_2 = r2.solucao;
        rg.solucao_3 = r3.solucao;

        if (rg.solucao_2) {
            rg.raiz_2 = r2.raiz;
            rg.num_iter_2 = r2.num_iter;
            rg.f_d_2 = r2.f_d;
        }

        if (rg.solucao_3) {
            rg.raiz_3 = r3.raiz;
            rg.num_iter_3 = r3.num_iter;
            rg.f_d_3 = r3.f_d;
        }
    }

    return rg;
}

resultado_individual newton_FL(int i_x, int i_a, double lambda_) {
    double x = iniciais[i_x];
    int maxIteracoes = 100;
    int numIteracoes;

    bool solucao = false;
    double xw;

    for (int i = 0; i < maxIteracoes; ++i) {
        double y = fpendulo(x);
        double dy = derivada_fpendulo(x);

        if (abs(y) < precisao) {
            solucao = true;
            numIteracoes = i + 1;
            break;
        }

        double FL;
        if (abs(dy) > lambda_) {
            xw = x;
            FL = dy;
        } else {
            FL = derivada_fpendulo(xw);
        }

        x = x - y / FL;
    }
    resultado_individual res;
    res.solucao = solucao;
    if (solucao) {
        res.raiz = x;
        res.f_d = fpendulo(x);
        res.num_iter = numIteracoes;
    }
    res.isolamento[0] = isolamentos[i_a];
    res.isolamento[1] = isolamentos[i_a + 1];
    return res;
}

resultado_geral newton_FL_geral(double _lambda) {

    resultado_individual r1 = newton_FL(0, 0, _lambda);

    resultado_geral rg;

    rg.lambda = lambda_list.front();

    rg.a3 = a3;
    rg.a2 = a2;
    rg.metodo = "FL";
    rg.isolamento_1[0] = r1.isolamento[0];
    rg.isolamento_1[1] = r1.isolamento[1];

    rg.solucao_1 = r1.solucao;

    if (rg.solucao_1) {
        rg.raiz_1 = r1.raiz;
        rg.num_iter_1 = r1.num_iter;
        rg.f_d_1 = r1.f_d;
    }

    if (qtd_raizes == 3) {
        resultado_individual r2 = newton_FL(1, 2, _lambda);
        resultado_individual r3 = newton_FL(2, 4, _lambda);

        rg.isolamento_2[0] = r2.isolamento[0];
        rg.isolamento_2[1] = r2.isolamento[1];
        rg.isolamento_3[0] = r3.isolamento[0];
        rg.isolamento_3[1] = r3.isolamento[1];

        rg.solucao_2 = r2.solucao;
        rg.solucao_3 = r3.solucao;

        if (rg.solucao_2) {
            rg.raiz_2 = r2.raiz;
            rg.num_iter_2 = r2.num_iter;
            rg.f_d_2 = r2.f_d;
        }

        if (rg.solucao_3) {
            rg.raiz_3 = r3.raiz;
            rg.num_iter_3 = r3.num_iter;
            rg.f_d_3 = r3.f_d;
        }
    }

    return rg;
}

void cabecalho_resultado() {
    ofstream myFile;
    myFile.open("resultados.csv");

    myFile << "id,a3,a2,metodo,lambda, qtd raizes reais,"
              "isolamento raiz 1,raiz 1,f(raiz1),num interacoes 1,"
              "isolamento raiz 2,raiz 2,f(raiz2),num interacoes 2,"
              "isolamento raiz 3,raiz 3,f(raiz3),num interacoes 3"
           << endl;
}

void registrar_resultado(const resultado_geral &rg) {
    ofstream myFile;
    myFile.open("resultados.csv", std::fstream::out | std::fstream::app);

    myFile << registro_id << ","
           << rg.a3 << ","
           << rg.a2 << ","
           << rg.metodo << ",";

    if (rg.metodo == "FL") {
        myFile << rg.lambda << ",";
    } else {
        myFile << "" << ",";
    }

    myFile << qtd_raizes << ",";

    myFile << rg.isolamento_1[0] << " ~ " << rg.isolamento_1[1] << ",";

    if (rg.solucao_1) {
        myFile << rg.raiz_1 << "," << rg.f_d_1 << "," << rg.num_iter_1 << ",";
    } else {
        myFile << "" << "," << "" << "," << "" << ",";
    }

    if (qtd_raizes == 3) {
        myFile << rg.isolamento_2[0] << " ~ " << rg.isolamento_2[1] << ",";

        if (rg.solucao_2) {
            myFile << rg.raiz_2 << "," << rg.f_d_2 << "," << rg.num_iter_2 << ",";
        } else {
            myFile << "" << "," << "" << "," << "" << ",";
        }

        myFile << rg.isolamento_3[0] << " ~ " << rg.isolamento_3[1] << ",";

        if (rg.solucao_3) {
            myFile << rg.raiz_3 << "," << rg.f_d_3 << "," << rg.num_iter_3;
        } else {
            myFile << "" << "," << "" << "," << "";
        }
    } else {
        myFile << "" << "," << "" << "," << "" << "," << "" << "," << "" << "," << "" << "," << "" << ",";
    }

    myFile << endl;
    myFile.close();
    registro_id++;
}

void encontrar_raizes() {
    for (int i = 0; i < n_parametros; ++i) {
        a3 = a3_list.front();
        a2 = a2_list.front();

        gerar_vetor_isolamentos(isolamentos);
        gerar_vetor_pontos_iniciais(iniciais, isolamentos);

        resultado_geral original = newton_original_geral();
        registrar_resultado(original);

        resultado_geral derivada_numerica = newton_derivada_geral();
        registrar_resultado(derivada_numerica);

        for (int j = 0; j < n_lambda; ++j) {
            double lambda_front = lambda_list.front();
            resultado_geral FL = newton_FL_geral(lambda_front);

            registrar_resultado(FL);

            lambda_list.pop_front();
            lambda_list.push_back(lambda_front);
        }
        a3_list.pop_front();
        a2_list.pop_front();
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
}

void output() {
    system("clear");
    print_cabecalho();

    cout << "A tabela resultados.csv contendo todos os resultados foi criada no mesmo diretorio "
            "deste programa\n";
    cout << "Pressione qualquer tecla para fechar";
}

int main() {

    input();
    cabecalho_resultado();
    encontrar_raizes();
    output();

    return 0;
}


