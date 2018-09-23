#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

const double EPS = 1e-5;

struct HeEx {
    int hot_num;
    int cold_num;

    double F;
    double k;
    double Q;

    double W1, W2;

    HeEx(int hot_num, int cold_num, double k, double F, double Q, double W1, double W2): hot_num(hot_num), cold_num(cold_num), k(k), F(F), Q(Q), W1(W1), W2(W2) { }
};

void balance(double hot_begin, double cold_end, HeEx &he_ex, double &hot_end, double &cold_begin) {
    hot_end = hot_begin - he_ex.W2 * he_ex.Q / (he_ex.k * he_ex.F);
    cold_begin = cold_end - he_ex.W1 * he_ex.Q / (he_ex.k * he_ex.F);

    double delta_t = ((hot_begin - hot_end) - (cold_end - cold_begin)) / (log((hot_begin - hot_end) / (cold_end - cold_begin)));

    he_ex.Q = he_ex.F * he_ex.k * delta_t;
}

void find_temperature(const vector<double> &cold_begin, const vector<double> &hot_begin,
                      const vector<HeEx> &he_exes, const vector<int> &order,
                      vector<double> &cold_end, vector<double> &hot_end)
{
    unsigned long n = cold_begin.size();
    unsigned long m = hot_begin.size();
    vector<double> delta_cold(n);
    vector<double> cold_now(n);
    cold_end.resize(n);
    hot_end.resize(n);
    do {
        for (int i = 0; i < n; ++i) {
            cold_end[i] = cold_begin[i] + delta_cold[i];
            cold_now[i] = cold_end[i];
        }
        for (int j = 0; j < m; ++j) {
            hot_end[j] = hot_begin[j];
        }

        for (HeEx he_ex : he_exes) {
            //HeEx he_ex = he_exes[i];

            int h_num = he_ex.hot_num;
            int c_num = he_ex.cold_num;

            balance(hot_end[h_num], cold_now[c_num], he_ex, hot_end[h_num], cold_now[c_num]);
        }

        for (int i = 0; i < n; ++i) {
            delta_cold[i] = cold_begin[i] - cold_now[i];
        }
    } while (*min_element(delta_cold.begin(), delta_cold.end()) > EPS);
}

int main() {
    vector<double> hot_ans, cold_ans;

    find_temperature({ 0, 0 }, { 100, 100 },
                     {
                             HeEx(1, 0, 1, 1, 30, 0.1, 0.2),
                             HeEx(1, 1, 1, 1, 30, 0.2, 0.1),
                             HeEx(0, 0, 1, 1, 30, 0.3, 0.2),
                             HeEx(0, 1, 1, 1, 30, 0.1, 0.1),
                             HeEx(1, 0, 1, 1, 30, 0.2, 0.4),
                             HeEx(1, 1, 1, 1, 30, 0.1, 0.4),
                             HeEx(0, 0, 1, 1, 30, 0.3, 0.1),
                             HeEx(0, 1, 1, 1, 30, 0.2, 0.2)
                     }, {7, 4, 5, 2, 0, 1},
    cold_ans, hot_ans);

    for (double cold_an : cold_ans) cout << cold_an << " ";
    cout << "\n";
    for (double hot_an : hot_ans) cout << hot_an << " ";
    cout << "\n";

    return 0;
}