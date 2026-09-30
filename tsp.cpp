#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "lib.cpp"
using namespace std;

#ifdef _WIN32
#define NULL_DEVICE "> NUL"
#else
#define NULL_DEVICE "> /dev/null"
#endif

class Task : public lib::Task {
public:
    int n = 0;
    vector<double> X, Y;

    void load(istream& in) override {
        in >> n;
        X.resize(n);
        Y.resize(n);
        for (int i = 0; i < n; i++)
            in >> X[i] >> Y[i];
    }

    double f(const vector<int>& p) override {
        double sum = 0;
        for (int i = 0; i < n; i++) {
            int a = p[i];
            int b = p[(i + 1) % n];
            sum += sqrt((X[a] - X[b]) * (X[a] - X[b]) + (Y[a] - Y[b]) * (Y[a] - Y[b]));
        }
        return sum;
    }

    vector<int> neighbour(const vector<int>& p0) override {
        vector<int> p = p0;
        int i = rand() % n;
        int j = rand() % n;
        if (i > j) swap(i, j);
        reverse(p.begin() + i, p.begin() + j + 1);
        return p;
    }

    vector<int> random() override {
        vector<int> p(n);
        for (int i = 0; i < n; i++) p[i] = i;
        for (int i = n - 1; i > 0; i--) {
            int j = rand() % (i + 1);
            swap(p[i], p[j]);
        }
        return p;
    }
};

class LocalSearch : public lib::Solver {
    Task& t;
public:
    explicit LocalSearch(Task& task) : t(task) {}

    vector<int> solve() override {
        vector<int> best = t.random();
        double bestF = t.f(best);
        for (int it = 0; it < 200000; it++) {
            vector<int> y = t.neighbour(best);
            double yF = t.f(y);
            if (yF < bestF) {
                best = y;
                bestF = yF;
            }
        }
        return best;
    }
};

class SatTsp : public lib::Solver {
    Task& t;
    const string cnfFile = "tsp.cnf";
    const string outFile = "tsp.out";
    static const int MAXD = 50;

    vector<vector<int>> D;

    int varX(int v, int p) { return v * t.n + p + 1; }
    int varC(int e, int k, int B) { return t.n * t.n + (e - 1) * B + k; }

    void buildDistances() {
        int n = t.n;
        vector<vector<double>> d(n, vector<double>(n, 0));
        double mx = 0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                d[i][j] = sqrt((t.X[i] - t.X[j]) * (t.X[i] - t.X[j]) +
                               (t.Y[i] - t.Y[j]) * (t.Y[i] - t.Y[j]));
                mx = max(mx, d[i][j]);
            }
        double scale = mx > 0 ? MAXD / mx : 1;
        D.assign(n, vector<int>(n, 0));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) D[i][j] = (int)llround(d[i][j] * scale);
    }

    int cost(const vector<int>& p) {
        int s = 0;
        for (int i = 0; i < t.n; i++) s += D[p[i]][p[(i + 1) % t.n]];
        return s;
    }

    vector<int> nearestNeighbour() {
        int n = t.n;
        vector<int> p{0};
        vector<bool> used(n, false);
        used[0] = true;
        for (int i = 1; i < n; i++) {
            int cur = p.back(), best = -1;
            for (int j = 0; j < n; j++)
                if (!used[j] && (best == -1 || D[cur][j] < D[cur][best])) best = j;
            used[best] = true;
            p.push_back(best);
        }
        return p;
    }

    void writeCnf(int B) {
        int n = t.n;
        stringstream body;
        long long clauses = 0;

        for (int v = 0; v < n; v++) {
            for (int p = 0; p < n; p++) body << varX(v, p) << " ";
            body << "0\n";
            clauses++;
            for (int p1 = 0; p1 < n; p1++)
                for (int p2 = p1 + 1; p2 < n; p2++) {
                    body << -varX(v, p1) << " " << -varX(v, p2) << " 0\n";
                    clauses++;
                }
        }
        for (int p = 0; p < n; p++) {
            for (int v = 0; v < n; v++) body << varX(v, p) << " ";
            body << "0\n";
            clauses++;
            for (int v1 = 0; v1 < n; v1++)
                for (int v2 = v1 + 1; v2 < n; v2++) {
                    body << -varX(v1, p) << " " << -varX(v2, p) << " 0\n";
                    clauses++;
                }
        }
        body << varX(0, 0) << " 0\n";
        clauses++;

        for (int e = 1; e <= n; e++) {
            int a = e - 1, b = e % n;
            for (int u = 0; u < n; u++)
                for (int v = 0; v < n; v++) {
                    if (u == v) continue;
                    int d = D[u][v];
                    for (int k = 0; k <= B; k++) {
                        if (e == 1 && k > 0) break;
                        if (k == 0 && d == 0) continue;
                        body << -varX(u, a) << " " << -varX(v, b) << " ";
                        if (k > 0) body << -varC(e - 1, k, B) << " ";
                        if (k + d > B) body << "0\n";
                        else body << varC(e, k + d, B) << " 0\n";
                        clauses++;
                    }
                }
        }

        ofstream out(cnfFile);
        out << "p cnf " << n * n + n * B << " " << clauses << "\n" << body.str();
    }

    bool tryTour(int B, vector<int>& tour) {
        writeCnf(B);
        remove(outFile.c_str());
        string cmd = "minisat " + cnfFile + " " + outFile + " " + string(NULL_DEVICE) + " 2>&1";
        system(cmd.c_str());

        ifstream in(outFile);
        if (!in.is_open()) exit(1);
        string status;
        in >> status;
        if (status != "SAT") return false;

        tour.assign(t.n, -1);
        int lit;
        while (in >> lit && lit != 0) {
            if (lit < 0 || lit > t.n * t.n) continue;
            int v = (lit - 1) / t.n, p = (lit - 1) % t.n;
            tour[p] = v;
        }
        return true;
    }

public:
    explicit SatTsp(Task& task) : t(task) {}

    vector<int> solve() override {
        int n = t.n;
        if (n <= 3) {
            vector<int> p(n);
            for (int i = 0; i < n; i++) p[i] = i;
            return p;
        }
        buildDistances();

        vector<int> best = nearestNeighbour();
        int B = cost(best) - 1;
        while (B >= 0) {
            vector<int> tour;
            if (!tryTour(B, tour)) break;
            best = tour;
            B = cost(best) - 1;
        }
        remove(cnfFile.c_str());
        remove(outFile.c_str());
        return best;
    }
};


int main() {
    srand(time(0));
    Task t;
    t.load(cin);

    LocalSearch solver(t);
    vector<int> ans = solver.solve();
    cout << t.f(ans) << endl;

    SatTsp sat_solver(t);
    ans = sat_solver.solve();
    cout << t.f(ans);
    
    return 0;
}