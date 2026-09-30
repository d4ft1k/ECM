#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <sstream>
#include <string>
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
    int n = 0, m = 0;
    vector<pair<int, int>> edges;

    void load(istream& in) override {
        in >> n >> m;
        edges.resize(m);
        bool zeroBased = false;
        for (auto& e : edges) {
            in >> e.first >> e.second;
            if (!in) exit(1);
            if (e.first == 0 || e.second == 0) zeroBased = true;
        }
        int shift = zeroBased ? 0 : 1;
        for (auto& e : edges) {
            e.first -= shift;
            e.second -= shift;
            if (e.first < 0 || e.second < 0 || e.first >= n || e.second >= n)
                exit(1);
        }
    }

    double f(const vector<int>& x) override {
        int conflicts = 0;
        for (auto& e : edges)
            if (x[e.first] == x[e.second]) conflicts++;
        set<int> colors(x.begin(), x.end());
        return (double)conflicts * n + colors.size();
    }

    vector<int> neighbour(const vector<int>& x0) override {
        vector<int> x = x0;
        x[rand() % n] = rand() % n;
        return x;
    }

    vector<int> random() override {
        vector<int> x(n);
        for (int i = 0; i < n; i++) x[i] = rand() % n;
        return x;
    }
};

class SatColoring : public lib::Solver {
    Task& t;
    const string cnfFile = "coloring.cnf";
    const string outFile = "coloring.out";

    int var(int v, int c, int k) { return v * k + c + 1; }

    void writeCnf(int k) {
        stringstream body;
        int clauses = 0;

        for (int v = 0; v < t.n; v++) {
            for (int c = 0; c < k; c++) body << var(v, c, k) << " ";
            body << "0\n";
            clauses++;
            for (int c1 = 0; c1 < k; c1++)
                for (int c2 = c1 + 1; c2 < k; c2++) {
                    body << -var(v, c1, k) << " " << -var(v, c2, k) << " 0\n";
                    clauses++;
                }
        }
        for (auto& e : t.edges)
            for (int c = 0; c < k; c++) {
                body << -var(e.first, c, k) << " " << -var(e.second, c, k) << " 0\n";
                clauses++;
            }
        if (!t.edges.empty() && k >= 2) {
            body << var(t.edges[0].first, 0, k) << " 0\n";
            body << var(t.edges[0].second, 1, k) << " 0\n";
            clauses += 2;
        }

        ofstream out(cnfFile);
        out << "p cnf " << t.n * k << " " << clauses << "\n" << body.str();
    }

    bool tryColor(int k, vector<int>& color) {
        writeCnf(k);
        remove(outFile.c_str());
        string cmd = "minisat " + cnfFile + " " + outFile + " " + string(NULL_DEVICE) + " 2>&1";
        system(cmd.c_str());

        ifstream in(outFile);
        if (!in.is_open()) exit(1);
        string status;
        in >> status;
        if (status != "SAT") return false;

        color.assign(t.n, -1);
        int lit;
        while (in >> lit && lit != 0) {
            if (lit < 0) continue;
            int v = (lit - 1) / k, c = (lit - 1) % k;
            color[v] = c;
        }
        return true;
    }

public:
    explicit SatColoring(Task& task) : t(task) {}

    vector<int> solve() override {
        vector<int> color;
        if (t.n == 0) return color;
        int start = t.edges.empty() ? 1 : 2;
        for (int k = start; k <= t.n; k++) {
            if (tryColor(k, color)) break;
        }
        remove(cnfFile.c_str());
        remove(outFile.c_str());
        return color;
    }
};

int main() {
    srand(time(0));
    Task t;
    t.load(cin);

    SatColoring solver(t);
    vector<int> color = solver.solve();

    set<int> used(color.begin(), color.end());
    cout << used.size();
    return 0;
}