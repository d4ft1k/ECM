#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include "lib.cpp"
using namespace std;

class Task : public lib::Task {
public:
    int n = 0, W = 0;
    vector<int> v, w;

    void load(istream& in) override {
        in >> n >> W;
        v.resize(n);
        w.resize(n);
        for (int i = 0; i < n; i++)
            in >> v[i] >> w[i];
    }

    double f(const vector<int>& x) override {
        int sumV = 0, sumW = 0;
        for (int i = 0; i < n; i++) {
            if (x[i] == 1) {
                sumV += v[i];
                sumW += w[i];
            }
        }
        if (sumW > W) return 0;
        return sumV;
    }

    vector<int> neighbour(const vector<int>& x0) override {
        vector<int> x = x0;
        int i = rand() % n;
        x[i] = 1 - x[i];
        return x;
    }

    vector<int> random() override {
        vector<int> x(n);
        int sumW = 0;
        for (int i = 0; i < n; i++) {
            x[i] = rand() % 2;
            if (x[i] == 1) sumW += w[i];
        }
        while (sumW > W) {
            int i = rand() % n;
            if (x[i] == 1) {
                x[i] = 0;
                sumW -= w[i];
            }
        }
        return x;
    }
};

static vector<int> sortedByDensity(Task& t) {
    vector<int> order(t.n);
    for (int i = 0; i < t.n; i++) order[i] = i;
    sort(order.begin(), order.end(), [&](int a, int b) {
        return (double)t.v[a] / t.w[a] > (double)t.v[b] / t.w[b];
    });
    return order;
}

class LocalSearch : public lib::Solver {
    Task& t;
public:
    explicit LocalSearch(Task& task) : t(task) {}

    vector<int> solve() override {
        vector<int> best = t.random();
        for (int it = 0; it < 100000; it++) {
            vector<int> y = t.neighbour(best);
            if (t.f(y) > t.f(best))
                best = y;
        }
        return best;
    }
};

class Greedy : public lib::Solver {
    Task& t;
public:
    explicit Greedy(Task& task) : t(task) {}

    vector<int> solve() override {
        vector<int> order = sortedByDensity(t);
        vector<int> x(t.n, 0);
        int sumW = 0;
        for (int i = 0; i < t.n; i++) {
            int id = order[i];
            if (sumW + t.w[id] <= t.W) {
                x[id] = 1;
                sumW += t.w[id];
            }
        }
        return x;
    }
};

class Branches : public lib::Solver {
    Task& t;
    vector<int> order, cur, best;
    int bestV = 0;

    double bound(int k, int curV, int curW) {
        double res = curV;
        int freeW = t.W - curW;
        for (int i = k; i < t.n; i++) {
            int id = order[i];
            if (t.w[id] <= freeW) {
                freeW -= t.w[id];
                res += t.v[id];
            } else {
                res += (double)t.v[id] * freeW / t.w[id];
                break;
            }
        }
        return res;
    }

    void branch(int k, int curV, int curW) {
        if (curV > bestV) {
            bestV = curV;
            best = cur;
        }
        if (k == t.n) return;
        if (bound(k, curV, curW) <= bestV) return;

        int id = order[k];
        if (curW + t.w[id] <= t.W) {
            cur[id] = 1;
            branch(k + 1, curV + t.v[id], curW + t.w[id]);
            cur[id] = 0;
        }
        branch(k + 1, curV, curW);
    }

public:
    explicit Branches(Task& task) : t(task) {}

    vector<int> solve() override {
        order = sortedByDensity(t);
        cur.assign(t.n, 0);
        best.assign(t.n, 0);
        bestV = 0;
        branch(0, 0, 0);
        return best;
    }
};

class Genetic : public lib::Solver {
    Task& t;

    vector<int> bestOf(const vector<vector<int>>& pop) {
        vector<int> best = pop[0];
        for (size_t i = 1; i < pop.size(); i++)
            if (t.f(pop[i]) > t.f(best))
                best = pop[i];
        return best;
    }

public:
    explicit Genetic(Task& task) : t(task) {}

    vector<int> solve() override {
        int size = 50;
        vector<vector<int>> pop;
        for (int i = 0; i < size; i++)
            pop.push_back(t.random());

        for (int gen = 0; gen < 300; gen++) {
            vector<vector<int>> next;
            next.push_back(bestOf(pop));
            while ((int)next.size() < size) {
                vector<int> p1 = pop[rand() % size];
                vector<int> p1b = pop[rand() % size];
                if (t.f(p1b) > t.f(p1)) p1 = p1b;
                vector<int> p2 = pop[rand() % size];
                vector<int> p2b = pop[rand() % size];
                if (t.f(p2b) > t.f(p2)) p2 = p2b;

                int cut = rand() % t.n;
                vector<int> child = p1;
                for (int i = cut; i < t.n; i++)
                    child[i] = p2[i];

                if (rand() % 100 < 20)
                    child = t.neighbour(child);

                next.push_back(child);
            }
            pop = next;
        }
        return bestOf(pop);
    }
};

int main() {
    srand(time(0));
    Task t;
    t.load(cin);
    cout << t.f(LocalSearch(t).solve()) << endl;
    cout << t.f(Greedy(t).solve()) << endl;
    cout << t.f(Branches(t).solve()) << endl;
    cout << t.f(Genetic(t).solve()) << endl;
    return 0;
}
