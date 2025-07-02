#include <bits/stdc++.h>
using namespace std;

static const int MAXN = 50000;
using BS = bitset<MAXN>;

// Count triangles
long long count_triangles(const vector<BS>& adjbit, int n) {
    long long cnt = 0;
    for (int x = 0; x < n; ++x) {
        for (int y = x + 1; y < n; ++y) {
            if (!adjbit[x].test(y)) continue;
            BS common = adjbit[x] & adjbit[y];
            for (int z = common._Find_first(); z < n; z = common._Find_next(z)) {
                if (z > y) ++cnt;
            }
        }
    }
    return cnt;
}

// Count 4‑cliques
long long count_4cliques(const vector<BS>& adjbit, int n) {
    long long cnt = 0;
    vector<int> cn;
    for (int x = 0; x < n; ++x) {
        for (int y = x + 1; y < n; ++y) {
            if (!adjbit[x].test(y)) continue;
            BS common = adjbit[x] & adjbit[y];
            cn.clear();
            for (int z = common._Find_first(); z < n; z = common._Find_next(z)) {
                if (z > y) cn.push_back(z);
            }
            int sz = cn.size();
            for (int i = 0; i < sz; ++i)
                for (int j = i + 1; j < sz; ++j)
                    if (adjbit[cn[i]].test(cn[j])) ++cnt;
        }
    }
    return cnt;
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    cout << "== Motif Counter Starting ==" << endl;

    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <hppi.edgelist> [num_randomizations]" << endl;
        return 0;
    }

    string fname = argv[1];
    int R = (argc >= 3 ? stoi(argv[2]) : 3);

    cout << "Reading from file: " << fname << endl;

    // Read edges (dedupe & self-loop filter)
    set<pair<int,int>> edge_set;
    int u, v, max_id = -1;
    ifstream fin(fname);
    if (!fin) {
        cout << "[Warning] Could not open \"" << fname << "\". Falling back to stdin." << endl;
    }
    istream &in = fin.is_open() ? fin : cin;

    int line = 0;
    while (in >> u >> v) {
        line++;
        u--; v--;
        if (u < 0 || v < 0) continue;
        if (u == v) continue;
        if (u > v) swap(u, v);
        edge_set.emplace(u, v);
        max_id = max(max_id, v);
        if (line <= 5) 
            cout << "  input edge (0-based): (" << u << "," << v << ")" << endl;
    }
    cout << "Total unique edges: " << edge_set.size() << endl;

    if (max_id < 0) {
        cout << "[Error] No valid edges read." << endl;
        return 1;
    }
    int n = max_id + 1;
    cout << "Detected node count: " << n << endl;

    // Build adjacency + bitsets
    vector<vector<int>> adj(n);
    vector<BS> adjbit(n);
    for (auto &e : edge_set) {
        tie(u, v) = e;
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    for (int i = 0; i < n; ++i) {
        for (int w : adj[i]) adjbit[i].set(w);
    }
    cout << "Built adjacency lists and bitsets." << endl;

    // Count real motifs
    cout << "Counting real triangles..." << endl;
    long long real_tri = count_triangles(adjbit, n);
    cout << "Real triangles   = " << real_tri << endl;

    cout << "Counting real 4-cliques..." << endl;
    long long real_4cl = count_4cliques(adjbit, n);
    cout << "Real 4-cliques   = " << real_4cl << endl;

    // Prepare configuration model
    vector<int> degree(n);
    for (int i = 0; i < n; ++i) degree[i] = adj[i].size();

    mt19937_64 rng(12345);
    vector<long long> tri_rand(R), cl4_rand(R);

    // Randomizations
    for (int r = 0; r < R; ++r) {
        cout << "\n-- Randomization " << (r+1) << " of " << R << " --" << endl;
        vector<int> stubs;
        for (int i = 0; i < n; ++i)
            for (int t = 0; t < degree[i]; ++t)
                stubs.push_back(i);
        shuffle(stubs.begin(), stubs.end(), rng);

        vector<BS> bits(n);
        int added = 0;
        for (int i = 0; i + 1 < (int)stubs.size() && added < (int)edge_set.size(); i += 2) {
            int a = stubs[i], b = stubs[i+1];
            if (a == b || bits[a].test(b)) continue;
            bits[a].set(b);
            bits[b].set(a);
            ++added;
        }

        adjbit = bits;
        cout << "Counting randomized triangles..." << endl;
        tri_rand[r] = count_triangles(adjbit, n);
        cout << "Random triangles = " << tri_rand[r] << endl;

        cout << "Counting randomized 4-cliques..." << endl;
        cl4_rand[r] = count_4cliques(adjbit, n);
        cout << "Random 4-cliques = " << cl4_rand[r] << endl;
    }

    // Compute Z‑scores
    auto compute_z = [&](const vector<long long>& vals, long long real) {
        double mu = accumulate(vals.begin(), vals.end(), 0.0) / R;
        double var = 0;
        for (auto t : vals) var += (t - mu)*(t - mu);
        double sd = sqrt(var / R);
        return (real - mu) / sd;
    };

    double z_tri = compute_z(tri_rand, real_tri);
    double z_4cl = compute_z(cl4_rand, real_4cl);

    cout << "\n=== Final Z-scores ===" << endl;
    cout << "Triangle Z-score = " << fixed << setprecision(2) << z_tri << endl;
    cout << "4-Clique Z-score = " << fixed << setprecision(2) << z_4cl << endl;
    cout << "== Done ==" << endl;

    return 0;
}
