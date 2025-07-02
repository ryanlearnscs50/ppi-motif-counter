#include <bits/stdc++.h>
using namespace std;

static const int MAXN = 50000;  // ensure this ≥ number of proteins
using BS = bitset<MAXN>;

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0]
             << " <hppi.edgelist> [num_randomizations]\n";
        return 1;
    }

    string fname = argv[1];
    int R = (argc >= 3 ? stoi(argv[2]) : 3);

    // Read edge list
    vector<pair<int,int>> edges;
    int u, v, max_id = 0;
    {
        ifstream fin(fname);
        while (fin >> u >> v) {
            max_id = max(max_id, max(u, v));
            edges.emplace_back(u, v);
        }
    }
    int n = max_id + 1;
    int m = edges.size();
    cout << "Loaded graph: " << n << " nodes, " << m << " edges.\n";

    // Build adjacency and bitsets for the real graph
    vector<vector<int>> adj(n);
    vector<BS> adjbit(n);
    for (auto &e : edges) {
        tie(u, v) = e;
        if (u == v) continue;
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    for (int i = 0; i < n; ++i) {
        for (int w : adj[i]) {
            adjbit[i].set(w);
        }
    }

    // Function to count triangles
    auto count_triangles = [&]() {
        long long cnt = 0;
        for (int x = 0; x < n; ++x) {
            for (int y : adj[x]) {
                if (y <= x) continue;
                BS common = adjbit[x] & adjbit[y];
                for (int z = common._Find_first(); z < n; z = common._Find_next(z)) {
                    if (z > y) ++cnt;
                }
            }
        }
        return cnt;
    };

    // Count triangles in the real network
    long long real_tri = count_triangles();
    cout << "Real triangles = " << real_tri << "\n";

    // Precompute degree sequence for configuration model
    vector<int> degree(n);
    for (int i = 0; i < n; ++i) {
        degree[i] = adj[i].size();
    }

    mt19937_64 rng(12345);
    vector<long long> tri_rand(R);

    // Perform R randomizations using configuration model
    for (int r = 0; r < R; ++r) {
        cout << "\nRandomization " << (r + 1) << " of " << R << "\n";

        // Build stub list
        vector<int> stubs;
        stubs.reserve(2 * m);
        for (int i = 0; i < n; ++i) {
            for (int t = 0; t < degree[i]; ++t) {
                stubs.push_back(i);
            }
        }

        // Shuffle stubs
        shuffle(stubs.begin(), stubs.end(), rng);

        // Pair stubs to form edges
        vector<BS> bits(n);
        int added = 0;
        for (int k = 0; k + 1 < (int)stubs.size() && added < m; k += 2) {
            int a = stubs[k];
            int b = stubs[k + 1];
            if (a != b && !bits[a].test(b)) {
                bits[a].set(b);
                bits[b].set(a);
                ++added;
            }
        }

        // Rebuild adjacency and bitsets from bits
        for (int i = 0; i < n; ++i) {
            adj[i].clear();
            for (int j = bits[i]._Find_first(); j < n; j = bits[i]._Find_next(j)) {
                adj[i].push_back(j);
            }
            adjbit[i] = bits[i];
        }

        // Count triangles in randomized graph
        tri_rand[r] = count_triangles();
        cout << "Random triangles = " << tri_rand[r] << "\n";
    }

    // Compute Z-score
    double sum = accumulate(tri_rand.begin(), tri_rand.end(), 0.0);
    double mu  = sum / R;
    double sq  = 0.0;
    for (auto t : tri_rand) {
        sq += (t - mu) * (t - mu);
    }
    double sd = sqrt(sq / R);

    cout << fixed << setprecision(2)
         << "\nTriangle Z-score = " << (real_tri - mu) / sd << "\n";

    return 0;
}
