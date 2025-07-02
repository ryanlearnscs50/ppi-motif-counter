#include <bits/stdc++.h>
using namespace std;

static const int MAXN = 50000;  // make sure MAXN >= number of proteins
using BS = bitset<MAXN>;

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <hppi.edgelist> [num_randomizations]\n";
        return 1;
    }

    const string fname = argv[1];
    const int R = (argc >= 3 ? stoi(argv[2]) : 3);

    // 1) Read edge list
    vector<pair<int,int>> edges;
    int u, v, max_id = 0;
    {
        ifstream fin(fname);
        while (fin >> u >> v) {
            max_id = max({max_id, u, v});
            edges.emplace_back(u, v);
        }
    }
    const int n = max_id + 1;
    const int m = edges.size();
    cout << "Loaded graph: " << n << " nodes, " << m << " edges.\n";

    // 2) Build adjacency & bitsets for real graph
    vector<vector<int>> adj(n);
    vector<BS> adjbit(n);
    for (auto &e : edges) {
        tie(u, v) = e;
        if (u == v) continue;
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    for (int i = 0; i < n; ++i)
        for (int w : adj[i])
            adjbit[i].set(w);

    // 3) Triangle counting function
    auto count_triangles = [&]() {
        long long cnt = 0;
        for (int x = 0; x < n; ++x) {
            for (int y : adj[x]) if (y > x) {
                BS common = adjbit[x] & adjbit[y];
                for (int z = common._Find_first(); z < n; z = common._Find_next(z))
                    if (z > y) ++cnt;
            }
        }
        return cnt;
    };

    // 4) Count real triangles
    long long real_tri = count_triangles();
    cout << "Real triangles = " << real_tri << "\n";

    // 5) Generate R random G(n,m) and count
    mt19937_64 rng(12345);
    vector<long long> tri_rand(R);

    for (int r = 0; r < R; ++r) {
        cout << "\n--- ER Randomization " << (r+1) << " of " << R << " ---\n";

        // build an empty bitset graph
        vector<BS> bits(n);
        int added = 0;
        while (added < m) {
            int a = rng() % n;
            int b = rng() % n;
            if (a != b && !bits[a].test(b)) {
                bits[a].set(b);
                bits[b].set(a);
                ++added;
            }
        }

        // rebuild adj & adjbit
        for (int i = 0; i < n; ++i) {
            adj[i].clear();
            for (int j = bits[i]._Find_first(); j < n; j = bits[i]._Find_next(j))
                adj[i].push_back(j);
            adjbit[i] = bits[i];
        }

        tri_rand[r] = count_triangles();
        cout << "Random triangles = " << tri_rand[r] << "\n";
    }

    // 6) Compute Z-score
    double sum = accumulate(tri_rand.begin(), tri_rand.end(), 0.0);
    double mu  = sum / R;
    double sq = 0;
    for (auto &t : tri_rand) sq += (t - mu)*(t - mu);
    double sd = sqrt(sq / R);

    cout << fixed << setprecision(2)
         << "\nTriangle Z-score = " << (real_tri - mu) / sd << "\n";

    return 0;
}
