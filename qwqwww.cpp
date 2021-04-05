// Skyqwq
#include <iostream>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cstring>
#define rint register int
#define pb push_back
using namespace std;

typedef long long LL;

const int N = 405, P = 998244353;

int n, m, K, bc[N][N], hg[N][N], fact[N], infact[N], wp[N];

bool ban[N][N];

int f[N][N], g[N][N];

void inline add(rint &x, rint y) {
    x += y;
    if (x >= P) x -= P;
}

struct E{
    int L, p;
    bool operator < (const E &b) const {
        return L < b.L;
    }
};

int inline power(int a, int b) {
    int res = 1;
    while (b) {
        if (b & 1) res = (LL)res * a % P;
        a = (LL)a * a % P;
        b >>= 1;
    }
    return res;
}
void inline factPrework(int n) {
    fact[0] = infact[0] = 1;
    for (int i = 1; i <= n; i++) fact[i] = (LL)fact[i - 1] * i % P;
    infact[n] = power(fact[n], P - 2);
    for (int i = n - 1; i; i--) infact[i] = infact[i + 1] * (i + 1ll) % P;
}


struct Te{
    int d, c;
    bool operator < (const Te b) const {
        return d < b.d;
    }
} e[N];

vector<E> h[N];

int inline get(int x) {
    return lower_bound(e + 1, e + 1 + n, (Te) { x, 0 }) - e;
}

int main() {
    scanf("%d%d%d", &n, &m, &K); factPrework(400);
    for (int i = 1; i <= n; i++) scanf("%d%d", &e[i].d, &e[i].c);
    sort(e + 1, e + 1 + n);
    for (int i = 1, l, r, p; i <= m; i++) {
        scanf("%d%d%d", &l, &r, &p);
        l = get(l), r = get(r + 1) - 1;
        h[r].pb((E) { l, p });
    }
    for (int i = 1; i <= n; i++) sort(h[i].begin(), h[i].end());
    for (int r = 1; r <= n; r++) {
        int s1 = 1;
        for (int l = 1, i = 0; l <= r; l++) {
            while (i < h[r].size() && h[r][i].L <= l) {
                s1 = (LL)s1 * (1ll + P - h[r][i].p) % P;
                i++;
            }
            bc[r][l] = s1;
            hg[r][l] = (1ll + P - s1) % P;
        }
        bc[r][0] = 1;
    }
    f[0][0] = 1;
    for (int i = 1; i <= n; i++) {
        memcpy(g, f, sizeof g);
        //memset(f, 0, sizeof f);
        for (int v = 1, s = e[i].c; v <= K; v++, s = (LL)s * e[i].c % P) {
            wp[v] = (LL)s * infact[v] % P;
        }
        for (rint j = 0; j <= n; j++) {
            rint t = j == 0 ? i : j;
            for (rint k = 0; k <= K; k++) {
                if (!g[j][k]) continue;
                //add(f[j][k], g[j][k]);
                rint o = K - k;
                for (rint v = 1; v <= o; v++) {
                    add(f[t][v + k], (LL)g[j][k] * wp[v] % P);
                }
            }
        }
        memcpy(g, f, sizeof g);
        memset(f, 0, sizeof f);
        for (int j = 0; j <= n; j++) {
            for (int k = 0; k <= K; k++) {
                if (!g[j][k]) continue;
                add(f[j][k], (LL)g[j][k] * bc[i][j] % P);
                add(f[0][k], (LL)g[j][k] * hg[i][j] % P);
            }
        }
    }
    printf("%lld\n", (LL)f[0][K] * fact[K] % P);
    return 0;
}