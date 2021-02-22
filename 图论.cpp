int head[N], numE = 1;

struct E{
    int next, v, w;
} e[M << 1];

void add(int u, int v, int w) {
    e[++numE] = (E) { head[u], v, w };
    head[u] = numE;
}

// Dijkstra

typedef pair<LL, int> PII;
 
priority_queue<PII, vector<PII>, greater<PII> > q;
 
void dijkstra() {
    memset(d, 0x3f, sizeof d);
    d[1] = 0; q.push(make_pair(0, 1));
    while(!q.empty()) {
        PII u = q.top(); q.pop();
        if(vis[u.second]) continue;
        vis[u.second] = true;
        for (int i = head[u.second]; i; i = e[i].next) {
            int v = e[i].v;
            if(d[u.second] + e[i].w < d[v]) {
                d[v] = d[u.second] + e[i].w;
                q.push(make_pair(d[v], v));
            }
        }
    }
}

// Spfa

// Prufer
void inline fToP() {
	for (int i = 1; i < n; i++) d[f[i]]++;
	for (int i = 1, j = 1; i <= n - 2; j++) {
		while (d[j]) j++;
		p[i++] = f[j];
		while (i <= n - 2 && --d[p[i - 1]] == 0 && p[i - 1] < j) p[i++] = f[p[i - 1]];
	}
}

void inline pToF() {
    for (int i = 1; i <= n - 2; i++) d[p[i]]++;
    p[n - 1] = n;
    for (int i = 1, j = 1; i < n; i++, j++) {
    	while (d[j]) j++;
    	f[j] = p[i];
    	while (i < n - 1 && --d[p[i]] == 0 && p[i] < j) f[p[i]] = p[i + 1], ++i;
    }
}

// 有向图 tarjan

void tarjan(int u) {
    dfn[u] = low[u] = ++dfncnt;
    s[++top] = u, ins[u] = true;
    for (int i = head[u]; i; i = e[i].next) {
        int v = e[i].v;
        if (!dfn[v]) {
            tarjan(v), low[u] = min(low[u], low[v]);
        } else if (ins[v]) low[u] = min(low[u], dfn[v]);
    }
    if (low[u] == dfn[u]) {
        int v; ++cnt;
        do {
            v = s[top--], ins[v] = false, col[v] = cnt;
        } while (v != u);
    }
}


// 最大流
bool inline bfs() {
    int hh = 0, tt = -1;
    memset(d, 0, sizeof d);
    q[++tt] = s, d[s] = 1, cur[s] = head[s];
    while (hh <= tt) {
        int u = q[hh++];
        for (int i = head[u]; i; i = e[i].next) {
            int v = e[i].v;
            if (!d[v] && e[i].w) {
                cur[v] = head[v];
                q[++tt] = v, d[v] = d[u] + 1;
                if (v == t) return true;
            }
        }
    }
    return false;
}

LL inline dinic(int u, LL flow) {
    if (u == t) return flow;
    LL rest = flow;
    for (int i = cur[u]; i && rest; i = e[i].next) {
        cur[u] = i;
        int v = e[i].v;
        if (e[i].w && d[v] == d[u] + 1) {
            int k = dinic(v, min((LL)e[i].w, rest));
            if (!k) d[v] = 0;
            rest -= k, e[i].w -= k, e[i ^ 1].w += k;
        }
    }
    return flow - rest;
}

// ----重链剖分----

int sz[SZ], fa[SZ], dep[SZ], top[SZ], hson[SZ];

void dfs1(int u) {
    sz[u] = 1;
    for (int i = head[u]; i; i = e[i].next) {
        int v = e[i].v;
        if (v == fa[u]) continue;
        dep[v] = dep[u] + 1, fa[v] = u;
        dfs1(v);
        sz[u] += sz[v];
        if (sz[v] > sz[hson[u]]) hson[u] = v;
    }
}

void dfs2(int u, int tp) {
    top[u] = tp;
    if (hson[u]) dfs2(hson[u], tp);
    for (int i = head[u]; i; i = e[i].next) {
        int v = e[i].v;
        if (v == fa[u] || v == hson[u]) continue;
        dfs2(v, v);
    }
}

int lca(int x, int y) {
    while (top[x] != top[y]) {
        if (dep[top[x]] < dep[top[y]]) swap(x, y);
        x = fa[top[x]];
    }
    if (dep[x] < dep[y]) swap(x, y);
    return y;
}


// ---- End 重链剖分----


// Start ：最小树形图

int rt = 1, col, in[N];

int vis[N], id[N], pre[N];

struct E{
    int u, v, w;
} e[M];

int inline edmonds() {
    int ans = 0;
    while (true) {
        for (int i = 1; i <= n; i++) in[i] = INF;
        memset(vis, 0, sizeof vis);
        memset(id, 0, sizeof id);
        for (int i = 1; i <= m; i++) 
            if (e[i].w < in[e[i].v]) in[e[i].v] = e[i].w, pre[e[i].v] = e[i].u;
        for (int i = 1; i <= n; i++)
            if (in[i] == INF && i != rt) return -1;
        col = 0;
        for (int i = 1; i <= n; i++) {
            if (i == rt) continue;
            ans += in[i];
            int v = i;
            while (!vis[v] && !id[v] && v != rt)
                vis[v] = i, v = pre[v];
            if (v != rt && vis[v] == i) {
                id[v] = ++col;
                for (int x = pre[v]; x != v; x = pre[x]) id[x] = col;
            }
        }
        if (!col) break;
        for (int i = 1; i <= n; i++) if (!id[i]) id[i] = ++col;
        int tot = 0;
        for (int i = 1; i <= m; i++) {
            int a = id[e[i].u], b = id[e[i].v];
            if (a == b) continue;
            e[++tot] = (E) { a, b, e[i].w - in[e[i].v] };
        }
        m = tot, n = col, rt = id[rt];
    }
    return ans;
}

// -- End

// Start : 长链剖分 + O(1) k 级祖先

int d[N], dep[N];
int g[N], son[N], fa[N][L], top[N];
LL res;

vector<int> U[N], D[N];


void dfs1(int u) {
    dep[u] = d[u] = d[fa[u][0]] + 1;
    for (int i = 1; fa[u][i - 1]; i++) fa[u][i] = fa[fa[u][i - 1]][i - 1];
    for (int i = head[u]; i; i = e[i].next) {
        int v = e[i].v;
        dfs1(v);
        if (dep[v] > dep[u]) dep[u] = dep[v], son[u] = v;
    }
} 

void dfs2(int u, int tp) {
    top[u] = tp;
    if (u == tp) {
        for (int x = u, i = 0; i <= dep[u] - d[u]; i++)
            U[u].push_back(x), x = fa[x][0];
        for (int x = u, i = 0; i <= dep[u] - d[u]; i++)
            D[u].push_back(x), x = son[x];
    }
    if (son[u]) dfs2(son[u], tp);
    for (int i = head[u]; i; i = e[i].next) {
        int v = e[i].v;
        if (v != son[u]) dfs2(v, v);
    }
}

int inline query(int x, int k) {
    if (!k) return x;
    x = fa[x][g[k]], k -= (1 << g[k]) + d[x] - d[top[x]], x = top[x];
    return k < 0 ? D[x][-k] : U[x][k];
}


// --End

// 最小费用最大流 EK

const int N = ?, M = ?;

const int INF = 0x3f3f3f3f;

int n, m, s, t, maxflow, cost, d[N], incf[N], pre[N];

int q[N];

int head[N], numE = 1;

bool vis[N];

struct E{
    int next, v, w, c;
} e[M];

void inline add(int u, int v, int w, int c) {
    e[++numE] = (E) { head[u], v, w, c };
    head[u] = numE;
}

// Spfa || 
bool spfa() {
    memset(vis, false, sizeof vis);
    memset(d, 0x3f, sizeof d);
    int hh = 0, tt = 1;
    q[0] = s; d[s] = 0; incf[s] = 2e9;
    while (hh != tt) {
        int u = q[hh++]; vis[u] = false;
        if (hh == N) hh = 0;
        for (int i = head[u]; i; i = e[i].next) {
            int v = e[i].v;
            if (e[i].w && d[u] + e[i].c < d[v]) {
                d[v] = d[u] + e[i].c;
                pre[v] = i;
                incf[v] = min(incf[u], e[i].w);
                if (!vis[v]) {
                    q[tt++] = v;
                    vis[v] = true;
                    if (tt == N) tt = 0;
                }
            }
        }
    } 
    return d[t] != INF;
}

void update() {
    int x = t;
    while (x != s) {
        int i = pre[x];
        e[i].w -= incf[t], e[i ^ 1].w += incf[t];
        x = e[i ^ 1].v;
    }
    maxflow += incf[t];
    cost += d[t] * incf[t];
}

// --End

// 匈牙利

int match[N];
bool vis[N];

bool find(int u) {
    for (int i = head[u]; i; i = e[i].next) {
        int v = e[i].v;
        if (vis[v]) continue;
        vis[v] = true;
        if (!match[v] || find(match[v])) {
            match[v] = u; return true;
        }
    }
    return false;
}
 
// End

// 点分治

void getRoot(int u, int last) {
    sz[u] = 1;
    int s = 0;
    for (int i = head[u]; i; i = e[i].next) {
        int v = e[i].v;
        if (vis[v] || v == last) continue;
        getRoot(v, u);
        sz[u] += sz[v];
        s = max(s, sz[v]);
    }
    s = max(s, S - sz[u]);
    if (s < maxPart)
        maxPart = s, rt = u;
}

void dfs(int u, int last, int col, int w, int dep) {
    maxDep = max(maxDep, dep);
    if (L <= dep && dep <= R)
        ans = max(ans, w);
    d[dep] = max(d[dep], w);
    for (int i = head[u]; i; i = e[i].next) {
        int v = e[i].v;
        if (v == last || vis[v])
            continue;
        dfs(v, u, e[i].w, w + (col == e[i].w ? 0 : c[e[i].w]), dep + 1);
    }
}

void dfs0(int u, int last, int dep) {
    maxDep = max(maxDep, dep);
    for (int i = head[u]; i; i = e[i].next) {
        int v = e[i].v;
        if (v == last || vis[v])
            continue;
        dfs0(v, u, dep + 1);
    }
}

int inline work(int a[], int len1, int b[], int len2) {
    int res = -INF;

    int hh = 0, tt = -1;
    int l = max(1, L - 1), r = R - 1;
    if (l > r)
        return res;
    len1 = min(len1, R - 1);
    for (int i = min(r, len2); i >= l; i--) {
        while (hh <= tt && b[q[tt]] < b[i]) tt--;
        q[++tt] = i;
    }
    if (hh <= tt)
        res = max(res, a[1] + b[q[hh]]);
    for (int i = 2; i <= len1; i++) {
        if (q[hh] == r)
            hh++;
        r--;
        if (l > 1) {
            --l;
            while (hh <= tt && b[q[tt]] < b[l]) tt--;
            q[++tt] = l;
        }
        if (hh <= tt)
            res = max(res, a[i] + b[q[hh]]);
    }

    return res;
}

void solve(int x) {
    if (S == 1)
        return;
    maxPart = 2e9, getRoot(x, 0), vis[rt] = true;
    tot = 0;
    for (int i = head[rt]; i; i = e[i].next) {
        int v = e[i].v;
        if (vis[v])
            continue;
        maxDep = 0;
        dfs0(v, rt, 1);
        mxDep[e[i].w] = max(mxDep[e[i].w], maxDep);
        sons[++tot] = (Son){ v, e[i].w, maxDep };
    }
    sort(sons + 1, sons + 1 + tot);
    zDep = 0;
    int nowDep = 0;
    for (int i = 1; i <= tot; i++) {
        maxDep = 0, now = c[sons[i].c];
        mxDep[e[i].w] = 0;
        dfs(sons[i].x, rt, sons[i].c, c[sons[i].c], 1);
        zDep = max(zDep, maxDep);
        nowDep = max(nowDep, maxDep);
        ans = max(ans, work(d, maxDep, nowVal, nowDep) - c[sons[i].c]);
        for (int j = 1; j <= maxDep; j++) {
            nowVal[j] = max(nowVal[j], d[j]);
            d[j] = -INF;
        }
        if (i == tot || sons[i].c != sons[i + 1].c) {
            ans = max(ans, work(nowVal, nowDep, val, zDep));
            for (int j = 1; j <= nowDep; j++) {
                val[j] = max(val[j], nowVal[j]);
                nowVal[j] = -INF;
            }
            nowDep = 0;
        }
    }
    for (int i = 1; i <= zDep; i++) val[i] = -INF;
    for (int i = head[rt]; i; i = e[i].next) {
        int v = e[i].v;
        if (vis[v])
            continue;
        S = sz[v], solve(v);
    }
}

// End