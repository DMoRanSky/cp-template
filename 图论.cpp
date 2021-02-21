int head[N], numE = 1;

struct E{
    int next, v, w;
} e[M << 1];

void add(int u, int v, int w) {
    e[++numE] = (E) { head[u], v, w };
    head[u] = numE;
}

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

bool spfa(int s) {
    memset(d, 0x3f, sizeof d);
    d[1] = 0;
    int hh = 0, tt = 0;
    q[0] = s;
    
    while(!q.empty()){
        int u = q[hh++];
        vis[u] = false;
        for(int i = head[u]; i; i = e[i].next){
            int v = e[i].v;
            if(d[u] + e[i].d < d[v]){
                d[v] = d[u] + e[i].d;
                cnt[v] = cnt[u] + 1;
                if(cnt[v] >= n) return false;
                if(!vis[v]) vis[v] = true, q.push(v);
            }
        }
    }
    return true;
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