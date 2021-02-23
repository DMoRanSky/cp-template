// ST 表
struct ST{
	void inline STPrework(int n) {
		g[0] = -1;
		for (int i = 1; i <= n; i++) 
			f[i][0] = a[i], g[i] = g[i >> 1] + 1;
		for (int j = 1; j <= g[n]; j++)
			for (int i = 1; i + (1 << j) - 1 <= n; i++)
				f[i][j] = max(f[i][j - 1], f[i + (1 << (j - 1))][j - 1]);
	}

	int inline query(int l, int r) {
		int k = g[r - l + 1];
		return max(f[l][k], f[r - (1 << k) + 1][k]);
	}
}

// 线性基
struct Linear{
	int idx = 1, tr[SZ * 32][2], cnt[SZ * 32];
	void insert(LL x) {
		int p = 1;
		for (int i = 31; ~i; i--) {
			int ch = x >> i & 1;
			if (!tr[p][ch]) tr[p][ch] = ++idx;
			p = tr[p][ch], cnt[p]++;
		}
	}

	LL query(LL x, int k) {
		int p = 1; LL res = 0;
		for (int i = 31; ~i; i--) {
			int ch = x >> i & 1;
			if (k <= cnt[tr[p][!ch]]) res |= 1ll << i, p = tr[p][!ch];
			else k -= cnt[tr[p][!ch]], p = tr[p][ch];
		}
		return res;
	}
}

// 普通线段树

struct Seg{
	void inline pushup(int p) {
		
	}

	void inline pushdown(int p, int l, int r, int mid) {
	
	}

	void build(int p, int l, int r) {
		if(l == r) { 
			return; 
		}
		int mid = (l + r) >> 1;
	    build(p << 1, l, mid);
	    build(p << 1 | 1, mid + 1, r);
	    pushup(p);
	}

	void change(int p, int l, int r, int x, int y, int k, int c) {
	    if(x <= l && r <= y) {
	        return ;
	    }
		int mid = (l + r) >> 1;
		pushdown(p, l, r, mid);
		if(x <= mid) change(p << 1, l, mid, x, y, k, c);
		if(mid + 1 <= y) change(p << 1 | 1, mid + 1, r, x, y, k, c);
		pushup(p);
	}

	int query(int p, int l, int r, int x, int y) {
		if(x <= l && r <= y) return ?;
		int mid = (l + r) >> 1, s = 0;
		pushdown(p, l, r, mid);
		if(x <= mid) s += query(p << 1, l, mid, x, y);
		if(mid + 1 <= y) s += query(p << 1 | 1, mid + 1, r, x, y);
		return s % P;
	}
}

// 用来动态开点的池
struct T{
	int l, r, val, rnd, sz;
} t[SZ];
int idx;


struct Fhq{
	int rt;
	void pushup(int p) {
		
	}
	// value(A) < value(B)
	int merge(int A, int B) {
		if (!A || !B) return A + B;
		else if(t[A].rnd > t[B].rnd) {
			t[A].r = merge(t[A].r, B);
			pushup(A);
			return A;
		} else {
			t[B].l = merge(A, t[B].l);
			pushup(B);
			return B;
		}
	}

	// 按值分裂
	void split(int p, int k, int &x, int &y) {
		if (!p) x = y = 0;
		else {
			if (t[p].val <= k) 
			x = p, split(t[p].r, k, t[p].r, y);
			else y = p, split(t[p].l, k, x, t[p].l);
			pushup(p);
		}
	}
	int getNode(int val) {
		t[++idx] = (T) { 0, 0, val, rand(), 1 };
		return idx;
	}

	void insert(int val) {
		int x, y;
		split(rt, val, x, y);
		rt = merge(merge(x, getNode(val)), y);
	}

	int get(int l, int r) {
		int x, y, z;
		split(rt, l - 1, x, y);
		split(y, r, y, z);
		int res = t[y].N;
		rt = merge(x, merge(y, z));
		return res;
	}

	void del(int val) {
		int x, y, z;
		split(rt, val - 1, x, y);
		split(y, val, y, z);
		y = merge(t[y].l, t[y].r);
		rt = merge(x, merge(y, z));
	}
}

struct LCT{
	int ch[N][2], fa[N], mx[N], w[N], rev[N];

	void inline pushup(int p) {
		
	}

	void inline pushdown(int p) {
		if (rev[p]) { swap(ls, rs), rev[ls] ^= 1, rev[rs] ^= 1, rev[p] = 0; }
	}

	void inline rotate(int x) {
		int y = fa[x], z = fa[y], k = get(x);
		if (!isRoot(y)) ch[z][get(y)] = x;
		ch[y][k] = ch[x][!k], fa[ch[y][k]] = y;
		ch[x][!k] = y, fa[y] = x, fa[x] = z;
		pushup(y); pushup(x);
	}

	void inline update(int p) {
		if (!isRoot(p)) update(fa[p]);
		pushdown(p);
	}

	void inline splay(int p) {
		update(p);
		for (int f = fa[p]; !isRoot(p); rotate(p), f = fa[p]) 
			if (!isRoot(f)) rotate(get(p) == get(f) ? f : p);
	}

	void inline access(int x) {
		for (int p = 0; x; p = x, x = fa[x]) {
			splay(x), ch[x][1] = p, pushup(x);
		}
	}

	int inline find(int p) {
		access(p), splay(p);
		while (ls) pushdown(p), p = ls;
		splay(p);
		return p;
	}

	void inline makeRoot(int x) {
		access(x), splay(x), rev[x] ^= 1;
	}

	void inline split(int x, int y) {
		makeRoot(x), access(y), splay(y);
	}

	void inline link(int x, int y) {
		makeRoot(x), fa[x] = y;
	}

	void inline cut(int x, int y) {
		split(x, y);
		ch[y][0] = 0, fa[x] = 0;
		pushup(y);
	}

}

// 主席树
struct PersisSeg{
	struct T{
		int l, r;
		LL v;
	} t[SZ];

	int rt[SZ], idx;

	void inline update(int &p, int q, int l, int r, int x, int k) {
		t[p = ++idx] = t[q];
		t[p].v += k;
		if (l == r) return;
		int mid = (l + r) >> 1;
		if (x <= mid) update(t[p].l, t[q].l, l, mid, x, k);
		else update(t[p].r, t[q].r, mid + 1, r, x, k);
	}

	LL inline query(int p, int l, int r, int x, int y) {
		if (!p || x > y) return 0;
		if (x <= l && r <= y) return t[p].v;
		int mid = (l + r) >> 1; LL res = 0;
		if (x <= mid) res += query(t[p].l, l, mid, x, y);
		if (mid < y) res += query(t[p].r, mid + 1, r, x, y);
		return res;
	}
}


// 并查集
struct DSU{
	int f[N], sz[N];
	void init(int n) { for (int i = 1; i <= n; i++) f[i] = i, sz[i] = 1; }
	int inline find(int x) { return f[x] == x ? x : f[x] = find(f[x]); }
	void inline merge(int x, int y) {
		x = find(x), y = find(y);
		if (x == y) return;
		if (sz[x] > sz[y]) swap(x, y);
		sz[y] += sz[x], f[x] = y;
	}
};

// 树状数组

struct BIT{
	int n;
	LL c[SZ];
	void inline init(int len, LL a[]) {
		n = len;
		for (int i = 1; i <= n; i++) {
			c[i] += a[i];
			if (i + (i & -i) <= n) c[i + (i & -i)] += c[i];
		}
	}
	void inline add(int x, LL k) {
		for (; x <= n; x += x & -x) c[x] += k;
	}
	LL inline ask(int x) {
		LL res = 0;
		for (; x; x -= x & -x) res += c[x];
		return res;
	}
} ;

// 区间加 区间查的树状数组
struct exBIT{
	BIT t1, t2;
	int n;
	void inline init(int len, int a[]) {
		n = len;
		for (int i = 1; i <= n; i++) 
			b[i] = a[i] - a[i - 1];
		t1.init(n, b);
		for (int i = 1; i <= n; i++) b[i] *= i;
		t2.init(n, b);
	}
	void inline add(int l, int r, LL c) {
		t1.add(l, c), t1.add(r + 1, -c);
		t2.add(l, c * l), t2.add(r + 1, -c * (r + 1));
	}
	LL inline ask(int x) {
		return (x + 1) * t1.ask(x) - t2.ask(x);
	}
	LL inline ask(int x, int y) { return ask(y) - ask(x - 1); }
};

// 左偏树

struct LeftistTree{
	struct T{
	    int l, r, v, d, f;
	    // l, r 表示左右儿子, v 表示值
	    // d 表示从当前节点到最近叶子节点的距离, f 表示当前节点的父亲
	} t[SZ];

	int find(int x) {
	    return t[x].f == x ? x : t[x].f = find(t[x].f);
	}

	int merge(int x, int y) { // 递归合并函数
	    if (!x || !y) return x + y;
	    if (t[x].v > t[y].v || (t[x].v == t[y].v && x > y)) swap(x, y);
	    rs = merge(rs, y);
	    if (t[ls].d < t[rs].d) swap(ls, rs);
	    t[x].d = t[rs].d + 1;
	    return x;
	}

	int work(int x, int y) { // 合并 x, y 两个堆。
	    if (x == y) return 0;
		if (!x || !y) return t[x + y].f = x + y;
	    if (t[x].v > t[y].v || (t[x].v == t[y].v && x > y)) swap(x, y);
	    t[x].f = t[y].f = x;
	    merge(x, y); return x;
	}

	void del(int x) {
	    t[x].f = work(ls, rs), t[x].v = -1;
	}
}

// ---
// 回文自动机
struct PAM{
	int n, ch[SZ][26], fail[SZ], len[SZ], sz[SZ], idx = -1, lastans, last;

	char s[SZ];

	int inline newNode(int x) {	len[++idx] = x; return idx; }
	int inline getFail(int x) {
		while (s[n - len[x] - 1] != s[n]) x = fail[x];
		return x;
	}

	int inline insert(char c) {
		int k = c - 'a';
		s[++n] = c;
		int p = getFail(last), x;
		if (!ch[p][k]) {
			x = newNode(len[p] + 2);
			fail[x] = ch[getFail(fail[p])][k];
			ch[p][k] = x, sz[x] = 1 + sz[fail[x]];
		} else x = ch[p][k];
		last = x;
		return sz[x];
	}

	void inline build() {
		newNode(0), newNode(-1);
		s[0] = '$', fail[0] = 1, last = 0;
	}
}

// 后缀自动机

struct SAM{
	struct T{
		int nxt[26], len, link;
	} t[SZ << 2];

	void inline extend(int c) {
		int x = ++idx, p = last; sz[x] = 1;
		t[x].len = t[last].len + 1;
		while (p && !t[p].nxt[c]) 
			t[p].nxt[c] = x, p = t[p].link;
		if (!p) t[x].link = 1;
		else {
			int q = t[p].nxt[c];
			if (t[p].len + 1 == t[q].len) t[x].link = q;
			else {
				int y = ++idx;
				t[y] = t[q], t[y].len = t[p].len + 1;
				while (p && t[p].nxt[c] == q)
					t[p].nxt[c] = y, p = t[p].link;
				t[q].link = t[x].link = y;
			}
		}
		last = x;
	}
}

// Cdq 分治

void cdq(int l, int r) {
	if (l == r) return;
	int mid = (l + r) >> 1;
	cdq(l, mid), cdq(mid + 1, r);
	// Do sth

}

//

// 莫队

int pos[N], L[N], R[N], t;

struct Q {
	int l, r, id;
	bool operator < (const Q &b) const {
		if (pos[l] != pos[b.l]) return pos[l] < pos[b.l];
		return r < b.r;
	}
} q[N];

t = sqrt(n);
for (int i = 1; i <= n; i++) {
	pos[i] = (i - 1) / t + 1;
	if (!L[pos[i]]) L[pos[i]] = i;
	R[pos[i]] = i;
}

sort(q + 1, q + 1 + m);

// 回滚

int l = 1, r = 0, last = -1;
for (int i = 1; i <= m; i++) {
	if (pos[q[i].l] == pos[q[i].r]) {
		// 块内暴力
		continue;
	}
	if (pos[q[i].l] != last) {
		// 新的左块
		res = 0, top = 0, r = R[pos[q[i].l]], l = r + 1;
		last = pos[q[i].l];
	}
	while (r < q[i].r) {
		++r;
		// insert r
	}
	int bl = l, tp = res; // 记录
	while (l > q[i].l) {
		--l;
		// insert l
	}
	// 恢复
	ans[q[i].id] = res; res = tp;
}

// End

// DLX 1精确覆盖

namespace DLX1{
	int n, m, U[N], D[N], L[N], R[N], idx, s[N], hh, tt, X[N], Y[N];

	int ans[M], top;

	void inline init() {
	    for (int i = 0; i <= m; i++)
	        L[i] = i - 1, R[i] = i + 1, U[i] = D[i] = i;
	    L[0] = m, R[m] = 0, idx = m;
	}

	void inline add(int x, int y) {
	    X[++idx] = x, Y[idx] = y, s[y]++; 
	    L[idx] = hh, R[idx] = tt, L[tt] = R[hh] = idx;
	    U[idx] = U[y], D[idx] = y, D[U[y]] = idx, U[y] = idx;
	    hh = idx;
	} 

	// 删除第 p 列

	void del(int p) {
	    L[R[p]] = L[p], R[L[p]] = R[p];
	    for (int i = D[p]; i != p; i = D[i]) {
	        for (int j = R[i]; j != i; j = R[j]) {
	            s[Y[j]]--, U[D[j]] = U[j], D[U[j]] = D[j];
	        }
	    }
	}

	void resume(int p) {
	    L[R[p]] = p, R[L[p]] = p;
	    for (int i = U[p]; i != p; i = U[i]) {
	        for (int j = L[i]; j != i; j = L[j]) {
	            s[Y[j]]++, U[D[j]] = j, D[U[j]] = j;
	        }
	    }
	}

	bool inline dfs() {
	    if (!R[0]) return true;
	    int p = R[0];
	    for (int i = R[0]; i; i = R[i])
	        if (s[i] < s[p]) p = i;
	    if (!s[p]) return false;
	    del(p);
	    for (int i = D[p]; i != p; i = D[i]) {
	        ans[++top] = X[i];
	        for (int j = R[i]; j != i; j = R[j]) del(Y[j]);
	        if (dfs()) return true;
	        for (int j = L[i]; j != i; j = L[j]) resume(Y[j]);
	        --top;
	    }
	    resume(p);
	    return false;
	}
	main::
	add(i, j) if (i, j) = 1
	dfs()
}

namespace DLX2{
	int n, m, U[N], D[N], L[N], R[N], idx, s[N], hh, tt, X[N], Y[N];

	int ans[M], top, dep, d[M];

	bool st[N];

	void inline init() {
	    for (int i = 0; i <= m; i++)
	        L[i] = i - 1, R[i] = i + 1, U[i] = D[i] = i;
	    L[0] = m, R[m] = 0, idx = m;
	}

	void inline add(int x, int y) {
	    X[++idx] = x, Y[idx] = y, s[y]++; 
	    L[idx] = hh, R[idx] = tt, L[tt] = R[hh] = idx;
	    U[idx] = U[y], D[idx] = y, D[U[y]] = idx, U[y] = idx;
	    hh = idx;
	} 

	// 删除第 p 列

	void inline del(int p) {
	    for (int i = D[p]; i != p; i = D[i])
	        L[R[i]] = L[i], R[L[i]] = R[i];
	}

	void resume(int p) {
	    for (int i = U[p]; i != p; i = U[i])
	        L[R[i]] = i, R[L[i]] = i;
	}

	int inline h() {
	    memset(st, false, sizeof st);
	    int cnt = 0;
	    for (int i = R[0]; i; i = R[i]) {
	        if (st[i]) continue;
	        cnt++;
	        for (int j = D[i]; j != i; j = D[j])
	            for (int k = R[j]; k != j; k = R[k]) st[Y[k]] = true;
	    }
	    return cnt;
	}

	bool inline dfs() {
	    if (top + h() > dep) return false;
	    if (!R[0]) return true;
	    int p = R[0];
	    for (int i = R[0]; i; i = R[i])
	        if (s[i] < s[p]) p = i;
	    if (!s[p]) return false;
	    for (int i = D[p]; i != p; i = D[i]) {
	        ans[++top] = X[i];
	        del(i);
	        for (int j = R[i]; j != i; j = R[j]) del(j);
	        if (dfs()) return true;
	        for (int j = L[i]; j != i; j = L[j]) resume(j);
	        resume(i);
	        --top;
	    }
	    return false;
	}

	main::
	add(i, j) if (i, j) = 1
	dep = 1;
    while(!dfs()) dep++;
    printf("%d\n", dep);
    for (int i = 1; i <= dep; i++) printf("%d ", ans[i]);

}