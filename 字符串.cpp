// 中间添加 #
char s[N], g[N];

void change() {
	n = strlen(s + 1) * 2;
	g[0] = 0;
	for (int i = 1; i <= n; i++) {
		if (i % 2) g[i] = 1;
		else g[i] = s[i >> 1];
	}
	g[++n] = 1, g[n + 1] = 2; 
	manacher();
}

void manacher() {
	int r = 0, mid = 0;
	for (int i = 1; i <= n; i++) {
		p[i] = i <= r ? min(r - i + 1, p[2 * mid - i]) : 1;
		while (g[i - p[i]] == g[i + p[i]]) ++p[i];
		if (i + p[i] - 1 > r) mid = i, r = i + p[i] - 1;
		ans = max(ans, p[i] - 1);
	}
}

struct KMP{
	int n, nxt[SZ];
	void inline build(char s[]) {
		n = strlen(s + 1);
		nxt[1] = 0;
		for (int i = 2, j = 0; i <= n; i++) {
			while (j && s[j + 1] != s[i]) j = nxt[j];
			if (s[j + 1] == s[i]) j++;
			nxt[i] = j;
		}
	}
	void inline match(char a[], int m) {
		for (int i = 1, j = 0; i <= m; i++) {
			while (j && s[j + 1] != a[i]) j = nxt[j];
			if (s[j + 1] == a[i]) j++;
			if (j == n) {
				j = nxt[j];
			}
		}
	}
} kmp;

struct SA{
	int rk[SZ], sa[SZ], cnt[SZ], oldrk[SZ], id[SZ], n, m, p, height[SZ];
	bool inline cmp(int i, int j, int k) {
		return oldrk[i] == oldrk[j] && oldrk[i + k] == oldrk[j + k];
	}
	void inline build(char s[]) {
		n = strlen(s + 1), m = 221;
		for (int i = 1; i <= n; i++) cnt[rk[i] = s[i]]++;
		for (int i = 1; i <= m; i++) cnt[i] += cnt[i - 1];
		for (int i = n; i; i--) sa[cnt[rk[i]]--] = i;
		for (int w = 1; w < n; w <<= 1, m = p) {
			p = 0;
			for (int i = n; i > n - w; i--) id[++p] = i;
			for (int i = 1; i <= n; i++)
				if (sa[i] > w) id[++p] = sa[i] - w;
			for (int i = 1; i <= m; i++) cnt[i] = 0;
			for (int i = 1; i <= n; i++) cnt[rk[i]]++, oldrk[i] = rk[i];
			for (int i = 1; i <= m; i++) cnt[i] += cnt[i - 1];
			for (int i = n; i; i--) sa[cnt[rk[id[i]]]--] = id[i];
			p = 0;
			for (int i = 1; i <= n; i++) {
				rk[sa[i]] = cmp(sa[i], sa[i - 1], w) ? p : ++p;
			}
			if (p == n) break;
		}
		for (int i = 1; i <= n; i++) {
			int j = sa[rk[i] - 1], k = max(0, height[rk[i - 1]] - 1);
			while (s[i + k] == s[j + k]) k++;
			height[rk[i]] = k;
		}
	}
};

// 切记复制一倍到后面, 最小表示法，返回开始下标
int inline minExp(int a[], int n) {
	int i = 1, j = 2;
	while (i <= n && j <= n) {
		int k;
		for (k = 0; k < n && a[i + k] == a[j + k]; k++);
		if (k == n) break;
		if (a[i + k] < a[j + k]) j += k + 1;
		else i += k + 1;
		if (i == j) i++;
	}
	return min(i, j);
}

typedef unsigned long long ULL;

// 哈希

struct Hash{
	int b, P, p[N], h[N];
	int inline get(int l, int r){
	    return (h[r] - (LL)h[l - 1] * p[r - l + 1] % P + P) % P;
	}
	void inline build(int n, int tb, int tp) {
		b = tb, P = tp;
		p[0] = 1;
	    for(int i = 1; i <= n; i++){
	        p[i] = (LL)p[i - 1] * b % P;
	        h[i] = ((LL)h[i - 1] * b + s[i]) % P;
	    }
	}
}

// Z 函数

z[1] = n;
for (int i = 2, r = 0, j = 0; i <= n; i++) {
	if (i <= r) z[i] = min(r - i + 1, z[i - j + 1]);
	while (i + z[i] <= n && a[i + z[i]] == a[1 + z[i]]) z[i]++;
	if (i + z[i] - 1 > r) r = i + z[i] - 1, j = i; 
}

for (int i = 1, r = 0, j = 0; i <= m; i++) {
	if (i <= r) p[i] = min(r - i + 1, z[i - j + 1]);
	while (i + p[i] <= m && b[i + p[i]] == a[1 + p[i]]) p[i]++;
	if (i + p[i] - 1 > r) r = i + p[i] - 1, j = i; 
}

// End

// AC 自动机

struct ACAutomation{
	int tr[SZ][26], nxt[SZ], idx, q[SZ];
	void inline insert(char s[]) {
		int p = 0;
		for (int j = 0; s[j]; j++) {
			int ch = s[j] - 'a';
			if(!tr[p][ch]) tr[p][ch] = ++idx;
			p = tr[p][ch];
		}
	}
	void build() {
		int hh = 0, tt = -1;
		for (int i = 0; i < 26; i++) 
			if (tr[0][i]) q[++tt] = tr[0][i];
		while (hh <= tt) {
			int u = q[hh++];
			for (int i = 0; i < 26; i++) {
				int v = tr[u][i];
				if (!v) tr[u][i] = tr[nxt[u]][i];
				else nxt[v] = tr[nxt[u]][i], q[++tt] = v;
			}
		}
	}
}

// Runs

struct Runs{
	typedef unsigned long long ULL;

	const int N = 1e6 + 5;

	int n, tot, b1[N], b2[N];

	ULL H[N], P[N], B = 31;

	char s[N];

	struct Node {
		int a, b, c;
		bool operator < (const Node &y) const {
			if (a != y.a) return a < y.a;
			if (b != y.b) return b < y.b;
			return c < y.c;
		}
		bool operator == (const Node &y) const {
			return a == y.a && b == y.b && c == y.c;
		}
	} ans[N << 1];

	ULL inline get(int l, int r) {
		return H[r] - H[l - 1] * P[r - l + 1];
	}

	int inline lcp(int i, int j) {
		if (s[i] != s[j]) return 0;
		int l = 1, r = min(n - i + 1, n - j + 1);
		while (l < r) {
			int mid = (l + r + 1) >> 1;
			if (get(i, i + mid - 1) == get(j, j + mid - 1)) l = mid;
			else r = mid - 1;
		}
		return r;
	}

	int inline lcs(int i, int j) {
		if (s[i] != s[j]) return 0;
		int l = 1, r = min(i, j);
		while (l < r) {
			int mid = (l + r + 1) >> 1;
			if (get(i - mid + 1, i) == get(j - mid + 1, j)) l = mid;
			else r = mid - 1;
		}
		return r;
	}

	int inline cmp(int i, int j) {
		int k = lcp(i, j);
		return s[i + k] < s[j + k];
	}

	int inline cmp2(int i, int j) {
		int k = lcp(i, j);
		if (i + k == n + 1) return false;
		return s[i + k] < s[j + k];
	}


	void inline add(int i, int j) {
		int L = lcs(i - 1, j), len = j - i + 1;
		if (true) {
			int R = lcp(i, j + 1);
			if (L + R >= len) ans[++tot] = (Node) { i - L, j + R, len };
		}
	}

	void inline prework() {
		P[0] = 1;
		for (int i = 1; i <= n; i++) {
			P[i] = P[i - 1] * B;
			H[i] = (H[i - 1] * B + s[i] - 'a');
		}
	}

	void inline Runs() {
		for (int i = n; i; i--) {
			b1[i] = b2[i] = i + 1;
			while (b1[i] < n && cmp(i, b1[i])) b1[i] = b1[b1[i]];
			while (b2[i] < n && cmp2(b2[i], i)) b2[i] = b2[b2[i]];
			add(i, b1[i] - 1);
			if (b1[i] != b2[i]) add(i, b2[i] - 1);
		}
		sort(ans + 1, ans + 1 + tot);
		tot = unique(ans + 1, ans + 1 + tot) - ans - 1;
	}

}

struct SAM{
	int idx, last;
	struct SAM_{
		int nxt[26], len, link;
	} t[N];
	void inline init() {
		last = idx = 1;
	}
	
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
} t;

struct GSAM{
	int idx, last;
	struct SAM{
		int ch[26], len, link;
	} t[N];
	void inline init() {
		last = idx = 1;
	}
	void inline insert(int c) {
		int p = last;
		if (t[p].ch[c]) {
			int q = t[p].ch[c];
			if (t[q].len == t[p].len + 1) last = q;
			else {
				int y = ++idx; t[y] = t[q];
				t[y].len = t[p].len + 1;
				while (p && t[p].ch[c] == q)
					t[p].ch[c] = y, p = t[p].link;
				t[q].link = y;
				last = y;	
			}
			return;
		}
		int x = ++idx; t[x].len = t[p].len + 1;
		while (p && !t[p].ch[c]) t[p].ch[c] = x, p = t[p].link;
		int q, y;
		if (!p) t[x].link = 1;
		else {
			q = t[p].ch[c];
			if (t[q].len == t[p].len + 1) t[x].link = q;
			else {
				int y = ++idx; t[y] = t[q];
				t[y].len = t[p].len + 1;
				while (p && t[p].ch[c] == q)
					t[p].ch[c] = y, p = t[p].link;
				t[q].link = t[x].link = y;
				last = y;	
			}
		}
		last = x;
	}
} t;

// 回文自动机
struct PAM{
    int n, ch[N][26], fail[N], len[N], sz[N], idx = -1, last;
    char s[N];
    void clr() {
        n = 0;
        for (int i = 0; i <= idx; i++) {
            sz[i] = len[i] = fail[i] = 0;
            for (int j = 0; j < 26; j++)
                ch[i][j] = 0;
        }
        idx = -1;
        last = 0;
    }
 
    int newNode(int x) { len[++idx] = x; return idx; }
    int getFail(int x) {
        while (s[n - len[x] - 1] != s[n]) x = fail[x];
        return x;
    }
    int insert(char c) {
        int k = c - 'a';
        s[++n] = c;
        int p = getFail(last), x;
        if (!ch[p][k]) {
            x = newNode(len[p] + 2);
            fail[x] = ch[getFail(fail[p])][k];
            ch[p][k] = x, sz[x] = 1 + sz[fail[x]];
        } else x = ch[p][k];
        last = x;
        return x;
    }
    void bd() {
        // -1:idx jigen
        newNode(0), newNode(-1);
        s[0] = '$', fail[0] = 1, last = 0;
    }
} pam;

