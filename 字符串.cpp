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



