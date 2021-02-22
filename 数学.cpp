LL inline exgcd(LL a, LL b, LL &x, LL &y) {
	if (b == 0) {
		x = 1, y = 0;
		return a;
	}
	LL d = exgcd(b, a % b, y, x);
	y -= a / b * x;
	return d;
}

LL inline mod(LL a, LL b) {
	return (a % b + b) % b;
}

LL inline mul(LL a, LL b, LL p) {
	return ((a * b - (LL)((long double)a * b / p) * p) % p + p) % p;
}

int inline CRT() {
	LL m = 1; int res = 0;
	for (int i = 1; i <= n; i++) M *= m[i];
	for (int i = 1; i <= n; i++) {
        LL x, y;
        LL d = exgcd(M / m[i], m[i], x, y);
        LL k = m[i] / d;
        x = (x % k + k) % k;
        res = (res + M / m[i] * a[i] * x) % M;
    }
    return res;
}

LL inline exCRT() {
	LL a1 = a[1], p1 = p[1];
	for (int i = 2; i <= n; i++) {
		LL x, y, t = a[i] - a1;
		LL d = exgcd(p1, p[i], x, y);
		if (t % d) return -1;
		x = mul(x, t / d, p[i] / d);
		LL k = p1 / d * p[i];
		a1 = mod(a1 + mul(x, p1, k), k);
		p1 = k;
	}
	LL t = max(0ll, (lim - a1 + p1 - 1) / p1);
	return a1 + t * p1;
}
int inline power(int a, int b, int P) {
	int res = 1;
	while (b) {
		if (b & 1) res = (LL)res * a % P;
		a = (LL)a * a % P;
		b >>= 1;
	}
	return res;
}

unordered_map<int, int> mp;

int BSGS(int a, int b, int P) {
	int t = sqrt(P) + 1; mp.clear(); b %= P;
	for (int j = 0, s = b; j < t; j++) 
		mp[s] = j, s = (LL)s * a % P;
	a = power(a, t, P);
	for (int i = 1, s = 1; i <= t; i++) {
		s = (LL)s * a % P;
		if (mp.count(s) && i * t - mp[s] >= 0)
			return i * t - mp[s];
	}
	return -1;
}

int gcd(int a, int b) {
	return b ? gcd(b, a % b) : a;
}

int exBSGS(int a, int b, int P) {
	int x, y, d, A = 1, k = 0;
	while ((d = gcd(a, P)) > 1) {
		if (b % d) return -1;
		b /= d, P /= d, k++, A = (LL)A * (a / d) % P;
		if (A == b) return k;
	}
	exgcd(A, P, x, y); x = (x % P + P) % P;
	int res = BSGS(a, (LL)b * x % P, P);
	return res == -1 ? -1 : res + k;
}

const int N = 5000005, S = 3000;
const LL INF = 9e18;

LL p1[N], p2[S], m1[N], m2[S];

int n, primes[N], tot;

bool vis[N];

// 杜教筛 phi
LL s1(int x) {
	if (x < N) return p1[x];
	else if (p2[n / x] != INF) return p2[n / x]; 
	LL res = x * (x + 1ll) / 2;
	for (LL l = 2, r; l <= x; l = r + 1) {
		r = x / (x / l);
		res -= (r - l + 1) * s1(x / l);
	}
	return p2[n / x] = res;
}

// 杜教筛 mu

LL s2(int x) {
	if (x < N) return m1[x];
	else if (m2[n / x] != INF) return m2[n / x]; 
	LL res = 1;
	for (LL l = 2, r; l <= x; l = r + 1) {
		r = x / (x / l);
		res -= (r - l + 1) * s2(x / l);
	}
	return m2[n / x] = res;
}

// 线性筛
void inline linear() {
	m1[1] = p1[1] = 1;
	for (int i = 2; i < N; i++) {
		if (!vis[i]) primes[++tot] = i, m1[i] = -1, p1[i] = i - 1;
		for (int j = 1; i * primes[j] < N; j++) {
			vis[i * primes[j]] = true;
			if (i % primes[j] == 0) {
				p1[i * primes[j]] = p1[i] * primes[j];
				break;
			}
			p1[i * primes[j]] = p1[i] * (primes[j] - 1);
			m1[i * primes[j]] = -m1[i];
		}
	}
}

// 矩阵

struct Mat{
	int n, m, w[N][N];
	Mat operator * (const Mat &b) const {
		Mat c; c.n = n, c.m = b.m;
		memset(c.w, 0, sizeof c.w);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < b.m; j++)
				for (int k = 0; k < m; k++)
					c.w[i][j] = (c.w[i][j] + (LL)w[i][k] * b.w[k][j]) % P;
		return c;
	}
} res;

bool inline gauss() {
	int r, c;
	for (r = 1, c = 1; c <= n; c++) {
		int u = r;
		for (int i = r + 1; i <= n; i++) if (fabs(a[i][c]) > fabs(a[u][c])) u = i;
		if (fabs(a[u][c]) < eps) break;
		for (int i = c; i <= n + 1; i++) swap(a[r][i], a[u][i]);
		for (int i = n + 1; i >= c; i--) a[r][i] /= a[r][c];
		for (int i = 1; i <= n; i++) {
			if (i != r) {
				for (int j = 1; j <= n + 1; j++)
					if (j != c) a[i][j] -= a[r][j] * a[i][c];
			}
		}
		r++;		
	}
	return r == n + 1;
}


// --- 行列式求值

int inline det() {
	int res = 1, v = 1;
	for (int i = 1; i <= n; i++) {
		int t = -1;
		for (int j = i; j <= n; j++)
			if (a[j][i] && (t == -1 || a[j][i] > a[t][i])) t = j;
		if (t == -1) return 0;
		if (i != t) swap(a[t], a[i]), v *= -1;
		for (int j = i + 1; j <= n; j++) {
			if (a[j][i] > a[i][i]) swap(a[j], a[i]), v *= -1;
			while (a[j][i]) {
				int t = a[i][i] / a[j][i];
				for (int k = i; k <= n; k++) {
					a[i][k] = (a[i][k] - 1ll * a[j][k] * t % P + P) % P;
					swap(a[j][k], a[i][k]);
				}
				v *= -1;
			}
		}
		res = (LL)res * a[i][i] % P;
	}
	return (LL)res * (v + P) % P;
}

// 拉格朗日插值

int inline Interpo(int k, int n, int x[], int y[]){
	int res = 0;
	for (int i = 1; i <= n; i++) {
		int v1 = y[i], v2 = 1;
		for (int j = 1; j <= n; j++) {
			if (i != j) {
				v1 = (LL) v1 * (K - x[j]) % P;
				v2 = (LL) v2 * (x[i] - x[j]) % P;
			}
		}
		res = ((res + (LL)v1 * power(v2, P - 2) % P) % P + P) % P;
	}
	return res;
}

// Min25

int inv2 = power(2, P - 2), inv6 = power(6, P - 2);

// 求 g_k 函数: <= x 的和
int inline getS(LL x, int k) {
	if (k == 1) return (x % P * (x % P + 1ll) % P * inv2 + P - 1ll)  % P;
	if (k == 2) return (P - 1ll + x % P * (x % P + 1ll) % P * (2ll * x % P + 1) % P * inv6) % P;
}

int inline getV(LL x, int k) {
	if (k == 1) return x % P;
	if (k == 2) return (LL)x % P * x % P;
}

bool vis[M];

int primes[M], tot;

void inline linear(int n) {
	for (int i = 2; i <= n; i++) {
		if (!vis[i]) primes[++tot] = i;
		for (int j = 1; primes[j] <= n / i; j++) {
			vis[i * primes[j]] = true;
			if (i % primes[j] == 0) break;
		}
	}
}

// 预处理 g_k 处所有 n / i 形式的质数前缀和

struct MP1{
	int m, g[M], pos1[M], pos2[M], len, id;
	LL n, d[M];
	int inline getPos(LL x) {
		return x <= m ? pos1[x] : pos2[n / x];
	}
	void inline add(LL v) {
		d[++len] = v;
		g[len] = getS(v, id);
		if (v <= m) pos1[v] = len;
		else pos2[n / v] = len; 
	}
	void build(LL sum, int t) {
		m = sqrt(n = sum); id = t;
		for (LL i = 1, j; i <= n; i = j + 1) {
			LL v = n / i; j = n / v;
			if (v <= m) break;
			add(v);
		}
		for (int i = m; i; i--) add(i);
		for (int i = 1; i <= tot && (LL)primes[i] * primes[i] <= n; i++) {
			LL pr = primes[i];
			for (int j = 1; j <= len && pr * pr <= d[j]; j++) {
				int k = getPos(d[j] / pr);
				g[j] = (g[j] - (LL)getV(pr, id) * (g[k] - g[getPos(primes[i - 1])] + P) % P + P) % P;
			}
		}
	}
	int inline s(LL x) { return g[getPos(x)]; }
} t1, t2;

int inline get(LL x) {
	return (t2.s(x) - t1.s(x) + P) % P;
}

int inline calc(LL x) {
	return x % P * (x % P - 1ll + P) % P;
}

void inline add(int &x, int y) {
	(x += y) %= P;
}

int inline s(LL n, int t) {
	if (primes[t] >= n) return 0;
	int ans = (get(n) - get(primes[t]) + P) % P;
	for (int i = t + 1; i <= tot && (LL)primes[i] * primes[i] <= n; i++) {
		int pr = primes[i];
		LL v = pr;
		for (int j = 1; v <= n; v = v * pr, j++) {
			add(ans, (LL)calc(v) * ((j != 1) + s(n / v, i)) % P);
		}
	}
	return ans;
}

//---