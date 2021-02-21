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
