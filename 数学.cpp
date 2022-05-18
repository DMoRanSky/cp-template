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
	LL M = 1; int res = 0;
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
    Mat operator + (const Mat &b) const {
        Mat c; c.n = n, c.m = b.m;
        memset(c.w, 0, sizeof c.w);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < b.m; j++)
                c.w[i][j] = (w[i][j] + b.w[i][j]) % P;
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

// FMT / FWT

void inline OR(int n, int a[], int o) {
	for (int w = 1; w < n; w <<= 1) 
		for (int i = 0; i < n; i += (w << 1)) 
			for (int j = 0; j < w; j++)
				add(a[i + j + w], o * a[i + j]);
}

void inline AND(int n, int a[], int o) {
	for (int w = 1; w < n; w <<= 1) 
		for (int i = 0; i < n; i += (w << 1)) 
			for (int j = 0; j < w; j++)
				add(a[i + j], o * a[i + j + w]);
}


// 反向传 1/2
void inline XOR(int n, int a[], int o) {
	for (int w = 1; w < n; w <<= 1) 
		for (int i = 0; i < n; i += (w << 1)) 
			for (int j = 0; j < w; j++) {
				int u = a[i + j], v = a[i + j + w];
				a[i + j] = ((LL)u + v + P) * o % P;
				a[i + j + w] = ((LL)u - v + P) * o % P;
			}
}

// 子集卷积


void inline SubConv(int n, int a[], int b[], int c[]) {
	for (int i = 0; i < (1 << n); i++) {
		f[get(i)][i] = a[i];
		g[get(i)][i] = b[i];
	}
	for (int i = 0; i <= n; i++)
		OR(1 << n, f[i], 1), OR(1 << n, g[i], 1);
	for (int i = 0; i <= n; i++)
		for (int j = 0; j <= i; j++)
			for (int k = 0; k < (1 << n); k++)
				add(h[i][k], (LL)f[j][k] * g[i - j][k] % P);
	for (int i = 0; i <= n; i++) OR(1 << n, h[i], -1);
	for (int i = 0; i < (1 << n); i++) c[i] = h[get(i)][i];
}

// 矩阵求逆


bool inline matInv(int n, int a[N][N], int b[N][N]) {
	int c[N][N << 1];
	memset(c, 0, sizeof c);
	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= n; j++) c[i][j] = a[i][j];
		c[i][i + n] = 1;
	}
	for (int i = 1; i <= n; i++) {
		int t = -1;
		for (int j = i; j <= n; j++) if (c[j][i]) t = j;
		if (t == -1) return false;
		if (t != i) swap(c[t], c[i]);
		int inv = power(c[i][i], P - 2);
		for (int j = 1; j <= 2 * n; j++) c[i][j] = (LL)c[i][j] * inv % P;
		for (int j = 1; j <= n; j++) {
			if (i == j) continue;
			int v = c[j][i];
			if (!v) continue;
			for (int k = 1; k <= 2 * n; k++)
				c[j][k] = (c[j][k] - (LL)c[i][k] * v % P + P) % P;
		}
	}
	for (int i = 1; i <= n; i++)
		for (int j = 1; j <= n; j++) b[i][j] = c[i][j + n];
	return true;
}

// St:

namespace Red{
	int MOD, I;

	struct Num{
		int x, y;
	};

	Num operator * (const Num a, const Num b) {
		return (Num) { (int)(((LL)a.x * b.x + (LL)I * a.y % P * b.y) % P) , (int)(((LL)a.x * b.y + (LL)a.y * b.x) % P) }; 
	}


	Num inline power(Num a, int b) {
		Num res = (Num) { 1, 0 };
		while (b) {
			if (b & 1) res = res * a;
			a = a * a;
			b >>= 1;
		}
		return res;
	}

	int inline EulerCri(int x) {
		return power(x, (MOD - 1) >> 1, MOD);
	}

	int inline Sqrt(int n, int o) {
		MOD = o;
		if (n == 0) return 0;
		if (EulerCri(n) == P - 1) return -1;
		int a = rand() % P, t;
		while (EulerCri(t = ((LL)a * a + MOD - n) % P) != MOD - 1) a = rand() % MOD;
		I = t;
		int x0 = power((Num) { a, 1 }, (P + 1) >> 1).x;
		return min(x0, P - x0);
	}	
}


// _end

// Exgcd 求逆

int inline inv(int a, int P) {
	LL x, y;
	LL d = exgcd(a, P, x, y);
    LL k = P / d;
    x = (x % k + k) % k;
    return x;
}

// ExLucas


namespace Cpk{
	int h[M], p, k, v, P;
	void init(int Po, int K, int V) {
		p = Po, k = K, P = v = V;
		h[0] = 1;
		for (int i = 1; i <= v; i++) {
			h[i] = h[i - 1];
			if (i % p) h[i] = (LL)h[i] * i % P;	
		} 
	}
	int inline f(LL n) {
		if (!n) return 1;
		return (LL)f(n / p) * power(h[v], n / v, v) % P * h[n % v] % P;
	}
	LL inline g(LL n) {
		LL res = 0;
		while (n) {
			res += n / p;
			n /= p;
		}
		return res;
	}
	int inline C(LL n, LL m) {
		return (LL)power(p, g(n) - g(m) - g(n - m), P) * f(n) % P * inv((LL)f(m) * f(n - m) % P, P) % P;
	}
}

int inline exLucas(LL n, LL m, int P) {
	for (int i = 2; i <= P; i++) {
		if (P % i == 0) {
			int k = 0, t = 1;
			while (P % i == 0) k++, P /= i, t *= i;
			++len, p[len] = t;
			Cpk::init(i, k, t);
			a[len] = Cpk::C(n, m);
		}
	}
	return CRT(len);
}

// 
LL inline mul(LL a, LL b, LL P) {
	a %= P, b %= P;
	return ((a * b - (LL)((long double)a * b / P) * P) % P + P) % P; 
}

LL inline power(LL a, LL b, LL Mod) {
	LL res = 1;
	while (b) {
		if (b & 1) res = mul(res, a, Mod);
		a = mul(a, a, Mod);
		b >>= 1;
	}
	return res;
}

class MillerRabin{
private:
	const int p[12] = { 2, 3, 5, 7, 11, 13, 17, 19, 61, 2333, 4567, 24251 };
	bool inline check(LL x, int P) {
		if (x % P == 0 || power(P % x, x - 1, x) == 0) return false;
		LL k = x - 1;
		while (k % 2 == 0) {
			LL v = power(P % x, k >>= 1, x);
			if (v != 1 && v != x - 1) return false;
			if (v == x - 1) return true;
		}
		return true;
	}
public:
	bool inline isPrime(LL x) {
		if (x < 2) return false;
		for (int i = 0; i < 12; i++) {
			if (p[i] == x) return true;
			if (!check(x, p[i])) return false;
		}
		return true;
	}
} mr;

// PollardRho 

LL inline Abs(LL x) {
	return x < 0 ? -x : x;
}

LL inline f(LL x, LL c, LL P) {
	return (mul(x, x, P) + c) % P;
}

LL gcd(LL a, LL b) {
	return b ? gcd(b, a % b) : a;
}

LL inline PollardRho(LL n) {
	LL c = (LL)rand() * rand() % (n - 1) + 1;
	LL p = 0, q = 0, v = 1;
	for (int i = 1; ; i <<= 1, p = q, v = 1) {
		for (int j = 0; j < i; j++) {
			q = f(q, c, n);
			v = mul(v, Abs(q - p), n);
			if (j % 100 == 0) {
				LL d = gcd(v, n);
				if (d > 1) return d;
			}
		}
		LL d = gcd(v, n);
		if (d > 1) return d;
	}
	return n;
}

LL n, ans;

// 找最大质因子

void fact(LL n) {
	if (n < 2 || n <= ans) return;
	if (mr.isPrime(n)) {
		if (n > ans) ans = n;
		return;
	}	
	LL p = n;
	while (p >= n) p = PollardRho(n);
	while (n % p == 0) n /= p;
	fact(n), fact(p);
}

// _end


// 线性基

LL a[N];

int insert(LL x) {
	for (int i = L - 1; ~i; i--) {
		if(!(x >> i & 1)) continue; 
		if (a[i]) x ^= a[i];
		else { a[i] = x; return i; }
	}
	return -1;
}

// 动态插入一行的高斯消元（可撤销，就是直接复制）

struct Gauss{
    int a[N][N], m;
    bool o;
    void inline ins(int b[]) {
        for (int i = 1; i <= n + 1; i++) c[i] = b[i];
        for (int i = 1; i <= n; i++) {
            if (c[i]) {
                if (!a[i][i]) {
                    int inv = power(c[i], P - 2);
                    for (int j = 1; j <= n + 1; j++)
                        a[i][j] = (LL)c[j] * inv % P;
                    for (int j = i + 1; j <= n; j++) {
                        if (a[j][j] && a[i][j]) {
                            int t = a[i][j];
                            for (int k = 1; k <= n + 1; k++)
                                a[i][k] = (a[i][k] - 1ll * a[j][k] * t % P + P) % P;
                        }
                    }
                    for (int j = 1; j < i; j++) {
                        if (a[j][i]) {
                            int t = a[j][i];
                            for (int k = 1; k <= n + 1; k++)
                                a[j][k] = (a[j][k] - 1ll * a[i][k] * t % P + P) % P;
                        }
                    }
                    ++m;
                    return;
                } else {
                    int t = c[i];
                    for (int j = 1; j <= n + 1; j++)
                        c[j] = (c[j] - 1ll * t * a[i][j] % P + P) % P;
                }
            } 
        }
        if (c[n + 1]) {
            o = 1;
        }
    }
    void inline sh() {
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n + 1; j++) cout << a[i][j] << " ";
            cout << endl;
        }
    }
    void inline out() {
        if (o) {
            puts("Impossible!");
        } else {
            for (int i = 1; i <= m; i++) {
                for (int j = 1; j <= n + 1; j++) z[j] = 0;
                z[i] = 1;
                ins(z);
            }
            //puts("Succc!!");
            memset(ans, 0, sizeof ans);
            for (int i = 1; i <= n; i++) {
                int p = 0;
                for (int j = 1; j <= n; j++)
                    if (a[i][j]) p = j;
                ans[p] = a[i][n + 1];
            }
            for (int i = 1; i <= n; i++)
                printf("%d ", ans[i]);
            puts("");
        }
    }
} t[15];

// 万能欧几里得

struct Data{
} IN;
// IN 无贡献时。


// 信息合并
Data operator + (const Data &a, const Data &b) {

}

// 信息叠加 b 次
Data operator * (const Data &a, int b) {

}

Data calc(LL n, LL p, LL q, LL r, Data A, Data B) {
	if (!n) return IN;
    r %= q;
	if (p >= q) return calc(n, p % q, q, r, B * (p / q) + A, B);
    LL x = ((I)p * n + r) / q;
    if (!x) return A * n;
    return A * ((q - r - 1) / p) + B + calc(x - 1, q, p, q - r - 1, B, A) + A * (n - ((I)q * x - r - 1 + p) / p + 1); 
}