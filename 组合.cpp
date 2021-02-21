int inline power(int a, int b) {
	int res = 1;
	while (b) {
		if (b & 1) res = (LL)res * a % P;
		a = (LL)a * a % P;
		b >>= 1;
	}
	return res;
}

void inline preInv(int n) {
	inv[1] = 1;
	for (int i = 2; i <= n; i++)
		inv[i] = ((LL)P - P / i) * inv[P % i] % P;
}

void inline factPrework(int n) {
	fact[0] = infact[0] = 1;
	for (int i = 1; i <= n; i++) fact[i] = (LL)fact[i - 1] * i % P;
	infact[n] = power(fact[n], P - 2);
	for (int i = n - 1; i; i--) infact[i] = infact[i + 1] * (i + 1ll) % P;
}

int inline C(int a, int b) {
	if (a < b) return 0;
	return (LL)fact[a] * infact[b] % P * infact[a - b] % P;
}

int lucas(int a, int b) {
	if (a < P) return C(a, b);
	return (LL)C(a % P, b % P) * lucas(a / P, b / P) % P;
}

int inline det() {
	int res = 1, v = 1;
	for (int i = 1; i <= n; i++) {
		for (int j = i + 1; j <= n; j++) {
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