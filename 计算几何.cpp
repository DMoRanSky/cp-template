

const double eps = 1e-4;
typedef pair<double, double> PDD;
struct Line{
    PDD s, t;
};

int inline cmp(double x, double y) {
    if (fabs(x - y) < eps) return 0;
    return x < y ? -1 : 1;
}

double inline cross(PDD a, PDD b) { return a.fi * b.se - a.se * b.fi; }
PDD operator - (const PDD &a, const PDD &b) { return make_pair(a.fi - b.fi, a.se - b.se); }
PDD operator + (const PDD &a, const PDD &b) { return make_pair(a.fi+ b.fi, a.se+ b.se); }
PDD operator / (const PDD &a, double b) { return make_pair(a.fi / b, a.se / b); }
PDD operator * (const PDD &a, double b) { return make_pair(a.fi * b, a.se * b); }
double inline area(PDD a, PDD b, PDD c) { return cross(b - a, c - a); }
double inline dot(PDD a, PDD b) { return a.fi * b.fi + a.se * b.se; }
double inline len(PDD a) { return sqrt(dot(a, a)); }
double inline project(PDD a, PDD b, PDD c) { return dot(b - a, c - a) / len(b - a); }
double inline dist(PDD a, PDD b) { return sqrt((a.fi - b.fi) * (a.fi - b.fi) + (a.se - b.se) * (a.se - b.se)); }
// 顺时针转 x
PDD inline rotate(PDD a, double x) { return make_pair ( cos(x) * a.fi + sin(x) * a.se, -sin(x) * a.fi + cos(x) * a.se ); }
PDD inline norm(PDD a) { return a / len(a); }
double angle(PDD a, PDD b) {
    return acos(dot(a, b) / len(a) / len(b));
}
int sign(double fi) {
    if (fabs(fi) < eps) return 0;
    if (fi < 0) return -1;
    return 1;
}

bool segInter(PDD a1, PDD a2, PDD b1, PDD b2) {
    double c1 = cross(a2 - a1, b1 - a1), c2 = cross(a2 - a1, b2 - a1);
    double c3 = cross(b2 - b1, a2 - b1), c4 = cross(b2 - b1, a1 - b1);
    return sign(c1) * sign(c2) <= 0 && sign(c3) * sign(c4) <= 0;
}

bool cmp2 (const Line &a, const Line &b) {
    double A = getAngle(a), B = getAngle(b);
    if (A != B) return A < B;
    else return area(a.s, a.t, b.t) < 0;
}

PDD getInter(PDD p, PDD v, PDD q, PDD w) {
    PDD u = p - q;
    double t = cross(w, u) / cross(v, w);
    return make_pair(p.fi + t * v.fi, p.se + t * v.se);
}

PDD getInter(Line a, Line b) { return getInter(a.s, a.t - a.s, b.s, b.t - b.s); }

bool inline Right(Line a, Line b, Line c) {
    PDD u = getInter(b, c);
    return area(a.s, a.t, u) <= 0;
}

// 凸包

void inline andrew() {
    sort(p + 1, p + 1 + n);
    for (int i = 1; i <= n; i++) {
        while (top > 1 && area(p[s[top - 1]], p[s[top]], p[i]) < 0) {
            if (area(p[s[top - 1]], p[s[top]], p[i]) <= 0) st[s[top--]] = false;
            else top--;
        }
        st[i] = true, s[++top] = i;
    }
    st[1] = false;
    for (int i = n; i; i--) {
        if (!st[i]) {
            while (top > 1 && area(p[s[top - 1]], p[s[top]], p[i]) <= 0) 
                st[s[top--]] = false;
            st[i] = true, s[++top] = i;
        }
    }
    for (int i = 0; i < top; i++) s[i] = s[i + 1];
    top--;
}

struct Line{
    PDD s, t;
    int id;
} e[N];

// 半平面交
double HPI() {
    sort(e + 1, e + 1 + n, cmp2);
    int hh = 0, tt = -1;
    for (int i = 1; i <= n; i++) {
        if (i && getAngle(e[i]) == getAngle(e[i - 1])) continue;
        while (hh < tt && Right(e[i], e[q[tt - 1]], e[q[tt]])) tt--;
        while (hh < tt && Right(e[i], e[q[hh]], e[q[hh + 1]])) hh++;
        q[++tt] = i;
    }
    while (hh < tt && Right(e[q[hh]], e[q[tt - 1]], e[q[tt]])) tt--;
    while (hh < tt && Right(e[q[tt]], e[q[hh]], e[q[hh + 1]])) hh++;
    q[++tt] = q[hh];
    tot = 0;
    for (int i = hh; i < tt; i++)
        p[++tot] = getInter(e[q[i]], e[q[i + 1]]);
    double res = 0;
    for (int i = 1; i < tot; i++)
        res += area(p[1], p[i], p[i + 1]);
    return res / 2;
}

Point inline getCircle(Point a, Point b, Point c) {
    return Inter((a + b) / 2, rotate(b - a, PI / 2), (a + c) / 2, rotate(c - a, PI / 2));
}

// 最小圆覆盖

void inline minCircle(PDD a[]) {
    random_shuffle(a + 1, a + 1 + n);
    double r = 0; Point u = a[1]; 
    for (int i = 2; i <= n; i++) {
        if (cmp(r, len(u - a[i])) == -1) {
            r = 0, u = a[i];
            for (int j = 1; j < i; j++) {
                if (cmp(r, len(u - a[j])) == -1) {
                    r = len(a[i] - a[j]) / 2, u = (a[i] + a[j]) / 2;
                    for (int k = 1; k < j; k++) {
                        if (cmp(r, len(u - a[k])) == -1) {
                            u = getCircle(a[i], a[j], a[k]), r = len(a[i] - u);
                        } 
                    }
                }
            }
        }
    }
}


// 自适应辛普森积分
double inline f(double fi) {
    return ?;
}
double inline s(double l, double r) {
    double mid = (l + r) / 2;
    return (r - l) * (f(l) + 4 * f(mid) + f(r)) / 6;
}

double inline asr(double l, double r) {
    double mid = (l + r) / 2, v = s(l, r);
    double a = s(l, mid), b = s(mid, r);
    if (fabs(a + b - v) < eps) return v;
    else return asr(l, mid) + asr(mid, r);
}

// https://codeforces.com/contest/1284/problem/E 的怨念 不丢精度的极角排序

LL inline cross(PII x, PII y) {
    return 1ll * x.fi * y.se - 1ll * x.se * y.fi;
}

int inline quad(PII x) {
    if (x.fi >= 0 && x.se >= 0) return 1;
    if (x.fi <= 0 && x.se >= 0) return 2;
    if (x.fi <= 0 && x.se <= 0) return 3;
    if (x.fi >= 0 && x.se <= 0) return 4;
    return 0;
}

// PII andrew + mincowf

LL operator * (PII a, PII b) {
    return (LL)a.fi * b.se - (LL)a.se * b.fi;
}

PII operator + (PII a, PII b) {
    return mp(a.fi + b.fi, a.se + b.se);
}

PII operator - (PII a, PII b) {
    return mp(a.fi - b.fi, a.se - b.se);
}

LL dot (PII a, PII b) {
    return (LL)a.fi * a.se + (LL)b.fi * b.se;
}

vector<PII> inline andrew(vector<PII> a) {
    int n = a.size();
    top = 0;
    sort(a.begin(), a.end());
    for (int i = 0; i < n; i++) {
        while (top > 1 && (a[i] - a[s[top - 1]]) * (a[s[top]] - a[s[top - 1]]) > 0) {
            vis[s[top--]] = 0;
        }
        vis[i] = 1, s[++top] = i;
    }
    vis[0] = 0;
    for (int i = n - 1; i >= 0; i--) {
        if (!vis[i]) {
            while (top > 1 && (a[i] - a[s[top - 1]]) * (a[s[top]] - a[s[top - 1]]) > 0)
                vis[s[top--]] = 0;
            vis[i] = 1, s[++top] = i;
        }
    }
    --top;
    vector<PII> ret;
    for (int i = 1; i <= top; i++) ret.pb(a[s[i]]);
    for (int i = 0; i < n; i++) vis[i] = 0;
    return ret;
}

// 有

vector<PII> calc(vector<PII> a, vector<PII> b) {
    vector<PII> c;
    c.pb(a[0] + b[0]);
    vector<PII> dx, dy;
    for (int i = 1; i < a.size(); i++) dx.pb(a[i] - a[i - 1]);
    dx.pb(a[0] - a.back());
    for (int i = 1; i < b.size(); i++) dy.pb(b[i] - b[i - 1]);
    dy.pb(b[0] - b.back());
    int i = 0, j = 0;
    while (i < dx.size() && j < dy.size()) {
        if (dx[i] * dy[j] > 0)
            c.pb(c.back() + dx[i++]);
        else if (dx[i] * dy[j] == 0 && c.size() > 1) {
            // 共线放一起不然是错的！！！！
            if (dot(c.back() - c[c.size() - 2], dx[i]) > 0)
                c.pb(c.back() + dx[i++]);
            else c.pb(c.back() + dy[j++]);
        } else {
            c.pb(c.back() + dy[j++]);
        }
    }
    while (i < dx.size()) c.pb(c.back() + dx[i++]);
    while (j < dy.size()) c.pb(c.back() + dy[j++]);
    assert(c.back() == c[0]);
    c.pop_back();
    return c;
}