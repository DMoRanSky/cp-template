typedef pair<double, double> PDD;
int inline cmp(double x, double y) {
    if (fabs(x - y) < eps) return 0;
    return x < y ? -1 : 1;
}



double inline cross(PDD a, PDD b) { return a.x * b.y - a.y * b.x; }
PDD operator - (const PDD &a, const PDD &b) { return make_pair(a.x - b.x, a.y - b.y); }
PDD operator + (const PDD &a, const PDD &b) { return make_pair(a.x+ b.x, a.y+ b.y); }
PDD operator / (const PDD &a, double b) { return make_pair(a.x / b, a.y / b); }
PDD operator * (const PDD &a, double b) { return make_pair(a.x * b, a.y * b); }
double inline area(PDD a, PDD b, PDD c) { return cross(b - a, c - a); }
double inline dot(PDD a, PDD b) { return a.x * b.x + a.y * b.y; }
double inline len(PDD a) { return sqrt(dot(a, a)); }
double inline project(PDD a, PDD b, PDD c) { return dot(b - a, c - a) / len(b - a); }
double inline dist(PDD a, PDD b) { return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y)); }
PDD inline rotate(PDD a, double x) { return make_pair ( cos(x) * a.x + sin(x) * a.y, -sin(x) * a.x + cos(x) * a.y ); }
PDD inline norm(PDD a) { return a / len(a); }

int sign(double x) {
    if (fabs(x) < eps) return 0;
    if (x < 0) return -1;
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
    return make_pair(p.x + t * v.x, p.y + t * v.y);
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
            if (area(p[s[top - 1]], p[s[top]], p[i]) < 0) st[s[top--]] = false;
            else top--;
        }
        st[i] = true, s[++top] = i;
    }
    st[1] = false;
    for (int i = n; i; i--) {
        if (!st[i]) {
            while (top > 1 && area(p[s[top - 1]], p[s[top]], p[i]) <s 0) 
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
double inline f(double x) {
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
