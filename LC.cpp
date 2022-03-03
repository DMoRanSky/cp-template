// Skyqwq
#define pb push_back
#define fi first
#define se second
#define mp make_pair

using namespace std;

typedef pair<int, int> PII;
typedef long long LL;

template <typename T> bool chkMax(T &x, T y) { return (y > x) ? x = y, 1 : 0; }
template <typename T> bool chkMin(T &x, T y) { return (y < x) ? x = y, 1 : 0; }

class Solution {
public:
    int c[100005], n, a[100005], b[100005], p1[100005], p2[100005];
    void inline add(int x, int k) {
        for (; x <= n; x += x & -x) c[x] += k;
    }
    void inline clr() {
        for (int i = 1; i <= n; i++) c[i] = 0;
    }
    int inline ask(int x) {
        int ret = 0;
        for (; x ; x -= x & -x) ret += c[x];
        return ret;
    }
    LL L[100005], R[100005];
    long long goodTriplets(vector<int>& A, vector<int>& B) {
        n = A.size();
        for (int i = 0; i < n; i++) a[i + 1] = A[i] + 1, b[i + 1] = B[i] + 1;
        LL c = 0;
        for (int i = 1; i <= n; i++) {
            p1[a[i]] = i;
            p2[b[i]] = i;
        }
        for (int i = 1; i <= n; i++) {
            L[i] = ask(p2[a[i]]);
            add(p2[a[i]], 1);
        }
        clr();
        for (int i = n; i; i--) {
            R[i] = ask(n) - ask(p2[a[i]]);
            add(p2[a[i]], 1);
            c += L[i] * R[i];
        }
        return c;
    }
};