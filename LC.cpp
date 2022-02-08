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
    vector<LL> a1, a0, b1, b0;
    LL chk(LL x) {
        int cnt = 0;
        LL s = 0;
        for (int i = (int)a1.size() - 1, j = 0; i >= 0; i--) {
            while (j < b1.size() && b1[j] * a1[i] <= x) j++, cnt++;
            s += cnt;
        }
        cnt = 0;
        for (int i = (int)a1.size() - 1, j = (int)b0.size() - 1; i < a0.size(); i--) {
            while (j >= 0 && b0[j] * a0[i] <= x) j--, cnt++;
            s += cnt;
        }
        cnt = 0;
        for (int i = 0, j = 0; i < a1.size(); i++) {
            while (j < b0.size() && b0[j] * a1[i] <= x) j++, cnt++;
            s += cnt;
        }
        cnt = 0;
        for (int i = 0, j = 0; i < b1.size(); i++) {
            while (j < a0.size() && a0[j] * b1[i] <= x) j++, cnt++;
            s += cnt;
        }
        return s;
    }
    LL kthSmallestProduct(vector<int>& a, vector<int>& b, long long k) {
        for (int v: a) {
            if (v >= 0) a1.pb(v);
            else a0.pb(v);
        }
        for (int v: b) {
            if (v >= 0) b1.pb(v);
            else b0.pb(v);
        }
        cout << chk(-54) << endl;
        LL l = -1e11, r = 1e11;
        while (l < r) {
            LL mid = (l + r) >> 1;
            if (chk(mid) >= k) r = mid;
            else l = mid + 1;
        }
        return r;
    }
};