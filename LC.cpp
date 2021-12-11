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
    int n, m, a[1005][1005], b[1005][1005], c[1005];
    int inline work(int d[1005]) {
        for (int i = 1; i <= m; i++) c[i] = d[i];
        for (int i = 2; i <= m; i++)
            chkMin(c[i], c[i - 1] + 1);
        for (int i = m; i >= 2; i--)
            chkMin(c[i - 1], c[i] + 1);
        int ret = 0;
        for (int i = 1; i <= m; i++) ret += c[i];
        return ret;
    }
    int inline work() {
        int ret = 0;
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= m; j++) {
                if (a[i][j]) b[i][j] = b[i - 1][j] + 1;
                else b[i][j] = 0;
            }
            ret += work(b[i]);
        }
        return ret;
    }
    int countPyramids(vector<vector<int>>& g) {
        n = g.size(), m = g[0].size();
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                a[i + 1][j + 1] = g[i][j];
        int ans = work();
        memcpy(b, a, sizeof b);
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= m; j++)
                a[i][j] = b[n - i + 1][j];
        }
        ans += work();
        return ans;
    }
};