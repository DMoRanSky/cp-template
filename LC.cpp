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

typedef unsigned long long ULL;

const int INF = 0x3f3f3f3f;

class Solution {
public:
    int s[1005];
    int longestSubsequence(string g, int k) {
        int ans = 0;
        for (char c: g) {
            if (c == '0') ans++;
        }
        for (int i = 0; i < g.size(); i++) {
            s[i + 1] = s[i] + (g[i] == '0');
        }
        int h  =0;
        for (int i = 0; i < 30 ;i++)
            if (k >> i & 1) h = i;
        // h + 1 
        int n = g.size();
        if (n < h + 1) {
            return n;
        }
        vector<int> t;
        for (int i = h; i >= 0; i--)
            t.pb(k >> i & 1);
        chkMax(ans, s[n - h] + h);
        for (int i = 0; i < g.size(); i++) {
            if (g[i] == '1') {
                int sum = s[i] + t.size(), now = 1;
                for (int j = i + 1; j < g.size(); j++) {
                    if (now == t.size()) break;
                    if (g.size() - j < t.size() - now) break;
                    if (t[now] == g[j] - '0') {
                        now++;
                    } else {
                        if (t[now] == 1 && g[j] == '0') {
                            now = t.size();
                            break;
                        }
                    }
                }
                if (now == t.size()) chkMax(ans, sum);
            }
        }
        return ans;
    }

};