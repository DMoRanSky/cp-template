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

class Solution {
public:
    set<ULL> M;
    int t[26][26];
    ULL hs(string u) {
        ULL w = 0;
        for (char v: u)
            w = w * 131 + v;
        return w;
    }
    long long distinctNames(vector<string>& w) {
        for (auto v: w) M.insert(hs(v));
        for (auto v: w) {
            string s = v;
            int A = s[0] - 'a';
            for (int i = 0; i < 26; i++) {
                if (i != A) {
                    s[0] = i + 'a';
                    if (!M.count(hs(s))) t[A][i]++;
                }
            }
        }LL ans = 0;
        for(int i = 0; i < 26; i++) {
            for (int j = 0; j < 26; j++) {
                ans += t[i][j] * t[j][i];
            }
        }
        return ans;
    }
};