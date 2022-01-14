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
	int n, m;
	vector<int> s[100005];
	vector<int> w[100005], c[100005];
	int inline get(int x1, int y1, int x2, int y2) {
		return s[x2][y2] - s[x1 - 1][y2] - s[x2][y1 - 1] + s[x1 - 1][y1 - 1];
	}
    bool possibleToStamp(vector<vector<int>>& g, int x, int y) {
    	n = g.size(), m = g[0].size();
    	s[0].resize(m + 2, 0);
    	c[0].resize(m + 2, 0);
    	for (int i = 1; i <= n; i++) {
    		w[i].resize(m + 2, 0);
    		s[i].resize(m + 2, 0);
    		c[i].resize(m + 2, 0);
    		for (int j = 1; j <= m; j++) {
    			w[i][j] = g[i - 1][j - 1];
    			s[i][j] = s[i - 1][j] + s[i][j - 1] - s[i - 1][j - 1] + w[i][j];
    		}
    	}
    	c[n + 1].resize(m + 2, 0);
    	for (int i = 1; i + x - 1 <= n; i++) {
    		for (int j = 1; j + y - 1 <= m; j++) {
    			if (!get(i, j, i + x - 1, j + y - 1)) {
    				c[i][j]++, c[i + x][j]--;
    				c[i][j + y]--, c[i + x][j + y]++;
    			}
    		}
    	}
    	for (int i = 1; i <= n; i++) {
    		for (int j = 1; j <= m; j++) {
    			int v = c[i][j] + w[i][j];
    			if (!v) return 0;
    		}
    	}
    	return 1;
    }
};