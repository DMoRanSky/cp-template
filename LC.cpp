// Skyqwq
#include <bits/stdc++.h>

#define pb push_back
#define fi first
#define se second
#define mp make_pair

using namespace std;

typedef pair<int, int> PII;
typedef long long LL;

template <typename T> bool chkMax(T &x, T y) { return (y > x) ? x = y, 1 : 0; }
template <typename T> bool chkMin(T &x, T y) { return (y < x) ? x = y, 1 : 0; }

template <typename T> void inline read(T &x) {
    int f = 1; x = 0; char s = getchar();
    while (s < '0' || s > '9') { if (s == '-') f = -1; s = getchar(); }
    while (s <= '9' && s >= '0') x = x * 10 + (s ^ 48), s = getchar();
    x *= f;
}
class Solution {
public:
	int largestVariance(string str) {
       int n = str.size(), ans = 0; 
    	for (int i = 0; i < 26; i++) {
    		for (int j = 0; j < 26; j++) {
    			if (i == j) continue;
    			int mn = 0, s = 0, t = 0, pre = -1;
    			for (int x = 0, u = 0; x < n; x++) {
    				char v = str[x];
    				if (v - 'a' == i) s++;
    				if (v - 'a' == j) pre = x, s--;
    				while (u < pre) {
    					if (str[u] - 'a' == i) t++;
    					if (str[u] - 'a' == j) t--;
    					chkMin(mn, t); 
    					u++;
    				}
    				if (pre != -1) chkMax(ans, s - mn);
    			}
    		}
    	}
    	return ans;
    }
};