/*
 * Used to perform a precise dynamic programming alignment on the ends of reads to verify containment
 * It uses a chain of kmer matches to seed the alignment and speed it up
 */
package hydroplane;

import java.util.ArrayList;

import hydroplane.ReadIndex.Hit;

public class DynamicProgrammingAligner {
	double dpContains(Read indexedRead, int k, Read r, ArrayList<Hit> sharedKmers, int[] matchChain, boolean checkMiddle)
	{
		int rs = sharedKmers.get(0).rp.rs;
		int strand = rs % 2;
		Hit leftmost = null, rightmost = null;
		for(int idx : matchChain)
		{
			Hit h = sharedKmers.get(idx);
			if(leftmost == null || h.myIndex < leftmost.myIndex) leftmost = h;
			if(rightmost == null || h.myIndex > rightmost.myIndex) rightmost = h;
		}
		int theirMin = Math.min(leftmost.rp.p, rightmost.rp.p);
		int theirMax = leftmost.rp.p + rightmost.rp.p - theirMin;
		
		int subStart = Math.max(0, theirMin - 2 * leftmost.myIndex);
		int subEnd = Math.min(indexedRead.s.length()-1, theirMax + 2 * (r.s.length() - rightmost.myIndex));
		
		String leftquery = r.s.substring(0, leftmost.myIndex);
		String rightquery = r.s.substring(rightmost.myIndex + k);
		
		String leftCand = indexedRead.s.substring(subStart, theirMin);
		String rightCand = indexedRead.s.substring(theirMax + k, subEnd);
		
		double res = 1.0;
		
		if(strand == 0)
		{
			res = Math.min(dp(revComp(leftquery), revComp(leftCand)), dp(rightquery, rightCand));
		}
		else
		{
			res = Math.min(dp(rightquery, revComp(leftCand)), dp(revComp(leftquery), rightCand));
		}
		
		if(checkMiddle)
		{
			String middlequery = r.s.substring(leftmost.myIndex, rightmost.myIndex);
			String middleCand = indexedRead.s.substring(theirMin, theirMax);
			if(strand != 0)
			{
				middleCand = revComp(middleCand);
			}
			res = Math.min(res, dp(middlequery, middleCand));
		}
		
		return res;
	}
	double dp(String query, String candidate)
	{
		int n = query.length();
		int m = candidate.length();
		if(n == 0)
		{
			return 1.0;
		}
		
		int[][] dp = new int[2][m+1];
		int maxScore = -n;
		for(int i = 0; i<=n; i++)
			for(int j = Math.max(0, i - 50); j <= Math.min(m, i + 50); j++)
			{
				dp[i&1][j] = -i;
				if(i > 0) dp[i&1][j] = Math.max(dp[i&1][j], dp[(i&1)^1][j] - 5);
				if(j > 0) dp[i&1][j] = Math.max(dp[i&1][j], dp[i&1][j-1] - 1);
				if(i > 0 && j > 0)
				{
					int score = (query.charAt(i-1) == candidate.charAt(j-1)) ? 1 : -1;
					dp[i&1][j] = Math.max(dp[i&1][j], dp[(i&1)^1][j-1] + score);
				}
				if(i == n)
				{
					maxScore = Math.max(maxScore, dp[i&1][j]);
				}
			}
		
		double prop = maxScore * 1.0 / n;
		
		return prop;
	}
	String revComp(String s)
	{
		char[] res = new char[s.length()];
		for(int i = 0; i<s.length(); i++)
		{
			char c = s.charAt(s.length()-1-i);
			if(c == 'A') res[i] = 'T';
			else if(c == 'C') res[i] = 'G';
			else if(c == 'G') res[i] = 'C';
			else res[i] = 'A';
		}
		return new String(res);
	}
}
