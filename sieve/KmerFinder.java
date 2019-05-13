package sieve;

import java.util.ArrayDeque;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Queue;

public class KmerFinder {
	static long[][] kmerize(String s, int k)
	{
		long[][] res = new long[2][s.length() - k + 1];
		long[] slidingKmer = new long[2];
		for(int i = 0; i<k; i++)
		{
			char ch = s.charAt(i);
			slidingKmer = updateSlidingKmer(slidingKmer, ch, k);
		}
		for(int i = k; i <= s.length(); i++)
		{
			res[0][i-k] = slidingKmer[0];
			res[1][i-k] = slidingKmer[1];
			// Update kmers
			if(i == s.length()) break;
			char ch = s.charAt(i);
			slidingKmer = updateSlidingKmer(slidingKmer, ch, k);
		}
		return res;
	}
	static long hash64(long key, int k)
	{
		long mask = (1L << (2*k)) - 1;
		key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
		key = key ^ key >> 24;
		key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
		key = key ^ key >> 14;
		key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
		key = key ^ key >> 28;
		key = (key + (key << 31)) & mask;
		return key;
	}
	static long[] minimizers(String s, int k, int w, int posStrandBits)
	{
		long[][] kmers = kmerize(s, k);
		MinQueue mq = new MinQueue(posStrandBits);
		HashSet<Long> minimizers = new HashSet<Long>();
		for(int i = 0; i<w; i++)
		{
			for(int strand = 0; strand < 2; strand++)
			{
				mq.add(strand + (i << 1) + (hash64(kmers[strand][i], k) << posStrandBits));
			}
		}
		for(int i = w; i<kmers[0].length; i++)
		{
			// Add minimizers to set
			long[] cur = mq.mins();
			for(long x : cur)
			{
				minimizers.add(x);
			}
			
			// Update queue
			for(int strand = 0; strand < 2; strand++)
			{
				mq.remove();
				mq.add(strand + (i << 1) + (hash64(kmers[strand][i], k) << posStrandBits));
			}
		}
		
		// Move minimizers to an array
		long[] res = new long[minimizers.size()];
		int idx = 0;
		for(long x : minimizers) res[idx++] = x;
		return res;
	}
	static long[] updateSlidingKmer(long[] cur, char ch, int k)
	{
		long kmer = cur[0];
		long revComp = cur[1];
		int c = charToInt(ch);
		
		// Update kmer
		long allButHighest = ((1L<< (2*k - 2)) - 1);
		kmer = ((kmer&allButHighest) << 2) | c;
		
		// Update reverse complement kmer
		long rc = c ^ 3;
		revComp >>= 2;
		revComp |= (rc << (2*k - 2));
		
		return new long[] {kmer, revComp};
	}
	static int charToInt(char c)
	{
		if(c >= 'a' && c <= 'z') c += 'A' - 'a';
		if(c == 'A') return 0;
		else if(c == 'C') return 1;
		else if(c == 'G') return 2;
		else return 3;
	}
	static class MinQueue
	{
		int posStrandBits;
		ArrayDeque<Long> d;
		Queue<Long> q;
		public MinQueue(int psb)
		{
			posStrandBits = psb;
			d = new ArrayDeque<Long>();
			q = new LinkedList<Long>();
		}
		long[] mins()
		{
			int count = 0;
			for(long x : d)
			{
				if((x>>posStrandBits) == (d.getFirst()>>posStrandBits)) count++;
				else break;
			}
			long[] res = new long[count];
			int idx = 0;
			for(long x : d)
			{
				res[idx++] = x;
				if(idx == res.length) break;
			}
			return res;
		}
		long min()
		{
			return d.getFirst();
		}
		void add(long x)
		{
			q.add(x);
			while(!d.isEmpty() && d.getLast() > x) d.removeLast();
			d.add(x);
		}
		void remove()
		{
			if(q.isEmpty()) return;
			q.poll();
			if(!d.isEmpty() && d.getFirst().equals(q.peek())) d.poll();
		}
	}
}
