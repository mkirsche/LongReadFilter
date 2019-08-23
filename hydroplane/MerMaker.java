package hydroplane;

import java.util.ArrayDeque;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Queue;

/*
 * This class handles the efficient generation of kmers and minimizers for a given string
 */
public class MerMaker {
	
	/*
	 * Generates a 2 by (n-k+1) array of all kmers of a string
	 */
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
	
	/*
	 * Generates an array of all kmers of a string with alternating strand
	 */
	static long[] kmers(String s, int k)
	{
		long[][] all = kmerize(s, k);
		int n = all[0].length;
		long[] res = new long[2*n];
		for(int i = 0; i<n; i++)
		{
			res[2*i] = all[0][i];
			res[2*i+1] = all[1][i];
		}
		return res;
	}
	
	/*
	 * Same hash function used in minimap to get a hash of a number with with 2*k bits 
	 */
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
	
	/*
	 * Computes the (w, k) minimizers of a String s.  Each minimizer is encoded as a 64-bit integer as follows:
	 *     Lowest order bit is the strand (0 for given strand, 1 for reverse complement)
	 *     The next <posStrandBits> bits are the position in the string where the kmer begins
	 *         Note for the reverse strand the position is the minimum (so for s[5..2] it would be 2)
	 *     The next 2*k bits are the contents of the kmer itself (after going through a random hash)
	 */
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
	
	/*
	 * Appends a character to the end of an encoded kmer (and beginning of its reverse complement)
	 *  and removes a character from the other end to keep the length as k.
	 */
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
	
	/*
	 * A/a = 0, C/c = 1, G/g = 2, T/t = 3
	 */
	static int charToInt(char c)
	{
		if(c >= 'a' && c <= 'z') c += 'A' - 'a';
		if(c == 'A') return 0;
		else if(c == 'C') return 1;
		else if(c == 'G') return 2;
		else return 3;
	}
	
	/*
	 * A minimum queue data structure, which supports the following operations:
	 * 
	 * 	mins() -> queries all elements with the minimum hashed kmer value
	 *  min() -> returns any one element with the minimum hashed kmer value
	 *  add(x) -> Adds x to the end of the queue
	 *  remove() -> Removes an element from the beginning of the queue
	 *  
	 *  All of these operations are amortized constant time
	 */
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
