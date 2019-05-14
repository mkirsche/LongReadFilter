package sieve;

import java.util.*;

public class ReadIndex {
	int n;
	int k;
	int w;
	String[] names;
	int[] countContaining;
	HashMap<Integer, LongList>[] kmerMap;
	HashSet<Long> badKmers;
	int logNumMaps = 16;
	int posStrandBits = 23;
	boolean minimizers;
	boolean verbose;
	CommandLineParser clp;
	
	@SuppressWarnings("unchecked")
	ReadIndex(ReadLengthSeparator re, CommandLineParser clp)
	{
		this.clp = clp;
		w = clp.w;
		minimizers = w > 1;
		verbose = clp.verbose;
		k = clp.k;
		n = re.data.length;
		countContaining = new int[n];
		kmerMap = new HashMap[1<<logNumMaps];
		for(int i = 0; i<(1<<logNumMaps); i++) kmerMap[i] = new HashMap<>();
		badKmers = new HashSet<>();
		long totalReadLength = 0;
		System.err.println("Building index with " + n + " reads");
		for(int i = 0; i<n; i++)
		{
			add(i, re.data[i]);
			if(clp.verbose && i > 0 && i%1000 == 0)
			{
				System.err.println("Added " + i + " reads");
			}
			totalReadLength += re.data[i].s.length();
			re.data[i] = null;
		}
		System.err.println("Index built with " + n + " reads having total length " + totalReadLength);
		int numKmers = 0;
		for(HashMap<Integer, LongList> hm : kmerMap) numKmers += hm.size();
		System.err.println("Number of kmers: " + numKmers);
		System.err.println("Bad kmers: " + badKmers.size());
	}
	
	void init()
	{
		
	}
	boolean contains(Read r)
	{
		double threshold = clp.p;
		int endLength = clp.el;
		int readLength = r.s.length();
		if(readLength < 5000)
		{
			return true;
		}
		
		HashMap<Integer, ArrayList<Hit>> hits = new HashMap<>();
		
		if(!minimizers)
		{
			long[][] kmers = KmerFinder.kmerize(r.s, k);
			
			// Loop through both strands since only forward strand in index
			for(int strand = 0; strand<2; strand++)
			{
				for(int i = 0; i<kmers[strand].length; i++)
				{
					long kmer = kmers[strand][i];
					if(badKmers.contains(kmer)) continue;
					LongList currentHits = kmerMap[(int)(kmer&((1<<logNumMaps)-1))].get((int)(kmer>>logNumMaps));
					if(currentHits == null) continue;
					for(int hitIndex = 0; hitIndex<currentHits.size; hitIndex++)
					{
						ReadPosition rp = ReadPosition.decode(currentHits.a[hitIndex]);
						if(rp.rs%2 != strand) rp.rs ^= 1;
						int readKey = rp.rs;
						if(!hits.containsKey(readKey))
						{
							hits.put(readKey, new ArrayList<>());
						}
						hits.get(readKey).add(new Hit(rp, i));
					}
				}
			}
		}
		else
		{
			readLength = (int)(readLength * 2.0 / (w+1));
			long[] miniKmers = KmerFinder.minimizers(r.s, k, w, posStrandBits);
			for(long miniKmer : miniKmers)
			{
				int strand = (int)(miniKmer&1);
				int i = (int) (miniKmer & ((1L << (posStrandBits)) - 1));
				long kmer = miniKmer >> posStrandBits;
				if(badKmers.contains(kmer)) continue;
				LongList currentHits = kmerMap[(int)(kmer&((1<<logNumMaps)-1))].get((int)(kmer>>logNumMaps));
				if(currentHits == null) continue;
				for(int hitIndex = 0; hitIndex<currentHits.size; hitIndex++)
				{
					ReadPosition rp = ReadPosition.decode(currentHits.a[hitIndex]);
					if(rp.rs%2 != strand) rp.rs |= 1;
					else if(rp.rs%2 == 1) rp.rs--;
					int readKey = rp.rs;
					if(!hits.containsKey(readKey))
					{
						hits.put(readKey, new ArrayList<>());
					}
					hits.get(readKey).add(new Hit(rp, i));
				}
			}
		}
		
		for(int readKey : hits.keySet())
		{
			ArrayList<Hit> sharedKmers = hits.get(readKey);
			if(sharedKmers.size() < threshold * readLength)
			{
				continue;
			}
			int readIndex = readKey / 2;
			int theirStrand = readKey % 2;
			
			if(theirStrand != 0)
			{
				Collections.reverse(sharedKmers);
			}
			
			int[] matchChain = lis(readLength, sharedKmers, theirStrand == 0);
			
			boolean chainContaining = chainContaining(sharedKmers, matchChain, readLength, threshold, endLength);
			if(chainContaining)
			{
				countContaining[readIndex]++;
				if(verbose)
				{
					System.err.println("Read " + r.i + " contained by long read " + readIndex + " on strand " + theirStrand);
				}
				return true;
			}
		}
		return false;
	}
	boolean chainContaining(ArrayList<Hit> sharedKmers, int[] matchChain, int readLength, double threshold, int endLength)
	{
		int countLeftEnd = 0, countRightEnd = 0;
		int totalMatches = matchChain.length;
		if(totalMatches < threshold * readLength)
		{
			return false;
		}
		for(int index : matchChain)
		{
			Hit cur = sharedKmers.get(index);
			if(cur.myIndex < endLength) countLeftEnd++;
			if(cur.myIndex + k + endLength > readLength) countRightEnd++;
		}
		if(countLeftEnd == 0 || countRightEnd == 0) return false;
		return true;
	}
	// Gets the sequence of indices in the longest increasing subsequence
	// Kmers near the end count as 5 matches
	int[] lis(int readLength, ArrayList<Hit> sharedKmers, boolean increasing)
	{
		int n = sharedKmers.size();
		int[] maxVal = new int[n];
		int[] backPointer = new int[n];
		int bestEnd = 0;
		for(int i = 0; i<n; i++)
		{
			Hit cur = sharedKmers.get(i);
			int currentVal = 1;
			if(cur.myIndex < clp.el || cur.myIndex + clp.el + k > readLength)
			{
				currentVal = 5;
			}
			maxVal[i] = currentVal;
			backPointer[i] = -1;
			for(int j = i-1; j>=0; j--)
			{
				Hit last = sharedKmers.get(j);
				boolean validTransition = (increasing && cur.rp.p > last.rp.p) || (!increasing && cur.rp.p < last.rp.p);
				int myJump = Math.abs(cur.myIndex - last.myIndex);
				int theirJump = Math.abs(cur.rp.p - last.rp.p);
				validTransition &= theirJump >= .9 * myJump && theirJump <= 1.1 * myJump;
				if(validTransition && maxVal[j] + currentVal > maxVal[i])
				{
					backPointer[i] = j;
					maxVal[i] = maxVal[j] + currentVal;
				}
			}
			if(maxVal[i] > maxVal[bestEnd])
			{
				bestEnd = i;
			}
		}
		ArrayList<Integer> indices = new ArrayList<>();
		while(bestEnd != -1)
		{
			indices.add(bestEnd);
			bestEnd = backPointer[bestEnd];
		}
		int[] res = new int[indices.size()];
		for(int i = 0; i<res.length; i++)
		{
			res[i] = indices.get(indices.size()-1-i);
		}
		return res;
	}
	void add(int index, Read cur)
	{
		if(!minimizers)
		{
			long[][] kmers = KmerFinder.kmerize(cur.s, k);
			for(int i = 0; i<kmers[0].length; i++)
			{
				ReadPosition forward = new ReadPosition(index, i);
				addKeyValue(kmers[0][i], forward);
			}
		}
		else
		{
			long[] miniKmers = KmerFinder.minimizers(cur.s, k, w, posStrandBits);
			for(long miniKmer : miniKmers)
			{
				int strand = (int)(miniKmer&1);
				int i = (int) (miniKmer & ((1L << (posStrandBits)) - 1));
				long kmer = miniKmer >> posStrandBits;
				ReadPosition toAdd = new ReadPosition(index, i, strand);
				addKeyValue(kmer, toAdd);
			}
		}
	}
	void addKeyValue(long fullKey, ReadPosition value)
	{
		if(badKmers.contains(fullKey))
		{
			return;
		}
		HashMap<Integer, LongList> addingTo = kmerMap[(int)(fullKey & ((1<<logNumMaps)-1))];
		int key = (int)(fullKey >> logNumMaps);
		if(!addingTo.containsKey(key)) addingTo.put((int)key, new LongList());
		if(addingTo.get(key).size >= 100)
		{
			addingTo.remove(key);
			badKmers.add(fullKey);
			return;
		}
		addingTo.get(key).add(value.encode());
	}
	class Hit
	{
		ReadPosition rp;
		int myIndex;
		Hit(ReadPosition readPos, int index)
		{
			rp = readPos;
			myIndex = index;
		}
	}
}
