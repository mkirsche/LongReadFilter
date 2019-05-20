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
	ArrayList<String> longReadNames;
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
		longReadNames = new ArrayList<>();
		names = new String[n];
		long totalReadLength = 0;
		System.err.println("Building index with " + n + " reads");
		for(int i = 0; i<n; i++)
		{
			names[i] = re.data[i].n;
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
	
	double[][] getParamInfo(Read[] rs)
	{
		int numReads = rs.length;
		double[][] res = new double[numReads][];
		for(int i = 0; i<numReads; i++)
		{
			Read r = rs[i];
			long[] kmers = minimizers ? KmerFinder.minimizers(r.s, k, w, posStrandBits) : KmerFinder.kmers(r.s, k);
			HashMap<Integer, ArrayList<Hit>> hits = getHits(r, kmers);
			
			int longestChain = 0;
			
			for(int readKey : hits.keySet())
			{
				ArrayList<Hit> sharedKmers = hits.get(readKey);
				
				Collections.sort(sharedKmers);
				
				if(readKey % 2 != 0)
				{
					Collections.reverse(sharedKmers);
				}
				
				int[] matchChain = lis(r.s.length(), sharedKmers, readKey % 2 == 0, 0);

				if(matchChain.length > longestChain)
				{
					longestChain = matchChain.length;
				}
			}
			res[i] = new double[] {longestChain * 1.0 / kmers.length};
		}
		return res;
	}
	
	
	boolean contains(Read r, Logger logger)
	{
		Logger.LogElement le = new Logger.LogElement();
		double threshold = clp.p;
		int endLength = clp.el;
		int readLength = r.s.length();
		if(readLength < clp.minLength)
		{
			return true;
		}
		
		long[] kmers = minimizers ? KmerFinder.minimizers(r.s, k, w, posStrandBits) : KmerFinder.kmers(r.s, k);
		int numMinimizers = kmers.length;
		
		HashMap<Integer, ArrayList<Hit>> hits = getHits(r, kmers);
		
		le.numCandidates = hits.keySet().size();
		le.readName = r.n;
		le.readLength = r.s.length();
		le.numMinimizers = kmers.length;
		le.longestChain = 0;
		le.leftEnd = le.readLength;
		le.rightEnd = le.readLength;
		le.contained = false;
		le.containingName = "none";
		le.chain = new int[] {};
		
		for(int readKey : hits.keySet())
		{
			ArrayList<Hit> sharedKmers = hits.get(readKey);
			if(sharedKmers.size() < threshold * numMinimizers)
			{
				continue;
			}
			int readIndex = readKey / 2;
			int theirStrand = readKey % 2;
			
			Collections.sort(sharedKmers);
			
			if(theirStrand != 0)
			{
				Collections.reverse(sharedKmers);
			}
			
			int[] matchChain = lis(r.s.length(), sharedKmers, theirStrand == 0, threshold);
			
			if(!le.contained && matchChain.length > le.longestChain)
			{
				le.longestChain = matchChain.length;
				int[] endLengths = getUnalignedEnds(sharedKmers, matchChain, r.s.length());
				le.leftEnd = endLengths[0];
				le.rightEnd = endLengths[1];
				le.chain = new int[matchChain.length];
				le.containingName = names[readKey/2];
				for(int i = 0; i<le.chain.length; i++)
				{
					le.chain[i] = sharedKmers.get(matchChain[i]).myIndex;
				}
			}
			boolean chainContaining = chainContaining(sharedKmers, matchChain, numMinimizers, r.s.length(), threshold, endLength);
			if(chainContaining)
			{
				countContaining[readIndex]++;
				if(verbose)
				{
					System.err.println("Read " + r.i + " contained by long read " + readIndex + " on strand " + theirStrand);
				}
				if(logger == null)
				{
					return true;
				}
				else
				{
					if(matchChain.length > le.longestChain || !le.contained)
					{
						le.longestChain = matchChain.length;
						int[] endLengths = getUnalignedEnds(sharedKmers, matchChain, r.s.length());
						le.leftEnd = endLengths[0];
						le.rightEnd = endLengths[1];
						le.chain = new int[matchChain.length];
						le.containingName = names[readKey/2];
						for(int i = 0; i<le.chain.length; i++)
						{
							le.chain[i] = sharedKmers.get(matchChain[i]).myIndex;
						}
						le.contained = true;
					}
				}
			}
		}
		if(logger != null)
		{
			logger.data.add(le);
			return le.contained;
		}
		return false;
	}
	HashMap<Integer, ArrayList<Hit>> getHits(Read r, long[] kmers)
	{
		HashMap<Integer, ArrayList<Hit>> hits = new HashMap<>();
		
		if(!minimizers)
		{
			// Loop through both strands since only forward strand in index
			for(int strand = 0; strand<2; strand++)
			{
				for(int i = strand; i<kmers.length; i+=2)
				{
					long kmer = kmers[i];
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
			for(long miniKmer : kmers)
			{
				int strand = (int)(miniKmer&1);
				int i = ((int) (miniKmer & ((1L << (posStrandBits)) - 1))) >> 1;
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
		return hits;
	}
	int[] getUnalignedEnds(ArrayList<Hit> sharedKmers, int[] matchChain, int readLength)
	{
		int leftEnd = readLength, rightEnd = readLength;
		for(int index : matchChain)
		{
			Hit h = sharedKmers.get(index);
			leftEnd = Math.min(leftEnd, h.myIndex);
			rightEnd = Math.min(rightEnd, readLength - k - h.myIndex);
		}
		return new int[] {leftEnd, rightEnd};
	}
	boolean chainContaining(ArrayList<Hit> sharedKmers, int[] matchChain, int numMinimizers, int readLength, double threshold, int endLength)
	{
		int countLeftEnd = 0, countRightEnd = 0;
		int totalMatches = matchChain.length;
		if(totalMatches < threshold * numMinimizers)
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
	// Kmers near the end count as 50 matches
	int[] lis(int readLength, ArrayList<Hit> sharedKmers, boolean increasing, double threshold)
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
				currentVal = 50;
			}
			maxVal[i] = currentVal;
			backPointer[i] = -1;
			for(int j = i-1; j>=0; j--)
			{
				Hit last = sharedKmers.get(j);
				boolean validTransition = ((increasing && cur.rp.p > last.rp.p) || (!increasing && cur.rp.p < last.rp.p)) && (cur.myIndex != last.myIndex);
				int myJump = Math.abs(cur.myIndex - last.myIndex);
				int theirJump = Math.abs(cur.rp.p - last.rp.p);
				validTransition &= theirJump >= .8 * myJump && theirJump <= 1.2 * myJump;
				validTransition &= myJump * threshold < 30.0;
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
		longReadNames.add(cur.n);
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
				int i = ((int) (miniKmer & ((1L << (posStrandBits)) - 1))) >> 1;
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
	class Hit implements Comparable<Hit>
	{
		ReadPosition rp;
		int myIndex;
		Hit(ReadPosition readPos, int index)
		{
			rp = readPos;
			myIndex = index;
		}
		@Override
		public int compareTo(Hit o) {
			return myIndex - o.myIndex;
		}
	}
}
