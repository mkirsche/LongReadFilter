package sieve;

import java.io.IOException;
import java.util.*;

public class ReadIndex {
	int n;
	int k;
	int w;
	int[] countContaining;
	HashMap<Integer, LongList>[] kmerMap;
	HashSet<Long> badKmers;
	ArrayList<String> longReadNames;
	
	boolean verbose;
	CommandLineParser clp;
	int lengthThreshold;
	Read[] data;
	DynamicProgrammingAligner dpa;
	
	@SuppressWarnings("unchecked")
	ReadIndex(ReadLengthSeparator re, CommandLineParser clp) throws IOException
	{
		dpa = new DynamicProgrammingAligner();
		lengthThreshold = re.lengthThreshold;
		this.data = re.data;
		this.clp = clp;
		w = clp.w;
		verbose = clp.verbose;
		k = clp.k;
		n = re.data.length;
		countContaining = new int[n];
		kmerMap = new HashMap[1<<clp.logNumMaps];
		for(int i = 0; i<(1<<clp.logNumMaps); i++) kmerMap[i] = new HashMap<>();
		badKmers = new HashSet<>();
		longReadNames = new ArrayList<>();
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
		}
		System.err.println("Index built with " + n + " reads having total length " + totalReadLength);
		int numKmers = 0;
		for(HashMap<Integer, LongList> hm : kmerMap) numKmers += hm.size();
		System.err.println("Number of kmers: " + numKmers);
		System.err.println("Bad kmers: " + badKmers.size());
	}
	
	// Gets parameter information based on querying a sample of reads - used by ParameterLearner
	double[][] getParamInfo(Read[] rs)
	{
		int numReads = rs.length;
		double[][] res = new double[numReads][];
		for(int i = 0; i<numReads; i++)
		{
			Read r = rs[i];
			long[] kmers = MerMaker.minimizers(r.s, k, w, clp.posStrandBits);
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
	
	
	boolean contains(Read r, Logger logger) throws InterruptedException, IOException
	{
		Logger.LogElement le = new Logger.LogElement();
		double threshold = clp.p;
		int endLength = clp.el;
		int readLength = r.s.length();
		if(readLength < clp.minLength)
		{
			return true;
		}
		
		// Kmerize the query read
		long[] kmers = MerMaker.minimizers(r.s, k, w, clp.posStrandBits);
		int numMinimizers = kmers.length;
		
		// Get all kmer matches against any database reads
		HashMap<Integer, ArrayList<Hit>> hits = getHits(r, kmers);
		
		// Initialize logging information
		le.numCandidates = hits.keySet().size();
		le.readName = r.n;
		le.readLength = r.s.length();
		le.numMinimizers = kmers.length;
		le.longestChain = -1;
		le.leftEnd = -le.readLength;
		le.rightEnd = -le.readLength;
		le.contained = false;
		le.containingName = "none";
		le.chain = new int[] {};
		le.theirChain = new int[] {};
		le.ctScore = 0.0;
		
		int bestKey = -1;
		
		TreeMap<Integer, ArrayList<Integer>> scoreToKey = new TreeMap<>();
		
		for(int readKey : hits.keySet())
		{
			ArrayList<Hit> sharedKmers = hits.get(readKey);
			
			// If there are too few hits to possibly form a long enough chain, ignore
			if(sharedKmers.size() < threshold * numMinimizers)
			{
				continue;
			}
			
			// Get the read index and strand
			int readIndex = readKey / 2;
			int theirStrand = readKey % 2;

			// If the database read is too short to possibly contain the query, ignore
			if(readLength >= data[readIndex].s.length())
			{
				continue;
			}
			
			// Sort matches by their position in query read
			Collections.sort(sharedKmers);
			if(theirStrand != 0)
			{
				Collections.reverse(sharedKmers);
			}
			
			// Get matches which form longest subsequence of increasing positions in database read
			int[] matchChain = lis(r.s.length(), sharedKmers, theirStrand == 0, threshold);
			
			// Store how close this chain gets to end of query - 
			//   to be used when selecting which alignments to investigate further
			int[] endLengths = getUnalignedEnds(sharedKmers, matchChain, r.s.length());
			int maxEnd = Math.max(endLengths[0], endLengths[1]);
			if(!scoreToKey.containsKey(maxEnd))
			{
				scoreToKey.put(maxEnd, new ArrayList<>());
			}
			scoreToKey.get(maxEnd).add(readKey);
			
			int oldMax = Math.max(le.leftEnd, le.rightEnd);
			if(oldMax < 0) oldMax = le.readLength + 1;
			
			// If this gets close to the ends, update best alignment for potential logging
			if(!le.contained && maxEnd < oldMax)
			{
				bestKey = readKey;
				le.longestChain = matchChain.length;
				le.leftEnd = endLengths[0];
				le.rightEnd = endLengths[1];
				le.chain = new int[matchChain.length];
				le.theirChain = new int[matchChain.length];
				le.containingName = longReadNames.get(readKey/2);
				for(int i = 0; i<le.chain.length; i++)
				{
					le.chain[i] = sharedKmers.get(matchChain[i]).myIndex;
					le.theirChain[i] = sharedKmers.get(matchChain[i]).rp.p;
				}
			}
			
			// Check if the chain length/proximity to ends is good enough to call this contained
			// Don't do this if querying a database read since we want to be more careful about removing longer reads
			boolean chainContaining = r.s.length() < lengthThreshold &&
					chainContaining(sharedKmers, matchChain, numMinimizers, r.s.length(), threshold, endLength);
			
			// If it is contained, return true (or mark and keep going if logging results)
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
					if(maxEnd < oldMax || !le.contained)
					{
						bestKey = readKey;
						le.longestChain = matchChain.length;
						le.leftEnd = endLengths[0];
						le.rightEnd = endLengths[1];
						le.chain = new int[matchChain.length];
						le.theirChain = new int[matchChain.length];
						le.containingName = longReadNames.get(readKey/2);
						for(int i = 0; i<le.chain.length; i++)
						{
							le.chain[i] = sharedKmers.get(matchChain[i]).myIndex;
							le.theirChain[i] = sharedKmers.get(matchChain[i]).rp.p;
						}
						le.contained = true;
					}
				}
			}
		}
		
		// Now, if close matches were found, use dynamic programming to investigate further
		// Iterate in order of how close the kmer chains got to the ends of the query read
		int attempts = 0;
		int last = -1;
		le.dpNames = new ArrayList<String>();
		if(bestKey != -1 && (logger == null || !le.contained))
		{
			boolean found = false;
			while(!found && attempts < clp.maxAttempts && (last == -1 || scoreToKey.higherKey(last) != null))
			{
				int curScore = last == -1 ? scoreToKey.firstKey() : scoreToKey.higherKey(last);
				last = curScore;
				ArrayList<Integer> keys = scoreToKey.get(curScore);
				for(int i = 0; i<keys.size() && attempts < clp.maxAttempts && !found; i++, attempts++)
				{
					int curKey = keys.get(i);
					ArrayList<Hit> sharedKmers = hits.get(curKey);
					int theirStrand = curKey % 2;
					le.dpNames.add(data[curKey/2].n);
					int[] matchChain = lis(r.s.length(), sharedKmers, theirStrand == 0, threshold);
					int[] endLengths = getUnalignedEnds(sharedKmers, matchChain, r.s.length());
					
					// If this is a database read, add extra end length criteria
					if(r.s.length() >= lengthThreshold && endLength * 5 < Math.max(endLengths[0], endLengths[1]))
					{
						continue;
					}
					
					// Get alignment score and return true if high enough
					double score = dpa.dpContains(data[sharedKmers.get(0).rp.rs/2], k, r, sharedKmers, matchChain, 
							matchChain.length < threshold * numMinimizers || r.s.length() >= lengthThreshold);
					le.ctScore = Math.max(le.ctScore, score);
					if(score > clp.dpCutoff)
					{
						countContaining[curKey/2]++;
						found = true;
						if(logger != null)
						{
							le.contained = true;
						}
						else return true;
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
		for(long miniKmer : kmers)
		{
			int strand = (int)(miniKmer&1);
			int i = ((int) (miniKmer & ((1L << (clp.posStrandBits)) - 1))) >> 1;
			long kmer = miniKmer >> clp.posStrandBits;
			if(badKmers.contains(kmer)) continue;
			LongList currentHits = kmerMap[(int)(kmer&((1<<clp.logNumMaps)-1))].get((int)(kmer>>clp.logNumMaps));
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
				currentVal = 5;
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
				validTransition &= myJump * threshold < 25.0;
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
		long[] miniKmers = MerMaker.minimizers(cur.s, k, w, clp.posStrandBits);
		for(long miniKmer : miniKmers)
		{
			int strand = (int)(miniKmer&1);
			int i = ((int) (miniKmer & ((1L << (clp.posStrandBits)) - 1))) >> 1;
			long kmer = miniKmer >> clp.posStrandBits;
			ReadPosition toAdd = new ReadPosition(index, i, strand);
			addKeyValue(kmer, toAdd);
		}
	}
	void addKeyValue(long fullKey, ReadPosition value)
	{
		if(badKmers.contains(fullKey))
		{
			return;
		}
		HashMap<Integer, LongList> addingTo = kmerMap[(int)(fullKey & ((1<<clp.logNumMaps)-1))];
		int key = (int)(fullKey >> clp.logNumMaps);
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
