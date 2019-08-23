/*
 * Used for testing - uses ground truth scoring file to check each read only against its most containing read
 * Useful for debugging reads which are being wrongly kept
 */

package sieve;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Scanner;
import java.util.TreeMap;

public class ScratchReadIndex extends ReadIndex {
	
	HashMap<String, Integer> map;
	int counter = 0;
	
	void initMap(String fn) throws IOException
	{
		HashMap<String, Integer> nameToIndex = new HashMap<>();
		for(int i = 0; i<n; i++)
		{
			nameToIndex.put(data[i].n, i);
		}
		map = new HashMap<>();
		@SuppressWarnings("resource")
		Scanner input = new Scanner(new FileInputStream(new File(fn)));
		while(input.hasNext())
		{
			String[] s = input.nextLine().split(" ");
			String key = s[0];
			String val = s[s.length-1];
			if(nameToIndex.containsKey(val))
			{
				map.put(key, nameToIndex.get(val));
			}
			else
			{
				map.put(key, -1);
			}
		}
	}
	
	ScratchReadIndex(ReadLengthSeparator re, CommandLineParser clp, String scorefn) throws IOException
	{
		super(re, clp);
		initMap(scorefn);
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
		
		int allowedIndex = -1;
		
		if(map.containsKey(r.n))
		{
			allowedIndex = map.get(r.n);
		}
		
		int count = 0;
		int att = 0;
		
		for(int readKey : hits.keySet())
		{
			// Get the read index and strand
			int readIndex = readKey / 2;
			int theirStrand = readKey % 2;
					
			if(readIndex != allowedIndex)
			{
				continue;
			}
			att++;
			
			ArrayList<Hit> sharedKmers = hits.get(readKey);
			
			// If there are too few hits to possibly form a long enough chain, ignore
			if(sharedKmers.size() < threshold * numMinimizers)
			{
				count++;
				continue;
			}
			
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
		
		if(count == att)
		{
			counter++;
			System.err.println(r.n+" "+numMinimizers+" "+count+" "+allowedIndex);
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
}
