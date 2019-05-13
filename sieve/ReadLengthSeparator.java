/*
 * Given a read file, sorts out the large reads for building an index
 * The length threshold used is calculated such that a given proportion of reads is used in the index
 * In addition, this class allows reading the file either one read at a time or skipping large reads
 * 
 */

package sieve;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;

public class ReadLengthSeparator {
	ReadReader rr;
	Read[] data;
	int lengthThreshold;
	int n;
	String fn;
	ReadLengthSeparator(String fn, double prop) throws IOException
	{
		this.fn = fn;
		n = 0;
		System.err.println("Calculating length threshold");
		lengthThreshold = getLengthThreshold(new ReadReader(fn), prop);
		System.err.println("Building index with reads having length at least: " + lengthThreshold);
		rr = new ReadReader(fn);
		ArrayList<Read> res = new ArrayList<>();
		while(rr.hasNext())
		{
			Read r = getNextRead();
			if(r == null)
			{
				break;
			}
			if(r.s.length() >= lengthThreshold)
			{
				res.add(r);
			}
		}
		Read[] toArray = new Read[res.size()];
		for(int i = 0; i<res.size(); i++) toArray[i] = res.get(i);
		data = toArray;
		
		n = rr.readCount;
		rr = new ReadReader(fn);
	}
	Read nextShortRead()
	{
		Read res = null;
		while(rr.hasNext() && res == null)
		{
			Read r = getNextRead();
			if(r.s.length() < lengthThreshold)
			{
				res = r;
			}
		}
		return res;
	}
	int getLengthThreshold(ReadReader rr, double p)
	{
		ArrayList<Integer> lengths = new ArrayList<>();
		while(rr.hasNext())
		{
			rr.getNextName();
			String read = rr.getNextRead();
			lengths.add(read.length());
		}
		
		Collections.sort(lengths);
		int keepReads = (int) (p * lengths.size());
		int cutoff = lengths.get(lengths.size() - keepReads) + 1;
		return cutoff;
	}
	Read getNextRead()
	{
		if(!rr.hasNext()) return null;
		String name = rr.getNextName();
		String read = rr.getNextRead();
		return new Read(name, read, rr.readCount-1);
	}
	void printUncontainedReads(boolean[] contained, String ofn) throws IOException
	{
		long totalLen = 0, uncontainedLen = 0;
		int numReads = 0, numUncontained = 0;
		PrintWriter out = new PrintWriter(new File(ofn));
		rr = new ReadReader(fn);
		int idx = 0;
		while(rr.hasNext())
		{
			Read r = getNextRead();
			totalLen += r.s.length();
			numReads++;
			if(!contained[idx])
			{
				uncontainedLen += r.s.length();
				numUncontained++;
				out.println(">"+r.n);
				out.println(r.s);
			}
			idx++;
		}
		out.close();
		
		System.err.printf("Kept %d out of %d reads (%.2f%%)\n", numUncontained, numReads, 100.0 * numUncontained / numReads);
		System.err.printf("Kept %d out of %d bases (%.2f%%)\n", uncontainedLen, totalLen, 100.0 * uncontainedLen / totalLen);
	}
}
