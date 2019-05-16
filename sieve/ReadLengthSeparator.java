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
import java.util.HashMap;
import java.util.Scanner;

public class ReadLengthSeparator {
	ReadReader rr;
	Read[] data;
	int lengthThreshold;
	int n;
	String fn;
	@SuppressWarnings("resource")
	ReadLengthSeparator(String fn, double prop, Timer timer, String readSplitScript) throws IOException, InterruptedException
	{
		this.fn = fn;
		n = 0;
		System.out.println(readSplitScript);
		if(readSplitScript == null || !(new File(readSplitScript)).exists())
		{
			System.err.println("Calculating length threshold");
			lengthThreshold = getLengthThreshold(new ReadReader(fn), prop);
			System.err.println("Building index with reads having length at least: " + lengthThreshold);
			System.err.println(timer.time());
			rr = new ReadReader(fn);
			fillDataFromReadReader();
		}
		else
		{
			try 
			{
				String command = readSplitScript + " " + fn + " " + prop;
				System.err.println("Running script to split reads: " + command);
				// TODO get output from command
				Runtime rt = Runtime.getRuntime();
				Process pr = rt.exec(command);
				Scanner prReader = new Scanner(pr.getInputStream()); 
				HashMap<String, String> prOutput = new HashMap<>();
				while (prReader.hasNext()) 
				{
					String line = prReader.nextLine();
					System.out.println(line);
					String[] tokens = line.split(" ");
					if(tokens.length != 2 || !tokens[0].endsWith(":")) continue;
					String key = tokens[0].substring(0, tokens[0].length()-1);
					String val = tokens[1];
					prOutput.put(key, val);
				}
				int retVal = pr.waitFor();
				if(retVal != 0)
				{
					System.err.println("Read splitting failed");
				}
				
				lengthThreshold = Integer.parseInt(prOutput.get("LengthCutoff"));
				
				System.err.println("Done splitting reads " + timer.time());
				rr = new ReadReader(fn + ".long");
				fillDataFromReadReader();
				n = Integer.parseInt(prOutput.get("NumReads"));
			} 
			catch(Exception e) 
			{
				System.err.println("Calculating length threshold");
				lengthThreshold = getLengthThreshold(new ReadReader(fn), prop);
				System.err.println("Building index with reads having length at least: " + lengthThreshold);
				System.err.println(timer.time());
				rr = new ReadReader(fn);
				fillDataFromReadReader();
			};
		}
		rr = new ReadReader(fn);
	}
	void fillDataFromReadReader()
	{
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
