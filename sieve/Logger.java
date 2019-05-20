package sieve;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.concurrent.ConcurrentLinkedQueue;

public class Logger {
	String ofn;
	ConcurrentLinkedQueue<LogElement> data;
	Logger(String ofn)
	{
		this.ofn = ofn;
		data = new ConcurrentLinkedQueue<LogElement>();
	}
	
	public void print() throws IOException
	{
		PrintWriter out = new PrintWriter(new File(ofn));
		for(LogElement le : data)
		{
			out.println(le);
		}
		out.close();
	}
	
	static class LogElement
	{
		String readName;
		int longestChain;
		int readLength;
		int numCandidates;
		int numMinimizers;
		int[] chain;
		int[] theirChain;
		boolean contained;
		String containingName;
		int leftEnd, rightEnd;
		LogElement()
		{
			
		}
		public String toString()
		{
			return readName + "\t" + longestChain + "\t" + Arrays.toString(chain) + "\t" + readLength + "\t" 
					+ numCandidates + "\t" + numMinimizers + "\t" + contained + "\t" + leftEnd + "\t" + rightEnd + "\t"
					+ containingName + "\t" + Arrays.toString(theirChain);
		}
	}
}
