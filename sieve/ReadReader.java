/*
 * A class which facilitates reading genomic read data from either a FASTA or FASTQ file
 */

package sieve;

import java.util.*;
import java.io.*;

public class ReadReader {
	Scanner input;
	String last;
	int readCount;
	boolean lastFasta = false;
	boolean lastFastq = false;
	ReadReader(String fn) throws IOException
	{
		readCount = 0;
		input = new Scanner(new FileInputStream(new File(fn)));
	}
	String getNextName()
	{
		lastFastq = false;
		lastFasta = false;
		while(input.hasNext() && !validNameLine(last))
		{
			last = input.nextLine();
		}
		if(validNameLine(last))
		{
			String res = last;
			last = null;
			if(res.charAt(0) == '@')
			{
				lastFastq = true;
			}
			else if(res.charAt(0) == '>')
			{
				lastFasta = true;
			}
			
			return res.substring(1);
		}
		else
		{
			return null;
		}
	}
	String getNextRead()
	{
		readCount++;
		if(lastFastq)
		{
			String res = input.nextLine();
			input.nextLine();
			input.nextLine();
			return res;
		}
		else if(lastFasta)
		{
			StringBuilder res = new StringBuilder("");
			while(input.hasNext())
			{
				String s = input.nextLine();
				if(validNameLine(s))
				{
					last = s;
					break;
				}
				else
				{
					res.append(s);
				}
			}
			return res.toString();
		}
		else
		{
			return null;
		}
	}
	boolean validNameLine(String s)
	{
		if(s == null || s.length() == 0) return false;
		if(s.charAt(0) == '>' || s.charAt(0) == '@') return true;
		return false;
	}
	boolean hasNext()
	{
		return last != null || input.hasNext();
	}
	
	Read[] getAllReads()
	{
		ArrayList<Read> list = new ArrayList<Read>();
		while(hasNext())
		{
			String nextName = getNextName();
			String nextRead = getNextRead();
			list.add(new Read(nextName, nextRead, readCount-1));
		}
		int n = list.size();
		Read[] res = new Read[n];
		for(int i = 0; i<n; i++) res[i] = list.get(i);
		return res;
	}
}
