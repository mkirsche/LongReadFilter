package sieve;

public class CommandLineParser
{
	String fn = "/home/mkirsche/github/ContainedReadRemoval/sim/ERR2173373.fastq";
	//String fn = "/home/mkirsche/github/ContainedReadRemoval/sim/simulatedreads.fa";
	String ofn = "/home/mkirsche/github/ContainedReadRemoval/sim/sievedreads.out";
	double indexSize = 0.01; // What proportion of the reads should be in the index
	int k = 15;
	int w = 21;
	double p = 0.005;
	int el = 500;
	int nt = 4;
	boolean verbose = false;
	CommandLineParser(String[] args)
	{
		for(String s : args)
		{
			int idx = s.indexOf('=');
			if(idx == -1)
			{
				if(s.equals("-v") || s.equals("-verbose"))
				{
					verbose = true;
				}
				continue;
			}
			String argName = s.substring(0, idx).toLowerCase();
			String val = s.substring(idx+1);
			if(argName.equals("fn"))
			{
				fn = val;
			}
			else if(argName.equals("ofn"))
			{
				ofn = val;
			}
			else if(argName.equals("indexsize"))
			{
				indexSize = Double.parseDouble(val);
				if(indexSize >= 1)
				{
					indexSize /= 100;
				}
			}
			else if(argName.equals("k"))
			{
				k = Integer.parseInt(val);
			}
			else if(argName.equals("w"))
			{
				w = Integer.parseInt(val);
			}
			else if(argName.equals("nt"))
			{
				nt = Integer.parseInt(val);
			}
			else if(argName.equals("el"))
			{
				el = Integer.parseInt(val);
			}
			else if(argName.equals("p"))
			{
				p = Double.parseDouble(val);
			}
		}
	}
}
