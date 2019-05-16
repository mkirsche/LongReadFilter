package sieve;

public class CommandLineParser
{
	String fn = "/home/mkirsche/github/ContainedReadRemoval/sim/ERR2173373.fastq";
	//String fn = "/home/mkirsche/github/ContainedReadRemoval/sim/simulatedreads.fa";
	String ofn = "/home/mkirsche/github/ContainedReadRemoval/sim/sievedreads.out";
	double indexSize = 0.02; // What proportion of the reads should be in the index
	String readSplitScript = "/home/mkirsche/github/LongReadFilter/split_reads.sh";
	String uncontainedReadFile = "uncontainedreadnames.txt";
	//String logfile = "log.txt";
	String logfile = null;
	int k = 15;
	int w = 9;
	double p = 0.01;
	int el = (1+w) * 50;
	int nt = 4;
	boolean verbose = false;
	int minLength = 10000;
	CommandLineParser(String[] args)
	{
		boolean setEl = false;
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
				setEl = true;
				el = Integer.parseInt(val);
			}
			else if(argName.equals("lf"))
			{
				logfile = val;
			}
			else if(argName.equals("p"))
			{
				p = Double.parseDouble(val);
			}
			else if(argName.equals("rss"))
			{
				readSplitScript = val;
			}
			else if(argName.equals("urf"))
			{
				uncontainedReadFile = val;
			}
			else if(argName.equals("ml"))
			{
				minLength = Integer.parseInt(val);
			}
		}
		if(!setEl)
		{
			el = (1+w) * 50;
		}
	}
}
