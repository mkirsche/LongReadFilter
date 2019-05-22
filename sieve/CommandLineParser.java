package sieve;

public class CommandLineParser
{
	static boolean localDebug = true;
	String fn = localDebug ? "/home/mkirsche/github/ContainedReadRemoval/sim/simulatedreads.fa"
			: "/home/mkirsche/github/ContainedReadRemoval/sim/ERR2173373.fastq";
	String ofn = localDebug ? "/home/mkirsche/github/ContainedReadRemoval/sim/sievedreads.out"
			: "/home/mkirsche/github/ContainedReadRemoval/sim/sievedarabareads.out";
	double indexSize = 0.02; // What proportion of the reads should be in the index
	String readSplitScript = "/home/mkirsche/github/LongReadFilter/split_reads.sh";
	String uncontainedReadFile = localDebug ? "uncontainedreadnames.txt"
			: "arabauncontainedreadnames.txt";
	String logfile = localDebug ? "log.txt" : null;
	int k = 15;
	int w = 11;
	
	// These two are overridden if learn is set to true, which it is by default
	double p = 0.01;
	int el = 500;
	int nt = 4;
	
	double dpCutoff = .5;
	
	boolean verbose = false;
	boolean learn = true;
	int minLength = 12000;
	
	// Estimated proportion of short reads we want to label as uncontained - used for learning threshold
	double propUncontained = .08;
	
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
				else if(s.equals("-mt") || s.equals("-manualthreshold"))
				{
					learn = false;
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
			else if(argName.equals("puc"))
			{
				propUncontained = Double.parseDouble(val);
			}
		}
		if(!setEl)
		{
			setEl();
		}
	}
	void setEl()
	{
		el = (int)(.6 * (1+w) / p);
	}
}
