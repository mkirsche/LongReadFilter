package hydroplane;

public class CommandLineParser
{
	// localDebug is whether or not the software is running on simulated data for debugging
	static boolean localDebug = true;
	
	// fn is the default input read file name
	String fn = localDebug ? "/home/mkirsche/git/ContainedReadRemoval/sim/simulatedreads.fa"
			: "/home/mkirsche/reads/ERR2173373.fastq";
	
	// ofn is the default output read file name
	String ofn = localDebug ? "/home/mkirsche/git/ContainedReadRemoval/sim/sievedreads.out"
			: "/home/mkirsche/git/ContainedReadRemoval/sim/sievedarabareads.out";
	
	// indexSize is the proportion of reads stored be in the index
	double indexSize = 0.02;
	
	// readSplitScript is the path to the awk script for splitting indexed vs. non-indexed reads
	String readSplitScript = "/home/mkirsche/git/LongReadFilter/split_reads.sh";
	
	// uncontainedReadFile is the output file where names of uncontained reads are written
	String uncontainedReadFile = localDebug ? "uncontainedreadnames.txt"
			: "arabauncontainedreadnames.txt";
	
	// logFile is the name of the file used for outputting logging information
	String logfile = localDebug ? "log.txt" : null;
	
	// Size of kmers
	int k = 15;
	
	// Window size for minimizers
	int w = 11;
	
	// These are overridden if learn is set to true, which it is by default
	// p is the shared kmer proportion cutoff
	// el is how close kmer chains must get to both ends of a read to call it contained
	// dpCutoff is the alignment score necessary in the ends of a read to call a 
	// read contained which didn't meet the above cutoffs
	double p = 0.01;
	int el = 500;
	double dpCutoff = .5;
	
	// nt is the number of threads to use
	int nt = 4;
	
	// verbose is whether or not to output additional logging info
	boolean verbose = false;
	
	// learn is whether or not to let the reads inform thresholds
	boolean learn = true;
	
	// minLength is the read length below which all reads are called as contained
	int minLength = 12000;
	
	// propUncontained is the estimated proportion of short reads we want to label as uncontained
	double propUncontained = .08;
	
	// These parameters are not (currently) customizable
	int logNumMaps = 16;
	int posStrandBits = 23;
	int maxAttempts = 5;
	
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
			else if(argName.equalsIgnoreCase("localDebug"))
			{
				if(val.equalsIgnoreCase("true") || val.equalsIgnoreCase("t"))
				{
					localDebug = true;
				}
				if(val.equalsIgnoreCase("false") || val.equalsIgnoreCase("f"))
				{
					localDebug = false;
				}
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
