package sieve;

import java.io.IOException;
import java.util.Arrays;

public class Sieve {
public static void main(String[] args) throws IOException, InterruptedException
{
	Timer timer = new Timer(Timer.Unit.SECOND);
	CommandLineParser clp = new CommandLineParser(args);
	ReadLengthSeparator re = new ReadLengthSeparator(clp.fn, clp.indexSize, timer, clp.readSplitScript);
	System.err.println(timer.time());
	ReadIndex index = new ReadIndex(re, clp);
	
	ConcurrentReadProcessor crp = new ConcurrentReadProcessor(re, index, clp, clp.nt, timer);
	crp.run();
	
	re.printUncontainedReads(crp.contained, clp.ofn);
	System.err.println(Arrays.toString(index.countContaining));
	System.err.println(timer.time());
}
}
