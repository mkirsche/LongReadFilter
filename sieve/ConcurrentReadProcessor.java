package sieve;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicInteger;

public class ConcurrentReadProcessor {
	
	ConcurrentLinkedQueue<Read> toProcess;
	ConcurrentLinkedQueue<String> toWriteName;
	ReadLengthSeparator re;
	ReadIndex index;
	CommandLineParser clp;
	AtomicInteger readsProcessed;
	AtomicInteger countContained;
	Timer timer;
	
	MyThread[] threads;
	WriterThread wt;
	boolean[] contained;
	
	Logger logger;
	
	ConcurrentReadProcessor(ReadLengthSeparator re, ReadIndex index, CommandLineParser clp, int numThreads, Timer timer) throws InterruptedException, FileNotFoundException
	{
		this.timer = timer;
		this.re = re;
		this.index = index;
		this.clp = clp;
		if(clp.logfile != null && clp.logfile.length() > 0)
		{
			logger = new Logger(clp.logfile);
		}
		toProcess = new ConcurrentLinkedQueue<Read>();
		toWriteName = new ConcurrentLinkedQueue<String>();
		threads = new MyThread[numThreads];
		readsProcessed = new AtomicInteger(0);
		countContained = new AtomicInteger(0);
		contained = new boolean[re.n];
		for(int i = 0; i<numThreads; i++)
		{
			threads[i] = new MyThread();
			threads[i].start();
		}
		wt = new WriterThread(clp.uncontainedReadFile);
		wt.start();
		System.err.println("All threads launched " + timer.time());
	}
	
	void run() throws InterruptedException
	{
		while(re.rr.hasNext())
		{
			while(toProcess.size() > 100)
			{
				Thread.sleep(1000);
			}
			Read cur = re.nextShortRead();
			if(cur == null)
			{
				break;
			}
			toProcess.add(cur);
		}
		for(String s : index.longReadNames)
		{
			toWriteName.add(s);
		}
		while(true)
		{
			if(readsProcessed.get() == re.n - index.n)
			{
				for(int i = 0; i<threads.length; i++)
				{
					threads[i].done = true;
				}
				break;
			}
			Thread.sleep(1000);
		}
		for(int i = 0; i<threads.length; i++) threads[i].join();
		while(!toWriteName.isEmpty())
		{
			Thread.sleep(1000);
		}
		wt.done = true;
		wt.join();
		System.err.println("All threads finished - outputting uncontained reads");
	}

	class MyThread extends Thread
	{
		boolean done;
		public MyThread()
		{
			done = false;
		}
		public void run() 
		{
			 try 
			 {
				 Read cur = null;
				 while(!done)
				 {
					 cur = toProcess.poll();
					 if(cur != null)
					 {
						 boolean c = index.contains(cur, logger);
						 int cc = countContained.get();
						 if(c)
						 {
							 countContained.incrementAndGet();
							 cc++;
							 contained[cur.i] = true;
						 }
						 else
						 {
							 toWriteName.add(cur.n);
						 }
						 int rp = readsProcessed.incrementAndGet();
						 if(rp%5000 == 0)
						 {
							 System.err.println("So far " + cc + " contained out of " + rp + " " + timer.time());
						 }
					 }
				 }
			 }
			 catch(Exception e) 
			 {
				 System.err.println("Error: " + e.getMessage());
				 e.printStackTrace();
			 }
		}
	}
	
	class WriterThread extends Thread
	{
		boolean done;
		PrintWriter out;
		public WriterThread(String fn) throws FileNotFoundException
		{
			out = new PrintWriter(new File(fn));
			done = false;
		}
		public void run() 
		{
			try 
			{
				String cur = null;
				while(!done)
				{
					cur = toWriteName.poll();
					if(cur != null)
					{
						String name = cur.split(" ")[0];
						out.println(name);
					}
					else
					{
						Thread.sleep(1000);
					}
				}
				out.close();
			}
			catch(Exception e) 
			{
				System.err.println("Error: " + e.getMessage());
				e.printStackTrace();
			}
		}
	}
}
