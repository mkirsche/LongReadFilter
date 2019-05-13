package sieve;

import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicInteger;

public class ConcurrentReadProcessor {
	
	ConcurrentLinkedQueue<Read> toProcess;
	ReadLengthSeparator re;
	ReadIndex index;
	CommandLineParser clp;
	AtomicInteger readsProcessed;
	AtomicInteger countContained;
	Timer timer;
	
	MyThread[] threads;
	boolean[] contained;
	
	ConcurrentReadProcessor(ReadLengthSeparator re, ReadIndex index, CommandLineParser clp, int numThreads, Timer timer) throws InterruptedException
	{
		this.timer = timer;
		this.re = re;
		this.index = index;
		this.clp = clp;
		toProcess = new ConcurrentLinkedQueue<Read>();
		threads = new MyThread[numThreads];
		readsProcessed = new AtomicInteger(0);
		countContained = new AtomicInteger(0);
		contained = new boolean[re.n];
		for(int i = 0; i<numThreads; i++)
		{
			threads[i] = new MyThread();
			threads[i].start();
		}
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
						 boolean c = index.contains(cur);
						 int cc = countContained.get();
						 if(c)
						 {
							 countContained.incrementAndGet();
							 cc++;
							 contained[cur.i] = true;
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
}
