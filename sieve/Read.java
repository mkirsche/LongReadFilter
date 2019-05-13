package sieve;

import java.io.IOException;
import java.util.ArrayList;

public class Read 
{
	String n, s;
	int i;
	Read(String name, String seq, int index)
	{
		n = name;
		s = seq;
		i = index;
	}
	Read[] getFromFile(String fn, int max) throws IOException
	{
		ReadReader rr = new ReadReader(fn);
		ArrayList<Read> res = new ArrayList<>();
		int index = 0;
		while(rr.hasNext())
		{
			String name = rr.getNextName();
			String read = rr.getNextRead();
			res.add(new Read(name, read, index));
			index++;
		}
		Read[] toArray = new Read[res.size()];
		for(int i = 0; i<res.size(); i++) toArray[i] = res.get(i);
		return toArray;
	}
}
