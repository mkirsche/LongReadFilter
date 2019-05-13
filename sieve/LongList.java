package sieve;

public class LongList {
	long[] a;
	int size;
	LongList()
	{
		size = 0;
		a = new long[2];
	}
	void add(long x)
	{
		if(a.length == size)
		{
			long[] na = new long[a.length*2];
			for(int i = 0; i<a.length; i++) na[i] = a[i];
			a = na;
		}
		a[size++] = x; 
	}
}
