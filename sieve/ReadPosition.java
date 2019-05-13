package sieve;

public class ReadPosition implements Comparable<ReadPosition> 
{
	int rs;
	int p;
	ReadPosition(int readIndex, int readPos)
	{
		rs = readIndex*2;
		p = readPos;
	}
	ReadPosition(int readIndex, int readPos, int strand)
	{
		rs = readIndex * 2 + strand;
		p = readPos;
	}
	public int compareTo(ReadPosition o)
	{
		if(rs != o.rs) return rs - o.rs;
		return p - o.p;
	}
	static ReadPosition decode(long x)
	{
		int rs = (int)(x & ((1<<30) - 1));
		int p = (int)(x>>30);
		ReadPosition res = new ReadPosition(rs/2, p);
		res.rs |= rs%2;
		return res;
	}
	long encode()
	{
		return ((long)p << 30) + rs;
	}
}
