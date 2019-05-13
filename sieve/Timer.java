package sieve;

public class Timer {
	static enum Unit {
		SECOND,
		MILLISECOND,
		MINUTE
	}
	long startTime;
	long elapsedTime;
	Unit u;
	Timer(Unit u)
	{
		elapsedTime = 0;
		this.u = u;
		startTime = System.currentTimeMillis();
	}
	void pause()
	{
		elapsedTime += System.currentTimeMillis() - startTime;
	}
	void resume()
	{
		startTime = System.currentTimeMillis();
	}
	String time()
	{
		return formatTime(elapsedTime + System.currentTimeMillis() - startTime);
	}
	String formatTime(long time)
	{
		if(u == Unit.MILLISECOND)
		{
			return String.format("(Time: %ldms)", time);
		}
		else if(u == Unit.SECOND)
		{
			return String.format("(Time: %.3fs)", time / 1000.0);
		}
		else if(u == Unit.MINUTE)
		{
			double seconds = time / 1000.0;
			return String.format("(Time %d:%02.3f)", (int)(seconds/60), seconds - 60*(int)(seconds/60));
		}
		else return null;
	}
}
