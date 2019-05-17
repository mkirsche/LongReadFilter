/*
 * Uses a sample of reads to learn what thresholds will work well for a given dataset
 * Assumes that k and w are fixed, and learns other parameters from there
 */
package sieve;

import java.util.Arrays;
import java.util.Comparator;

public class ParameterLearner {
	
	double sharedCutoff;

	ParameterLearner(ReadIndex ri, Read[] sample, double propUncontained)
	{
		double[][] data = ri.getParamInfo(sample);
		Integer[] sortedByPropShared = sortByIndex(data, 0);
		sharedCutoff = data[sortedByPropShared[(int)(sample.length * propUncontained)]][0];
	}
	
	static Integer[] sortByIndex(double[][] data, int index)
	{
		int n = data.length;
		Integer[] is = new Integer[n];
		for(int i = 0; i<n; i++) is[i] = i;
		Arrays.sort(is, new Comparator<Integer>() {

			public int compare(Integer a, Integer b) {
				return Double.compare(data[a][index], data[b][index]);
			}});
		return is;
	}
}
