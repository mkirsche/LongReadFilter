/*
 * Uses a sample of reads to learn what thresholds will work well for a given dataset
 * Assumes that k and w are fixed, and learns other parameters from there
 */
package sieve;

import java.util.Arrays;
import java.util.Comparator;

public class ParameterLearner {
	
	double sharedCutoff;
	double dpCutoff;

	ParameterLearner(ReadIndex ri, Read[] sample, double propUncontained)
	{
		double[][] data = ri.getParamInfo(sample);
		Integer[] sortedByPropShared = sortByIndex(data, 0);
		sharedCutoff = data[sortedByPropShared[(int)(sample.length * propUncontained)]][0];
		double error = 1 - Math.pow(data[sortedByPropShared[(int)(sample.length * propUncontained)]][0], .5 / ri.k);
		double expectedDifferentBases = error*2 - error*error;
		dpCutoff = 1 - 2 * expectedDifferentBases;
		dpCutoff = dpCutoff * dpCutoff;
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
