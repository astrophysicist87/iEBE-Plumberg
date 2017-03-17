#ifndef STATS_H
#define STATS_H

#include <vector>
#include <algorithm>

using namespace std;

int place_in_range(int idx, int minidx, int maxidx);
double get_median(double * array, int length);
double get_median_absolute_deviation(double * array, int length, double * median);
void find_outliers_Hampel(double * array, int length, bool * results, double * med, double Hampel_factor);
void find_outliers_window_Hampel(double * array, int length, bool * results, double * meds, double Hampel_factor, int width);

#endif

//End of file
