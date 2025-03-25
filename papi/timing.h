#ifndef TIMING_H
#define TIMING_H
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <papi.h>
#include <string.h>


// PAPI counters management
#define CYCLES 0
#define CACHE_MISS_L1 1
#define CACHE_MISS_L2 2
#define N_CTR 3

static const char *EVENT_NAMES[] = {
#ifdef CYCLES
	"CPU_CLK_UNHALTED:excl=0",
#endif
#ifdef CACHE_MISS_L1
	"PAPI_L1_DCM",
#endif
#ifdef CACHE_MISS_L2
	"PAPI_L2_DCM",
#endif
};

void print_counters(long long* ctrs);


// Timing function
typedef struct timing {
	int64_t cycles;
	int64_t ns;
} timing_t;

long get_ns();
long get_ms();

timing_t get_timestamp();
timing_t delta(timing_t t_before, timing_t t_after);
timing_t timing_add(timing_t t1, timing_t t2);
timing_t timing_div(timing_t t, int rep);
char* from_timing(timing_t);
char* from_timing_with_peak(timing_t, int64_t );

// papi functions
typedef struct papi_info {
	int eventset;
	size_t num_events; 
	char ** event_names;
} papi_info_t;

int compare_timing(const void* t1, const void* t2);

// Papi Timing function
timing_t papi_get_timestamp();
void init_papi();
papi_info_t build_papi_info();
void record_events(papi_info_t info);
void retrieve_results(papi_info_t papi_info, long long* results);
void show_all_records(papi_info_t papi_info, long long* papi_values);

#endif
