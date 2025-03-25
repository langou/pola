#include "timing.h"
#include <assert.h>



void print_counters(long long* ctrs) {
#ifdef CYCLES
  printf("%lld\n", ctrs[CYCLES]);
#endif
#ifdef CACHE_MISS_L1
  printf("%lld\n", ctrs[CACHE_MISS_L1]);
#endif
#ifdef CACHE_MISS_L2
  printf("%lld\n", ctrs[CACHE_MISS_L2]);
#endif
}

long get_ns() {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (long)ts.tv_sec * 1000000000 + ts.tv_nsec;
}

long get_ms() {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (long)ts.tv_sec * 1000000 + ts.tv_nsec / 1000;
}


unsigned long get_cycle() {
  unsigned long result ;
  unsigned long result2;
  asm volatile("rdtsc":"=a" (result), "=d" (result2));
  return ((unsigned long long)result2 << 32) + result;
}

timing_t get_timestamp() {
	timing_t res;
	unsigned long long dtime, dtime_ns;
	dtime_ns = get_ns();
	dtime = get_cycle();
	res.ns = dtime_ns;
	res.cycles = dtime;
	return res;
}

char* from_timing(timing_t timing) {
	char* str = (char *) malloc(1000);
	float time = (float)timing.ns / 1000000000;
	if (sprintf(str, "{time : %e sec, %ld (%e) cycles}", time, timing.cycles, (float)timing.cycles) < 0) exit(0);
	return str;
}


char* from_timing_with_peak(timing_t timing, int64_t peak) {
	char* str = (char *) malloc(1000);
	float time = (float)timing.ns / 1000000000;
	if (sprintf(str, "{time : %e sec, %ld (%e) cycles, ratio peak/perf: %f %%}", time, timing.cycles, (float)timing.cycles, 100 * (float) peak / timing.cycles) < 0) exit(0);
	return str;
}

timing_t timing_add(timing_t t1, timing_t t2) {
  timing_t res;
  res.ns = t1.ns + t2.ns;
  res.cycles = t1.cycles + t2.cycles;
  return res;
}

timing_t delta(timing_t t_before, timing_t t_after) {
	timing_t res;
	unsigned long long dcycles, dns;
	dns = t_after.ns - t_before.ns;
	dcycles = t_after.cycles - t_before.cycles;
	res.cycles = dcycles;
	res.ns = dns;
	return res;
}

timing_t timing_div(timing_t t, int rep) {
	timing_t res;
	res.cycles = t.cycles / rep;
	res.ns = t.ns / rep;
	return res;
}

int compare_timing(const void* t1, const void* t2) {
  unsigned long t1_ns, t2_ns;
  t1_ns =((timing_t *)t1)->ns; 
  t2_ns =((timing_t *)t2)->ns; 
  if ( t1_ns < t2_ns) {
    return -1;
  }
  else if (t1_ns == t2_ns) {
    return 0;
  }
  else {
    return 1;
  }
}

/// statement that GCC (currently) cannot reorder around
#define COMPILER_BARRIER() asm volatile("" ::: "memory");
void init_papi() {
	const int retval = PAPI_library_init(PAPI_VER_CURRENT);
	if (retval != PAPI_VER_CURRENT) {
		const char *err = PAPI_strerror(retval);
		printf("%s:%d::PAPI_library_init failed: (%d) %s\n", __FILE__, __LINE__, retval, err);
		exit(1);
	}
}


void free_papi_info(papi_info_t info) {
  for (size_t i = 0; i < info.num_events; i++)
    free(info.event_names[i]);
  free(info.event_names);
}

papi_info_t build_papi_info() {
	papi_info_t papi_info;
  /*
	const char *EVENT_NAMES[] = {
#ifdef CYCLES
		"CPU_CLK_UNHALTED:excl=0",
#endif
#ifdef CACHE_MISS_L1
		"PAPI_L1_DCM",
#endif
#ifdef CACHE_MISS_L2
		"PAPI_L2_DCM",
#endif
#ifdef CACHE_MISS_L3
		"MEM_LOAD_RETIRED:L3_MISS",
#endif
#ifdef BR_PRED
		"PAPI_BR_PRC",
#endif
#ifdef BR_MISPRED
		"PAPI_BR_MSP",
#endif
#ifdef INST
		"PAPI_TOT_INS",
#endif
	};
  */
	const size_t num_events = sizeof(EVENT_NAMES) / sizeof(EVENT_NAMES[0]);
  assert(num_events == N_CTR);

	int papi_eventset = PAPI_NULL;

  papi_info.num_events = num_events;
  papi_info.event_names = malloc(sizeof(char*) * num_events);

	COMPILER_BARRIER();
	{
		// create event set
		{
			const int retval = PAPI_create_eventset(&papi_eventset);
			if (retval != PAPI_OK) {
				const char *err = PAPI_strerror(retval);
				printf("%s:%d::PAPI_create_eventset failed: (%d) %s\n", __FILE__, __LINE__, retval, err);
				exit(1);
			}
		}
		// add events
    for (size_t i = 0; i < num_events; i++) {
      const char *event_name = EVENT_NAMES[i];

      int event_code = 0;
      {
        const int retval = PAPI_event_name_to_code(event_name, &event_code);
        if (retval != PAPI_OK) {
          const char *err = PAPI_strerror(retval);
          printf("%s:%d::PAPI_event_name_to_code(\"%s\") failed: (%d) %s\n", __FILE__, __LINE__, event_name, retval, err);
          exit(1);
        }
      }
      {
        const int retval = PAPI_add_event(papi_eventset, event_code);
        if (retval != PAPI_OK) {
          const char *err = PAPI_strerror(retval);
          printf("%s:%d::PAPI_add_event(%d) (\"%s\") failed: (%d) %s\n", __FILE__, __LINE__, event_code, event_name, retval, err);
          exit(1);
        }
      }
      int len = strlen(event_name);
      char* event_n = malloc((len + 1) * sizeof(char));
      strcpy(event_n, event_name);
      papi_info.event_names[i] = event_n;
    }
	}
  papi_info.eventset = papi_eventset;
	return papi_info;
}

void record_events(papi_info_t info) {
  {
    const int retval = PAPI_start(info.eventset);
    if (retval != PAPI_OK) {
      const char *err = PAPI_strerror(retval);
      printf("%s:%d::PAPI_start failed: (%d) %s\n", __FILE__, __LINE__, retval, err);
      exit(1);
    }
  }
  COMPILER_BARRIER();
}

void retrieve_results(papi_info_t papi_info, long long* results) {
  COMPILER_BARRIER();

  long long int papi_values[64];
  {
    const int retval = PAPI_stop(papi_info.eventset, papi_values);
    if (retval != PAPI_OK) {
      const char *err = PAPI_strerror(retval);
      printf("%s:%d::PAPI_stop failed: (%d) %s\n", __FILE__, __LINE__, retval, err);
      exit(1);
    }
  }
  COMPILER_BARRIER();
  for (unsigned i = 0; i < papi_info.num_events; i++) {
   results[i] = papi_values[i];
  }
}

void show_all_records(papi_info_t papi_info, long long* papi_values) {
  for (unsigned i = 0; i < papi_info.num_events; i++) {
    printf("%-30s %15lli\n", papi_info.event_names[i], papi_values[i]);
  }
}

timing_t papi_get_timestamp() {
	timing_t res;
	unsigned long long cycles, ns;
	//ns = PAPI_get_real_nsec();
  ns = 100;
	cycles = PAPI_get_real_cyc();
	res.ns = ns;
	res.cycles = cycles;
	return res;
}

