/* Function declaration for Pyrex with C linkage 
 * 
 * These have equivalents in adda.h where they are
 * declared as 'extern "C"'
 * 
 * */

int radar_run_from_files(const char * seq,
			 const char * mali,
			 const char * lfasta1,
			 const char * lfasta2,
			 int random_seed);

void radar_setLogLevel(int f);
