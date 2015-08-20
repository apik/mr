/* Names of the various functions, including equivalent permutations. */

#define NUM_U_PERMS 2
#define NUM_T_PERMS 2
#define NUM_S_PERMS 6
#define NUM_B_PERMS 2
#define NUM_V_PERMS 2

const char *uname[4][2] = {{"Uzxyv","Uzxvy"},
			   {"Uuyxv","Uuyvx"},
			   {"Uxzuv","Uxzvu"},
			   {"Uyuzv","Uyuvz"}};

const char *tname[6][2] = {{"Tvyz","Tvzy"},
			   {"Tuxv","Tuvx"},
			   {"Tyzv","Tyvz"},
			   {"Txuv","Txvu"},
			   {"Tzyv","Tzvy"},
			   {"Tvxu","Tvux"}};

const char *sname[2][6] = {{"Svyz","Szvy","Syzv","Svzy","Syvz","Szyv"},
			   {"Suxv","Svux","Sxvu","Suvx","Sxuv","Svxu"}};

const char *bname[2][2] = {{"Bxz","Bzx"}, 
                           {"Byu","Buy"}};

const char *vname[4][2] = {{"Vzxyv","Vzxvy"},
			   {"Vuyxv","Vuyvx"},
			   {"Vxzuv","Vxzvu"},
			   {"Vyuzv","Vyuvz"}};

const char *tbarname[6][2] = {{"TBARvyz","TBARvzy"},
			      {"TBARuxv","TBARuvx"},
			      {"TBARyzv","TBARyvz"},
			      {"TBARxuv","TBARxvu"},
			      {"TBARzyv","TBARzvy"},
			      {"TBARvxu","TBARvux"}};

/* These are only used for display purposes and by the test program,
   so we really don't need permutations for them. */
const char *uuname[4][3] = {{"UUzxyv0", "UUzxyv1", "UUzxyv2"},
			    {"UUuyxv0", "UUuyxv1", "UUuyxv2"},
			    {"UUxzuv0", "UUxzuv1", "UUxzuv2"},
			    {"UUyuzv0", "UUyuzv1", "UUyuzv2"}};
const char *vvname[4][3] = {{"VVzxyv0", "VVzxyv1", "VVzxyv2"},
			    {"VVuyxv0", "VVuyxv1", "VVuyxv2"},
			    {"VVxzuv0", "VVxzuv1", "VVxzuv2"},
			    {"VVyuzv0", "VVyuzv1", "VVyuzv2"}};
const char *ttname[6][3] = {{"TTvyz0", "TTvyz1", "TTvyz2"},
			    {"TTuxv0", "TTuxv1", "TTuxv2"},
			    {"TTyzv0", "TTyzv1", "TTyzv2"},
			    {"TTxuv0", "TTxuv1", "TTxuv2"},
			    {"TTzyv0", "TTzyv1", "TTzyv2"},
			    {"TTvxu0", "TTvxu1", "TTvxu2"}};
const char *ssname[2][3] = {{"SSvyz0", "SSvyz1", "SSvyz2"},
			    {"SSuxv0", "SSuxv1", "SSuxv2"}};

/* But here they are in case we want them some day: */
/* const char *uuname0[4][2] = {{"UUzxyv0","UUzxvy0"}, */
/* 			     {"UUuyxv0","UUuyvx0"}, */
/* 			     {"UUxzuv0","UUxzvu0"}, */
/* 			     {"UUyuzv0","UUyuvz0"}}; */

/* const char *uuname1[4][2] = {{"UUzxyv1","UUzxvy1"}, */
/* 			     {"UUuyxv1","UUuyvx1"}, */
/* 			     {"UUxzuv1","UUxzvu1"}, */
/* 			     {"UUyuzv1","UUyuvz1"}}; */

/* const char *uuname2[4][2] = {{"UUzxyv2","UUzxvy2"}, */
/* 			     {"UUuyxv2","UUuyvx2"}, */
/* 			     {"UUxzuv2","UUxzvu2"}, */
/* 			     {"UUyuzv2","UUyuvz2"}}; */

/* const char *vvname0[4][2] = {{"VVzxyv0","VVzxvy0"}, */
/* 			     {"VVuyxv0","VVuyvx0"}, */
/* 			     {"VVxzuv0","VVxzvu0"}, */
/* 			     {"VVyuzv0","VVyuvz0"}}; */

/* const char *vvname1[4][2] = {{"VVzxyv1","VVzxvy1"}, */
/* 			     {"VVuyxv1","VVuyvx1"}, */
/* 			     {"VVxzuv1","VVxzvu1"}, */
/* 			     {"VVyuzv1","VVyuvz1"}}; */

/* const char *vvname2[4][2] = {{"VVzxyv2","VVzxvy2"}, */
/* 			     {"VVuyxv2","VVuyvx2"}, */
/* 			     {"VVxzuv2","VVxzvu2"}, */
/* 			     {"VVyuzv2","VVyuvz2"}}; */

/* const char *ttname0[6][2] = {{"TTvyz0","TTvzy0"}, */
/* 			     {"TTuxv0","TTuvx0"}, */
/* 			     {"TTyzv0","TTyvz0"}, */
/* 			     {"TTxuv0","TTxvu0"}, */
/* 			     {"TTzyv0","TTzvy0"}, */
/* 			     {"TTvxu0","TTvux0"}}; */

/* const char *ttname1[6][2] = {{"TTvyz1","TTvzy1"}, */
/* 			     {"TTuxv1","TTuvx1"}, */
/* 			     {"TTyzv1","TTyvz1"}, */
/* 			     {"TTxuv1","TTxvu1"}, */
/* 			     {"TTzyv1","TTzvy1"}, */
/* 			     {"TTvxu1","TTvux1"}}; */

/* const char *ttname2[6][2] = {{"TTvyz2","TTvzy2"}, */
/* 			     {"TTuxv2","TTuvx2"}, */
/* 			     {"TTyzv2","TTyvz2"}, */
/* 			     {"TTxuv2","TTxvu2"}, */
/* 			     {"TTzyv2","TTzvy2"}, */
/* 			     {"TTvxu2","TTvux2"}}; */

/* const char *ssname0[2][6] = {{"SSvyz0","SSzvy0","SSyzv0","SSvzy0","SSyvz0","SSzyv0"}, */
/* 			     {"SSuxv0","SSvux0","SSxvu0","SSuvx0","SSxuv0","SSvxu0"}}; */

/* const char *ssname1[2][6] = {{"SSvyz1","SSzvy1","SSyzv1","SSvzy1","SSyvz1","SSzyv1"}, */
/* 			     {"SSuxv1","SSvux1","SSxvu1","SSuvx1","SSxuv1","SSvxu1"}}; */

/* const char *ssname2[2][6] = {{"SSvyz2","SSzvy2","SSyzv2","SSvzy2","SSyvz2","SSzyv2"}, */
/* 			     {"SSuxv2","SSvux2","SSxvu2","SSuvx2","SSxuv2","SSvxu2"}}; */
