
typedef struct GEBMERGEPAR
{
  int addT0;
  int vetocubes;
  char vetocubefn[512];
  int nevents;
  int xyz_smear;
  float smear_x;
  float smear_y;
  float smear_z;

  int GT2AGG4;
  char GT2AGG4_fn[256];

  int waitusec;

} GEBMERGEPAR;
