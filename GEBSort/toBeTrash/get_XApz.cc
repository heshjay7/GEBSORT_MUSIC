/* This function loops through GS detectors and gets ehiraw and writes to spe.*/

{

char p = 0x22;
char spefil [50];
int i,n;

for (i=1;i<50;i++)
  {

  if(i<10) 
    n=sprintf(spefil,"pz00%i.spe",i);
  else if (i<100)
    n=sprintf(spefil,"pz0%i.spe",i);
  else
    n=sprintf(spefil,"pz%i.spe",i);
  
  // spefil=sprintf("%c%i%c",fil1c,i,fil2);
  
  printf("%i %s\n",i,spefil);
  
  pjy("XApzraw","p",i,i);
  //d1("p",0,5000);
  wrspe("p",spefil);
  }
}
