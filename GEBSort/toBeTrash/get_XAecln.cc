/* This function loops through GS detectors and gets ehiraw and writes to spe.*/

{

char p = 0x22;
char spefil [50];
int i,n;

for (i=0;i<50;i++)
  {

  if(i<10) 
    n=sprintf(spefil,"ehi00%i.spe",i);
  else if (i<100)
    n=sprintf(spefil,"ehi0%i.spe",i);
  else
    n=sprintf(spefil,"ehi%i.spe",i);
  
  // spefil=sprintf("%c%i%c",fil1c,i,fil2);
  
  printf("%i %s\n",i,spefil);
  
  pjx("XAEhiCln","p",i,i);
  //d1("p",0,5000);
  wrspe("p",spefil);

  pjx("XAEhiCln","p",1,110);
  wrspe("p","ref.spe");

  }
}
