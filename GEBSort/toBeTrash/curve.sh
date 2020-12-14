

findPT keep/agata-code-data-check-nosing 10 5 50 12846181  | awk '{print $5, $6}' > x.2
echo "0.0623 0.6437" > x.2
echo ".061512846181 .5667919245" > x.2
echo "0.0627 0.6474" > x.2

# remember that simulated data has twice the number of gamma rays

findPT fom01.spe 10 5 50 12846181 | awk '{print $5, $6}' > x.3
findPT fom02.spe 10 5 50 12846181 | awk '{print $5, $6}' >> x.3
findPT fom03.spe 10 5 50 12846181 | awk '{print $5, $6}' >> x.3
findPT fom04.spe 10 5 50 12846181 | awk '{print $5, $6}' >> x.3
findPT fom05.spe 10 5 50 12846181 | awk '{print $5, $6}' >> x.3
findPT fom06.spe 10 5 50 12846181 | awk '{print $5, $6}' >> x.3
findPT fom07.spe 10 5 50 12846181 | awk '{print $5, $6}' >> x.3
findPT fom08.spe 10 5 50 12846181 | awk '{print $5, $6}' >> x.3
findPT fom09.spe 10 5 50 12846181 | awk '{print $5, $6}' >> x.3
findPT fom10.spe 10 5 50 12846181 | awk '{print $5, $6}' >> x.3
findPT fom11.spe 10 5 50 12846181 | awk '{print $5, $6}' >> x.3
findPT fom12.spe 10 5 50 12846181 | awk '{print $5, $6}' >> x.3
findPT fom13.spe 10 5 50 12846181 | awk '{print $5, $6}' >> x.3
findPT fom14.spe 10 5 50 12846181 | awk '{print $5, $6}' >> x.3
findPT fom15.spe 10 5 50 12846181 | awk '{print $5, $6}' >> x.3
findPT fom16.spe 10 5 50 12846181 | awk '{print $5, $6}' >> x.3
findPT fom17.spe 10 5 50 12846181 | awk '{print $5, $6}' >> x.3
findPT fom18.spe 10 5 50 12846181 | awk '{print $5, $6}' >> x.3
findPT fom19.spe 10 5 50 12846181 | awk '{print $5, $6}' >> x.3
findPT fom20.spe 10 5 50 12846181 | awk '{print $5, $6}' >> x.3

xmgrace  -noask -timestamp -free -nosafe -nosigcatch x.3 x.2 -par p.par

exit

findPT fom02 10 5 50 12846181
findPT fom04 10 5 50 12846181
findPT fom06 10 5 50 12846181
findPT fom08 10 5 50 12846181
findPT fom10 10 5 50 12846181
findPT fom12 10 5 50 12846181
findPT fom14 10 5 50 12846181
findPT fom16 10 5 50 12846181
findPT fom18 10 5 50 12846181
findPT fom20 10 5 50 12846181 

findPT keep/agata-code-data-check-nosing 10 5 50 12846181  | awk '{print $10, $11}'> x.4

findPT fom02 10 5 50 12846181 | awk '{print $10, $11}' > x.5
findPT fom04 10 5 50 12846181 | awk '{print $10, $11}' >> x.5
findPT fom06 10 5 50 12846181 | awk '{print $10, $11}' >> x.5
findPT fom08 10 5 50 12846181 | awk '{print $10, $11}' >> x.5
findPT fom10 10 5 50 12846181 | awk '{print $10, $11}' >> x.5
findPT fom12 10 5 50 12846181 | awk '{print $10, $11}' >> x.5
findPT fom14 10 5 50 12846181 | awk '{print $10, $11}' >> x.5
findPT fom16 10 5 50 12846181 | awk '{print $10, $11}' >> x.5
findPT fom18 10 5 50 12846181 | awk '{print $10, $11}' >> x.5
findPT fom20 10 5 50 12846181 | awk '{print $10, $11}' >> x.5

xmgrace x.5 x.4 -par p.par
exit
pjx("fomXe","x",0,0.2); wrspe("x","fom02.spe");
pjx("fomXe","x",0,0.4); wrspe("x","fom04.spe");
pjx("fomXe","x",0,0.6); wrspe("x","fom06.spe");
pjx("fomXe","x",0,0.8); wrspe("x","fom08.spe");
pjx("fomXe","x",0,1.0); wrspe("x","fom10.spe");
pjx("fomXe","x",0,1.2); wrspe("x","fom12.spe");
pjx("fomXe","x",0,1.4); wrspe("x","fom14.spe");
pjx("fomXe","x",0,1.6); wrspe("x","fom16.spe");
pjx("fomXe","x",0,1.8); wrspe("x","fom18.spe");
pjx("fomXe","x",0,2.0); wrspe("x","fom20.spe");

