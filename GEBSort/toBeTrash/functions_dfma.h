
//**************************\\
// TRACE ANALYSIS FUNCTIONS \\
//**************************\\


#include <iostream>
#include <fstream>
#include <cmath>
#include <array>







// smooth3ch filter
void smoothFilter3ch(short int trace[], int trLen, short int smoothTrace2[]) {

 int i;

 for (i=5;i<trLen;i++) {
	smoothTrace2[i]= (trace[i-1] + 2*trace[i] + trace[i+1])/4;
 }
}


// smooth5ch filter
void smoothFilter5ch(short int trace[], int trLen, short int smoothTrace2[]) {

 int i;

 for (i=5;i<trLen;i++) {
	smoothTrace2[i]= (trace[i-2] + 2*trace[i-1] + 4*trace[i] + 2*trace[i+1] + trace[i+2])/10;
 }
}

// smooth7ch filter
void smoothFilter7ch(short int trace[], int trLen, short int smoothTrace2[]) {

 int i;

 for (i=5;i<trLen;i++) {
	smoothTrace2[i]= (trace[i-3] + 2*trace[i-2] + 4*trace[i-1] + 8*trace[i] + 4*trace[i+1] +2*trace[i+2] + trace[i+3])/22;
 }
}



// D filter 
void dFilter(short int trace[], int trLen, int dTrace[]) {

 int i;

 for (i=1; i<trLen; i++) {
   dTrace[i]=trace[i]-trace[i-1];
 }
}

// D filter 
void dFilterNEW(short int trace[], int trLen, short int dTrace[]) {

 int i;

 for (i=1; i<trLen; i++) {
   dTrace[i]=trace[i]-trace[i-5];
 }
}


void averaging(short int trace[], int trLen, short int dTrace[]){

  int i;

  for(i=1;i<trLen;i++){
    dTrace[i] = (trace[i-2] + trace[i-1] + trace[i] + trace[i+1] + trace[i+2])/5;
  }



}


// DD filter 
void ddFilter(int trace[], int trLen, int ddTrace[]){

 int i;

 for (i=10; i<trLen-4; i++) {
   ddTrace[i]=trace[i]-2*trace[i-1]+trace[i-2];
 }

}

// DD filter small signals
void ddFilterNEW2(int trace[], int trLen, int k, int d, short int ddTrace[]){

 int i;

 for (i=k+d+1; i<trLen-4; i++) {
   ddTrace[i]=trace[i]-trace[i-d]-(trace[i-k]-trace[i-k-d]);
 }  

}


// DD filter small signals
void ddFilterNEW(int trace[], int trLen, short int ddTrace[]){

 int i;

 for (i=10; i<trLen-4; i++) {
   ddTrace[i]=trace[i]-2*trace[i-5]+trace[i-10];
 }  

}


// LE filter  
void leFilter(short int trace[], int trLen, int k, int leTrace[], int& t0, int& t1, int leThr) {

 int i;
 t0 = 0;
 t1 = 0;

  

 for (i=k;i<trLen;i++) {  
   leTrace[i]=trace[i]-trace[i-k];
   //if(i<100)printf("trace[%i]: %i, trace[%i]: %i, leTrace[%i]: %i\n",i,trace[i],i-k,trace[i-k],leTrace[i]);
   
   if( (leTrace[i]>leThr) && (leTrace[i-1]<leThr) && (t0>0) && (t1==0) ){
     t1 = i;
   }
   if( (leTrace[i]>leThr) && (leTrace[i-1]< leThr) && (t0==0)){
     t0 = i;
   }

 }

//printf("\n#2a trace[10]: %i\n",trace[10]);
//printf("trLen: %i, k: %i\n",trLen,k);
//printf("\n#2b letrace[50]: %i\n",leTrace[50]);

}



// trigger 1 on trace
void trigTrace(int trace[], int trLen, int trStart, int trStop, int thr, int &t0) {

  int i;

      for (i=trStart;i<trStop;i++) {  
         if ((trace[i]>=thr)&&(trace[i+1]<thr)&&(t0==0)) {
	    t0=i;
            break;
         }
      }
        
}


void negd2(int trace[], int trLen, int trStart, int trStop, int thr, int &numtrigs) {

  int i;
   
  for(i=trStart;i<trStop;i++){
    if((trace[i]>=thr) && (trace[i+1] < thr)){
      numtrigs++;
    }
  } 


}
 

// Trapezoidal filter
void trapFilter(short int trace[], int trLen, int m, int k, int trapTrace[]){

   int i;
   for(i=m+k+m;i<trLen;i++){
      trapTrace[i] = trapTrace[i-1] + (trace[i]-trace[i-m])-(trace[i-m-k]-trace[i-m-k-m]);
   }

}


// extract energy from trace
void getEn2(int t0, short int trace[], int trlen, float &e0, int m, int k, int n, float trise) {

  int i;
  float ebl;
  float s1,s2;
  float corr;

  
  ebl=0;
  for (i=t0-m; i<t0; i++) {
    ebl=ebl+trace[i];
  }
  for (i=t0+k; i<t0+m+k; i++) {
    e0=e0+trace[i];
  }
    e0=e0-ebl;
 
  // PU correction

  /*
  corr=0;
  if ((n>0)&&(t0>2*n)) {
  s1=0; s2=0;
  for (i=1; i<n+1; i++) {
    s1=s1+trace[i];
  }
  for (i=t0-n; i<t0; i++) {
    s2=s2+trace[i];
  }
  corr=(s1-s2)/float(n)*(float(m)+float(k))/(float(t0)-float(n))*float(m);
  }
  e0=e0+corr;*/


}




// extract energies from PU traces
void GetPUEn(int t1, int t0, int trace[], int trlen, float &e0, float &e1, int m, int k, int n, float trise) {

   //int tr[1000];
   float ebl;
   float s1, s2;
   float corr;
   float dt;
   int iw0, iw1;
   int i;

   // globals: m, kk, trise

   //tr=trace;

   /*
   for (i=0; i<10; i++) {
     fprintf(ftest,"%d  %d \n", i, trace[i]);
   }
   */

   if (t1-t0>k) { // time separation longer than k 



        // integration windows
        if (t1-t0-k>m) iw0=m; else iw0=t1-t0-k; // 11-11 correction
        //if (iw0<0) iw0=1;
        if (trlen-t1-k>m) iw1=m; else iw1=trlen-t1-k; //11-11 correction
        //if (iw1<0) iw1=1;
        //if (iw1>iw0) iw1=iw0;




        // front energies
        e0=0; e1=0;
        for (i=t0+k;i<t0+k+iw0;i++) {
	  e0=e0+trace[i]-trace[i-iw0-k];
        }
        for (i=t1+k;i<t1+iw1+k;i++) { 
	  e1=e1+trace[i]-trace[i-k-t1+t0-iw1];
        }



        // PZ correction
        /*if (iw0>2*n) { // gate wide enough to get slope   
          s1=0; s2=0;
          for (i=t0+k;i<t0+k+n;i++) {
	    s1=s1+trace[i];
          }
          for (i=t0+k+iw0-n;i<t0+k+iw0;i++) {
	    s2=s2+trace[i];
        }
	  corr=float(s1-s2)/float(n)*float(iw0+k)/float(iw0-n)*float(iw1); // 11-11 correction
        //corr=0;
        e1=e1+corr;
        }*/


   
        e0=e0*float(m)/float(iw0); // gain correction // 11-11 correction
        e1=e1*float(m)/float(iw1); //  gain correction // 11-11 corrrection
        e1=e1-e0;



   } else { // t1-t0 < k



     if (t1-t0<=trise) {
   
 
        e0=0; e1=0;
        dt=t1-t0;
        // integration of the begining of 1st signal
        for (i=t0; i<t1;i++) {
	  e0=e0+trace[i];
        }
        //printf("integral %f\n", e0);
        ebl=0;    
        // base line inegration
        for (i=t0-m; i<t0;i++) {
	  ebl=ebl+trace[i];
        }
        e0=e0-ebl*dt/float(m); // energy of 1st signal
        // integration of 1st and 2nd signal
        for (i=t1+k;i<t1+m+k;i++) { 
	  e1=e1+trace[i];
        }           

        e1=e1-ebl; // energy sum of 1st and 2nd signal
        e0=2.0*trise*float(m)*e0/dt/dt; // energy of 1st signal
        e1=e1-e0; // enery of 2nd signal

     } else {

        e0=0; e1=0;
        dt=t1-t0-trise;
        // integration of the begining of 1st signal
        for (i=t0+trise; i<t1;i++) {
	  e0=e0+trace[i];
        }
        //printf("integral %f\n", e0);
        ebl=0;    
        // base line inegration
        for (i=t0-m; i<t0;i++) {
	  ebl=ebl+trace[i];
        }
        e0=e0-ebl*dt/float(m); // energy of 1st signal
        // integration of 1st and 2nd signal
        for (i=t1+k;i<t1+m+k;i++) { 
	  e1=e1+trace[i];
        }           

        e1=e1-ebl; // energy sum of 1st and 2nd signal
        e0=float(m)*e0/dt; // energy of 1st signal
        e1=e1-e0; // enery of 2nd signal

     }

   }

  }


/*

//load the range data; ranges is global, be careful!
std::array<std::array<double, 2>, 100> ranges;
void load_ranges();

//single escape:
//returns the initial alpha energy, energy loss to dssd dead layer, energy loss to box dead layer and the angle between dssd normal and alpha trajectory in rad, respectively.
std::array<double, 4> get_dead_layer_corrections(double dssd_E, int dssd_strip_x, int dssd_strip_y, double box_E, int box_wall, int box_strip, int box_detector);

//double escape:
std::array<double, 2> double_escape(double dssd_E, int dssd_strip_x, int dssd_strip_y, double box1_E, int box1_wall, int box1_strip, int box1_detector, double box2_E, int box2_wall, int box2_strip, int box2_detector);


//finds the implantation depth in the case of double alpha escape, which both are detected in box. Based on the implantation depth this function calculates the original alpha
//energy for both escapes.
std::array<double, 2> double_escape(double dssd_E, int dssd_strip_x, int dssd_strip_y, double box1_E, int box1_wall, int box1_strip, int box1_detector, double box2_E, int box2_wall, int box2_strip, int box2_detector) {

	//process the two escapes individually:
	std::array<double, 4> corrections_box1 = get_dead_layer_corrections(dssd_E, dssd_strip_x, dssd_strip_y, box1_E, box1_wall, box1_strip, box1_detector);
	std::array<double, 4> corrections_box2 = get_dead_layer_corrections(dssd_E, dssd_strip_x, dssd_strip_y, box2_E, box2_wall, box2_strip, box2_detector);

	//energies when alpha enters the dssd dead layer (i.e. energy not observed in dssd);
	double energy_rest_alpha1 = box1_E + corrections_box1[2] + corrections_box1[1];
	double energy_rest_alpha2 = box2_E + corrections_box2[2] + corrections_box2[1];

	//initial implantation depth in um: (excluding dead layer)
	double implantation_depth = 0.1;

	//initialize the energy deposit to active layer of dssd:
	double to_dssd_alpha1 = 0;
	double to_dssd_alpha2 = 0;


	//increase the implantation depth until the calculated (trajectory based) event energy in dssd agrees with the observed one.
	while (to_dssd_alpha1 + to_dssd_alpha2 < dssd_E) {

		//path length of alpha in active dssd layer. corrections_boxX[3] is the angle between dssd normal and alpha trajectory in radians:
		double alpha1_path_in_active_dssd = implantation_depth / std::cos(corrections_box1[3]);
		double alpha2_path_in_active_dssd = implantation_depth / std::cos(corrections_box2[3]);

		//1) find the range range_energy_rest_alpha1 of alpha particle with energy energy_rest_alpha1:
		double lowerenergy = 0;
		double higherenergy = 0;
		double lowerrange = 0;
		double higherrange = 0;

		int j = 1;
		while (energy_rest_alpha1 > ranges[j][0]) {
			j++;
		}

		lowerenergy = ranges[j - 1][0];
		higherenergy = ranges[j][0];
		lowerrange = ranges[j - 1][1];
		higherrange = ranges[j][1];

		//linear approximation:
		double slope = (higherrange - lowerrange) / (higherenergy - lowerenergy);
		double range_energy_rest_alpha1 = lowerrange + slope*(energy_rest_alpha1 - lowerenergy);

		//2) find the "initial" energy of the alpha1 particle:
		double range_total_alpha1 = range_energy_rest_alpha1 + alpha1_path_in_active_dssd;

		lowerenergy = 0;
		higherenergy = 0;
		lowerrange = 0;
		higherrange = 0;

		j = 1;
		while (range_total_alpha1 > ranges[j][1]) {
			j++;
		}
		lowerenergy = ranges[j - 1][0];
		higherenergy = ranges[j][0];
		lowerrange = ranges[j - 1][1];
		higherrange = ranges[j][1];

		//linear approximation:
		slope = (higherenergy - lowerenergy) / (higherrange - lowerrange);
		double energy_total_alpha1 = lowerenergy + slope*(range_total_alpha1 - lowerrange);

		//energy to active dssd (alpha1):
		to_dssd_alpha1 = energy_total_alpha1 - energy_rest_alpha1;

		// repeat for alpha2:
		//3) find the range range_energy_rest_alpha2 of alpha particle with energy energy_rest_alpha2:
		lowerenergy = 0;
		higherenergy = 0;
		lowerrange = 0;
		higherrange = 0;

		j = 1;
		while (energy_rest_alpha2 > ranges[j][0]) {
			j++;
		}

		lowerenergy = ranges[j - 1][0];
		higherenergy = ranges[j][0];
		lowerrange = ranges[j - 1][1];
		higherrange = ranges[j][1];

		//linear approximation:
		slope = (higherrange - lowerrange) / (higherenergy - lowerenergy);
		double range_energy_rest_alpha2 = lowerrange + slope*(energy_rest_alpha2 - lowerenergy);

		//2) find the "initial" energy of the alpha particle:
		double range_total_alpha2 = range_energy_rest_alpha2 + alpha2_path_in_active_dssd;

		lowerenergy = 0;
		higherenergy = 0;
		lowerrange = 0;
		higherrange = 0;

		j = 1;
		while (range_total_alpha2 > ranges[j][1]) {
			j++;
		}
		lowerenergy = ranges[j - 1][0];
		higherenergy = ranges[j][0];
		lowerrange = ranges[j - 1][1];
		higherrange = ranges[j][1];

		//linear approximation:
		slope = (higherenergy - lowerenergy) / (higherrange - lowerrange);
		double energy_total_alpha2 = lowerenergy + slope*(range_total_alpha2 - lowerrange);

		//energy to active dssd (alpha1):
		to_dssd_alpha2 = energy_total_alpha2 - energy_rest_alpha2;

		//increase implantation depth for the next round:
		implantation_depth = implantation_depth + 0.001;

	}


	std::cout << "to_dssd_alpha1 : " << to_dssd_alpha1 << std::endl;
	std::cout << "to_dssd_alpha2 : " << to_dssd_alpha2 << std::endl;
	std::cout << "implantation_depth : " << implantation_depth << std::endl;
	std::cout << "dssd_E : " << dssd_E << std::endl;
	std::cout << "DELTA_E : " << dssd_E - to_dssd_alpha1 - to_dssd_alpha2 << std::endl;
	std::cin.get();


	std::array<double, 2> results = { 0,0 };
	return  results;
}


//finds the actual location of a dssd event and a box event based on the strip/detector/wall numbers. After this estimates the energy loss of an alpha particle in dssd dead layer and box dead layer.
//returns the energy losses in 1x4 array, dead_layer_corrections[0] is the initial alpha energy. dead_layer_corrections[1] is the energy loss to dssd dead layer and corrections[2] in box dead layer.
//corrections[3] is the angle between alpha trajectory and dssd normal in radians.
std::array<double, 4> get_dead_layer_corrections(double dssd_E, int dssd_strip_x, int dssd_strip_y, double box_E, int box_wall, int box_strip, int box_detector) {

	//the strip widths in mm
	double strip_width_dssd = 0.5;
	double strip_width_box = 5.0;

	//the detector geometry in  mm
	double dssd_width = 80.0;
	double box_lenght = 40.0;

	//thicness of dead layers in um
	double dssd_dead_layer = 0.8;
	double box_dead_layer = 0.8;

	//position of dssd event (x,y,z)
	std::array<double, 3> dssd_event_place = { 0,0,0 };
	dssd_event_place[0] = (dssd_strip_x - 79 - 0.5)*strip_width_dssd; //x-coordinate
	dssd_event_place[1] = (dssd_strip_y - 79 - 0.5)*strip_width_dssd; //y-coordinate


	//position of box event (x,y,z)
	std::array<double, 3> box_event_place = { 0,0,0 };

	//wall#1, detector#1
	if (box_wall == 1 && box_detector == 1) {
		box_event_place[0] = (box_strip - 4)*strip_width_box;
		box_event_place[1] = 0.5*dssd_width;
		box_event_place[2] = 0.5*box_lenght;
	}
	//wall#1, detector#2
	else if (box_wall == 1 && box_detector == 2) {
		box_event_place[0] = (box_strip - 4)*strip_width_box;
		box_event_place[1] = 0.5*dssd_width;
		box_event_place[2] = 1.5*box_lenght;
	}
	//wall#2, detector#1
	else if (box_wall == 2 && box_detector == 1) {
		box_event_place[0] = 0.5*dssd_width;
		box_event_place[1] = -(box_strip - 4)*strip_width_box;
		box_event_place[2] = 0.5*box_lenght;
	}
	//wall#2, detector#2
	else if (box_wall == 2 && box_detector == 2) {
		box_event_place[0] = 0.5*dssd_width;
		box_event_place[1] = -(box_strip - 4)*strip_width_box;
		box_event_place[2] = 1.5*box_lenght;
	}
	//wall#3, detector#1
	else if (box_wall == 3 && box_detector == 1) {
		box_event_place[0] = -(box_strip - 4)*strip_width_box;
		box_event_place[1] = -0.5*dssd_width;
		box_event_place[2] = 0.5*box_lenght;
	}
	//wall#3, detector#2
	else if (box_wall == 3 && box_detector == 2) {
		box_event_place[0] = -(box_strip - 4)*strip_width_box;
		box_event_place[1] = -0.5*dssd_width;
		box_event_place[2] = 1.5*box_lenght;
	}
	//wall#4, detector#1
	else if (box_wall == 4 && box_detector == 1) {
		box_event_place[0] = -0.5*dssd_width;
		box_event_place[1] = (box_strip - 4)*strip_width_box;
		box_event_place[2] = 0.5*box_lenght;
	}
	//wall#4, detector#2
	else if (box_wall == 4 && box_detector == 2) {
		box_event_place[0] = -0.5*dssd_width;
		box_event_place[1] = (box_strip - 4)*strip_width_box;
		box_event_place[2] = 1.5*box_lenght;
	}


	//std::cout << "dssd_x: " << dssd_event_place[0] << std::endl;
	//std::cout << "dssd_y: " << dssd_event_place[1] << std::endl;
	//std::cout << "dssd_z: " << dssd_event_place[2] << std::endl;
	//std::cout << "box_x: " << box_event_place[0] << std::endl;
	//std::cout << "box_y: " << box_event_place[1] << std::endl;
	//std::cout << "box_z: " << box_event_place[2] << std::endl;
	//std::cin.get();

	//new cordinates: dssd_event_place=(0,0,0), then
	std::array<double, 3> alpha_vector = { 0,0,0 };
	alpha_vector[0] = box_event_place[0] - dssd_event_place[0];
	alpha_vector[1] = box_event_place[1] - dssd_event_place[1];
	alpha_vector[2] = box_event_place[2] - dssd_event_place[2];

	//length of the alpha_vector:
	double alpha_vector_length = std::sqrt(std::pow(alpha_vector[0], 2) + std::pow(alpha_vector[1], 2) + std::pow(alpha_vector[2], 2));

	//path lenght of the alpha particle in the dead layer of the dssd, dssd normal vector (0,0,1):
	double dssd_theta_angle = std::acos(alpha_vector[2] / alpha_vector_length);
	double path_alpha_dssd = dssd_dead_layer / (alpha_vector[2] / alpha_vector_length);

	//path lenght of alpha particle in the dead layer of the box. Box normal vector is (0,-1,0), (-1,0,0), (0,1,0), (1,0,0) for walls 1,2,3,4, respectively.
	//normal is always pointing inside the box. Note that the alpha_vector is pointing inside the box wall.:
	double path_alpha_box = 0;
	if (box_wall == 1) {
		path_alpha_box = box_dead_layer / (alpha_vector[1] / alpha_vector_length);
	}
	else if (box_wall == 2) {
		path_alpha_box = box_dead_layer / (alpha_vector[0] / alpha_vector_length);
	}
	else if (box_wall == 3) {
		path_alpha_box = box_dead_layer / (-alpha_vector[1] / alpha_vector_length);
	}
	else if (box_wall == 4) {
		path_alpha_box = box_dead_layer / (-alpha_vector[0] / alpha_vector_length);
	}

	//std::cout << "path_alpha_dssd: " << path_alpha_dssd << std::endl;
	//std::cout << "path_alpha_box: " << path_alpha_box << std::endl;
	//std::cout << "dssd_theta_angle: " << dssd_theta_angle << std::endl;
	//std::cin.get();



	//energy loss to the dead layer of the box:
	//1) find the range range_box_E of alpha particle with energy E_box:
	int j = 1;
	while (box_E > ranges[j][0]) {
		j++;
	}
	double lowerenergy = ranges[j - 1][0];
	double higherenergy = ranges[j][0];
	double lowerrange = ranges[j - 1][1];
	double higherrange = ranges[j][1];

	//linear approximation:
	double slope = (higherrange - lowerrange) / (higherenergy - lowerenergy);
	double range_box_E = lowerrange + slope*(box_E - lowerenergy);
	//std::cout << "range_box_E: " << range_box_E << std::endl;

	//2) find the energy of the alpha particle when entering box:
	double range_box_total = range_box_E + path_alpha_box;

	lowerenergy = 0;
	higherenergy = 0;
	lowerrange = 0;
	higherrange = 0;

	j = 1;
	while (range_box_total > ranges[j][1]) {
		j++;
	}
	lowerenergy = ranges[j - 1][0];
	higherenergy = ranges[j][0];
	lowerrange = ranges[j - 1][1];
	higherrange = ranges[j][1];

	//linear approximation:
	slope = (higherenergy - lowerenergy) / (higherrange - lowerrange);
	double energy_box_total = lowerenergy + slope*(range_box_total - lowerrange);

	//3) find the energy of the alpha particle when entering dssd dead layer:
	double range_total = range_box_E + path_alpha_box + path_alpha_dssd;

	lowerenergy = 0;
	higherenergy = 0;
	lowerrange = 0;
	higherrange = 0;

	j = 1;
	while (range_total > ranges[j][1]) {
		j++;
	}
	lowerenergy = ranges[j - 1][0];
	higherenergy = ranges[j][0];
	lowerrange = ranges[j - 1][1];
	higherrange = ranges[j][1];

	//linear approximation:
	slope = (higherenergy - lowerenergy) / (higherrange - lowerrange);
	double energy_total_not_observed_in_dssd = lowerenergy + slope*(range_total - lowerrange);


	//4) prepare the results
	double to_dssd_dead_layer = energy_total_not_observed_in_dssd - energy_box_total;
	double to_box_dead_layer = energy_box_total - box_E;
	double alpha_energy_total = energy_total_not_observed_in_dssd + dssd_E;



	//std::cout << "range_box_total: " << range_box_total << std::endl;
	//std::cout << "energy_box_total: " << energy_box_total << std::endl;
	//std::cout << "range_total: " << range_total << std::endl;
	//std::cout << "energy_total_not_observed_in_dssd: " << energy_total_not_observed_in_dssd << std::endl;
	//std::cin.get();

	//initialize the results array
	std::array<double, 4>  results = { alpha_energy_total,to_dssd_dead_layer,to_box_dead_layer, dssd_theta_angle };

	return results;
}

*/
