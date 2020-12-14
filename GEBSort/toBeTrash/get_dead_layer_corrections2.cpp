#include <iostream>
#include <fstream>
#include <cmath>
#include <array>


//load the range data; ranges is global, be careful!
std::array<std::array<double, 2>, 100> ranges;
void load_ranges();

//single escape:
//returns:
//[0]the initial alpha energy,
//[1]energy loss to dssd dead layer, 
//[2]energy loss to box dead layer, 
//[3]angle between dssd normal and alpha trajectory in rad
std::array<double, 4> get_dead_layer_corrections(double dssd_E, int dssd_strip_x, int dssd_strip_y, double box_E, int box_wall, int box_strip, int box_detector);

//single escape, one non escape:
//Estimates the energy loss of the escaped alpha particle to the active layer of the dssd based on the implantation depth. 
//returns: 
//[0]energy of the alpha which did not escape, 
//[1]energy loss of the escaped alpha to active layer of the dssd, 
//[2]energy loss of the escaped alpha to dead layer of the dssd,
//[3]energy loss of the escaped alpha to dead layer of the box.
std::array<double, 4> single_escape_one_nonescape(double dssd_E, int dssd_strip_x, int dssd_strip_y, double box1_E, int box1_wall, int box1_strip, int box1_detector);


//double escape:
//Returns:
//[0]energy loss to the active layer of the dssd for alpha1
//[1]energy loss to the active layer of the dssd for alpha2
//[2]energy loss to the dead layer of the dssd for alpha1
//[3]energy loss to the dead layer of the box for alpha1
//[4]energy loss to the dead layer of the dssd for alpha2
//[5]energy loss to the dead layer of the box for alpha2
std::array<double, 6> double_escape(double dssd_E, int dssd_strip_x, int dssd_strip_y, double box1_E, int box1_wall, int box1_strip, int box1_detector, double box2_E, int box2_wall, int box2_strip, int box2_detector);

/*
int main() {

	load_ranges();

	//measured energies
	double dssd_E, box1_E, box2_E;
	//coordinates
	int dssd_strip_x, dssd_strip_y, box1_wall, box1_strip, box1_detector, box2_wall, box2_strip, box2_detector;

	//dssd test event
	dssd_E = 3;
	dssd_strip_x = 80;
	dssd_strip_y = 80;

	//box1 test event
	box1_E = 2;
	box1_wall = 2;
	box1_strip = 4;
	box1_detector = 1;

	//box2 test event
	box2_E = 2;
	box2_wall = 2;
	box2_strip = 4;
	box2_detector = 2;


	std::array<double, 4> dead_layer_corrections = get_dead_layer_corrections(dssd_E, dssd_strip_x, dssd_strip_y, box1_E, box1_wall, box1_strip, box1_detector);
	//std::cout << "to_dssd_dead_layer: " << dead_layer_corrections[1] << std::endl;
	//std::cout << "to_box_dead_layer: " << dead_layer_corrections[2] << std::endl;
	//std::cout << "alpha_energy_total : " << dead_layer_corrections[0] << std::endl;
	//std::cout << "dssd_theta_angle : " << dead_layer_corrections[3] << std::endl;
	//std::cin.get();

	
	std::array<double, 4> single_escape_corrections = single_escape_one_nonescape(dssd_E, dssd_strip_x, dssd_strip_y, box1_E, box1_wall, box1_strip, box1_detector);
	//std::cout << "to_dssd_alpha1_nonescape : " << single_escape_corrections[0] << std::endl;
	//std::cout << "to_dssd_alpha2_escape : " << single_escape_corrections[1] << std::endl;
	//std::cout << "to_box_dead_layer_alpha1_escape : " << single_escape_corrections[2] << std::endl;
	//std::cout << "to_box_dead_layer_alpha1_escape : " << single_escape_corrections[3] << std::endl;
	//std::cin.get();

	std::array<double, 6> double_escape_corrections = double_escape(dssd_E, dssd_strip_x, dssd_strip_y, box1_E, box1_wall, box1_strip, box1_detector, box2_E, box2_wall, box2_strip, box2_detector);
	//std::cout << "to_dssd_alpha1 : " << double_escape_corrections[0] << std::endl;
	//std::cout << "to_dssd_alpha2 : " << double_escape_corrections[1] << std::endl;
	//std::cout << "to_dssd_dead_layer_alpha1 : " << double_escape_corrections[2] << std::endl;
	//std::cout << "to_box_dead_layer_alpha1 : " << double_escape_corrections[3] << std::endl;
	//std::cout << "to_dssd_dead_layer_alpha2 : " << double_escape_corrections[4] << std::endl;
	//std::cout << "to_box_dead_layer_alpha2 : " << double_escape_corrections[5] << std::endl;
	//std::cin.get();

	return 0;
}
*/

//load the range data of alphas in silicon from a file alpha_ranges.txt. 1st column is the energy of the alpha particle in MeV and 2nd is the range in um.
void load_ranges() {

	std::ifstream file;
	file.open("alpha_ranges.txt");

	for (int x = 0; x < 99; x++) {
		for (int y = 0; y < 1; y++) {
			file >> ranges[x][0];
			file >> ranges[x][1];
		}
	}

	file.close();

}

//single escape, one non escape:
//Estimates the energy loss of the escaped alpha particle to the active layer of the dssd based on the implantation depth. 
//returns: energy of the alpha which did not escape, energy loss of the escaped alpha to active layer of the dssd, energy loss of the escaped alpha to dead layer of the dssd, energy loss of the escaped alpha to dead layer of the box.
std::array<double, 4> single_escape_one_nonescape(double dssd_E, int dssd_strip_x, int dssd_strip_y, double box1_E, int box1_wall, int box1_strip, int box1_detector) {
  

	//process the escaped alpha:
	std::array<double, 4> corrections_box1 = get_dead_layer_corrections(dssd_E, dssd_strip_x, dssd_strip_y, box1_E, box1_wall, box1_strip, box1_detector);

	//implantation depth in um. This does not include dead layer. WARNING: your results will rely on this value.
	double implantation_depth = 5.0;

	//energy when alpha enters the dssd dead layer (i.e. energy not observed in dssd);
	double energy_rest_alpha_escape = box1_E + corrections_box1[2] + corrections_box1[1];

	//path length of alpha in active dssd layer. corrections_box1[3] is the angle between dssd normal and alpha trajectory in radians:
	double alpha_escape_path_in_active_dssd = implantation_depth / std::cos(corrections_box1[3]);

	//1) find the range range_energy_rest_alpha1 of alpha particle with energy energy_rest_alpha1:
	double lowerenergy = 0;
	double higherenergy = 0;
	double lowerrange = 0;
	double higherrange = 0;

	int j = 1;
	while (energy_rest_alpha_escape > ranges[j][0]) {
		j++;
	}

	lowerenergy = ranges[j - 1][0];
	higherenergy = ranges[j][0];
	lowerrange = ranges[j - 1][1];
	higherrange = ranges[j][1];

	//linear approximation:
	double slope = (higherrange - lowerrange) / (higherenergy - lowerenergy);
	double range_energy_rest_alpha1 = lowerrange + slope*(energy_rest_alpha_escape - lowerenergy);

	//2) find the "initial" energy of the alpha1 particle:
	double range_total_alpha1 = range_energy_rest_alpha1 + alpha_escape_path_in_active_dssd;

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

	//energy to active dssd (the escaped one):
	double to_dssd_alpha_escape = energy_total_alpha1 - energy_rest_alpha_escape;

	//energy to active dssd (the nonescaped one):
	double to_dssd_alpha_nonescape = dssd_E - to_dssd_alpha_escape;

	//energy of the alpha which did not escape,
	//energy loss of the escaped alpha to active layer of the dssd, 
	//energy loss of the escaped alpha to dead layer of the dssd, 
	//energy loss of the escaped alpha to dead layer of the box.
	std::array<double, 4> results = { to_dssd_alpha_nonescape, to_dssd_alpha_escape, corrections_box1[1], corrections_box1[2] };

	return results;
	
	
}





//finds the implantation depth in the case of double alpha escape, which both are detected in box. Based on the implantation depth this function calculates the original alpha
//energy for both escapes.
std::array<double, 6> double_escape(double dssd_E, int dssd_strip_x, int dssd_strip_y, double box1_E, int box1_wall, int box1_strip, int box1_detector, double box2_E, int box2_wall, int box2_strip, int box2_detector) {

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


	//std::cout << "implantation_depth : " << implantation_depth << std::endl;
	//std::cout << "dssd_E : " << dssd_E << std::endl;
	//std::cout << "DELTA_E : " << dssd_E - to_dssd_alpha1 - to_dssd_alpha2 << std::endl;

	//energy loss to the active layer of the dssd for alpha1
	//energy loss to the active layer of the dssd for alpha2
	//energy loss to the dead layer of the dssd for alpha1
	//energy loss to the dead layer of the box for alpha1
	//energy loss to the dead layer of the dssd for alpha2
	//energy loss to the dead layer of the box for alpha2
	std::array<double, 6> results = {to_dssd_alpha1, to_dssd_alpha2, corrections_box1[1], corrections_box1[2], corrections_box2[1], corrections_box2[2]};
	return  results;
}


//finds the actual location of a dssd event and a box event based on the strip/detector/wall numbers. After this estimates the energy loss of an alpha particle in dssd dead layer and box dead layer.
//returns the energy losses in 1x4 array, dead_layer_corrections[0] is the initial alpha energy. dead_layer_corrections[1] is the energy loss to dssd dead layer and corrections[2] in box dead layer.
//corrections[3] is the angle between alpha trajectory and dssd normal in radians.
std::array<double, 4> get_dead_layer_corrections(double dssd_E, int dssd_strip_x, int dssd_strip_y, double box_E, int box_wall, int box_strip, int box_detector) {

  //std::cout << "1 " << dssd_E << std::endl;
  //std::cout << "2 " << dssd_strip_x << std::endl;
  //std::cout << "3 " << dssd_strip_y << std::endl;
  //std::cout << "4 " << box_E << std::endl;
  //std::cout << "5 " << box_wall << std::endl;
  //std::cout << "6 " << box_strip << std::endl;
  //std::cout << "7 " << box_detector << std::endl;

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
