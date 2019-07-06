#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>
#include <stdbool.h>
#include <inttypes.h>
#include <immintrin.h>
#include <GL/glew.h>
#include <GL/glut.h>

#define G 6.67e-11
#define RAND_64_MAX ~(0ULL)


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

void get_arg(double* val)
{
	char* t = strtok(NULL, ",");
	*val = atof(t);
}

void get_argf(double *val)
{
	char* t = strtok(NULL, ",");
	*val = atof(t);
}

void get_args(char **val)
{
	char* t = strtok(NULL, ",");
	*val = t;
}

void get_argi(int *val)
{
	char* t = strtok(NULL, ",");
	*val = atoi(t);
}

void get_argzu(size_t *val)
{
	char* t = strtok(NULL, ",");
	*val = atoi(t);
}

uint64_t rand_64()
{
    unsigned long long val;
    while(!_rdrand64_step(&val));
    return (uint64_t)val;
}

void readConfig(char* filename, size_t* input_bodies, char** input_filename, int* prec, char** output_filename, int* time_steps, double* max_time, size_t* ring_assigned, double* assigned_angle, size_t* ring_bodies, 
				double* ring_radius, double* ring_veloc, size_t* target)
{

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
	//
	// Reads data from input file.
	//
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

	FILE* f = fopen(filename, "r");

	if (f == NULL)
		{
			printf("Could not find config file! \"%s\"\n", filename);
			return;
		}

	char* buffer = NULL;

	size_t lines_read = 0;
	size_t size;

	getline(&buffer, &size, f);

	size_t i = 0;

	while (getline(&buffer, &size, f) != EOF)
	{
		char* t = strtok(buffer, ",");

		get_argzu(input_bodies);
		get_args(input_filename);

		get_argi(prec);
		get_args(output_filename);
		
		get_argi(time_steps);
		get_argf(max_time);

		get_argzu(ring_assigned);
		get_argf(assigned_angle);
		get_argzu(ring_bodies);
		get_argf(ring_radius);
		get_argf(ring_veloc);
		
		get_argzu(target);

		fflush(stdout);

		if (i >= 1)
		{
			exit(1);
		}

		i++;

	}

	/*

	if (input_filename != NULL)
	{
		
		char temp_filename[512];
   		sprintf(temp_filename, "%s.csv", input_filename);
   		*input_filename = temp_filename;

	}
	*/

	fclose(f);
}

void readInputArgs(int argc, char** argv, size_t* input_bodies, char** input_filename, int* prec, char** output_filename, int* time_steps, double* max_time, 
					size_t* ring_assigned, double* assigned_angle, size_t* ring_n_bodies, double* ring_radius, double* ring_veloc, size_t* target)
{
	int c;
	opterr = 0;

	while ((c = getopt(argc, argv, "n:i:p:o:t:m:a:q:e:r:v:")) != -1)
	{
		switch (c)
		{
			case 'n': *input_bodies = atoi(optarg); break;
			case 'i': *input_filename = optarg; break;

			case 'p': *prec = atoi(optarg); break;
			case 'o': *output_filename = optarg; break;

			case 't': *time_steps = atoi(optarg); break;
			case 'm': *max_time = atof(optarg); break;

			case 'a': *ring_assigned = atoi(optarg); break;
			case 'q': *assigned_angle = atof(optarg);break;
			case 'e': *ring_n_bodies = atoi(optarg); break;
			case 'r': *ring_radius = atof(optarg); break;	
			case 'v': *ring_veloc = atof(optarg); break;

			case 'x': *target = atoi(optarg); break;

		}
	}	

	/*

	if (input_filename != NULL)
	{
		
		char temp_filename[512];
   		sprintf(temp_filename, "%s.csv", input_filename);
   		*input_filename = temp_filename;

	}

	*/
}

void readFile(char* filename, double* mass, double* radius, double* posit, double* veloc, size_t num_bodies, size_t* bodies_read)
{

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
	//
	// Reads data from input file.
	//
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

	FILE* f = fopen(filename, "r");

	if (f == NULL)
		{
			printf("Could not find input file! \"%s\"\n", filename);
			return;
		}

	char* buffer = NULL;

	size_t lines_read = 0;
	size_t size;

	getline(&buffer, &size, f);

	size_t i = 0;

	while (getline(&buffer, &size, f) != EOF)
	{
		char* t = strtok(buffer, ",");

		get_arg(&mass[i]);
		get_arg(&radius[i]);

		get_arg(&posit[3*i]);
		get_arg(&posit[(3*i)+1]);
		get_arg(&posit[(3*i)+2]);

		get_arg(&veloc[3*i]);
		get_arg(&veloc[(3*i)+1]);
		get_arg(&veloc[(3*i)+2]);

		if (i >= num_bodies)
		{
			fprintf(stderr, "File bigger than num bodies!\n");
			exit(1);
		}

		i++;
	}

	*bodies_read = i;
	fclose(f);
}

double calcAbs(double* vector)
{
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
	//
	// Returns absoloute magnitude of 3d vector.
	//
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

	return sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
}

void Arange(size_t length, double* array)
{
	for (size_t i = 0; i < (length); ++i) { array[i] = i; }
}

void calcMax(double* array, size_t array_height, size_t array_width, size_t array_max_width, double* max)
{

	for (size_t i = 0; i < array_width; ++i)
	{
		double temp_max = array[0]; 
		for (size_t t = 0; t < array_height; ++t)
		{
			if (array[t*array_max_width + i] > temp_max)
			{
				temp_max = array[t*array_max_width + i];
			}

		}

		max[i] = temp_max;

	}
}

void calcMin(double* array, size_t array_height, size_t array_width, size_t array_max_width, double* max)
{

	for (size_t i = 0; i < array_width; ++i)
	{
		double temp_max = array[0]; 
		for (size_t t = 0; t < array_height; ++t)
		{
			if (array[t*array_max_width + i] < temp_max)
			{
				temp_max = array[t*array_max_width + i];
			}

		}

		max[i] = temp_max;

	}
}

void calcDifference(double* final, double* initial, double* vector)
{
	for (size_t d = 0; d < 3; ++d) { vector[d] = final[d] - initial[d]; }
}

void calcDisplcaement(double* posit, double* n_posit, size_t target, size_t num_n_bodies, size_t max_n_bodies, size_t max_bodies, size_t start_time_steps, size_t end_time_steps, double* displace)
{
	for (size_t t = start_time_steps; t < (end_time_steps + start_time_steps); ++t)
	{
		for (size_t i = 0; i < num_n_bodies ; ++i)
		{
			double difference[3];
			calcDifference(&n_posit[t*max_n_bodies*3 + (i*3)], &posit[t*max_bodies*3 + (target*3)], difference);
			displace[t*max_n_bodies + i] = calcAbs(&difference[0]);
		}

	}		
}

void calcGravAccel(double* affected_posit, double* effector_posit, double* effector_mass, int affected_no, int num_effector, double* accel)
{

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
	//
	// Returns (accel) gravitational acceleration at point (affected_posit) given number (num_effector) of graviational attractors at an array of positions (effector_postit) with masses (effector_mass).
	// If affected object is also and effector than its number in the position array (effector_posit), should be given as (affected_no), otherwise set affected_no to -1. 
	//
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

	// Assigning current position to variables:

	for (size_t j = 0; j < num_effector; ++j)	
	{
		if (affected_no == j) continue;  // <-- If comparing body to itself iteration skips to prevent devide by 0 errors.

		// Calculating diplacement and assigning to variables:

		double displace[3];
		calcDifference(&effector_posit[(j*3)], affected_posit, displace );

		// Calculating magnitude of gravitational acceleration:

		double grav_mag = effector_mass[j]/(calcAbs(displace)*calcAbs(displace)*calcAbs(displace));

		// Calculating gravitational acelleration vector:

		for (size_t d = 0; d < 3; ++d) { accel[d] += grav_mag*displace[d]; }

				
	}
}

void calcAllGravAccel(double* affected_posit, double* effector_posit, double* effector_mass, bool same, int num_affected, int num_effector, double* affected_accel)
{
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
	//
	// Returns (accel) gravitational acceleration at array of points (affected_posit) given number (num_effector) of graviational attractors at an array of positions (effector_postit) with masses (effector_mass).
	// If affected objects are also an effector than (same) should be set to 1, elsewize 0. 
	//
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

	for (size_t i = 0; i < num_affected; ++i) { calcGravAccel(&affected_posit[(3*i) + 0], effector_posit, effector_mass, i*same - (1 - same) , num_effector, &affected_accel[(3*i)]); } 
}
		
void calcGravVerlet(size_t current_time, double* posit, double* veloc, double* accel, double* effector_mass, double* effector_posit, double dt, int num_affected, int num_effector, int max_affected, int max_effector, bool same)
{

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
	//
	// Returns (posit), acceleration (veloc), velocity (acell) of (num_affected) objects at (posit) after one time interval duration dt, given graviitational attraction between objects and num_effector) objects with mass (mass_effector) 
	// at positions (pefffector_posit).  
	//
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

	for (size_t i = 0; i < num_affected; ++i)
	{

		// ~~ Saving old acceleration and resetting to 0:

		double prev_accel[3];
		for (size_t d = 0; d < 3; ++d) { prev_accel[d] = accel[(3*i) + d]; }

		accel[(3*i) + 0] = 0, accel[(3*i) + 1] = 0, accel[(3*i) + 2] = 0; 

		// ~~~~ Calculating new position:

		for (size_t d = 0; d < 3; ++d) { posit[(current_time*max_affected*3) + (i*3) + d] = posit[((current_time - 1)*max_affected*3) + (i*3) + d] + veloc[(3*i) + d]*dt + prev_accel[d]*dt*dt*0.5; }

		// ~~~~ Calculating new acceleration:

		calcGravAccel(&posit[((current_time*max_affected*3)) + (i*3)], &effector_posit[(current_time - 1)*max_effector*3], effector_mass, i*same - (1 - same), num_effector, &accel[(i*3)]); //

		// ~~~~ Calcualting new velocity:

		for (size_t d = 0; d < 3; ++d) { veloc[(3*i) + d] += (prev_accel[d] + accel[(3*i) + d])*0.5*dt; }

	}
}

void calcRotatedVector(double* init_vector, double* centre, double* rotation, double* fin_vector)
{
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
	//
	// Returns (fin_vector) the postion of a 3d vector after rotation (rotation) around a point (centre). 
	//
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

	//Assigning values to variables for ease of use:

	double x = init_vector[0], y = init_vector[1], z = init_vector[2];
	double a = centre[0],      b = centre[1],      c = centre[2];
	double u = rotation[0],    v = rotation[1],    w = rotation[2], angle = rotation[3];

	double cos_angle = cos((2.*M_PI/360.) * angle); double sin_angle = sin((2.*M_PI/360.) * angle);

	//Calculating rotation:

	fin_vector[0] = (a*(v*v + w*w) - u*(b*v + c*w - u*x - v*y - w*z))*(1-cos_angle) + x*cos_angle + (- c*v + b*w - w*y + v*z)*sin_angle;
	fin_vector[1] = (b*(u*u + w*w) - v*(a*u + c*w - u*x - v*y - w*z))*(1-cos_angle) + y*cos_angle + (  c*u - a*w + w*x - u*z)*sin_angle;
	fin_vector[2] = (c*(u*u + v*v) - w*(a*u + b*v - u*x - v*y - w*z))*(1-cos_angle) + z*cos_angle + (- b*u + a*v - v*x + u*y)*sin_angle;
}

void rotateBody(double* posit, double* veloc, double* accel, double* effector_posit, double* effector_mass, double* centre, double* rotation, size_t body_no, size_t max_affected, size_t max_effector, size_t num_effector, 
				size_t current_time, bool same)
{
	calcRotatedVector(&posit[(current_time*max_affected*3 + body_no*3)], centre, rotation, &posit[(current_time*max_affected*3 + body_no*3)]);
	calcRotatedVector(&veloc[(body_no*3)], centre, rotation, &veloc[(body_no*3)]);

	// Calculating initial acceleration:

	calcGravAccel(&posit[3*max_affected*current_time + 3*body_no], &effector_posit[3*max_effector*current_time], effector_mass, body_no*same - (1 - same), num_effector, &accel[3*body_no]);
}

void calcPolygonVertices(double* centre, double* rotation, double radius, int no_of_vertices, double* vertices_positions) 
{
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
	//
	// Returns (vertices_position) 3d vector positions of the vertices of a regular polygon with radius (radius) and number of sides (no_of_vertices),
	// centered around a point (centre) with rotation (rotation).
	//
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

	for (size_t i = 0; i < (no_of_vertices+1); ++i)
	{
				
		vertices_positions[(i*3) + 0] = centre[0] + cos((2*M_PI)/no_of_vertices*i)*radius;
		vertices_positions[(i*3) + 1] = centre[1] + sin((2*M_PI)/no_of_vertices*i)*radius;
		vertices_positions[(i*3) + 2] = centre[2];
				
		calcRotatedVector(&vertices_positions[i*3], centre, rotation, &vertices_positions[(i*3)]);
			
	}
}

void calcSphereRandPoints(double* centre, double radius, int no_of_vertices, double* vertices_positions) 
{
	for (int i = 0; i < no_of_vertices; ++i)
	{
		for (size_t d = 0; d < 3; ++d) { vertices_positions[3*i + d] = (double)rand_64()/(double)RAND_64_MAX;}
			
		double mag = calcAbs(&vertices_positions[3*i]);

		for (size_t d = 0; d < 3; ++d) { vertices_positions[3*i + d] /= mag; vertices_positions[3*i + d] *= radius; vertices_positions[3*i + d] += centre[d]; }
	}
}

void calcVelocityFromPoint(double* posit, double* centre, double* init_veloc, double extra_veloc, int num_n_bodies, double* veloc)
{
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
	//
	// Returns (veloc) velocities of an array of vectors (posit) of number (num_n_bodies) given a constant velocity (extra_veloc) away from a point (centre) plus an initial
	// velocity (init_veloc) 
	//
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //


	for (size_t i = 0; i < (num_n_bodies+1); ++i)
	{

	 	// Calculating diplacement from point and assigning to variables:

		double displace[3];
		calcDifference(&posit[(i*3)], centre, displace);

		// Calculating magnetude of diplacement:

		double mag = calcAbs(displace);

		// Assigning new velocity:======================
		for (size_t d = 0; d < 3; ++d) { veloc[(i*3) + d] = init_veloc[d] + extra_veloc*(displace[d]/mag); }

	}
}

void addBodiesFromInputFile(double input_bodies, char* input_filename, size_t* current_bodies, size_t* effector_bodies, size_t current_time, size_t max_bodies, size_t max_n_bodies, double* mass, double* radius, double* posit, 
							double* veloc, double* accel, double* effector_posit, double* effector_mass, bool same)
{
	size_t val_current_bodies = *current_bodies;
	
	if (input_bodies <= max_bodies - val_current_bodies)
	{

		size_t bodies_read = 0;

		if (input_filename != NULL)
		{
			readFile(input_filename, &mass[val_current_bodies], &radius[val_current_bodies], &posit[3*max_bodies*current_time + 3*val_current_bodies], &veloc[3*val_current_bodies], input_bodies, &bodies_read);
			printf("Read %d bodies.\n", bodies_read);
			*current_bodies += input_bodies;

			// Putting mass in units of G:

			for (size_t i = 0; i < input_bodies; ++i) { mass[i] = mass[i]*G; }

			// Calculating initial acceleration:

			calcAllGravAccel(&posit[3*max_bodies*current_time + 3*val_current_bodies], &effector_posit[3*max_n_bodies*current_time], effector_mass, same, *current_bodies, *effector_bodies, &accel[3*val_current_bodies]);


		}

		else { printf("Could not find input file!\n"); }



	}

	else { printf("Max bodies reached, unable to add! Total bodies: %d, Max Bodies: %d, Trying to add: %d. \n", current_bodies, max_bodies, input_bodies); }
}

void addBodiesInRingAroundObject(int ring_bodies, double ring_radius, double* rotation, double ring_veloc, double* object_posit, double* inital_veloc, size_t* current_bodies, size_t* effector_bodies, size_t current_time,
								size_t max_bodies, size_t max_n_bodies, double* posit, double* veloc, double* accel, double* effector_posit, double* effector_mass, bool same)
{
	size_t val_current_bodies = *current_bodies;
	size_t val_effector_bodies = *effector_bodies;


	if (ring_bodies <= max_n_bodies - *current_bodies)
	{
		calcPolygonVertices(object_posit, rotation, ring_radius, ring_bodies, &posit[3*max_n_bodies*current_time + 3*val_current_bodies]);
		calcVelocityFromPoint(posit, object_posit, inital_veloc, ring_veloc, ring_bodies, &veloc[3*val_current_bodies]);
		*current_bodies += ring_bodies;

		// Calculating initial acceleration:

		calcAllGravAccel(&posit[3*max_bodies*current_time + 3*val_current_bodies], &effector_posit[3*max_n_bodies*current_time], effector_mass, same, *current_bodies, *effector_bodies, &accel[3*val_current_bodies]);


	}

	else { printf("Max bodies reached, unable to add! Total bodies: %d, Max Bodies: %d, Trying to add: %d. \n", current_bodies, max_bodies, ring_bodies); }
}

void addBodiesInRandSphereAroundObject(int sphere_bodies, double sphere_radius, double sphere_veloc, double* object_posit, double* inital_veloc, size_t* current_bodies, size_t* effector_bodies, size_t current_time,
								size_t max_bodies, size_t max_n_bodies, double* posit, double* veloc, double* accel, double* effector_posit, double* effector_mass, bool same)

{
	size_t val_current_bodies = *current_bodies;
	size_t val_effector_bodies = *effector_bodies;

	if (sphere_bodies <= max_n_bodies - *current_bodies)
	{
		calcSphereRandPoints(object_posit, sphere_radius, sphere_bodies, &posit[3*max_n_bodies*current_time + 3*val_current_bodies]);
		calcVelocityFromPoint(posit, object_posit, inital_veloc, sphere_veloc, sphere_bodies, &veloc[3*val_current_bodies]);
		*current_bodies += sphere_bodies;

		//Randomise Velocities

		for (size_t i = 0; i < sphere_bodies; ++i) 
		{
			double mag = calcAbs(&veloc[3*i]);
			for (size_t d = 0; d < 3; ++d) { veloc[3*i + d] += 3*(veloc[3*i + d]/mag)*sphere_veloc*(double)rand_64()/(double)RAND_64_MAX;}
		}

		// Calculating initial acceleration:

		calcAllGravAccel(&posit[3*max_bodies*current_time + 3*val_current_bodies], &effector_posit[3*max_n_bodies*current_time], effector_mass, same, *current_bodies, *effector_bodies, &accel[3*val_current_bodies]);

	}

	else { printf("Max bodies reached, unable to add! Total bodies: %d, Max Bodies: %d, Trying to add: %d. \n", current_bodies, max_bodies, sphere_bodies); }

}

void runSim(size_t time_steps, double run_time, size_t max_time_steps, size_t* current_time, double* posit, double* n_posit, double* veloc, double* n_veloc, double* accel, double* n_accel, 
			double* effector_mass, double* effector_posit, int current_bodies, int current_n_bodies, int max_bodies, int max_n_bodies)
{
	//Looping over timesteps using velocity verlet intergation:

	if ((time_steps + *current_time) > max_time_steps)
	{
		printf("Timesteps exceeds max time steps! Current Timesteps: \n", current_time);
		fflush(stdout);

		exit(1);
 
	}	

	else 
	{
		printf("Running from 0 - %fs, over %d timesteps, with %d bodies and %d ejecta...\n", run_time, time_steps, current_bodies, current_n_bodies);
		fflush(stdout);
	}

	//Calculate timestep:

	double dt = ((run_time)/(double)time_steps);

	*current_time += 1;

  	for (int t = *current_time; t < *current_time + time_steps; ++t)
	{
		
		//Looping over gravitational bodies:

		calcGravVerlet(t, posit, veloc, accel, effector_mass, posit, dt, current_bodies, current_bodies, max_bodies, max_bodies, true);

		//Looping over gravitationally negligable bodies:

		calcGravVerlet(t, n_posit, n_veloc, n_accel, effector_mass, posit, dt, current_n_bodies, current_bodies, max_n_bodies, max_bodies, false);

	}

	*current_time += time_steps;

	//Ensuring loop isn't removed:

	printf("%f %f %f\n", posit[(current_n_bodies*3) + (0*3) + 0], posit[(current_n_bodies*3) + (0*3) + 1], posit[(current_n_bodies*3) + (0*3) + 2]);
	fflush(stdout);

	printf("\nSimulation completed.\n");
}

void resetSim(size_t max_time_steps, size_t max_bodies, size_t max_n_bodies, size_t* current_bodies, size_t* current_n_bodies, size_t* current_time, double* mass, double* radius, double* posit, 
				double* n_posit, double* veloc, double* n_veloc, double* accel, double* n_accel)
{

	printf("Reseting simulation...\n");
	
	*current_time = 0;
	*current_bodies = 0;
	*current_n_bodies = 0;

	for (int i = 0; i < max_bodies*max_time_steps*3; ++i) { posit[i] = 0.0; }
	for (int i = 0; i < max_n_bodies*max_time_steps*3; ++i) { n_posit[i] = 0.0; }
	for (int i = 0; i < max_bodies*3; ++i) { veloc[i] = 0.0; accel[i] = 0.0; }
	for (int i = 0; i < max_n_bodies*3; ++i) { n_veloc[i] = 0.0; n_accel[i] = 0.0; }
	for (int i = 0; i < max_bodies; ++i) { mass[i] = 0.0; radius[i] = 0.0;}

	printf("Reset complete.\n");
}

void clearSim(size_t max_bodies, size_t max_n_bodies, size_t current_time, size_t* current_bodies, size_t* current_n_bodies, double* mass, double* radius, double* posit, double* n_posit, double* veloc, 
				double* n_veloc, double* accel, double* n_accel)
{
	printf("Removing all bodies from simulation...\n");

	*current_bodies = 0;
	*current_n_bodies = 0;

	for (int i = 0; i < max_bodies*3; ++i) { posit[(current_time*max_bodies*3) + i] = 0.0; }
	for (int i = 0; i < max_n_bodies*3; ++i) { n_posit[(current_time*max_n_bodies*3) + i] = 0.0; }
	for (int i = 0; i < max_bodies*3; ++i) { veloc[i] = 0.0; accel[i] = 0.0; }
	for (int i = 0; i < max_n_bodies*3; ++i) { n_veloc[i] = 0.0; n_accel[i] = 0.0; }
	for (int i = 0; i < max_bodies; ++i) { mass[i] = 0.0; radius[i] = 0.0;}

	printf("Bodies removed.\n");
}

void printPosit(char* output_filename, size_t detail, int start_time_steps, int end_time_steps, int prec, double* posit, double* n_posit, size_t max_bodies, size_t max_n_bodies, size_t current_bodies, size_t current_n_bodies, int print_number)
{
	int no_lines = (end_time_steps - start_time_steps)/detail;
	int start_lines = start_time_steps/detail;

	if (output_filename != NULL)
	{

		char filename_positions[512];
   		sprintf(filename_positions, "%s_posit_%d.csv", output_filename, print_number);

		printf("Printing positions to file: %s\n", filename_positions);
		FILE* f = fopen(filename_positions, "w");
		for (size_t t = start_lines; t < (start_lines + no_lines); ++t)
		{
			for (size_t i = 0; i < current_bodies; ++i)
			{
				fprintf(f, "%.*e,%.*e,%.*e,", prec, posit[(detail*t*max_bodies*3) + (i*3) + 0], prec, posit[(detail*t*max_bodies*3) + (i*3) + 1], prec, posit[(detail*t*max_bodies*3) + (i*3) + 2]);
			}
			for (size_t i = 0; i < current_n_bodies; ++i)
			{
				fprintf(f, "%.*e,%.*e,%.*e,", prec, n_posit[(detail*t*max_n_bodies*3) + (i*3) + 0], prec, n_posit[(detail*t*max_n_bodies*3) + (i*3) + 1], prec, n_posit[(detail*t*max_n_bodies*3) + (i*3) + 2] );
			}
			fprintf(f, "%s\n");

		}

		fclose(f);

		printf("Printing positions completed. \n");
	}
	else printf("No output selected. Canceling printing. \n" );
}

void printSpecial(char* output_filename, size_t detail, int start_time_steps, int end_time_steps, int prec, double* origin_displace, double* target_displace, size_t max_n_bodies, size_t current_n_bodies, int print_number)
{
	int no_lines = (end_time_steps - start_time_steps)/detail;
	int start_lines = start_time_steps/detail;

	if (output_filename != NULL)
	{

		char filename_displace[512];
   		sprintf(filename_displace, "%s_origin_displace_%d.csv", output_filename, print_number);

		printf("Printing origin displacement to file: %s\n", filename_displace);

		FILE* f = fopen(filename_displace, "w");
		for (size_t t = start_lines; t < (start_lines + no_lines); ++t)
		{
			for (size_t i = 0; i < current_n_bodies; ++i)
			{
				fprintf(f, "%.*e,", prec, origin_displace[(detail*t*max_n_bodies) + i]);
			}
			fprintf(f, "%s\n");

		}

		fclose(f);

		char filename_displace_2[512];
   		sprintf(filename_displace_2, "%s_target_displace_%d.csv", output_filename, print_number);

		printf("Printing target displacement to file: %s\n", filename_displace_2);

		FILE* f2 = fopen(filename_displace_2, "w");
		for (size_t t = start_lines; t < (start_lines + no_lines); ++t)
		{
			for (size_t i = 0; i < current_n_bodies; ++i)
			{
				fprintf(f2, "%.*e,", prec, target_displace[(detail*t*max_n_bodies) + i]);
			}
			fprintf(f2, "%s\n");

		}

		fclose(f2);

		printf("Printing special completed. \n");
	}
	else printf("No output selected. Canceling printing. \n" );
}

void printMax(char* output_filename, char* filename_suffix, int start_time_steps, int end_time_steps, int prec, double* array, size_t max_n_bodies, size_t current_n_bodies, int print_number)
{
	
	if (output_filename != NULL)
	{
		char filename_displace[512];
   		sprintf(filename_displace, "%s_%s_%d.csv", output_filename, filename_suffix, print_number);

		printf("Printing to file: %s\n", filename_displace);

		FILE* f = fopen(filename_displace, "w");
		for (size_t t = start_time_steps; t < (start_time_steps + end_time_steps); ++t)
		{
			for (size_t i = 0; i < current_n_bodies; ++i)
			{
				fprintf(f, "%.*e,", prec, array[(t*max_n_bodies) + i]);
			}
			fprintf(f, "%s\n");
		}

	}
	else printf("No output selected. Canceling printing. \n" );
}

int main(int argc, char** argv)
{

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

	size_t max_bodies = 6;
	size_t current_bodies = 0;

	size_t max_n_bodies = 36;
	size_t current_n_bodies = 0;

	size_t input_bodies = 1;
	char* input_filename = NULL;

	int prec = 7;
	int detail = 10;
	char* output_filename = NULL;
	
	size_t current_time = 0;
	int max_time_steps = 4000000;

	double run_time = 100*365.25*24*3600;
	int time_steps = 300000;

	size_t repeat = 40;

	size_t ring_assigned = 1;
	double assigned_angle = 0;
	size_t ring_n_bodies = 0;
	
	double ring_radius = 621000000;
	double ring_veloc = 1000;

	size_t target = 3;

	double veloc_increase = 200;

	//Reading config file:

	readConfig("config.csv", &input_bodies, &input_filename, &prec, &output_filename, &time_steps, &run_time, &ring_assigned, &assigned_angle, &ring_n_bodies, &ring_radius, 
				&ring_veloc, &target);

	//Reading console arguments

	readInputArgs(argc, argv, &input_bodies, &input_filename, &prec, &output_filename, &time_steps, &run_time, &ring_assigned, &assigned_angle, &ring_n_bodies, &ring_radius, 
				&ring_veloc, &target);

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

	//Creating z_rotation and origin vector

	double z_rotation[4] = {0, 0, 1, assigned_angle};
	double origin[3] = {0,0,0};

	//Allocating Memory for gravitational bodies

	double* mass = (double*)malloc(sizeof(double)*max_bodies);
	double* radius = (double*)malloc(sizeof(double)*max_bodies);

	double* posit = (double*)malloc(sizeof(double)*max_bodies*max_time_steps*3);
	double* veloc = (double*)calloc(max_bodies*3, sizeof(double));
	double* accel = (double*)calloc(max_bodies*3, sizeof(double));
	//Allocating Memory for non gravitation bodies

	double* n_posit = (double*)malloc(sizeof(double)*max_n_bodies*max_time_steps*3);
	double* n_veloc = (double*)calloc(max_n_bodies*3, sizeof(double));
	double* n_accel = (double*)calloc(max_n_bodies*3, sizeof(double));

	//Allocating Memory for Special arrays

	double* origin_displace = (double*)malloc(sizeof(double)*max_n_bodies*max_time_steps);
	double* target_displace = (double*)malloc(sizeof(double)*max_n_bodies*max_time_steps);

	double* max_origin_displace = (double*)malloc(sizeof(double)*max_n_bodies*repeat);
	double* max_target_displace = (double*)malloc(sizeof(double)*max_n_bodies*repeat);

	//printf("Total memory allocated = %dMB\n", sizeof(double)*(num_bodies + num_n_bodies)*((max_time_steps*3 + 5)/1000000));

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SIMULATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

	for (size_t i = 0; i < repeat; ++i)
	{
		//Adding bodies read from input_file:

		z_rotation[3] = (360*(double)rand_64()/(double)RAND_64_MAX); //< --- Rotaty buisness

		addBodiesFromInputFile(input_bodies, input_filename, &current_bodies, &current_bodies, current_time, max_bodies, max_n_bodies, mass, radius, posit, veloc, accel, posit, mass, true);
		
		//Rotate assigned planet and velocity around origin:

		for (size_t j = 0; j < current_bodies; ++j)
		{
			z_rotation[3] = (360*(double)rand_64()/(double)RAND_64_MAX);
			rotateBody(posit, veloc, accel, posit, mass, origin, z_rotation, j, max_bodies, max_bodies, current_bodies, current_time, true);
		}
	
		//Assigning ring positions and velocities:

		addBodiesInRandSphereAroundObject(ring_n_bodies, ring_radius, ring_veloc, &posit[(ring_assigned*3)], &veloc[(ring_assigned*3)], &current_n_bodies, &current_bodies,
									 current_time, max_bodies, max_n_bodies, n_posit, n_veloc, n_accel, posit, mass, false);

		printf("Angle = %f \n", z_rotation[3]);

		runSim(time_steps, run_time, max_time_steps, &current_time, posit, n_posit, veloc, n_veloc, accel, n_accel, mass, posit, current_bodies, current_n_bodies, max_bodies, max_n_bodies);

		//calcDisplcaement(posit, n_posit, 0, current_n_bodies, max_n_bodies, max_bodies, 0, time_steps, origin_displace);
		calcDisplcaement(posit, n_posit, target, current_n_bodies, max_n_bodies, max_bodies, 0, time_steps, target_displace);

		//calcMax(origin_displace, time_steps, current_n_bodies, max_n_bodies, &max_origin_displace[i*max_n_bodies]);
		calcMin(target_displace, time_steps, current_n_bodies, max_n_bodies, &max_target_displace[i*max_n_bodies]);
		
		if (i != (repeat-1)) { resetSim(max_time_steps, max_bodies, max_n_bodies, &current_bodies, &current_n_bodies, &current_time, mass, radius, posit, n_posit, veloc, n_veloc, accel, n_accel); }

	}

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PRINTING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

	// ~~~~~~~~~~ Printing to file ~~~~~~~~~~ //

	//printMax(output_filename, "_max_origin-veloc0-50", 0, repeat, prec, max_origin_displace, max_n_bodies, current_n_bodies, 1);
	printMax(output_filename, "_max_target-displace-sphere", 0, repeat, prec, max_target_displace, max_n_bodies, current_n_bodies, 106);

	//printPosit(output_filename, detail, 0, time_steps, prec, posit, n_posit, max_bodies, max_n_bodies, current_bodies, current_n_bodies, 1);
	//printSpecial(output_filename, detail, 0, time_steps, prec, origin_displace, target_displace, max_n_bodies, current_n_bodies, 1);

	printf("\nProgram terminated.\n");

	return 0;
}

