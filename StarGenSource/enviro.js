#include	<stdio.h>
#include	<math.h>

/*************************** #include	"const.h" *********************************/
#include <math.h>

#define PI						(3.1415926536)
#define RADIANS_PER_ROTATION	(2.0 * PI)

#ifndef TRUE
#define TRUE					(1)
#define FALSE					(0)
#endif

#define ECCENTRICITY_COEFF		(0.077)			/* Dole's was 0.077			*/
#define PROTOPLANET_MASS		(1.0E-15)		/* Units of solar masses	*/
#define CHANGE_IN_EARTH_ANG_VEL (-1.3E-15)		/* Units of radians/sec/year*/
#define SOLAR_MASS_IN_GRAMS		(1.989E33)		/* Units of grams			*/
#define SOLAR_MASS_IN_KILOGRAMS	(1.989E30)		/* Units of kg				*/
#define EARTH_MASS_IN_GRAMS		(5.977E27)		/* Units of grams			*/
#define EARTH_RADIUS			(6.378E8)		/* Units of cm				*/
#define EARTH_DENSITY			(5.52)			/* Units of g/cc			*/
#define KM_EARTH_RADIUS			(6378.0)		/* Units of km				*/
//      EARTH_ACCELERATION		(981.0)			/* Units of cm/sec2			*/
#define EARTH_ACCELERATION		(980.7)			/* Units of cm/sec2			*/
#define EARTH_AXIAL_TILT		(23.4)			/* Units of degrees			*/
#define EARTH_EXOSPHERE_TEMP	(1273.0)		/* Units of degrees Kelvin	*/
#define SUN_MASS_IN_EARTH_MASSES (332775.64)
#define ASTEROID_MASS_LIMIT		(0.001)			/* Units of Earth Masses	*/
#define EARTH_EFFECTIVE_TEMP	(250.0)			/* Units of degrees Kelvin (was 255)	*/
#define CLOUD_COVERAGE_FACTOR	(1.839E-8)		/* Km2/kg					*/
#define EARTH_WATER_MASS_PER_AREA	 (3.83E15)	/* grams per square km		*/
#define EARTH_SURF_PRES_IN_MILLIBARS (1013.25)
#define EARTH_SURF_PRES_IN_MMHG	(760.)			/* Dole p. 15				*/
#define EARTH_SURF_PRES_IN_PSI	(14.696)		/* Pounds per square inch	*/
#define MMHG_TO_MILLIBARS (EARTH_SURF_PRES_IN_MILLIBARS / EARTH_SURF_PRES_IN_MMHG)
#define PSI_TO_MILLIBARS (EARTH_SURF_PRES_IN_MILLIBARS / EARTH_SURF_PRES_IN_PSI)
#define H20_ASSUMED_PRESSURE	(47. * MMHG_TO_MILLIBARS) /* Dole p. 15      */
#define MIN_O2_IPP	(72. * MMHG_TO_MILLIBARS)	/* Dole, p. 15				*/
#define MAX_O2_IPP	(400. * MMHG_TO_MILLIBARS)	/* Dole, p. 15				*/
#define MAX_HE_IPP	(61000. * MMHG_TO_MILLIBARS)	/* Dole, p. 16			*/
#define MAX_NE_IPP	(3900. * MMHG_TO_MILLIBARS)	/* Dole, p. 16				*/
#define MAX_N2_IPP	(2330. * MMHG_TO_MILLIBARS)	/* Dole, p. 16				*/
#define MAX_AR_IPP	(1220. * MMHG_TO_MILLIBARS)	/* Dole, p. 16				*/
#define MAX_KR_IPP	(350. * MMHG_TO_MILLIBARS)	/* Dole, p. 16				*/
#define MAX_XE_IPP	(160. * MMHG_TO_MILLIBARS)	/* Dole, p. 16				*/
#define MAX_CO2_IPP (7. * MMHG_TO_MILLIBARS)	/* Dole, p. 16				*/
#define MAX_HABITABLE_PRESSURE (118 * PSI_TO_MILLIBARS)	/* Dole, p. 16		*/
// The next gases are listed as poisonous in parts per million by volume at 1 atm:
#define PPM_PRSSURE (EARTH_SURF_PRES_IN_MILLIBARS / 1000000.)
#define MAX_F_IPP	(0.1 * PPM_PRSSURE)			/* Dole, p. 18				*/
#define MAX_CL_IPP	(1.0 * PPM_PRSSURE)			/* Dole, p. 18				*/
#define MAX_NH3_IPP	(100. * PPM_PRSSURE)		/* Dole, p. 18				*/
#define MAX_O3_IPP	(0.1 * PPM_PRSSURE)			/* Dole, p. 18				*/
#define MAX_CH4_IPP	(50000. * PPM_PRSSURE)		/* Dole, p. 18				*/



#define EARTH_CONVECTION_FACTOR (0.43)			/* from Hart, eq.20			*/
//      FREEZING_POINT_OF_WATER (273.0)			/* Units of degrees Kelvin	*/
#define FREEZING_POINT_OF_WATER (273.15)		/* Units of degrees Kelvin	*/
//      EARTH_AVERAGE_CELSIUS   (15.5)			/* Average Earth Temperature */
#define EARTH_AVERAGE_CELSIUS   (14.0)			/* Average Earth Temperature */
#define EARTH_AVERAGE_KELVIN    (EARTH_AVERAGE_CELSIUS + FREEZING_POINT_OF_WATER)
#define DAYS_IN_A_YEAR			(365.256)		/* Earth days per Earth year*/
//		gas_retention_threshold = 5.0;  		/* ratio of esc vel to RMS vel */
#define GAS_RETENTION_THRESHOLD (6.0)			/* ratio of esc vel to RMS vel */

#define ICE_ALBEDO				(0.7)
#define CLOUD_ALBEDO			(0.52)
#define GAS_GIANT_ALBEDO		(0.5)			/* albedo of a gas giant	*/
#define AIRLESS_ICE_ALBEDO		(0.5)
#define EARTH_ALBEDO			(0.3)			/* was .33 for a while */
#define GREENHOUSE_TRIGGER_ALBEDO (0.20)
#define ROCKY_ALBEDO			(0.15)
#define ROCKY_AIRLESS_ALBEDO	(0.07)
#define WATER_ALBEDO			(0.04)

#define SECONDS_PER_HOUR		(3600.0)
#define CM_PER_AU				(1.495978707E13)/* number of cm in an AU	*/
#define CM_PER_KM				(1.0E5)			/* number of cm in a km		*/
#define KM_PER_AU				(CM_PER_AU / CM_PER_KM)
#define CM_PER_METER			(100.0)
//#define MILLIBARS_PER_BAR		(1013.25)
#define MILLIBARS_PER_BAR		(1000.00)

#define GRAV_CONSTANT			(6.672E-8)		/* units of dyne cm2/gram2	*/
#define MOLAR_GAS_CONST			(8314.41)		/* units: g*m2/(sec2*K*mol) */
#define K						(50.0)			/* K = gas/dust ratio		*/
#define B						(1.2E-5)		/* Used in Crit_mass calc	*/
#define DUST_DENSITY_COEFF		(2.0E-3)		/* A in Dole's paper		*/
#define ALPHA					(5.0)			/* Used in density calcs	*/
#define N						(3.0)			/* Used in density calcs	*/
#define J						(1.46E-19)		/* Used in day-length calcs (cm2/sec2 g) */
#ifdef HUGE_VAL
#define INCREDIBLY_LARGE_NUMBER HUGE_VAL
#else
#define INCREDIBLY_LARGE_NUMBER (9.9999E37)
#endif

/*	Now for a few molecular weights (used for RMS velocity calcs):	   */
/*	This table is from Dole's book "Habitable Planets for Man", p. 38  */

#define ATOMIC_HYDROGEN			(1.0)	/* H   */
#define MOL_HYDROGEN			(2.0)	/* H2  */
#define HELIUM					(4.0)	/* He  */
#define ATOMIC_NITROGEN			(14.0)	/* N   */
#define ATOMIC_OXYGEN			(16.0)	/* O   */
#define METHANE					(16.0)	/* CH4 */
#define AMMONIA					(17.0)	/* NH3 */
#define WATER_VAPOR				(18.0)	/* H2O */
#define NEON					(20.2)	/* Ne  */
#define MOL_NITROGEN			(28.0)	/* N2  */
#define CARBON_MONOXIDE			(28.0)	/* CO  */
#define NITRIC_OXIDE			(30.0)	/* NO  */
#define MOL_OXYGEN				(32.0)	/* O2  */
#define HYDROGEN_SULPHIDE		(34.1)	/* H2S */
#define ARGON					(39.9)	/* Ar  */
#define CARBON_DIOXIDE			(44.0)	/* CO2 */
#define NITROUS_OXIDE			(44.0)	/* N2O */
#define NITROGEN_DIOXIDE		(46.0)	/* NO2 */
#define OZONE					(48.0)	/* O3  */
#define SULPH_DIOXIDE			(64.1)	/* SO2 */
#define SULPH_TRIOXIDE			(80.1)	/* SO3 */
#define KRYPTON					(83.8)	/* Kr  */
#define XENON					(131.3) /* Xe  */

//	And atomic numbers, for use in ChemTable indexes
#define AN_H	1
#define AN_HE	2
#define AN_N	7
#define AN_O	8
#define AN_F	9
#define AN_NE	10
#define AN_P	15
#define AN_CL	17
#define AN_AR	18
#define AN_BR	35
#define AN_KR	36
#define AN_I	53
#define AN_XE	54
#define AN_HG	80
#define AN_AT	85
#define AN_RN	86
#define AN_FR	87

#define AN_NH3	900
#define AN_H2O	901
#define AN_CO2	902
#define AN_O3	903
#define AN_CH4	904
#define AN_CH3CH2OH	905

/*	The following defines are used in the kothari_radius function in	*/
/*	file enviro.c.														*/
#define A1_20					(6.485E12)		/* All units are in cgs system.	 */
#define A2_20					(4.0032E-8)		/*	 ie: cm, g, dynes, etc.		 */
#define BETA_20					(5.71E12)

#define JIMS_FUDGE				(1.004)

/*	 The following defines are used in determining the fraction of a planet	 */
/*	covered with clouds in function cloud_fraction in file enviro.c.		 */
#define Q1_36					(1.258E19)		/* grams	*/
#define Q2_36					(0.0698)		/* 1/Kelvin */

/* macros: */
#define pow2(a) ((a) * (a))
#define pow3(a) ((a) * (a) * (a))
#define pow4(a) ((a) * (a) * (a) * (a))
#define pow1_4(a)		sqrt(sqrt(a))
#define pow1_3(a)		pow(a,(1.0/3.0))




/*************************** #include	"structs.h" *********************************/
#pragma once

typedef struct dust_record	*dust_pointer;
typedef struct planets_record  *planet_pointer;
typedef struct gen *gen_pointer;

typedef enum planet_type {
	tUnknown,
	tRock,
	tVenusian,
	tTerrestrial,
	tGasGiant,
	tMartian,
	tWater,
	tIce,
	tSubGasGiant,
	tSubSubGasGiant,
	tAsteroids,
	t1Face
} planet_type;

typedef struct gas {
  	int         num;
	long double	surf_pressure;		/* units of millibars (mb)			 */
	} gas;

typedef struct sun {
	long double	luminosity;
	long double	mass;
	long double life;
	long double age;
	long double r_ecosphere;
	char		*name;
	} sun;

typedef struct planets_record {
	int			planet_no;
	long double a;					/* semi-major axis of solar orbit (in AU)*/
	long double e;					/* eccentricity of solar orbit		 */
	long double	axial_tilt;			/* units of degrees					 */
	long double mass;				/* mass (in solar masses)			 */
	int 		gas_giant;			/* TRUE if the planet is a gas giant */
	long double	dust_mass;			/* mass, ignoring gas				 */
	long double	gas_mass;			/* mass, ignoring dust				 */
									/*   ZEROES start here               */
	long double moon_a;				/* semi-major axis of lunar orbit (in AU)*/
	long double moon_e;				/* eccentricity of lunar orbit		 */
	long double	core_radius;		/* radius of the rocky core (in km)	 */
	long double radius;				/* equatorial radius (in km)		 */
	int 		orbit_zone;			/* the 'zone' of the planet			 */
	long double density;			/* density (in g/cc)				 */
	long double orb_period;			/* length of the local year (days)	 */
	long double day;				/* length of the local day (hours)	 */
	int 		resonant_period;	/* TRUE if in resonant rotation		 */
	long double	esc_velocity;		/* units of cm/sec					 */
	long double	surf_accel;			/* units of cm/sec2					 */
	long double	surf_grav;			/* units of Earth gravities			 */
	long double	rms_velocity;		/* units of cm/sec					 */
	long double molec_weight;		/* smallest molecular weight retained*/
	long double	volatile_gas_inventory;
	long double	surf_pressure;		/* units of millibars (mb)			 */
	int		 	greenhouse_effect;	/* runaway greenhouse effect?		 */
	long double	boil_point;			/* the boiling point of water (Kelvin)*/
	long double	albedo;				/* albedo of the planet				 */
	long double	exospheric_temp;	/* units of degrees Kelvin			 */
	long double estimated_temp;     /* quick non-iterative estimate (K)  */
	long double estimated_terr_temp;/* for terrestrial moons and the like*/
	long double	surf_temp;			/* surface temperature in Kelvin	 */
	long double	greenhs_rise;		/* Temperature rise due to greenhouse */
	long double high_temp;			/* Day-time temperature              */
	long double low_temp;			/* Night-time temperature			 */
	long double max_temp;			/* Summer/Day						 */
	long double min_temp;			/* Winter/Night						 */
	long double	hydrosphere;		/* fraction of surface covered		 */
	long double	cloud_cover;		/* fraction of surface covered		 */
	long double	ice_cover;			/* fraction of surface covered		 */
	sun*		sun;
	int			gases;				/* Count of gases in the atmosphere: */
	gas*		atmosphere;
	planet_type type;				/* Type code						 */
	int			minor_moons;
	planet_pointer first_moon;
									/*   ZEROES end here               */
	planet_pointer next_planet;
	} planets;

/*	Define the solar system for comparisons, etc. */
#define ZEROES 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,NULL,0,NULL,tUnknown

typedef struct dust_record {
	long double inner_edge;
	long double outer_edge;
	int 		dust_present;
	int 		gas_present;
	dust_pointer next_band;
	 } dust;

typedef struct star {
	long double		luminosity;
	long double		mass;
	long double		m2;
	long double		e;
	long double		a;
	planet_pointer	known_planets;
	char			*desig;
	int				in_celestia;
	char			*name;
	} star;

typedef struct catalog {
	int				count;
	char*			arg;
	star			(*stars)[];
	} catalog;

typedef	struct gen
{
	dust_pointer	dusts;
	planet_pointer	planets;
	gen_pointer		next;
} generation;

// From Keris

typedef struct ChemTableS
{
  int         	num;
  char       	*symbol;
  char			*html_symbol;
  char       	*name;
  long double	weight;
  long double	melt;
  long double	boil;
  long double	density;
  long double	abunde;
  long double	abunds;
  long double	reactivity;
  long double	max_ipp;	// Max inspired partial pressure im millibars
} ChemTable;




/*************************** #include	"enviro.h" *********************************/
long double luminosity(long double);
int orb_zone(long double, long double);
long double volume_radius(long double, long double);
long double kothari_radius(long double, int, int);
long double empirical_density(long double, long double, long double, int);
long double volume_density(long double, long double);
long double period(long double, long double, long double);
long double day_length(planet_pointer);
int inclination(long double);
long double escape_vel(long double, long double);
long double rms_vel(long double, long double);
long double molecule_limit(long double, long double, long double);
long double min_molec_weight (planet_pointer);
long double acceleration(long double, long double);
long double gravity(long double);
long double vol_inventory(long double, long double, long double, long double, int, int, int);
long double pressure(long double, long double, long double);
long double boiling_point(long double);
long double hydro_fraction(long double, long double);
long double cloud_fraction(long double, long double, long double, long double);
long double ice_fraction(long double, long double);
long double eff_temp(long double, long double, long double);
long double est_temp(long double, long double, long double);
int grnhouse(long double r_ecosphere, long double);
long double green_rise(long double, long double, long double);
long double planet_albedo(long double, long double, long double, long double);
long double opacity(long double, long double);
long double gas_life(long double, planet_pointer);
void calculate_surface_temp(planet_pointer, int, long double, long double,
							long double, long double, long double);
void iterate_surface_temp(planet_pointer);

long double inspired_partial_pressure (long double, long double);

unsigned int breathability (planet_pointer);

#define	NONE			0
#define	BREATHABLE		1
#define	UNBREATHABLE	2
#define	POISONOUS		3

extern char* breathability_phrase[4];

extern long double lim(long double x);
long double soft(long double v, long double max, long double min);
extern void set_temp_range(planet_pointer planet);



/*************************** #include	"stargen.h" *********************************/
/*
 *	StarGen main API
 */

typedef	enum actions {						// Callable StarGen can:
	aGenerate,								//	- Generate randon system(s)
	aListGases,								//	- List the gas table
	aListCatalog,							//	- List the stars in a catalog
	aListCatalogAsHTML,						//  - For creating a <FORM>
	aSizeCheck,								//  - List sizes of various types
	aListVerbosity,							//  - List values of the -v option
} actions;

int stargen (actions		action,			// One of the above
			 char			flag_char,
			 char *			path,			// OS path to where to write files
			 char *			url_path_arg,	// HTML path to parent of both the
			 								//  directory named in 'path' and
			 								//  the ref directory with images
			 char *			filename_arg,	// Output file name (optional)
			 char *			sys_name_arg,	// Human readble System name (opt.)
			 
			 FILE *			sgOut,			// Main stream to write to 
			 								//	Thumbnails will be written there
			 								//  for HTML format
			 FILE *			sgErr,			// Stream to write errors to (opt.)
			 char *			prognam,		// Name of program (opt.)
			 long double	mass_arg,		// Mass of star (not used with catalog)
			 long			seed_arg,		// Random number seed
			 int			count_arg,		// Number of systems (or cats) to do
			 int			incr_arg,		// Amount to increment seed by
			 catalog *		cat_arg,		// A star catalog (see below)
			 int			sys_no_arg,		// Star within a catalog (0 = all)
			 
			 long double	ratio_arg,		// Change dust density (experimental)
			 
			 int			flags_arg,		// Options (see below)
			 int			out_format,		// Output file formats (see below)
			 int			graphic_format	// Graphic file formats (see below)
			 );

										// Values of flags_arg:
#define	fUseSolarsystem			0x0001
#define	fReuseSolarsystem		0x0002
#define	fUseKnownPlanets		0x0004
#define fNoGenerate				0x0008

#define	fDoGases				0x0010
#define	fDoMoons				0x0020

#define fOnlyHabitable			0x0100
#define fOnlyMultiHabitable		0x0200
#define fOnlyJovianHabitable	0x0400
#define fOnlyEarthlike			0x0800

										// Values of out_format
#define	ffHTML				'HTML'
#define	ffTEXT				'TEXT'
#define	ffCELESTIA			'.SSC'
#define ffCSV				'.CSV'
#define ffCSVdl				'+CSV'
#define ffSVG				'.SVG'

										// Values of graphic_format
#define	gfGIF				'.GIF'
#define gfSVG				'.SVG'

										// The two predefined star catalogs.
extern catalog	solstation;
extern catalog	dole;
extern catalog  jimb;
										// You can roll your own (see main.c)

extern planets mercury;					// For building private catalogs


extern int          flag_verbose;		// Likely to move into stargen() args.

										// Various statistics that are kept:
extern int 		    total_earthlike;
extern int 		    total_habitable;

extern long double	min_breathable_terrestrial_g;
extern long double	min_breathable_g;
extern long double	max_breathable_terrestrial_g;
extern long double	max_breathable_g;
extern long double	min_breathable_terrestrial_l;
extern long double	min_breathable_l;
extern long double	max_breathable_terrestrial_l;
extern long double	max_breathable_l;
extern long double	min_breathable_temp;
extern long double	max_breathable_temp;
extern long double	min_breathable_p;
extern long double	max_breathable_p;

										// Experimental gas model variables
										//  Likely to go away or be changed
extern ChemTable    gases[];
extern int max_gas;

										// OS-specific constants for finding
										// the default output directory and
										// other dirs:
#ifdef WIN32
#define	SUBDIR	"html\\"
#define DIRSEP	"\\"
#else
#define	SUBDIR	"html/"
#define DIRSEP	"/"
#endif

extern char *	stargen_revision;		// RCS revision of stargen.c



/*************************** #include	"utils.h" *********************************/
long double random_number(long double, long double);
long double about(long double, long double);
long double random_eccentricity(void);

/************************************************************/


char* breathability_phrase[4] =
	{
		"none",
		"breathable",
		"unbreathable",
		"poisonous"
	};

long double luminosity(mass_ratio)
long double mass_ratio;
{
	long double n; 
	
	if (mass_ratio < 1.0)
		n = 1.75 * (mass_ratio - 0.1) + 3.325;
	else 
		n = 0.5 * (2.0 - mass_ratio) + 4.4;
	return(pow(mass_ratio,n));
}


/*--------------------------------------------------------------------------*/
/*	 This function, given the orbital radius of a planet in AU, returns		*/
/*	 the orbital 'zone' of the particle.									*/
/*--------------------------------------------------------------------------*/

int orb_zone(long double luminosity, long double orb_radius)
{
	if (orb_radius < (4.0 * sqrt(luminosity)))
		return(1);
	else if (orb_radius < (15.0 * sqrt(luminosity)))
		return(2);
	else 
		return(3);
}


/*--------------------------------------------------------------------------*/
/*	 The mass is in units of solar masses, and the density is in units		*/
/*	 of grams/cc.  The radius returned is in units of km.					*/
/*--------------------------------------------------------------------------*/

long double volume_radius(mass, density)
long double mass, density;
{
	long double volume; 
	
	mass = mass * SOLAR_MASS_IN_GRAMS;
	volume = mass / density;
	return(pow((3.0 * volume) / (4.0 * PI),(1.0 / 3.0)) / CM_PER_KM);
}

/*--------------------------------------------------------------------------*/
/*	 Returns the radius of the planet in kilometers.						*/
/*	 The mass passed in is in units of solar masses.						*/
/*	 This formula is listed as eq.9 in Fogg's article, although some typos	*/
/*	 crop up in that eq.  See "The Internal Constitution of Planets", by	*/
/*	 Dr. D. S. Kothari, Mon. Not. of the Royal Astronomical Society, vol 96 */
/*	 pp.833-843, 1936 for the derivation.  Specifically, this is Kothari's	*/
/*	 eq.23, which appears on page 840.										*/
/*--------------------------------------------------------------------------*/

long double kothari_radius(mass, giant, zone)
long double mass;
int giant, zone;
{
	volatile long double temp1;
	long double temp, temp2, atomic_weight, atomic_num;
	
	if (zone == 1)
	{
		if (giant)
		{
			atomic_weight = 9.5;
			atomic_num = 4.5;
		}
		else 
		{
			atomic_weight = 15.0;
			atomic_num = 8.0;
		}
	}
	else 
		if (zone == 2)
		{
			if (giant)
			{
				atomic_weight = 2.47;
				atomic_num = 2.0;
			}
			else 
			{
				atomic_weight = 10.0;
				atomic_num = 5.0;
			}
		}
		else 
		{
			if (giant)
			{
				atomic_weight = 7.0;
				atomic_num = 4.0;
			}
			else 
			{
				atomic_weight = 10.0;
				atomic_num = 5.0;
			}
		}
	
	temp1 = atomic_weight * atomic_num;
	
	temp = (2.0 * BETA_20 * pow(SOLAR_MASS_IN_GRAMS,(1.0 / 3.0)))
	     / (A1_20 * pow(temp1, (1.0 / 3.0)));
	
	temp2 = A2_20 * pow(atomic_weight,(4.0 / 3.0)) * pow(SOLAR_MASS_IN_GRAMS,(2.0 / 3.0));
	temp2 = temp2 * pow(mass,(2.0 / 3.0));
	temp2 = temp2 / (A1_20 * pow2(atomic_num));
	temp2 = 1.0 + temp2;
	temp = temp / temp2;
	temp = (temp * pow(mass,(1.0 / 3.0))) / CM_PER_KM;
	
	temp /= JIMS_FUDGE;			/* Make Earth = actual earth */
	
	return(temp);
}


/*--------------------------------------------------------------------------*/
/*	The mass passed in is in units of solar masses, and the orbital radius	*/
/*	is in units of AU.	The density is returned in units of grams/cc.		*/
/*--------------------------------------------------------------------------*/

long double empirical_density(long double mass, long double orb_radius, 
							  long double r_ecosphere, int gas_giant)
{
	long double temp; 
	
	temp = pow(mass * SUN_MASS_IN_EARTH_MASSES,(1.0 / 8.0));
	temp = temp * pow1_4(r_ecosphere / orb_radius);
	if (gas_giant)
		return(temp * 1.2);
	else 
		return(temp * 5.5);
}


/*--------------------------------------------------------------------------*/
/*	The mass passed in is in units of solar masses, and the equatorial		*/
/*	radius is in km.  The density is returned in units of grams/cc.			*/
/*--------------------------------------------------------------------------*/

long double volume_density(mass, equat_radius)
long double mass, equat_radius;
{
	long double volume; 
	
	mass = mass * SOLAR_MASS_IN_GRAMS;
	equat_radius = equat_radius * CM_PER_KM;
	volume = (4.0 * PI * pow3(equat_radius)) / 3.0;
	return(mass / volume);
}


/*--------------------------------------------------------------------------*/
/*	The separation is in units of AU, and both masses are in units of solar */
/*	masses.	 The period returned is in terms of Earth days.					*/
/*--------------------------------------------------------------------------*/

long double period(separation, small_mass, large_mass)
long double separation, small_mass, large_mass;
{
	long double period_in_years; 
	
	period_in_years = sqrt(pow3(separation) / (small_mass + large_mass));
	return(period_in_years * DAYS_IN_A_YEAR);
}


/*--------------------------------------------------------------------------*/
/*	 Fogg's information for this routine came from Dole "Habitable Planets	*/
/* for Man", Blaisdell Publishing Company, NY, 1964.  From this, he came	*/
/* up with his eq.12, which is the equation for the 'base_angular_velocity' */
/* below.  He then used an equation for the change in angular velocity per	*/
/* time (dw/dt) from P. Goldreich and S. Soter's paper "Q in the Solar		*/
/* System" in Icarus, vol 5, pp.375-389 (1966).	 Using as a comparison the	*/
/* change in angular velocity for the Earth, Fogg has come up with an		*/
/* approximation for our new planet (his eq.13) and take that into account. */
/* This is used to find 'change_in_angular_velocity' below.					*/
/*																			*/
/*	 Input parameters are mass (in solar masses), radius (in Km), orbital	*/
/* period (in days), orbital radius (in AU), density (in g/cc),				*/
/* eccentricity, and whether it is a gas giant or not.						*/
/*	 The length of the day is returned in units of hours.					*/
/*--------------------------------------------------------------------------*/

long double day_length(planet_pointer	planet)
{
	long double planetary_mass_in_grams = planet->mass * SOLAR_MASS_IN_GRAMS;
	long double	equatorial_radius_in_cm = planet->radius * CM_PER_KM;
	long double	year_in_hours			= planet->orb_period * 24.0;
	int giant = (planet->type == tGasGiant ||
				 planet->type == tSubGasGiant || 
				 planet->type == tSubSubGasGiant);
	long double	k2;
	long double	base_angular_velocity;
	long double	change_in_angular_velocity;
	long double	ang_velocity;
	long double	spin_resonance_factor;
	long double	day_in_hours;

	int stopped = FALSE;

	planet->resonant_period = FALSE;	/* Warning: Modify the planet */

	if (giant)
		k2 = 0.24;
	else 
		k2 = 0.33;

	base_angular_velocity = sqrt(2.0 * J * (planetary_mass_in_grams) /
								 (k2 * pow2(equatorial_radius_in_cm)));

/*	This next calculation determines how much the planet's rotation is	 */
/*	slowed by the presence of the star.								 */

	change_in_angular_velocity = CHANGE_IN_EARTH_ANG_VEL *
								 (planet->density / EARTH_DENSITY) *
								 (equatorial_radius_in_cm / EARTH_RADIUS) *
								 (EARTH_MASS_IN_GRAMS / planetary_mass_in_grams) *
								 pow(planet->sun->mass, 2.0) *
								 (1.0 / pow(planet->a, 6.0));
	ang_velocity = base_angular_velocity + (change_in_angular_velocity * 
											planet->sun->age);

/* Now we change from rad/sec to hours/rotation.						 */

	if (ang_velocity <= 0.0)
	{
	   stopped = TRUE;
	   day_in_hours = INCREDIBLY_LARGE_NUMBER ;
	}
	else 
		day_in_hours = RADIANS_PER_ROTATION / (SECONDS_PER_HOUR * ang_velocity);

	if ((day_in_hours >= year_in_hours) || stopped)
	{
		if (planet->e > 0.1)
		{
		  spin_resonance_factor 	= (1.0 - planet->e) / (1.0 + planet->e);
		  planet->resonant_period 	= TRUE;
		  return(spin_resonance_factor * year_in_hours);
		}
		else 
		  return(year_in_hours);
	}

	return(day_in_hours);
}


/*--------------------------------------------------------------------------*/
/*	 The orbital radius is expected in units of Astronomical Units (AU).	*/
/*	 Inclination is returned in units of degrees.							*/
/*--------------------------------------------------------------------------*/

int inclination(orb_radius)
long double orb_radius; 
{
	int temp; 
	
	temp = (int)(pow(orb_radius,0.2) * about(EARTH_AXIAL_TILT,0.4));
	return(temp % 360);
}


/*--------------------------------------------------------------------------*/
/*	 This function implements the escape velocity calculation.	Note that	*/
/*	it appears that Fogg's eq.15 is incorrect.								*/
/*	The mass is in units of solar mass, the radius in kilometers, and the	*/
/*	velocity returned is in cm/sec.											*/
/*--------------------------------------------------------------------------*/

long double escape_vel(mass, radius)
long double mass, radius;
{
	long double mass_in_grams, radius_in_cm;
	
	mass_in_grams = mass * SOLAR_MASS_IN_GRAMS;
	radius_in_cm = radius * CM_PER_KM;
	return(sqrt(2.0 * GRAV_CONSTANT * mass_in_grams / radius_in_cm));
}


/*--------------------------------------------------------------------------*/
/*	This is Fogg's eq.16.  The molecular weight (usually assumed to be N2)	*/
/*	is used as the basis of the Root Mean Square (RMS) velocity of the		*/
/*	molecule or atom.  The velocity returned is in cm/sec.					*/
/*	Orbital radius is in A.U.(ie: in units of the earth's orbital radius).	*/
/*--------------------------------------------------------------------------*/

long double rms_vel(long double molecular_weight, long double exospheric_temp)
{
	return(sqrt((3.0 * MOLAR_GAS_CONST * exospheric_temp) / molecular_weight)
		   * CM_PER_METER);
}


/*--------------------------------------------------------------------------*/
/*	 This function returns the smallest molecular weight retained by the	*/
/*	body, which is useful for determining the atmosphere composition.		*/
/*	Mass is in units of solar masses, and equatorial radius is in units of	*/
/*	kilometers.																*/
/*--------------------------------------------------------------------------*/

long double molecule_limit(mass, equat_radius, exospheric_temp)
long double mass, equat_radius, exospheric_temp;
{
	long double esc_velocity = escape_vel(mass,equat_radius);
	
	return ((3.0 * MOLAR_GAS_CONST * exospheric_temp) /
			(pow2((esc_velocity/ GAS_RETENTION_THRESHOLD) / CM_PER_METER)));

}

/*--------------------------------------------------------------------------*/
/*	 This function calculates the surface acceleration of a planet.	 The	*/
/*	mass is in units of solar masses, the radius in terms of km, and the	*/
/*	acceleration is returned in units of cm/sec2.							*/
/*--------------------------------------------------------------------------*/

long double acceleration(mass, radius)
long double mass, radius;
{
	return(GRAV_CONSTANT * (mass * SOLAR_MASS_IN_GRAMS) /
					   pow2(radius * CM_PER_KM));
}


/*--------------------------------------------------------------------------*/
/*	 This function calculates the surface gravity of a planet.	The			*/
/*	acceleration is in units of cm/sec2, and the gravity is returned in		*/
/*	units of Earth gravities.												*/
/*--------------------------------------------------------------------------*/

long double gravity(acceleration)
long double acceleration; 
{
	return(acceleration / EARTH_ACCELERATION);
}

/*--------------------------------------------------------------------------*/
/*	This implements Fogg's eq.17.  The 'inventory' returned is unitless.	*/
/*--------------------------------------------------------------------------*/

long double vol_inventory(mass, escape_vel, rms_vel, stellar_mass, zone, 
					 greenhouse_effect, accreted_gas)
long double mass, escape_vel, rms_vel, stellar_mass;
int zone, greenhouse_effect, accreted_gas;
{
	long double velocity_ratio, proportion_const, temp1, temp2, earth_units;
	
	velocity_ratio = escape_vel / rms_vel;
	if (velocity_ratio >= GAS_RETENTION_THRESHOLD)
	{
		switch (zone) {
			case 1:
				proportion_const = 140000.0;	/* 100 -> 140 JLB */
				break;
			case 2:
				proportion_const = 75000.0;
				break;
			case 3:
				proportion_const = 250.0;
				break;
			default:
				proportion_const = 0.0;
				printf("Error: orbital zone not initialized correctly!\n");
				break;
		}
		earth_units = mass * SUN_MASS_IN_EARTH_MASSES;
		temp1 = (proportion_const * earth_units) / stellar_mass;
		temp2 = about(temp1,0.2);
		temp2 = temp1;
		if (greenhouse_effect || accreted_gas)
			return(temp2);
		else 
			return(temp2 / 140.0);	/* 100 -> 140 JLB */
	}
	else 
		return(0.0);
}


/*--------------------------------------------------------------------------*/
/*	This implements Fogg's eq.18.  The pressure returned is in units of		*/
/*	millibars (mb).	 The gravity is in units of Earth gravities, the radius */
/*	in units of kilometers.													*/
/*																			*/
/*  JLB: Aparently this assumed that earth pressure = 1000mb. I've added a	*/
/*	fudge factor (EARTH_SURF_PRES_IN_MILLIBARS / 1000.) to correct for that	*/
/*--------------------------------------------------------------------------*/

long double pressure(volatile_gas_inventory, equat_radius, gravity)
long double volatile_gas_inventory, equat_radius, gravity;
{
	equat_radius = KM_EARTH_RADIUS / equat_radius;
	return(volatile_gas_inventory * gravity * 
			(EARTH_SURF_PRES_IN_MILLIBARS / 1000.) / 
			pow2(equat_radius));
}

/*--------------------------------------------------------------------------*/
/*	 This function returns the boiling point of water in an atmosphere of	*/
/*	 pressure 'surf_pressure', given in millibars.	The boiling point is	*/
/*	 returned in units of Kelvin.  This is Fogg's eq.21.					*/
/*--------------------------------------------------------------------------*/

long double boiling_point(surf_pressure)
long double surf_pressure; 
{
	long double surface_pressure_in_bars; 
	
	surface_pressure_in_bars = surf_pressure / MILLIBARS_PER_BAR;
	return (1.0 / ((log(surface_pressure_in_bars) / -5050.5) + 
				   (1.0 / 373.0) ));
	
}


/*--------------------------------------------------------------------------*/
/*	 This function is Fogg's eq.22.	 Given the volatile gas inventory and	*/
/*	 planetary radius of a planet (in Km), this function returns the		*/
/*	 fraction of the planet covered with water.								*/
/*	 I have changed the function very slightly:	 the fraction of Earth's	*/
/*	 surface covered by water is 71%, not 75% as Fogg used.					*/
/*--------------------------------------------------------------------------*/

long double hydro_fraction(volatile_gas_inventory, planet_radius)
long double volatile_gas_inventory, planet_radius;
{
	long double temp; 
	
	temp = (0.71 * volatile_gas_inventory / 1000.0)
			 * pow2(KM_EARTH_RADIUS / planet_radius);
	if (temp >= 1.0)
		return(1.0);
	else 
		return(temp);
}


/*--------------------------------------------------------------------------*/
/*	 Given the surface temperature of a planet (in Kelvin), this function	*/
/*	 returns the fraction of cloud cover available.	 This is Fogg's eq.23.	*/
/*	 See Hart in "Icarus" (vol 33, pp23 - 39, 1978) for an explanation.		*/
/*	 This equation is Hart's eq.3.											*/
/*	 I have modified it slightly using constants and relationships from		*/
/*	 Glass's book "Introduction to Planetary Geology", p.46.				*/
/*	 The 'CLOUD_COVERAGE_FACTOR' is the amount of surface area on Earth		*/
/*	 covered by one Kg. of cloud.											*/
/*--------------------------------------------------------------------------*/

long double cloud_fraction(surf_temp, smallest_MW_retained, equat_radius, hydro_fraction)
long double surf_temp, smallest_MW_retained, equat_radius,
	hydro_fraction;
{
	long double water_vapor_in_kg, fraction, surf_area, hydro_mass;
	
	if (smallest_MW_retained > WATER_VAPOR)
		return(0.0);
	else 
	{
		surf_area = 4.0 * PI * pow2(equat_radius);
		hydro_mass = hydro_fraction * surf_area * EARTH_WATER_MASS_PER_AREA;
		water_vapor_in_kg = (0.00000001 * hydro_mass) * 
							exp(Q2_36 * (surf_temp - EARTH_AVERAGE_KELVIN));
		fraction = CLOUD_COVERAGE_FACTOR * water_vapor_in_kg / surf_area;
		if (fraction >= 1.0)
			return(1.0);
		else 
			return(fraction);
	}
}


/*--------------------------------------------------------------------------*/
/*	 Given the surface temperature of a planet (in Kelvin), this function	*/
/*	 returns the fraction of the planet's surface covered by ice.  This is	*/
/*	 Fogg's eq.24.	See Hart[24] in Icarus vol.33, p.28 for an explanation. */
/*	 I have changed a constant from 70 to 90 in order to bring it more in	*/
/*	 line with the fraction of the Earth's surface covered with ice, which	*/
/*	 is approximatly .016 (=1.6%).											*/
/*--------------------------------------------------------------------------*/

long double ice_fraction(hydro_fraction, surf_temp)
long double hydro_fraction, surf_temp;
{
	long double temp; 
	
	if (surf_temp > 328.0) 
		surf_temp = 328.0;
	temp = pow(((328.0 - surf_temp) / 90.0), 5.0);
	if (temp > (1.5 * hydro_fraction))
		temp = (1.5 * hydro_fraction);
	if (temp >= 1.0)
		return(1.0);
	else 
		return(temp);
}


/*--------------------------------------------------------------------------*/
/*	This is Fogg's eq.19.  The ecosphere radius is given in AU, the orbital */
/*	radius in AU, and the temperature returned is in Kelvin.				*/
/*--------------------------------------------------------------------------*/

long double eff_temp(ecosphere_radius, orb_radius, albedo)
long double ecosphere_radius, orb_radius, albedo;
{
	return(sqrt(ecosphere_radius / orb_radius)
		  * pow1_4((1.0 - albedo) / (1.0 - EARTH_ALBEDO))
		  * EARTH_EFFECTIVE_TEMP);
}


long double est_temp(ecosphere_radius, orb_radius, albedo)
long double ecosphere_radius, orb_radius, albedo;
{
	return(sqrt(ecosphere_radius / orb_radius)
		  * pow1_4((1.0 - albedo) / (1.0 - EARTH_ALBEDO))
		  * EARTH_AVERAGE_KELVIN);
}


/*--------------------------------------------------------------------------*/
/* Old grnhouse:                                                            */
/*	Note that if the orbital radius of the planet is greater than or equal	*/
/*	to R_inner, 99% of it's volatiles are assumed to have been deposited in */
/*	surface reservoirs (otherwise, it suffers from the greenhouse effect).	*/
/*--------------------------------------------------------------------------*/
/*	if ((orb_radius < r_greenhouse) && (zone == 1)) */

/*--------------------------------------------------------------------------*/
/*	The new definition is based on the inital surface temperature and what	*/
/*	state water is in. If it's too hot, the water will never condense out	*/
/*	of the atmosphere, rain down and form an ocean. The albedo used here	*/
/*	was chosen so that the boundary is about the same as the old method		*/
/*	Neither zone, nor r_greenhouse are used in this version				JLB	*/
/*--------------------------------------------------------------------------*/

int grnhouse(long double r_ecosphere, long double orb_radius)
{
	long double	temp = eff_temp(r_ecosphere, orb_radius, GREENHOUSE_TRIGGER_ALBEDO);
	
	if (temp > FREEZING_POINT_OF_WATER)
		return(TRUE);
	else 
		return(FALSE);
}


/*--------------------------------------------------------------------------*/
/*	This is Fogg's eq.20, and is also Hart's eq.20 in his "Evolution of		*/
/*	Earth's Atmosphere" article.  The effective temperature given is in		*/
/*	units of Kelvin, as is the rise in temperature produced by the			*/
/*	greenhouse effect, which is returned.									*/
/*	I tuned this by changing a pow(x,.25) to pow(x,.4) to match Venus - JLB	*/
/*--------------------------------------------------------------------------*/

long double green_rise(optical_depth, effective_temp, surf_pressure)
long double optical_depth, effective_temp, surf_pressure;
{
	long double convection_factor = EARTH_CONVECTION_FACTOR * 
									pow(surf_pressure / 
										EARTH_SURF_PRES_IN_MILLIBARS, 0.4);
	long double rise = (pow1_4(1.0 + 0.75 * optical_depth) - 1.0) * 
					   effective_temp * convection_factor;
	
	if (rise < 0.0) rise = 0.0;
	
	return rise;
}


/*--------------------------------------------------------------------------*/
/*	 The surface temperature passed in is in units of Kelvin.				*/
/*	 The cloud adjustment is the fraction of cloud cover obscuring each		*/
/*	 of the three major components of albedo that lie below the clouds.		*/
/*--------------------------------------------------------------------------*/

long double planet_albedo(water_fraction, cloud_fraction, ice_fraction, surf_pressure)
long double water_fraction, cloud_fraction, ice_fraction, surf_pressure;
{
	long double rock_fraction, cloud_adjustment, components, cloud_part,
	rock_part, water_part, ice_part;
	
	rock_fraction = 1.0 - water_fraction - ice_fraction;
	components = 0.0;
	if (water_fraction > 0.0)
		components = components + 1.0;
	if (ice_fraction > 0.0)
		components = components + 1.0;
	if (rock_fraction > 0.0)
		components = components + 1.0;
	
	cloud_adjustment = cloud_fraction / components;
	
	if (rock_fraction >= cloud_adjustment)
		rock_fraction = rock_fraction - cloud_adjustment;
	else 
		rock_fraction = 0.0;
	
	if (water_fraction > cloud_adjustment)
		water_fraction = water_fraction - cloud_adjustment;
	else 
		water_fraction = 0.0;
		
	if (ice_fraction > cloud_adjustment)
		ice_fraction = ice_fraction - cloud_adjustment;
	else 
		ice_fraction = 0.0;
		
	cloud_part = cloud_fraction * CLOUD_ALBEDO;		/* about(...,0.2); */
	
	if (surf_pressure == 0.0)
	{
		rock_part = rock_fraction * ROCKY_AIRLESS_ALBEDO;	/* about(...,0.3); */
		ice_part = ice_fraction * AIRLESS_ICE_ALBEDO;		/* about(...,0.4); */
		water_part = 0;
	}
	else 
	{
		rock_part = rock_fraction * ROCKY_ALBEDO;	/* about(...,0.1); */
		water_part = water_fraction * WATER_ALBEDO;	/* about(...,0.2); */
		ice_part = ice_fraction * ICE_ALBEDO;		/* about(...,0.1); */
	}

	return(cloud_part + rock_part + water_part + ice_part);
}


/*--------------------------------------------------------------------------*/
/*	 This function returns the dimensionless quantity of optical depth,		*/
/*	 which is useful in determining the amount of greenhouse effect on a	*/
/*	 planet.																*/
/*--------------------------------------------------------------------------*/

long double opacity(molecular_weight, surf_pressure)
long double molecular_weight, surf_pressure;
{
	long double optical_depth; 
	
	optical_depth = 0.0;
	if ((molecular_weight >= 0.0) && (molecular_weight < 10.0))
		optical_depth = optical_depth + 3.0;
	if ((molecular_weight >= 10.0) && (molecular_weight < 20.0))
		optical_depth = optical_depth + 2.34;
	if ((molecular_weight >= 20.0) && (molecular_weight < 30.0))
		optical_depth = optical_depth + 1.0;
	if ((molecular_weight >= 30.0) && (molecular_weight < 45.0))
		optical_depth = optical_depth + 0.15;
	if ((molecular_weight >= 45.0) && (molecular_weight < 100.0))
		optical_depth = optical_depth + 0.05;

	if (surf_pressure >= (70.0 * EARTH_SURF_PRES_IN_MILLIBARS))
		optical_depth = optical_depth * 8.333;
	else 
		if (surf_pressure >= (50.0 * EARTH_SURF_PRES_IN_MILLIBARS))
			optical_depth = optical_depth * 6.666;
		else 
			if (surf_pressure >= (30.0 * EARTH_SURF_PRES_IN_MILLIBARS))
				optical_depth = optical_depth * 3.333;
			else 
				if (surf_pressure >= (10.0 * EARTH_SURF_PRES_IN_MILLIBARS))
					optical_depth = optical_depth * 2.0;
				else 
					if (surf_pressure >= (5.0 * EARTH_SURF_PRES_IN_MILLIBARS))
						optical_depth = optical_depth * 1.5;

	return(optical_depth);
}


/*
 *	calculates the number of years it takes for 1/e of a gas to escape
 *	from a planet's atmosphere. 
 *	Taken from Dole p. 34. He cites Jeans (1916) & Jones (1923)
 */
long double gas_life(long double molecular_weight, 
					 planet_pointer planet)
{
	long double v = rms_vel(molecular_weight, planet->exospheric_temp);
	long double g = planet->surf_grav * EARTH_ACCELERATION;
	long double r = (planet->radius * CM_PER_KM);
	long double t = (pow3(v) / (2.0 * pow2(g) * r)) * exp((3.0 * g * r) / pow2(v));
	long double years = t / (SECONDS_PER_HOUR * 24.0 * DAYS_IN_A_YEAR);
	
//	long double ve = planet->esc_velocity;
//	long double k = 2;
//	long double t2 = ((k * pow3(v) * r) / pow4(ve)) * exp((3.0 * pow2(ve)) / (2.0 * pow2(v)));
//	long double years2 = t2 / (SECONDS_PER_HOUR * 24.0 * DAYS_IN_A_YEAR);
		
//	if (flag_verbose & 0x0040)
//		fprintf (stderr, "gas_life: %LGs, V ratio: %Lf\n", 
//				years, ve / v);

	if (years > 2.0E10)
		years = INCREDIBLY_LARGE_NUMBER;
		
	return years;
}

long double min_molec_weight (planet_pointer planet)
{
	long double mass    = planet->mass;
	long double radius  = planet->radius;
	long double temp    = planet->exospheric_temp;
	long double target  = 5.0E9;
	
	long double guess_1 = molecule_limit (mass, radius, temp);
	long double	guess_2 = guess_1;
	
	long double life = gas_life(guess_1, planet);
	
	int	loops = 0;
	
	if (NULL != planet->sun)
	{
		target = planet->sun->age;
	}

	if (life > target)
	{
		while ((life > target) && (loops++ < 25))
		{
			guess_1 = guess_1 / 2.0;
			life 	= gas_life(guess_1, planet);
		}
	}
	else
	{
		while ((life < target) && (loops++ < 25))
		{
			guess_2 = guess_2 * 2.0;
			life 	= gas_life(guess_2, planet);
		}
	}

	loops = 0;

	while (((guess_2 - guess_1) > 0.1) && (loops++ < 25))
	{
		long double guess_3 = (guess_1 + guess_2) / 2.0;
		life			 	= gas_life(guess_3, planet);

		if (life < target)
			guess_1 = guess_3;
		else
			guess_2 = guess_3;
	}
	
	life = gas_life(guess_2, planet);

	return (guess_2);
}


/*--------------------------------------------------------------------------*/
/*	 The temperature calculated is in degrees Kelvin.						*/
/*	 Quantities already known which are used in these calculations:			*/
/*		 planet->molec_weight												*/
/*		 planet->surf_pressure												*/
/*		 R_ecosphere														*/
/*		 planet->a															*/
/*		 planet->volatile_gas_inventory										*/
/*		 planet->radius														*/
/*		 planet->boil_point													*/
/*--------------------------------------------------------------------------*/

void calculate_surface_temp(planet_pointer 	planet,
							int				first, 
							long double		last_water,
							long double 	last_clouds,
							long double 	last_ice,
							long double 	last_temp,
							long double 	last_albedo)
{
	long double effective_temp;
	long double water_raw;
	long double clouds_raw;
	long double greenhouse_temp;
	int			boil_off = FALSE;

	if (first)
	{
		planet->albedo = EARTH_ALBEDO;
	
		effective_temp 		= eff_temp(planet->sun->r_ecosphere, planet->a, planet->albedo);
		greenhouse_temp     = green_rise(opacity(planet->molec_weight, 
												 planet->surf_pressure), 
										 effective_temp, 
										 planet->surf_pressure);
		planet->surf_temp   = effective_temp + greenhouse_temp;

		set_temp_range(planet);
	}
	
	if (planet->greenhouse_effect
	 && planet->max_temp < planet->boil_point)
	{
		if (flag_verbose & 0x0010)
			fprintf (stderr, "Deluge: %s %d max (%Lf) < boil (%Lf)\n",
					planet->sun->name,
					planet->planet_no,
					planet->max_temp,
					planet->boil_point);

		planet->greenhouse_effect = 0;
		
		planet->volatile_gas_inventory 	= vol_inventory(planet->mass,
													    planet->esc_velocity,
													    planet->rms_velocity,
													    planet->sun->mass,
													    planet->orbit_zone,
													    planet->greenhouse_effect,
													    (planet->gas_mass
													     / planet->mass) > 0.000001);
		planet->surf_pressure 			= pressure(planet->volatile_gas_inventory,
												   planet->radius,
												   planet->surf_grav);

		planet->boil_point 			= boiling_point(planet->surf_pressure);
	}	

	water_raw     			=
	planet->hydrosphere		= hydro_fraction(planet->volatile_gas_inventory, 
											 planet->radius);
	clouds_raw     			=
	planet->cloud_cover 	= cloud_fraction(planet->surf_temp, 
											 planet->molec_weight, 
											 planet->radius, 
											 planet->hydrosphere);
	planet->ice_cover   	= ice_fraction(planet->hydrosphere, 
										   planet->surf_temp);
	
	if ((planet->greenhouse_effect)
	 && (planet->surf_pressure > 0.0))
		planet->cloud_cover	= 1.0;
	
	if ((planet->high_temp >= planet->boil_point)
	 && (!first)
	 && !((int)planet->day == (int)(planet->orb_period * 24.0) ||
		  (planet->resonant_period)))
	{
		planet->hydrosphere	= 0.0;
		boil_off = TRUE;
		
		if (planet->molec_weight > WATER_VAPOR)
			planet->cloud_cover = 0.0;
		else
			planet->cloud_cover = 1.0;
	}

	if (planet->surf_temp < (FREEZING_POINT_OF_WATER - 3.0))
		planet->hydrosphere	= 0.0;
	
	planet->albedo			= planet_albedo(planet->hydrosphere, 
											planet->cloud_cover, 
											planet->ice_cover, 
											planet->surf_pressure);
	
	effective_temp 			= eff_temp(planet->sun->r_ecosphere, planet->a, planet->albedo);
	greenhouse_temp     	= green_rise(opacity(planet->molec_weight,
												 planet->surf_pressure), 
										 effective_temp, 
										 planet->surf_pressure);
	planet->surf_temp   	= effective_temp + greenhouse_temp;

	if (!first)
	{
		if (!boil_off)
			planet->hydrosphere	= (planet->hydrosphere + (last_water * 2))  / 3;
		planet->cloud_cover	    = (planet->cloud_cover + (last_clouds * 2)) / 3;
		planet->ice_cover	    = (planet->ice_cover   + (last_ice * 2))    / 3;
		planet->albedo		    = (planet->albedo      + (last_albedo * 2)) / 3;
		planet->surf_temp	    = (planet->surf_temp   + (last_temp * 2))   / 3;
	}

	set_temp_range(planet);

	if (flag_verbose & 0x0020)
		fprintf (stderr, "%5.1Lf AU: %5.1Lf = %5.1Lf ef + %5.1Lf gh%c "
				"(W: %4.2Lf (%4.2Lf) C: %4.2Lf (%4.2Lf) I: %4.2Lf A: (%4.2Lf))\n", 
				planet->a,
				planet->surf_temp - FREEZING_POINT_OF_WATER,
				effective_temp - FREEZING_POINT_OF_WATER,
				greenhouse_temp,
				(planet->greenhouse_effect) ? '*' :' ',
				planet->hydrosphere, water_raw,
				planet->cloud_cover, clouds_raw,
				planet->ice_cover,
				planet->albedo);
}

void iterate_surface_temp(planet)
planet_pointer planet; 
{
	int			count = 0;
	long double initial_temp = est_temp(planet->sun->r_ecosphere, planet->a, planet->albedo);

	long double h2_life  = gas_life (MOL_HYDROGEN,    planet);
	long double h2o_life = gas_life (WATER_VAPOR,     planet);
	long double n2_life  = gas_life (MOL_NITROGEN,    planet);
	long double n_life   = gas_life (ATOMIC_NITROGEN, planet);
	
	if (flag_verbose & 0x20000)
		fprintf (stderr, "%d:                     %5.1Lf it [%5.1Lf re %5.1Lf a %5.1Lf alb]\n",
				planet->planet_no,
				initial_temp,
				planet->sun->r_ecosphere, planet->a, planet->albedo
				);

	if (flag_verbose & 0x0040)
		fprintf (stderr, "\nGas lifetimes: H2 - %Lf, H2O - %Lf, N - %Lf, N2 - %Lf\n",
				h2_life, h2o_life, n_life, n2_life);

	calculate_surface_temp(planet, TRUE, 0, 0, 0, 0, 0);

	for (count = 0;
		 count <= 25;
		 count++)
	{
		long double	last_water	= planet->hydrosphere;
		long double last_clouds	= planet->cloud_cover;
		long double last_ice	= planet->ice_cover;
		long double last_temp	= planet->surf_temp;
		long double last_albedo	= planet->albedo;
		
		calculate_surface_temp(planet, FALSE, 
							   last_water, last_clouds, last_ice, 
							   last_temp, last_albedo);
		
		if (fabs(planet->surf_temp - last_temp) < 0.25)
			break;
	}

	planet->greenhs_rise = planet->surf_temp - initial_temp;
	
	if (flag_verbose & 0x20000)
		fprintf (stderr, "%d: %5.1Lf gh = %5.1Lf (%5.1Lf C) st - %5.1Lf it [%5.1Lf re %5.1Lf a %5.1Lf alb]\n",
				planet->planet_no,
				planet->greenhs_rise,
				planet->surf_temp,
				planet->surf_temp - FREEZING_POINT_OF_WATER,
				initial_temp,
				planet->sun->r_ecosphere, planet->a, planet->albedo
				);
}

/*--------------------------------------------------------------------------*/
/*	 Inspired partial pressure, taking into account humidification of the	*/
/*	 air in the nasal passage and throat This formula is on Dole's p. 14	*/
/*--------------------------------------------------------------------------*/

long double inspired_partial_pressure (long double surf_pressure,
									   long double gas_pressure)
{
	long double pH2O = (H20_ASSUMED_PRESSURE);
	long double fraction = gas_pressure / surf_pressure;
	
	return	(surf_pressure - pH2O) * fraction;
}


/*--------------------------------------------------------------------------*/
/*	 This function uses figures on the maximum inspired partial pressures   */
/*   of Oxygen, other atmospheric and traces gases as laid out on pages 15, */
/*   16 and 18 of Dole's Habitable Planets for Man to derive breathability  */
/*   of the planet's atmosphere.                                       JLB  */
/*--------------------------------------------------------------------------*/

unsigned int breathability (planet_pointer planet)
{
	int	oxygen_ok	= FALSE;
	int index;

	if (planet->gases == 0)
		return NONE;
	
	for (index = 0; index < planet->gases; index++)
	{
		int	n;
		int	gas_no = 0;
		
		long double ipp = inspired_partial_pressure (planet->surf_pressure,
													 planet->atmosphere[index].surf_pressure);
		
		for (n = 0; n < max_gas; n++)
		{
			if (gases[n].num == planet->atmosphere[index].num)
				gas_no = n;
		}

		if (ipp > gases[gas_no].max_ipp)
			return POISONOUS;
			
		if (planet->atmosphere[index].num == AN_O)
			oxygen_ok = ((ipp >= MIN_O2_IPP) && (ipp <= MAX_O2_IPP));
	}
	
	if (oxygen_ok)
		return BREATHABLE;
	else
		return UNBREATHABLE;
}

/* function for 'soft limiting' temperatures */

long double lim(long double x)
{
  return x / sqrt(sqrt(1 + x*x*x*x));
}

long double soft(long double v, long double max, long double min)
{
  long double dv = v - min;
  long double dm = max - min;
  return (lim(2*dv/dm-1)+1)/2 * dm + min;
}

void set_temp_range(planet_pointer planet)
{
  long double pressmod = 1 / sqrt(1 + 20 * planet->surf_pressure/1000.0);
  long double ppmod    = 1 / sqrt(10 + 5 * planet->surf_pressure/1000.0);
  long double tiltmod  = fabs(cos(planet->axial_tilt * PI/180) * pow(1 + planet->e, 2));
  long double daymod   = 1 / (200/planet->day + 1);
  long double mh = pow(1 + daymod, pressmod);
  long double ml = pow(1 - daymod, pressmod);
  long double hi = mh * planet->surf_temp;
  long double lo = ml * planet->surf_temp;
  long double sh = hi + pow((100+hi) * tiltmod, sqrt(ppmod));
  long double wl = lo - pow((150+lo) * tiltmod, sqrt(ppmod));
  long double max = planet->surf_temp + sqrt(planet->surf_temp) * 10;
  long double min = planet->surf_temp / sqrt(planet->day + 24);

  if (lo < min) lo = min;
  if (wl < 0)   wl = 0;

  planet->high_temp = soft(hi, max, min);
  planet->low_temp  = soft(lo, max, min);
  planet->max_temp  = soft(sh, max, min);
  planet->min_temp  = soft(wl, max, min);
}

