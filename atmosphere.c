#include <windows.h>
#include <math.h>
#include <string.h>
#include "usersurf.h"



// Constants from the original Python code
const double k1 = 1167.05214528;
const double k2 = -724213.167032;
const double k3 = -17.0738469401;
const double k4 = 12020.8247025;
const double k5 = -3232555.03223;
const double k6 = 14.9151086135;
const double k7 = -4823.26573616;
const double k8 = 405113.405421;
const double k9 = -0.238555575678;
const double k10 = 650.175348448;


// Function implementations
double omega(double t);
double A(double t);
double B(double t);
double C(double t);
double X(double t);
double y(double t);
double f(double t, double RH);
double nvaporf(double x, double Tf, double RH);
double nseaf(double x, double Tf, double RH);
double naltitudef(double x, double Tf, double RH, double Pf);
double refractive_index(double x, double Tf, double RH, double Pf, double zf) ;
double dispersion_arcsec(double x, double xref, double Tf, double RH, double Pf, double zf) ;



int __declspec(dllexport) APIENTRY UserDefinedSurface3(USER_DATA *UD, FIXED_DATA3 *FD);

/* a generic Snells law refraction routine */
int Refract(double thisn, double nextn, double *l, double *m, double *n, double ln, double mn, double nn);

BOOL WINAPI DllMain (HANDLE hInst, ULONG ul_reason_for_call, LPVOID lpReserved)
	{
   return TRUE;
   }

/* this DLL models a standard ZEMAX surface type, either plane, sphere, or conic */

int __declspec(dllexport) APIENTRY UserDefinedSurface3(USER_DATA *UD, FIXED_DATA3 *FD)
	{
   int i;
   double p2, alpha, power, a, b, c, rad, casp, t, zc;
   double zen_angle, height, Temp, Pres, RelH, latitude, wav,n_atm,wavp;
   double nx,ny,nz,ux,uy,uz;
   double delta,x,y,z,dpdx,dpdy;
   switch(FD->type)
   	{
      case 0:
         switch(FD->numb)
         	{
            case 0:
            	/* ZEMAX wants to know the name of the surface type; do not exceed 12 characters */
		         strcpy(UD->string,"Filippenko_Atmosphere");
               break;
            case 1:
            	/* ZEMAX wants to know if this surface is rotationally symmetric;  it is, so return any character in the string; otherwise, return a null string */
            	strcpy(UD->string, "1");
               break;
            case 2:
            	/* ZEMAX wants to know if this surface is a gradient index media */
               /* it is not, so return a null string */
            	UD->string[0] = '\0';
            	break;
            }
         break;
      case 1:
      	/* ZEMAX is requesting the names of the parameter columns */
         /* the value FD->numb will indicate which value ZEMAX wants. */
         /* they are all "Unused" for this surface type */
         /* returning a null string indicates that the parameter is unused. */
         switch(FD->numb)
         	{
            case 1:
                strcpy(UD->string, "Zenith");
               break;
            case 2:
                strcpy(UD->string, "Temprature");
                break;
            case 3:
                strcpy(UD->string, "Pressure");
                break;
            case 4:
                strcpy(UD->string, "Humidity");
                break;
            default:
            	UD->string[0] = '\0';
            	break;
            }
      	break;
      case 2:
      	/* ZEMAX is requesting the names of the extra data columns */
         /* the value FD->numb will indicate which value ZEMAX wants. */
         /* they are all "Unused" for this surface type */
         /* returning a null string indicates that the extradata value is unused. */
         switch(FD->numb)
         	{
            default:
            	UD->string[0] = '\0';
            	break;
            }
      	break;
      case 3:
      	/* ZEMAX wants to know the sag of the surface */
         /* if there is an alternate sag, return it as well */
         /* otherwise, set the alternate sag identical to the sag */
         /* The sag is sag1, alternate is sag2. */

         UD->sag1 = 0.0;
         UD->sag2 = 0.0;

			/* if a plane, just return */
			if (FD->cv == 0) return(0);
         p2 = UD->x * UD->x + UD->y * UD->y;
         alpha = 1 - (1+FD->k)*FD->cv*FD->cv*p2;

		 // If the absolute value of alpha is smaller than 1e-13, we're going to assume this instead means zero
		 // This assumption is based on the fact that floating point numbers on computers cannot be
		 // represented 100% accurately
		 if (fabs(alpha) < 1e-13)
			 alpha = 0;

         if (alpha < 0) return(-1);
         UD->sag1 = (FD->cv*p2)/(1 + sqrt(alpha));
         if (alpha != 1.0) UD->sag2 = (FD->cv*p2)/(1 - sqrt(alpha));
      	break;
      case 4:
      	/* ZEMAX wants a paraxial ray trace to this surface */
         /* x, y, z, and the path are unaffected, at least for this surface type */
         /* for paraxial ray tracing, the return z coordinate should always be zero. */
         /* paraxial surfaces are always planes with the following normals */
         UD->ln =  0.0;
         UD->mn =  0.0;
         UD->nn = -1.0;
         power = (FD->n2 - FD->n1)*FD->cv;
         if ((UD->n) != 0.0)
         	{
            (UD->l) = (UD->l)/(UD->n);
            (UD->m) = (UD->m)/(UD->n);

            (UD->l) = (FD->n1*(UD->l) - (UD->x)*power)/(FD->n2);
            (UD->m) = (FD->n1*(UD->m) - (UD->y)*power)/(FD->n2);

            /* normalize */
            (UD->n) = sqrt(1/(1 + (UD->l)*(UD->l) + (UD->m)*(UD->m) ) );
            /* de-paraxialize */
            (UD->l) = (UD->l)*(UD->n);
            (UD->m) = (UD->m)*(UD->n);
            }
         break;
      case 5:
        /* ZEMAX wants a real ray trace to this surface */
         if (FD->cv == 0.0)
         	{
	         UD->ln =  0.0;
   	      UD->mn =  0.0;
      	   UD->nn = -1.0;
			   if (Refract(FD->n1, FD->n2, &UD->l, &UD->m, &UD->n, UD->ln, UD->mn, UD->nn)) return(-FD->surf);
            }
         else
         	{
	         /* okay, not a plane. */
				a = (UD->n) * (UD->n) * FD->k + 1;
				b = ((UD->n)/FD->cv) - (UD->x) * (UD->l) - (UD->y) * (UD->m);
				c = (UD->x) * (UD->x) + (UD->y) * (UD->y);
				rad = b * b - a * c;
				if (rad < 0) return(FD->surf);  /* ray missed this surface */
				if (FD->cv > 0) t = c / (b + sqrt(rad));
				else           t = c / (b - sqrt(rad));
				(UD->x) = (UD->l) * t + (UD->x);
				(UD->y) = (UD->m) * t + (UD->y);
				(UD->z) = (UD->n) * t + (UD->z);
				UD->path = t;
				zc = (UD->z) * FD->cv;
				rad = zc * FD->k * (zc * (FD->k + 1) - 2) + 1;
				casp = FD->cv / sqrt(rad);
				UD->ln = (UD->x) * casp;
				UD->mn = (UD->y) * casp;
				UD->nn = ((UD->z) - ((1/FD->cv) - (UD->z) * FD->k)) * casp;
	         if (Refract(FD->n1, FD->n2, &UD->l, &UD->m, &UD->n, UD->ln, UD->mn, UD->nn)) return(-FD->surf);
            }
            /* Account for atmospheric dispersion */
            /* Atmospheric dispersion affects both the refraction angle and the optical path (OPD) */

            zen_angle = FD->param[1]; // Zenith angle in degrees   
            Temp = FD->param[2];      // Temperature in Kelvin
            Pres = FD->param[3];      // Pressure in hPa
            RelH = FD->param[4];      // Relative humidity in percentage
            wav = FD->wavelength;     // Wavelength in micrometers
            wavp = FD->pwavelength;

            /* Dispersion parameters */
            delta = dispersion_arcsec(wav, wavp, Temp, RelH, Pres, zen_angle); // Result in radians
            delta /= fabs(FD->n2); // Adjust for possible changes in refractive index
            if (FD->n1 * FD->n2 < 0) delta = -delta; // Correct sign for index mismatch

            /* Update optical path with dispersion */
            UD->path += delta * wav / FD->n1; // Include phase term in optical path

            /* Initialize direction cosines */
            nx = -UD->ln;
            ny = -UD->mn;
            nz = -UD->nn;

            /* Compute new ray direction cosines based on atmospheric effects */
            dpdx = 0.0;       // Placeholder for potential azimuthal dispersion
            dpdy = delta;     // Vertical dispersion due to atmospheric effects

            /* Adjust for direction cosines */
            ux = UD->l - dpdx;
            uy = UD->m - dpdy;
            uz = UD->n;

            /* Normalize the resulting direction cosines */
            rad = nx * ux + ny * uy + nz * uz;
            rad = 1.0 - (ux * ux + uy * uy + uz * uz) + rad * rad;
            if (rad <= 0.0) rad = 0.0;
            else rad = sqrt(rad);

            /* Update the unit direction cosines */
            UD->l = ux - (nx * ux + ny * uy + nz * uz) * nx + nx * rad;
            UD->m = uy - (nx * ux + ny * uy + nz * uz) * ny + ny * rad;
            UD->n = uz - (nx * ux + ny * uy + nz * uz) * nz + nz * rad;

            break;

      case 6: {

          UD->index = FD->n2;
          UD->dndx = 0.0;
          UD->dndy = 0.0;
          UD->dndz = 0.0;
          break;
      }

      case 7:
      	/* ZEMAX wants the "safe" data. */
         /* this is used by ZEMAX to set the initial values for all parameters and extra data */
         /* when the user first changes to this surface type. */
         /* this is the only time the DLL should modify the data in the FIXED_DATA FD structure */
         //param 0 is zenith angle
         FD->param[1]=60;
         //param 1 is height
         FD->param[2]=275;
         //param 2 is temperature
         FD->param[3]=620;
         //param 3 is pressure
         FD->param[4]=0.3;
         //param 4 is humidity
         break;
      case 8:
      	/* ZEMAX is calling the DLL for the first time, do any memory or data initialization here. */
         break;
      case 9:
      	/* ZEMAX is calling the DLL for the last time, do any memory release here. */
         break;
      }
   return 0;
   }


int Refract(double thisn, double nextn, double *l, double *m, double *n, double ln, double mn, double nn)
{
double nr, cosi, cosi2, rad, cosr, gamma;
if (thisn != nextn)
	{
	nr = thisn / nextn;
	cosi = fabs((*l) * ln + (*m) * mn + (*n) * nn);
	cosi2 = cosi * cosi;
	if (cosi2 > 1) cosi2 = 1;
	rad = 1 - ((1 - cosi2) * (nr * nr));
	if (rad < 0) return(-1);
	cosr = sqrt(rad);
	gamma = nr * cosi - cosr;
	(*l) = (nr * (*l)) + (gamma * ln);
	(*m) = (nr * (*m)) + (gamma * mn);
	(*n) = (nr * (*n)) + (gamma * nn);
	}
return 0;
}


// Function implementations

/**
 * Calculate the omega value for a given temperature.
 *
 * @param t  Temperature value.
 * @return   Computed omega value.
 */
double omega(double t) {
    // Calculate and return the omega value using the formula (t + k9) / (t - k10)
    return (t + k9) / (t - k10);
}

/**
 * Calculate the A value for a given temperature.
 *
 * @param t  Temperature value.
 * @return   Computed A value.
 */
double A(double t) {
    // Calculate and return the A value using the formula w^2 + k1 * w + k2
    return pow(omega(t), 2) + k1 * omega(t) + k2;
}


/**
 * Calculate the B value for a given temperature.
 *
 * @param t  Temperature value.
 * @return   Computed B value.
 */
double B(double t) {
    // Calculate and return the B value using the formula k3 * w^2 + k4 * w + k5
    return k3 * pow(omega(t), 2) + k4 * omega(t) + k5;
}


/**
 * Calculate the C value for a given temperature.
 *
 * @param t  Temperature value.
 * @return   Computed C value.
 *
 * The C value is calculated using the formula k6 * w^2 + k7 * w + k8, where
 * w is the omega value calculated by the function omega(t).
 */
double C(double t) {
    return k6 * pow(omega(t), 2) + k7 * omega(t) + k8;
}


/**
 * Calculate the quantity X(t) = -B(t) + sqrt(B(t)^2 - 4 * A(t) * C(t)),
 * where A(t), B(t), and C(t) are defined in the paper by Edlen 1966.
 *
 * @param t     Temperature in Kelvin.
 *
 * @return      The quantity X(t).
 */
double X(double t) {
    return -B(t) + sqrt(pow(B(t), 2) - 4.0 * A(t) * C(t));
}

/**
 * Calculate the quantity y(t) = 10^6 * (2 * C(t) / X(t))^4, where C(t) and X(t)
 * are defined in the paper by Edlen 1966.
 *
 * @param t     Temperature in Kelvin.
 *
 * @return      The quantity y(t).
 */
double y(double t) {
    return pow(10, 6) * pow((2 * C(t) / X(t)), 4);
}

/**
 * Calculate the refractive index of water vapor at the given conditions.
 *
 * @param x      Wavelength in micrometers.
 * @param Tf     Temperature in Kelvin.
 * @param RH     Relative Humidity in percentage.
 *
 * @return       The refractive index of water vapor at the given conditions.
 *
 * This function returns the refractive index of water vapor at the given
 * conditions. The refractive index is calculated using the formula given in
 * the paper by Edlen 1966.
 */
double f(double t, double RH) {
    // Conversion factor from Pascal to Torr
    return (RH / 100.0) * y(t) * 0.00750062;
}


/**
 * Calculate the refractive index of water vapor at the given conditions.
 *
 * @param x      Wavelength in micrometers.
 * @param Tf     Temperature in Kelvin.
 * @param RH     Relative Humidity in percentage.
 *
 * @return       The refractive index of water vapor at the given conditions.
 */
double nvaporf(double x, double Tf, double RH) {
    /**
     * Formula for refractive index of water vapor from the 1995 paper by
     * Ciddor.  The formula is
     *
     * (0.0624 - 0.680 / lambda^2) / (1 + 0.003661 * (T - 273.15))
     *
     * where lambda is the wavelength in micrometers and T is the temperature
     * in Kelvin.
     */
    return ((0.0624 - (0.000680 / pow(x, 2))) / (1.0 + 0.003661 * (Tf - 273.15))) * f(Tf, RH);
}


/**
 * Calculate the refractive index of air at sea level.
 *
 * @param x      Wavelength in micrometers.
 * @param Tf     Temperature in Kelvin.
 * @param RH     Relative Humidity in percentage.
 *
 * @return       The refractive index of air at sea level.
 */
double nseaf(double x, double Tf, double RH) {
    // Calculate sea level refractive index minus 1.0
    return (((64.328 + (29498.1 / (146.0 - pow(1.0 / x, 2))) + (255.4 / (41.0 - pow(1.0 / x, 2))))
        - nvaporf(x, Tf, RH)) / pow(10, 6)) + 1.0;
}

/**
 * Calculate the refractive index of air at a given altitude.
 *
 * @param x   Current altitude (wavelength in micrometers).
 * @param Tf  Temperature in Kelvin.
 * @param RH  Relative Humidity in percentage.
 * @param Pf  Pressure in hPa.
 *
 * @return    The refractive index of air at the specified altitude.
 */
double naltitudef(double x, double Tf, double RH, double Pf) {
    // Calculate sea level refractive index minus 1.0
    double n_seaf_level = nseaf(x, Tf, RH) - 1.0;

    // Calculate pressure-related term with conversion factor
    double pressure_term = (Pf / 1.3332239) * 
                           (1.0 + (1.049 - 0.0157 * (Tf - 273.15)) * 
                           (Pf / 1.3332239) / pow(10, 6));

    // Calculate temperature-related term for dry air
    double temperature_term = 720.883 * (1.0 + 0.003661 * (Tf - 273.15));

    // Combine all terms to calculate and return the refractive index
    return (n_seaf_level * pressure_term / temperature_term) + 1.0;
}

/**
 * Calculate the dispersion in arcseconds.
 *
 * @param x      Wavelenth.
 * @param xref   reference wavelength.
 * @param Tf     Temperature in Kelvin.
 * @param RH     Relative Humidity in percentage.
 * @param Pf     Pressure in hPa.
 * @param zf     Zenith angle in degrees.
 * @return       Dispersion in arcseconds.
 */
double dispersion_arcsec(double x, double xref, double Tf, double RH, double Pf, double zf) {
    // Calculate the difference in refractive indices at the given and reference altitudes
    double delta_n = naltitudef(x, Tf, RH, Pf) - naltitudef(xref, Tf, RH, Pf);
    
    // Convert zenith angle from degrees to radians and calculate dispersion in arcseconds
    return delta_n * tan(zf * 3.14159 / 180.0);
}
