#include <stdio.h>
#include <stdlib.h>
#include "mathlink.h"


#define return_null()  {\
    MLPutSymbol(stdlink, "Null"); \
    MLEndPacket(stdlink); \
    MLFlush(stdlink); \
    return;\
}

/* Halofit+ */
/* `real` was defined as double in smith2.h */

#ifdef HALOFIT

#include "halofit/smith2.h"

void HFset_parameters(real OMEGAM, real OMEGAV, real GAMMA, real SIGMA8,
		   real NSPEC, real BETAP, real Z0, int NONLINEAR){
    setparameters_(&OMEGAM, &OMEGAV, &GAMMA, &SIGMA8, &NSPEC, &BETAP, &Z0, &NONLINEAR);
     return_null();
}

#else

#define real double
void HFset_parameters(real OMEGAM, real OMEGAV, real GAMMA, real SIGMA8,
		   real NSPEC, real BETAP, real Z0, int NONLINEAR){
    MLPutSymbol(stdlink, "Null");
    MLEndPacket(stdlink);
    MLFlush(stdlink);
}
extern real P_NL(real a, real k);
extern real Pkappa(real ell);

#endif /*HALOFIT*/


#ifdef TRANSFER

/* TF Transfer function from Eisenstein&hu */

/* high baryon version */
extern void TFset_parameters(float omega0hh, float f_baryon, float Tcmb);
extern float TFfit_onek(float k, float *tf_baryon, float *tf_cdm); 
extern float TFnowiggles(float omega0, float f_baryon, float hubble, 
                                float Tcmb, float k_hmpc);
extern float TFzerobaryon(float omega0, float hubble, float Tcmb, float k_hmpc);
extern float TFsound_horizon_fit(float omega0, float f_baryon, float hubble);
extern float TFk_peak(float omega0, float f_baryon, float hubble);

/* other version */
extern int TFmdm_set_cosm(float omega_matter, float omega_baryon, float omega_hdm,
	int degen_hdm, float omega_lambda, float hubble, float redshift);
extern float TFmdm_onek_hmpc(float kk);


void TFfit_onek_wrap(float k){
    float tf_baryon, tf_cdm, tf_full;
    tf_full = TFfit_onek(k, &tf_baryon, &tf_cdm);

        /* return tf_full, tf_baryon, tf_cdm */
    MLPutFunction(stdlink, "List", 3);
    MLPutReal32(stdlink, tf_full);
    MLPutReal32(stdlink, tf_baryon);
    MLPutReal32(stdlink, tf_cdm);
    MLEndPacket(stdlink);
    MLFlush(stdlink);
}

void TFset_parameters_wrap(float omega0hh, float f_baryon, float Tcmb){
    TFset_parameters(omega0hh, f_baryon, Tcmb);
     /* We need to return *something* */
    MLPutSymbol(stdlink, "Null");
    MLEndPacket(stdlink);
    MLFlush(stdlink);
}

#else
void TFfit_onek_wrap(float k){
    return_null();
}
void TFset_parameters_wrap(float omega0hh, float f_baryon, float Tcmb){
    return_null();
}
extern float TFnowiggles(float omega0, float f_baryon, float hubble, float Tcmb, float k_hmpc);
extern float TFzerobaryon(float omega0, float hubble, float Tcmb, float k_hmpc);
extern float TFsound_horizon_fit(float omega0, float f_baryon, float hubble);
extern float TFk_peak(float omega0, float f_baryon, float hubble);
extern int TFmdm_set_cosm(float omega_matter, float omega_baryon, float omega_hdm,
	int degen_hdm, float omega_lambda, float hubble, float redshift);
extern float TFmdm_onek_hmpc(float kk);

#endif /*TRANSFER*/


/* CosmicEmulator version 1.1 */

#ifdef COSMICEMU
extern void emu(double *xstar, double *ystar, int *outtype);
extern void getH0fromCMB(double *xstar, double *stuff);

void CEget_PkNL(double omegaM, double omegaB, double ns, double sigma8, double w, double z ){
    const int output_length = 2*1995;
    double input[6], output[output_length], more_output[4];
    int type=2; // Output: P(k)

    input[0] = omegaM;
    input[1] = omegaB;
    input[2] = ns;
    input[3] = sigma8;
    input[4] = w;
    input[5] = z;

    getH0fromCMB(input, more_output);
    emu(input, output, &type);


    MLPutFunction(stdlink, "List", 2);
    MLPutReal64List(stdlink, (double*)output, output_length);
    MLPutReal64List(stdlink, (double*)more_output, 4);
    MLEndPacket(stdlink);
    MLFlush(stdlink);
}

#else
void CEget_PkNL(double omegaM, double omegaB, double ns, double sigma8, double w, double z ){
    return_null();
}
#endif /*COSMICEMU*/


/* FrankenEmu */

#ifdef FRANKENEMU

extern void emu(double *xstar, double *ystar, int *outtype);
extern void emu_noh(double *xstar, double *ystar, int *outtype);
extern void getH0fromCMB(double *xstar, double *stuff);

void franken_CEget_PkNL(double omegaM, double omegaB, double h, double ns, double sigma8, double w, double z ){
    const int output_length = 2*582;
    double input[7], output[output_length], more_output[4];
    int type=2; // Output: P(k)


    if (h>0) {
        input[0] = omegaB;
        input[1] = omegaM;
        input[2] = ns;
        input[3] = h*100;
        input[4] = w;
        input[5] = sigma8;
        input[6] = z;
        emu(input, output, &type);
    } else {
        input[0] = omegaB;
        input[1] = omegaM;
        input[2] = ns;
        input[3] = w;
        input[4] = sigma8;
        input[5] = z;
        emu_noh(input, output, &type);
        getH0fromCMB(input, more_output);
    }


    MLPutFunction(stdlink, "List", 2);
    MLPutReal64List(stdlink, (double*)output, output_length);
    MLPutReal64List(stdlink, (double*)more_output, 4);
    MLEndPacket(stdlink);
    MLFlush(stdlink);
}

/* Dummy functions need to be defined, but will never be called anyway
 * because FrankenEmu is alone in this module */
double P_NL(double a, double k){return 0;}
double Pkappa(double ell){return 0;}
float TFnowiggles(float omega0, float f_baryon, float hubble, float Tcmb, float k_hmpc){return 0.;}
float TFzerobaryon(float omega0, float hubble, float Tcmb, float k_hmpc){return 0.;}
float TFsound_horizon_fit(float omega0, float f_baryon, float hubble){return 0.;}
float TFk_peak(float omega0, float f_baryon, float hubble){return 0.;}

int TFmdm_set_cosm(float omega_matter, float omega_baryon, float omega_hdm,
	int degen_hdm, float omega_lambda, float hubble, float redshift){return 0;}
float TFmdm_onek_hmpc(float kk){return 0.;}

#else

void franken_CEget_PkNL(double omegaM, double omegaB, double ns, double sigma8, double w, double z ){
    return_null();
}

#endif /*FrankenEmu*/



/* CAMB */

#ifdef CAMB
extern void runcamb();

void CAMBrun(double *floats, long floats_len, int *ints, long ints_len){
    /* TODO Do longs need to be converted to ints here? */

    const int floats_out_len = 2000000;
    double *floats_out = malloc(sizeof(*floats_out) * floats_out_len);
    const int ints_out_len = 100;
    double ints_out[ints_out_len];

    runcamb(floats, &floats_len, ints, &ints_len, (double*)floats_out, &floats_out_len, (int*)ints_out, &ints_out_len);

    MLPutFunction(stdlink, "List", 2);
    MLPutReal64List(stdlink, (double*)floats_out, floats_out_len);
    MLPutIntegerList(stdlink, (int*)ints_out, ints_out_len);
    MLEndPacket(stdlink);
    MLFlush(stdlink);

    free(floats_out);
}

#else

// void runcamb(){}
void CAMBrun(double *floats, long floats_len, int *ints, long ints_len){
    return_null();
}

#endif /*CAMB*/


#ifdef CLASS

extern int class_init(int, char**, char*, char*);
extern double* class_get(int **);
extern int class_free(char*, char*);
    
void MLClass(char *inifile){
    double *result;
    int *limits;
    int i, size = 0;
    char e1[2048], e2[2048];
    char (*argv[2])[1024];
    int fail;
    

    fail = 0;
    argv[0] = "class";
    argv[1] = inifile;
    printf("%s\n", inifile);

    
    if (class_init(2, (char **)argv, e1, e2)==0){
        result = class_get(&limits);
        for (i=1; i<1024; i+=2) size += limits[i];
        MLPutFunction(stdlink,"List", 2);
        MLPutIntegerList(stdlink, limits, 1024);
        MLPutReal64List(stdlink, result, size);
        MLEndPacket(stdlink);
        MLFlush(stdlink);
        free(result);
        free(limits);
    }
    else fail = 1;

    if (fail==1 || class_free(e1, e2)!=0) fail = 1;

    if (fail==1) {
        printf("%s\n%s\n", e1, e2);
        MLPutFunction(stdlink,"List", 2);
        MLPutString(stdlink, e1);
        MLPutString(stdlink, e2);
        MLEndPacket(stdlink);
        MLFlush(stdlink);
    }

}

#else

void MLClass(int argc, char *argv){
    return_null();
}

#endif /*CLASS*/

int main(int argc, char* argv[]) {
    return MLMain(argc, argv);
}
