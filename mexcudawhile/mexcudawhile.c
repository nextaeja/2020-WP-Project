#include <mex.h>
#include <matrix.h>
#include <math.h>
#include <cufft.h>
#include "cuda.h"

#include "../MEX_helpers/complex.h"
#include "../MEX_helpers/cuda_helper.h"
#include "../Setup/cuda_setup_dynamic_potential.h"

//#include "unistd.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int it;
float A=1E-10;
float ps=1E-12;
float eV=1.60217662E-19;
int flleng = 10; //max length of floats in specification files inc space
float dx;
float dy;
float dz;
float dt;
int nx;
int ny;
int nz;
float eps;
float tStart;
float tFinish;
int savingsimulationrunning;
int propogationmethod;
int decaytype;
float zOffset;
float alpha;
char potentialfile[256];
float xsigma;
float ysigma;
float gausspeakval;
float welldepth;
int numAds;
float* psisavets;
float* ts;
float*** Adspos;
float*custpot;
float**** psi;


void setup() {
	
	/*char cwdtmp[256];

	//getcwd(cwd, sizeof(cwd));
	char* cwd = "C:\\Users\\jackl\\Documents\\MATLAB\\summer 2020\\2020-WP-Project-6_optimizing_setup_dynamic_gaussian_potential\\Beta4_2\\mexcudawhile";//make general
	char* l= strrchr(cwd,'\\');
	cwd[strlen(cwd)-strlen(l)+1]='\0';

	strcpy(cwdtmp,cwd);

	FILE *setup=fopen(strcat(cwdtmp,"setup.txt"),"r");//navigates up to setup file


	char buff[1024];

	fgets(buff,sizeof(buff),setup);
	dx=strtof(buff,NULL)*A;//reads setup options

	

	fgets(buff,sizeof(buff),setup);
	dy=strtof(buff,NULL)*A;

	fgets(buff,sizeof(buff),setup);
	dz=strtof(buff,NULL)*A;

	fgets(buff,sizeof(buff),setup);
	dt=strtof(buff,NULL)*ps;

	fgets(buff,sizeof(buff),setup);
	nx=atoi(buff);

	fgets(buff,sizeof(buff),setup);
	ny=atoi(buff);

	fgets(buff,sizeof(buff),setup);
	nz=atoi(buff);

	fgets(buff,sizeof(buff),setup);
	eps=strtof(buff,NULL);

	fgets(buff,sizeof(buff),setup);
	tStart=strtof(buff,NULL)*ps;

	fgets(buff,sizeof(buff),setup);
	tFinish=strtof(buff,NULL)*ps;

	fgets(buff,sizeof(buff),setup);
	savingsimulationrunning=atoi(buff);

	fgets(buff,sizeof(buff),setup);
	propogationmethod =atoi(buff);

	fgets(buff,sizeof(buff),setup);
	decaytype=atoi(buff);

	fgets(buff,sizeof(buff),setup);
	zOffset=strtof(buff,NULL)*A;

	fgets(buff,sizeof(buff),setup);
	alpha=strtof(buff,NULL);

	fgets(buff,sizeof(buff),setup);
	strcpy(potentialfile,buff);
	potentialfile[strlen(potentialfile)-1]='\0';

	fgets(buff,sizeof(buff),setup);
	xsigma=strtof(buff,NULL)*A;

	fgets(buff,sizeof(buff),setup);
	ysigma=strtof(buff,NULL)*A;

	fgets(buff,sizeof(buff),setup);
	gausspeakval=strtof(buff,NULL)*A;

	fgets(buff,sizeof(buff),setup);
	welldepth=strtof(buff,NULL)*eV;

	fgets(buff,sizeof(buff),setup);
	numAds=atoi(buff);


	fgets(buff,sizeof(buff),setup);
	int numsavetimes=atoi(buff);

	fgets(buff,sizeof(buff),setup);
	char* token=strtok(buff," ");

	float* psisavets= malloc(sizeof(float) * numsavetimes);


	for(int i=0;i<numsavetimes;i++){
		psisavets[i]=strtof(token,NULL);
		token = strtok(NULL, " ");
	}

	fgets(buff,sizeof(buff),setup);
	it=atoi(buff);

	fclose(setup);

	strcpy(cwdtmp,cwd);

	FILE *paths=fopen(strcat(cwdtmp,"brownianpaths.txt"),"r");

	char* buff2= malloc(sizeof(char) * flleng*(2*numAds+1));

	//fscanf(paths,"%f",buff2);

	float* ts= malloc(sizeof(float) * (1+(int)round(tFinish/dt)));

	//float* Adspos=[sizeof(ts)][numAds][2];

	Adspos = (float***)malloc(sizeof(ts) * sizeof(float**));
	for (int t = 0;t < sizeof(ts);t++) {
		Adspos[t] = (float**)malloc(numAds * sizeof(float*));
		for (int n = 0;n < numAds;n++) {
			Adspos[t][n] = (float*) malloc(2 * sizeof(float));
		}
	}

	int i=0;


	while(fgets(buff2,sizeof(buff2),paths) != NULL){

		ts[i]=strtof(strtok(buff2," "),NULL);

		for(int j=0;j < numAds;j++){
			Adspos[i][j][0]=strtof(strtok(NULL," "),NULL)*A;
			Adspos[i][j][1]=strtof(strtok(NULL," "),NULL)*A;
		}
		i++;
	}
	fclose(paths);

	if(decaytype==4){

		float* custpot = malloc(sizeof(float) * nz);

		strcpy(cwdtmp,cwd);

		FILE *pot=fopen(strcat(cwdtmp,potentialfile),"r");


		char buff3[256];
		i=0;
		while(fgets(buff3,sizeof(buff3),pot) != NULL){

			custpot[i]=strtof(buff3,NULL);
			i++;
		}

		fclose(pot);
	}

	strcpy(cwdtmp,cwd);

	/*FILE *psifile=fopen(strcat(cwdtmp,"psi.txt"),"r");
	char buff4 [256];

	psi= (float****)malloc(nx*sizeof(char***));//need to make this go bigger

	for(int x=0;x<nx;x++){
		psi[x]=(float***)malloc(ny*sizeof(char**));
		for(int y=0;y<ny;y++){
			psi[x][y]=(float**)malloc(nz*sizeof(char*));
			for(int z=0;z<nz;z++){
				psi[x][y][z]=(float*)malloc(2*sizeof(char));
				for(int im=0;im<2;im++){
					fgets(buff4,sizeof(buff4),psifile);
					psi[x][y][z][im]=strtof(buff4,NULL);
				}
			}
		}
	}
	*/
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	setup();

	long long expv_ptr = mxGetScalar(prhs[0]);
	long long z_offset_ptr = mxGetScalar(prhs[1]);
	long long gauss_pos_ptr = mxGetScalar(prhs[2]);
	long long x0_ptr = mxGetScalar(prhs[3]);
	long long y0_ptr = mxGetScalar(prhs[4]);
	long long expk_ptr = mxGetScalar(prhs[5]);
	long long psi_ptr = mxGetScalar(prhs[6]);
	int nx = mxGetScalar(prhs[7]);
	int ny = mxGetScalar(prhs[8]);
	int nz = mxGetScalar(prhs[9]);
	int decay_type = mxGetScalar(prhs[10]);
	double A = mxGetScalar(prhs[11]);
	double eV = mxGetScalar(prhs[12]);
	double h_bar = mxGetScalar(prhs[13]);
	double dt = mxGetScalar(prhs[14]);
	double dx = mxGetScalar(prhs[15]);
	double dy = mxGetScalar(prhs[16]);
	double dz = mxGetScalar(prhs[17]);
	int it = 0;
	int gfxSteps= mxGetScalar(prhs[18]);
	double t = mxGetScalar(prhs[19]);

	myComplex* dev_expv = reinterpret_cast<myComplex*>(expv_ptr);
	double* dev_z_offset = reinterpret_cast<double*>(z_offset_ptr);
	double* dev_gauss_pos = reinterpret_cast<double*>(gauss_pos_ptr);
	double* dev_x0 = reinterpret_cast<double*>(x0_ptr);
	double* dev_y0 = reinterpret_cast<double*>(y0_ptr);
	myComplex* dev_expk = reinterpret_cast<myComplex*>(expk_ptr);
	myComplex* dev_psi = reinterpret_cast<myComplex*>(psi_ptr);

	const mwSize* gauss_dims = mxGetDimensions(prhs[15]);

	while (it <= gfxSteps)) {
		/*double totprob = 0;
		for (int x = 0;x < nx;x++) {
			for (int y = 0;y < ny;y++) {
				for (int z = 0;z < nz;z++) {
					for (int im = 0;im < 2;im++) {
						totprob = totprob + psi[x][y][z][im] * psi[x][y][z][im];
					}
				}
			}
		}
		

		if (!(fabs(totprob - 1) < eps)) {
			nrhs = 1;
			prhs[0] = 5;//make actual error
		}*/
		mex_split_operator_step_3rd_vsplit_time_dependent(t, dev_expv, dev_z_offset, dev_gauss_pos, dev_x0, dev_y0, dev_expk, dev_psi, nx, ny, nz, decay_type, A, eV, h_bar, dt, dx, dy, dz, it, gauss_dims);
		it++;
		t = t + dt;
	}
	
}

