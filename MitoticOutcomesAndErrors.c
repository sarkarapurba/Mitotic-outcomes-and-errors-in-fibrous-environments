/**********************************************************************************/
/*                  Mitotic outcomes and errors in fibrous environments           */
/**********************************************************************************/
/*			THIS CODE USES                   		       */
/*      	      MONTE CARLO ALGORITHM           			       */ 
/*      	        TO CHARACTERIZE      			               */
/*      	     MITOTIC SPINDLE DEFECTS   				       */ 
/*                    IN FIBROUS ENVIRONMENTS                                     */
/**********************************************************************************/

/*   School of Mathematical and Computational Sciences, Indian Association for the Cultivation of Science, Jadavpur, Kolkata 700032, India    */
/*   Department of Mechanical Engineering, Virginia Tech, Blacksburg, VA 24061                                                                */
/*   Department of Biomedical Engineering and Mechanics, Virginia Tech, Blacksburg, VA 24061                                                  */
/*   Department of Chemical and Biological Physics, Weizmann Institute of Science, Rehovot 7610001, Israel                                    */
/*   Department of Biochemistry and Molecular Biology, Colorado State University, Fort Collins, CO 80523                                      */
/*   (c) 2023                                                                                                                                 */
  
/*   DESIGN & DEVELOPMENTS                                                                                                                    */
  
/*   Apurba Sarkar, Raja Paul                                                                                                                 */
/*   School of Mathematical and Computational Sciences, Indian Association for the Cultivation of Science, Jadavpur, Kolkata 700032, India    */

/*   Aniket Jana, Haonan Zhang, Atharva Agashe, Amrinder S. Nain                                                                              */
/*   Department of Mechanical Engineering, Virginia Tech, Blacksburg, VA 24061                                                                */
   
/*   Ji Wang, Amrinder S. Nain                                                                                                                */
/*   Department of Biomedical Engineering and Mechanics, Virginia Tech, Blacksburg, VA 24061                                                  */

/*   Nir S. Gov                                                                                                                               */
/*   Department of Chemical and Biological Physics, Weizmann Institute of Science, Rehovot 7610001, Israel                                    */

/*   Jennifer G DeLuca                                                                                                                        */
/*   Department of Biochemistry and Molecular Biology, Colorado State University, Fort Collins, CO 80523                                      */
/* *********************************************************************************/                   



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
const double pi=3.141592653589793238;
int idum=-12347;

/**********dyanamic memory allocation********/  
double   *vector ( int nrl, int nrh);
double  **matrix ( int nrl, int nrh, int ncl , int nch );
double ***tensor3 ( int nxl, int nxh, int nyl, int nyh, int nzl, int nzh );

int   *ivector ( int nrl, int nrh );
int  **imatrix ( int nrl, int nrh, int ncl, int nch );
int ***itensor3 ( int nxl, int nxh, int nyl, int nyh, int nzl, int nzh );


double *vector( int nl, int nh)
{
        double *v;

        v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
        if (!v) printf("allocation failure in vector()\n\n");
        return v-nl;
}

double **matrix ( int nrl, int nrh, int ncl, int nch )
{
        int i;
        double **m;

        m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
        if (!m) printf("allocation failure 1 in matrix()\n\n");
        m -= nrl;

        for(i=nrl;i<=nrh;i++) {
                m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
                if (!m[i]) printf("allocation failure 2 in matrix()\n\n");
                m[i] -= ncl;
        }
        return m;
}

double ***tensor3 ( int nxl, int nxh, int nyl, int nyh, int nzl, int nzh )
{
  int i, j;
  double ***m;
  
  m=(double ***) malloc((unsigned) (nxh-nxl+1)*sizeof(double*));
  if (!m) printf("allocation failure 1 in tensor3()\n\n");
  m -= nxl;
  
  for(i=nxl;i<=nxh;i++) {
    m[i]=(double **) malloc((unsigned) (nyh-nyl+1)*sizeof(double*));
    if (!m[i]) printf("allocation failure 2 in tensor3()\n\n");
    m[i] -= nyl;
  };
  
  for(i=nxl;i<=nxh;i++) {
    for(j=nyl;j<=nyh;j++) {
      m[i][j]=(double *) malloc((unsigned) (nzh-nzl+1)*sizeof(double));
      if (!m[i][j]) printf("allocation failure 3 in tensor3()\n\n");
      m[i][j] -= nzl;
    }
  };
    
  return m;

}
  

int *ivector( int nl, int nh)
{
        int *v;

        v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
        if (!v) printf("allocation failure in ivector()\n\n");
        return v-nl;
}

int **imatrix ( int nrl, int nrh, int ncl, int nch )
{
        int i;
        int **m;

        m=(int **) malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
        if (!m) printf("allocation failure 1 in matrix()\n\n");
        m -= nrl;

        for(i=nrl;i<=nrh;i++) {
                m[i]=(int *) malloc((unsigned) (nch-ncl+1)*sizeof(int));
                if (!m[i]) printf("allocation failure 2 in matrix()\n\n");
                m[i] -= ncl;
        }
        return m;
}

int ***itensor3 ( int nxl, int nxh, int nyl, int nyh, int nzl, int nzh )
{
  int i, j;
  int ***m;
  
  m=(int ***) malloc((unsigned) (nxh-nxl+1)*sizeof(int*));
  if (!m) printf("allocation failure 1 in tensor3()\n\n");
  m -= nxl;
  
  for(i=nxl;i<=nxh;i++) {
    m[i]=(int **) malloc((unsigned) (nyh-nyl+1)*sizeof(int*));
    if (!m[i]) printf("allocation failure 2 in tensor3()\n\n");
    m[i] -= nyl;
  };
  
  for(i=nxl;i<=nxh;i++) {
    for(j=nyl;j<=nyh;j++) {
      m[i][j]=(int *) malloc((unsigned) (nzh-nzl+1)*sizeof(int));
      if (!m[i][j]) printf("allocation failure 3 in tensor3()\n\n");
      m[i][j] -= nzl;
    }
  };
    
  return m;

}
/************************************************************/


/************************************************************/
/*       RANDOM NUMBER GENERATOR AND                        */
/*        LIST OF SUBROUTINES                               */
/************************************************************/

#define N_BITS 31
#define IRND1  250
#define IRND2  103
#define MAX_RN 0X7FFFFFFF
#define MAX_R2 0X40000000

int ir[256], irndx1[256], irndx2[256];


int     Pw(double n);
int     n_t1_t2( int n );

double  Time_elapsed(int *idum);
double  ran2(int *idum);

void    OccConfig();

void    NearestNeighbAllocation();
void    Monte_Carlo_Update();
void Particle_Exchange(int S, int X_S, int Y_S, int Z_S,
			   int NB, int X_NB, int Y_NB, int Z_NB);

int counter = 0, t, Lx = 0, Ly, Lz, Lxy, N, MAX_POWER = 0, t1,t2, Nclust,N_total,N_cent,N_ch, label_RF=0; 
int MAX_ENSEMBLE = 1, MIN_ENSEMBLE=1, Npoles=0,Npole_inst=0;
int Ax, Ay, Az, X0, Y0, Z0,Ax_ini,Ay_ini,Az_ini;
int *particle, *occ, *label_particle, *particle_dummy;
int  **neighb;
double Temp = 0.5, Lav = 10.0, Lav_ast = 5.0, fcs_csd=0.0, fcs_csk=0.0, fcs_kt=0.0, fcs_ch=0.0, fch_ch=0.0, fkt_kt=0.0, 
  r_cut, fch_wal=0.0, fcs_wal=0.0, fcs_surf=0.0, dist_merge=1.0, fcs_RetFib=0.0;
double xs,ys,zs,delang,avg_dist,d_fibfromOrigin;
int NoOfFiber,flag_stp;




/******random number generator rndini *******/
void rndini ( int seed ) {
    int irseed, imult=105, i2h23m1=8388607, warmup=10000;
    int i,j,irt,irhalf,sys_size = 256;

    irseed=920331+2*seed;
    irhalf=i2h23m1/2;
    irt=irseed;
    for ( i=0; i<=sys_size-1; i++ ) ir[i]=0;
    for ( i=1; i<=warmup; i++ ) irt = imult*irt % i2h23m1;
    for ( i=0; i<=sys_size-1; i++) {
      for (j=1;j<=105;j++) irt = imult*irt % i2h23m1;
      for (j=1;j<=N_BITS;j++) {
	irt = imult*irt % i2h23m1;
	ir[i] = (ir[i]<<1);
	if (irt > irhalf) ir[i] |= 1;
      };
    };
    for (i=0; i<=sys_size-1; i++) {
      irndx1[i] = (6+i)   & sys_size-1;
      irndx2[i] = (153+i) & sys_size-1;
    };
}

/********************************************************/

/*****Random number generator****/
int Rand(void){
  int j;
  j = counter++ & 255;
  ir[j] = ir[irndx1[j]] ^ ir[irndx2[j]];
  return(ir[j]);
}
/*******************************************************/  
 
void nrerror0(const char error_text[])
{
        printf("Numerical Recipes run-time error...\n");
        printf("%s\n",error_text);
        printf("...now exiting to system...\n");
        exit(1);
}

/**********************************************************/ 

/***** random number generator ran2 ***********************/
#include <math.h>

#define M 714025
#define IA 1366
#define IC 150889

double ran2(int *idum)

{
	static long iy,irn[98];
	static int iff=0;
	int j;
	void nrerror();

	if (*idum < 0 || iff == 0) {
		iff=1;
		if ((*idum=(IC-(*idum)) % M) < 0) *idum = -(*idum);
		for (j=1;j<=97;j++) {
			*idum=(IA*(*idum)+IC) % M;
			irn[j]=(*idum);
		}
		*idum=(IA*(*idum)+IC) % M;
		iy=(*idum);
	}
	j=1 +(int)( 97.0*iy/M);
	if (j > 97 || j < 1) nrerror0("RAN2: This cannot happen.");
	iy=irn[j];
	*idum=(IA*(*idum)+IC) % M;
	irn[j]=(*idum);
	return (double)iy/M;
}

#undef M
#undef IA
#undef IC
/****************************************************************/

/************* 2^n**********************************/ 
int Pw(double n) {
  return((int)(pow(2.0,n)));
}  

//// function for finding the square root of x, y and z
double squareroot(double x, double y, double z){
  return(sqrt((((x)*(x))+((y)*(y))+((z)*(z)))));
}

//// function for finding the dot product of two vectors
double dot(double x1,double y1,double z1,double x2,double y2,double z2){
  return(((x1)*(x2))+((y1)*(y2))+((z1)*(z2)));
}

///Subroutine for intial placement of the centrosomes and chromosomes within the spheroidal volume having sami-axes: Ax_ini, Ay_ini, Az_ini
int OutOfSpheroid_ini(int x, int y, int z){
  //Here, we assume the spheroid is centered at X0,Y0,Z0 and
  //Ax_ini,Ay_ini,Az_ini are the radii in 3 directions
  
  x = x-X0;
  y = y-Y0;
  z = z-Z0;
  
  double x2 = (double)(x*x);
  double y2 = (double)(y*y);
  double z2 = (double)(z*z);
  double Ax2 = (double)(Ax_ini*Ax_ini);
  double Ay2 = (double)(Ay_ini*Ay_ini);
  double Az2 = (double)(Az_ini*Az_ini);
  
  
  if(x2/Ax2+y2/Ay2+z2/Az2<1.0){
    return(0);//inside the spheroid
  }else return(1);//outside the spheroid
}


///Subroutines to confine the particle's movement within the spheroidal volume having semi-axes: Ax, Ay, Az
int OutOfSpheroid(int x, int y, int z){
  //Here, we assume the spheroid is centered at X0,Y0,Z0 and
  //Ax,Ay,Az are the radii in 3 directions
  
  x = x-X0;
  y = y-Y0;
  z = z-Z0;
  
  double x2 = (double)(x*x);
  double y2 = (double)(y*y);
  double z2 = (double)(z*z);
  double Ax2 = (double)(Ax*Ax);
  double Ay2 = (double)(Ay*Ay);
  double Az2 = (double)(Az*Az);
  
  
  if(x2/Ax2+y2/Ay2+z2/Az2<1.0){
    return(0);//inside the spheroid
  }else return(1);//outside the spheroid
}

/// subroutine to define the lattice sites inside the cellular volume 
void DefineCellVol(){
  int i,j,k,site;
  for(k=1;k<=Lz;k++){
    for(i=1;i<=Lx;i++){
      for(j=1;j<=Ly;j++){
	site = j + (i-1)*Lx + (k-1)*Lxy;
	occ[site] = 0;
	particle[site] = 0;
	if(OutOfSpheroid(i,j,k)==1){
	  occ[site] = 2;
	};
      };
    };
  };  
}

/// Subroutine to define the surface nodes of the spheroidal cell///
void DefineSurface(){
  int i,j,k,l,site,nb;
  for(k=1;k<=Lz;k++){
    for(i=1;i<=Lx;i++){
      for(j=1;j<=Ly;j++){
	site = j + (i-1)*Lx + (k-1)*Lxy;
	if(occ[site] == 2){
	  for(l=1;l<=6;l++){
	    nb=neighb[site][l];
	    if(occ[nb] !=2)particle[site]=3;
	  };
	};
      };
    };
  };
}


////Finding the retraction fiber spots depending on the inter-fiber spacing.
////Finding the azimuthal angle for the point where the fiber line hits the cell surface
double RetfractionFiberSpotDependondist(double d){                 
  double Ax2,Ay2,Az2,x0,y0,z0,x1,y1,z1,L,M,N,A,A1,A2,A3,A4,A5,t1,t2,x1s,y1s,z1s,x2s,y2s,z2s,xsu,ysu,zsu,cos_phi,phi;
  
  Ax2 = (double)(Ax*Ax);
  Ay2 = (double)(Ay*Ay);
  Az2 = (double)(Az*Az);

  x0=d;          //if the fiber is placed along the y axis  //here d is the separation between y axis and one of the fibers paraller to the y axis.
  y0=0.0;
  z0=0.0;

  x1=d;
  y1=1.0;
  z1=0.0;
    
  L=x1-x0;                                      
  M=y1-y0;
  N=z1-z0;
  
  A= (Ax2*Ay2*N*N)+(Ax2*Az2*M*M)+(Ay2*Az2*L*L);     

  A1=(Ax2*Ay2*N*z0)+(Ax2*Az2*M*y0)+(Ay2*Az2*L*x0);
 
  
  A2= Ax2*Ay2*Az2;
  A3=A;
  A4=(Ax2*M*M*z0*z0)+(Ay2*N*N*x0*x0)+(Ax2*N*N*y0*y0)+(Ay2*L*L*z0*z0)+(Az2*M*M*x0*x0)+(Az2*L*L*y0*y0);
  A5=(2.0*Ax2*M*N*y0*z0)+(2.0*Ay2*L*N*x0*z0)+(2.0*Az2*L*M*x0*y0);


  t1=(1.0/A)*(-A1+sqrt(A2*(A3-A4+A5)));
  t2=-(1.0/A)*(A1+sqrt(A2*(A3-A4+A5)));

  x1s=(double)x0+t1*L;
  y1s=(double)y0+t1*M;
  z1s=(double)z0+t1*N;

  x2s=(double)x0+t2*L;
  y2s=(double)y0+t2*M;
  z2s=(double)z0+t2*N;
  
  if(dot((double)(x1-x0),double(y1-y0),double(z1-z0),(x1s-(double)x0),(y1s-(double)y0),(z1s-(double)z0)) >=0){
    xsu=x1s;
    ysu=y1s;
    zsu=z1s;
  }
  if(dot((double)(z1-x0),double(y1-y0),double(z1-z0),(x2s-(double)x0),(y2s-(double)y0),(z2s-(double)z0)) >=0){
    xsu=x2s;
    ysu=y2s;
    zsu=z2s;
  }

  cos_phi= dot((double)(x1-x0),double(y1-y0),double(z1-z0),(xsu-(double)x0),(ysu-(double)y0),(zsu-(double)z0))/(squareroot(x0,y0,z0)*squareroot(xsu,ysu,zsu));
  cos_phi= dot(x0,y0,z0,xsu,ysu,zsu)/(squareroot(x0,y0,z0)*squareroot(xsu,ysu,zsu));

  phi=acos(cos_phi);
  if(dot(x0,y0,z0,xsu,ysu,zsu)==0){phi=pi/2.0;};
  return(phi);
}

//To find the actual point on the cell surface where the line connecting the cell center and simulated cortical node with particle[site]=3 hits the cell surface,
//since particle[site]=3 does not necessarily mean they are exactly on the cell surface, they may reside slightly inside the actual cell surface.
void EllipsoidToLineIntersectionPoint(int i, int j, int k){
  double Ax2,Ay2,Az2,L,M,N,A,A1,A2,A3,t1,t2,x1s,y1s,z1s,x2s,y2s,z2s;
  
  Ax2 = (double)(Ax*Ax);
  Ay2 = (double)(Ay*Ay);
  Az2 = (double)(Az*Az);
  
  L=(double)(i-X0);
  M=(double)(j-Y0);
  N=(double)(k-Z0);
  
  A= (Ax2*Ay2*N*N)+(Ax2*Az2*M*M)+(Ay2*Az2*L*L);
  
  A1= Ax2*Ay2*Az2;
  A2=A;
  
  t1=(1.0/A)*sqrt(A1*A2);
  t2=-(1.0/A)*sqrt(A1*A2);

  x1s=(double)X0+t1*L;
  y1s=(double)Y0+t1*M;
  z1s=(double)Z0+t1*N;

  x2s=(double)X0-t1*L;
  y2s=(double)Y0-t1*M;
  z2s=(double)Z0-t1*N;
  
  if(dot((double)(i-X0),double(j-Y0),double(k-Z0),(x1s-(double)X0),(y1s-(double)Y0),(z1s-(double)Z0)) >=0){
    xs=x1s;
    ys=y1s;
    zs=z1s;
  }
  if(dot((double)(i-X0),double(j-Y0),double(k-Z0),(x2s-(double)X0),(y2s-(double)Y0),(z2s-(double)Z0)) >=0){
    xs=x2s;
    ys=y2s;
    zs=z2s;
  }
}


////Finding the suitable RF spots (for fiber cases)/ bands (for 2D-rounded case) on the cell surface
void RetractionFiber(){   
  int i,j,k,site, Nl, nl;
  double tan_phi,cos_theta,sinep,cosinep,phi,theta, dfo, phi_RF;  
  for(k=1;k<=Lz;k++){
    for(i=1;i<=Lx;i++){
      for(j=1;j<=Ly;j++){
	site = j + (i-1)*Lx + (k-1)*Lxy;	
	if(particle[site]==3){
	  EllipsoidToLineIntersectionPoint(i,j,k);	  
	  sinep=(double)Ay*(ys-(double)Y0);
	  cosinep=(double)Ax*(xs-(double)X0);
	  phi=atan2(sinep,cosinep);
	  phi+=(2.0*pi);
	  phi=fmod(phi,(2.0*pi));
	  cos_theta=(zs-(double)Z0)/(double)Az;
	  theta=acos(cos_theta);
	  if(label_RF==1){
	    if (theta >= pi/2.0 && theta <= (pi/2.0 + 2.0*delang)){   //RF spots are extending from the equatorial plane towards below	 
	      if(NoOfFiber == 1){
		// Construction of two RF-spots for 1F-elonagted case
		if (phi >= (pi/2.0-delang) && phi <= (pi/2.0 + delang)){
		  particle[site]=4;
		}
		if (phi >= (pi+pi/2.0-delang) && phi <= (pi+pi/2.0 + delang)){	      
		  particle[site]=4;
		}	    
	      }
	    }
	  }
	  
	  if(label_RF==2){
	    if (theta >= pi/2.0 && theta <= (pi/2.0 + 2.0*delang)){   //RF spots are extending from the equatorial plane towards below	 
	      if(NoOfFiber % 2 == 0){  //For two fiber cases	      	
		Nl=NoOfFiber/2;
		for(nl=1;nl<=Nl;nl++){
		  dfo=d_fibfromOrigin*nl;
		  phi_RF= RetfractionFiberSpotDependondist(dfo);
		                                                // Construction of four RF-spots for 2F-elonagted case
		  if (phi >= (phi_RF-delang) && phi <= (phi_RF + delang)){
		    particle[site]=4;
		  }
		  if (phi >= (2.0*pi-phi_RF-delang) && phi <= (2.0*pi-phi_RF + delang)){	      
		    particle[site]=4;
		  }
		  
		  if (phi >= (pi-phi_RF-delang) && phi <= (pi-phi_RF + delang)){
		    particle[site]=4;
		  }
		  if (phi >= ((2.0*pi)-(pi-phi_RF)-delang) && phi <= ((2.0*pi)-(pi-phi_RF) + delang)){	      
		    particle[site]=4;
		  }		
		}
	      }
      
	      else{                                             // Construction of multiple RF-spots (2 times the number of external fibers) for multi-fiber scenario
		if (phi >= (pi/2.0-delang) && phi <= (pi/2.0 + delang)){
		  particle[site]=4;
		}
		if (phi >= (pi+pi/2.0-delang) && phi <= (pi+pi/2.0 + delang)){	      
		  particle[site]=4;
		}
		
		Nl=(NoOfFiber-1)/2;
		for(nl=1;nl<=Nl;nl++){
		  dfo=d_fibfromOrigin*nl;
		  phi_RF= RetfractionFiberSpotDependondist(dfo);
		  if (phi >= (phi_RF-delang) && phi <= (phi_RF + delang)){
		    particle[site]=4;
		  }
		  if (phi >= (2.0*pi-phi_RF-delang) && phi <= (2.0*pi-phi_RF + delang)){	      
		    particle[site]=4;
		  }
		  
		  if (phi >= (pi-phi_RF-delang) && phi <= (pi-phi_RF + delang)){
		    particle[site]=4;
		  }
		  if (phi >= ((2.0*pi)-(pi-phi_RF)-delang) && phi <= ((2.0*pi)-(pi-phi_RF) + delang)){	      
		    particle[site]=4;
		  }
	  
		}
	      }
      
	    }
	  }
  
	  if(label_RF==3){
	    if (theta >= pi/2.0 && theta <= (pi/2.0 + 2.0*delang)){   //RF spots are extending from the equatorial plane towards below	 
	      particle[site]=4;                                       //Construction of band-like RF distributions
	    }  
	  }
	  
	  if(label_RF==4){
	    if (theta >= pi/2.0 && theta <= (pi/2.0 + 2.0*delang)){
	                                                             // Construction of 3 spots for 2F-kite scenario
	      if (phi >= (pi/2.0-delang) && phi <= (pi/2.0 + delang)){   
		particle[site]=4;
	      }
	      if(phi >= ((2.0*pi)-delang) && phi <= delang){	      
		particle[site]=4;
	      }
	      if(phi >= (pi+ (pi/4.0) -delang) && phi <= (pi+ (pi/4.0) +delang)){	      
		particle[site]=4;
	      }
	    }
	  }
	}
      }		
    }
  }
}  


/**************Particle (centrosomes/chromosomes) Occupancy in the lattice*****************/  
void OccConfig(){
  int i,count;
  
  //Define cell volume
  DefineCellVol();
  
  //Define cell surface
  DefineSurface();
  
  //Define RF-regions on cell surface
  RetractionFiber();
  
  
  // distribute all the centrosomes and chromosomes within the volume having semi-axes Ax_ini, Ay_ini, Az_ini
  count = 0;
  if(N_cent>0){
    do{
      int x = X0 + 2*(int)((double)Ax_ini*(0.5-ran2(&idum)));
      int y = Y0 + 2*(int)((double)Ay_ini*(0.5-ran2(&idum)));
      int z = Z0 + 2*(int)((double)Az_ini*(0.5-ran2(&idum)));
      
      if(x==0)x=1;
      if(y==0)y=1;
      if(z==0)z=1;
      
      if(OutOfSpheroid_ini(x,y,z)==0){
	
	int site = x+(y-1)*Lx+(z-1)*Lxy;
	
	if(occ[site]==0){
	  occ[site]=1;
	  particle[site]= 1;
	  count++;
	};
      };
    }while(count<N_cent);
  };
  
  count = 0;
  if(N_ch>0){
    do{
      int x = X0 + 2*(int)((double)Ax_ini*(0.5-ran2(&idum)));
      int y = Y0 + 2*(int)((double)Ay_ini*(0.5-ran2(&idum)));
      int z = Z0 + 2*(int)((double)Az_ini*(0.5-ran2(&idum)));
      if(x==0)x=1;
      if(y==0)y=1;
      if(z==0)z=1;
      if(OutOfSpheroid_ini(x,y,z)==0){
	int site = x+(y-1)*Lx+(z-1)*Lxy;
	if(occ[site]==0){
	  occ[site]=1;
	  particle[site]= 2;
	  count++;
	};
      };
    }while(count<N_ch);
  };
}

/***********neighbouring array with pbc******************/    
void NearestNeighbAllocation(){
  int i,j;
  
  /* neigbor array*/
  for( i=1; i<=N; i++ ){  /*                             4  5  */
    neighb[i][1]=i+1;     /* right      neighbor         |/    */
    neighb[i][2]=i+Lx;     /* bottom     neighbor      3-i-1   */
    neighb[i][3]=i-1;     /* left       neighbor        /|     */
    neighb[i][4]=i-Lx;     /* top        neighbor     6  2     */
    neighb[i][5]=i-Lxy;   /* previous   neighbor               */
    neighb[i][6]=i+Lxy;   /* next       neighbor               */     
  };
  
  /* boundary condition */
  for(i=1;i<=Lx;i++){
    for ( j=1; j<=Ly; j++ ){
      neighb [j+(Lxy*i-Lx)]        [2] = j + (i-1)*Lxy;
      neighb [j+(i-1)*Lxy]         [4] = j + (Lxy*i -Lx);
      
      neighb [j*Lx+(i-1)*Lxy]      [1] = 1+(j-1)*Lx+(i-1)*Lxy;     
      neighb [1+(j-1)*Lx+(i-1)*Lxy] [3] = j*Lx+(i-1)*Lxy;
      
      neighb [j + (i-1)*Lz]         [5] = j + (i-1)*Lz  + N - Lxy;
      neighb [j + (i-1)*Lz + N -Lxy][6] = j + (i-1)*Lz;
    };
  };
}

/**************************************************************/
/******* Monte Carlo Update of particle's position *******/   
void Monte_Carlo_Update(){
  int i,j,k,l,m,n,x,y,z,nb,s;
  int excess;
  double flag;
  
  /* main monte carlo routine */
  for(k=1;k<=Lz;k++){
    for(i=1;i<=Lx;i++){
      for(j=1;j<=Ly;j++){
	s = j + (i-1)*Lx + (k-1)*Lxy;
	x = (int)(double)(s%Lxy)/((double)Lx+0.1)+1;
	if (s%Lxy==0) x=Lx;
	y = (s-1)%Lx+1;
	z = (int)((double)s/((double)Lxy+0.1))+1;	
	if(occ[s] == 1){
	  flag = ran2(&idum);
	  nb = (int)(6.0*flag)+1;
	  n = neighb[s][nb];	  
	  x = (int)(double)(n%Lxy)/((double)Lx+0.1)+1;
	  if (n%Lxy==0) x=Lx;
	  y = (n-1)%Lx+1;
	  z = (int)((double)n/((double)Lxy+0.1))+1;
	  if(occ[n] == 0) {  
	    Particle_Exchange(s,i,j,k,n,x,y,z);
	  }; 
	};
      };
    };
  };
}


/**************************************************************/
/***Subrountine for the movement of centrosomes and chromosomes to it's randomly chosen neighbor***/
void Particle_Exchange(int S, int X_S, int Y_S, int Z_S,
		     int NB, int X_NB, int Y_NB, int Z_NB){
  
  
  int i,j,k,l,sum,site;
  double DeltaE,flag,energy1,energy2,prob;
  int sumleft = 0;
  int nb_S,nb_NB;

  DeltaE = 0.0;
  for(k=1;k<=Lz;k++){
    for(i=1;i<=Lx;i++){
      for(j=1;j<=Ly;j++){
	site = j + (i-1)*Lx + (k-1)*Lxy;
	
	if((site != S) && (occ[site]!=0)){
	  double r1 = sqrt((double)((X_S-i)*(X_S-i)+ (Y_S-j)*(Y_S-j) + (Z_S-k)*(Z_S-k))); 
	  double r2 = sqrt((double)((X_NB-i)*(X_NB-i)+ (Y_NB-j)*(Y_NB-j) + (Z_NB-k)*(Z_NB-k)));
	  if(particle[S]==1){//particle[S]=1 is centrosome 
	    if(particle[site]==1){//particle[site]=1 is centrosome
	      DeltaE += fcs_csk*Lav*(((r2*exp(-r2/Lav))-(r1*exp(-r1/Lav)))+(Lav*(exp(-r2/Lav)-exp(-r1/Lav))));
	      DeltaE += fcs_csd*Lav*(exp(-r2/Lav)-exp(-r1/Lav));
	    };
	    if(particle[site]==2){//particle[S]=2 is kinetochore and chromosome	      
	      DeltaE += -fcs_kt*(r2-r1);
	      DeltaE += fcs_ch*Lav*(exp(-r2/Lav)-exp(-r1/Lav));
	      //printf("2 NB=%d S=%d site=%d i=%d j=%d DeltaE=%f r1=%f r2=%f\n",NB,S,site,i,j,DeltaE,r1,r2);
	    };
	    
	    if(particle[site]==3){//particle[S]=3 is the boundary; interaction with centrosome
	      DeltaE += fcs_surf*Lav_ast*(exp(-r2/ Lav_ast)-exp(-r1/ Lav_ast));
	    };

	    if(particle[site]==4){//particle[S]=4 is the Rf-spots; interaction with centrosome
	      DeltaE += fcs_RetFib*Lav_ast*(exp(-r2/ Lav_ast)-exp(-r1/ Lav_ast));
	    };
	    
	  };
	
	  if(particle[S]==2){//particle[S]=2 is kinetochore and chromosome
	    if(particle[site]==1){//particle[site]=1 is centrosome
	      DeltaE += -fcs_kt*(r2-r1);
	      DeltaE += fcs_ch*Lav*(exp(-r2/Lav)-exp(-r1/Lav));
	    };
	    
	    //centrosome-boundary steric repulsion
	    energy1=energy2=0.0;
	    for(l=1;l<=6;l++){
	      nb_S = neighb[S][l];
	      if(occ[nb_S]==2)
		energy1 += fcs_wal;
	      nb_NB = neighb[NB][l];
	      if(occ[nb_NB]==2)
		energy2 += fcs_wal;
	    };
	    DeltaE += energy2-energy1;
	      
	  };
	  	  
	  //chromosome-chromosome steric repulsion
	  if((particle[S]==2) && (particle[site]==2)){//particle[S]=2 is kinetochore and chromosome
	    energy1=energy2=0.0;
	    if(r1<=r_cut) 
	      energy1 = fch_ch*(1.0/r1);
	    if(r2<=r_cut) 
	      energy2 = fch_ch*(1.0/r2);
	    DeltaE += energy2-energy1;
	  };

	  //chromosome-boundary steric repulsion
	  if(particle[S]==2){
	    energy1=energy2=0.0;
	    for(l=1;l<=6;l++){
	      nb_S = neighb[S][l];
	      if(occ[nb_S]==2)
		energy1 += fch_wal;
	      nb_NB = neighb[NB][l];
	      if(occ[nb_NB]==2)
		energy2 += fch_wal;
	    };
	    DeltaE += energy2-energy1;
	  };
	  
	  //kt-kt long distance attraction
	  if((particle[S]==2) && (particle[site]==2)){//particle[S]=2 is kinetochore and chromosome
	    if((r1>4.*r_cut) && (r2>4.*r_cut))
	      DeltaE += fkt_kt*Lav*(exp(-r2/Lav)-exp(-r1/Lav));
	  };
	};
      };
    };
  };

  if(DeltaE<=0.0){
    particle[NB] = particle[S];
    particle[S] = 0;
    occ[NB]=1;
    occ[S]=0;
  }else{ 
    flag = ran2(&idum);
    prob = exp(-DeltaE/Temp) ;
    if (flag  < prob) {
      particle[NB] = particle[S];
      particle[S] = 0;
      occ[NB]=1;
      occ[S]=0;
    };
  };
}

///Subroutine to visualize the system using POV-Ray///
void Visualization_povray(int tt, int ensemble){    
    char fn[100];
    FILE *fp;    
    int i,j,k,s,ncnt;
    int x,y,z;
    
    sprintf(fn,"%dsnp%dcnt_%dNch%dAx%dEnsemble.pov",100000+tt,N_cent, N_ch, Ax,ensemble);
    fp = fopen(fn,"w");
    
    fprintf(fp,"background{color rgb <1, 1, 1>}\n");
    fprintf(fp,"\n");
    
   /*********** Cell ***********/
    fprintf(fp,"sphere{<%d, %d, %d>, 1 scale<%d, %d, %d>\n",0,0,0,Ax,Ay,Az);
    //fprintf(fp,"translate<%d, %d, %d>\n",X0,Y0,Z0);
    fprintf(fp,"finish {phong 0.3}\n");
    fprintf(fp,"pigment{color rgbt<0.85, 0.7, 0.9, 0.8>}}\n");
    fprintf(fp,"\n");
    /**************************/
    for(k=1;k<=Lz;k++){
      for(i=1;i<=Lx;i++){
	for(j=1;j<=Ly;j++){
	  s = j + (i-1)*Lx + (k-1)*Lxy;
	  
	  //Translate the center to 0,0,0
	  x = i-X0; y = j-Y0; z = k-Z0; 
	  if(particle[s] == 1){
	    fprintf(fp,"sphere{<%d, %d, %d>, %f\n",x,y,z,0.5);
	    fprintf(fp,"pigment{color rgbt <0.9, 0.5, 0.5, 0.0>}\n");
	    fprintf(fp,"finish{ambient 0.5 specular 0.5}}\n");
	  };
	  if(particle[s] == 2){
	    fprintf(fp,"sphere{<%d, %d, %d>, %f\n",x,y,z,0.5);
	    fprintf(fp,"pigment{color rgbt <0.5, 0.9, 0.5, 0.0>}\n");
	    fprintf(fp,"finish{ambient 0.5 specular 0.5}}\n");
	  };
	  
	  /*if(particle[s] == 3){  //(particle[s] == 3) corresponds to the cortex nodes.
	    fprintf(fp,"sphere{<%d, %d, %d>, %f\n",x,y,z,0.1);
	    fprintf(fp,"pigment{color rgbt <0.8, 0.0, 0.0, 0.2>}\n");
	    fprintf(fp,"finish{ambient 0.5 specular 0.5}}\n");
	    };*/
	  
	  /*if(particle[s] == 4){    //(particle[s] == 3) corresponds to the cortical regions enriched with RFs.
	    fprintf(fp,"sphere{<%d, %d, %d>, %f\n",x,y,z,0.1);
	    fprintf(fp,"pigment{color rgbt <0.0, 0.0, 0.8, 0.0>}\n");
	    fprintf(fp,"finish{ambient 0.5 specular 0.5}}\n");
	    };*/
			  
	};
	//fprintf(fp,"\n");
      };    
    };
   

  //light source    
  fprintf(fp,"\n");
  fprintf(fp,"light_source\n");
  fprintf(fp,"{\n");
  
  fprintf(fp,"<%d, %d, %d> color rgb <1, 1, 1>\n",-(Lx+100),(Ly+100),(Lz+100));
  fprintf(fp,"area_light <-150, 0, 0>, <0, 150, 0>, 25, 25\n");
  fprintf(fp,"jitter\n"); 
  fprintf(fp,"adaptive 1\n");
  fprintf(fp,"rotate <60, 60, 60>\n");
  fprintf(fp,"}\n");
  
  fprintf(fp,"\n");
  fprintf(fp,"light_source\n");
  fprintf(fp,"{\n");
  fprintf(fp,"<%d, %d, %d> color rgb <1, 1, 1>\n",(Lx+100),-(Ly+100),(Lz+100));
  fprintf(fp,"area_light <-150, 0, 0>, <0, 150, 0>, 25, 25\n");
  fprintf(fp,"jitter\n"); 
  fprintf(fp,"adaptive 1\n");
  fprintf(fp,"rotate <60, 60, 60>\n");
  fprintf(fp,"}\n");
  
  fprintf(fp,"\n");
  fprintf(fp,"light_source\n");
  fprintf(fp,"{\n");
  fprintf(fp,"<%d, %d, %d> color rgb <1, 1, 1>\n",(Lx+100),(Ly+100),-(Lz+100));
  fprintf(fp,"area_light <-150, 0, 0>, <0, 150, 0>, 25, 25\n");
  fprintf(fp,"jitter\n"); 
  fprintf(fp,"adaptive 1\n");
  fprintf(fp,"rotate <60, 60, 60>\n");
  fprintf(fp,"}\n");
  
  // camera position
  fprintf(fp,"\n");
  fprintf(fp,"camera\n");
  fprintf(fp,"{\n");
  fprintf(fp,"sky <0,0,1>           //Don't change this\n");
  fprintf(fp,"direction <-1,0,0>    //Don't change this\n");  
  fprintf(fp,"right <-4/3,0,0>      //Don't change this\n");
  fprintf(fp,"location <%d, %d, %d>\n",0,60,0);
  fprintf(fp,"look_at <%d, %d, %d>\n",0,0,0);
  fprintf(fp,"rotate <0,0,0>");
  
  fprintf(fp,"\n");
  fprintf(fp,"}\n");
  fclose(fp);    

}

///final positions of centrosomes and chromosomes////
void Final_configuration(int ensemble,int mc_max){
  FILE *fp1, *fp2;
  char fn1[200],fn2[200];
  int i,j,k,s,x,y,z;
  
  sprintf(fn1,"Centrosome_config_time_power%d_Ncent%d_Nch%d_Fcs_csk%1.1f_Fcs_csd%1.1f_Fcs_kt%1.1f_Fcs_ch%1.1f_Fcs_surf%1.2f_Fcs_RetFib%1.2f_ensemble%d.dat",mc_max,N_cent,N_ch,fcs_csk,fcs_csd,fcs_kt,fcs_ch,fcs_surf,fcs_RetFib,ensemble);
  fp1=fopen(fn1,"w"); 
  sprintf(fn2,"Chromosome_config_time_power%d_Ncent%d_Nch%d_Fcs_csk%1.1f_Fcs_csd%1.1f_Fcs_kt%1.1f_Fcs_ch%1.1f_Fcs_surf%1.2f_Fcs_RetFib%1.2f_ensemble%d.dat",mc_max,N_cent,N_ch,fcs_csk,fcs_csd,fcs_kt,fcs_ch,fcs_surf,fcs_RetFib,ensemble);
  fp2=fopen(fn2,"w"); 

  for(k=1;k<=Lz;k++){
    for(i=1;i<=Lx;i++){
      for(j=1;j<=Ly;j++){
	s = j + (i-1)*Lx + (k-1)*Lxy;
        //Translate the center to 0,0,0
	x = i-X0; y = j-Y0; z = k-Z0; 
	if(particle[s]==1){
	  fprintf(fp1,"%d %d %d\n",x,y,z);
	}
	if(particle[s]==2){
	  fprintf(fp2,"%d %d %d\n",x,y,z);
	}	
      }
    }
  }
  fclose(fp1);
  fclose(fp2);
}


/***************Subroutine to identify the number of CSs assemblies*********/
void IdentifyPoles(){
  int i,j,k,l,p,q,r,sum,site,site1,site2;
  double dist;
  
  for(site=1;site<=N;site++){
    label_particle[site]=0;
    particle_dummy[site]=0;
  };

  for(k=1;k<=Lz;k++){
    for(i=1;i<=Lx;i++){
      for(j=1;j<=Ly;j++){
	site1 = j + (i-1)*Lx + (k-1)*Lxy;
	if((particle[site1]==1) && (label_particle[site1]==0)){//particle[site]=1 is centrosome
	  label_particle[site1]=1;
	  particle_dummy[site1]=1;//put a centrosome in the dummy lattice 
	  for(p=1;p<=Lz;p++){
	    for(q=1;q<=Lx;q++){
	      for(r=1;r<=Ly;r++){
		site2 = r + (q-1)*Lx + (p-1)*Lxy;
		if((particle[site2]==1) && (site2!=site1)){// && (label_particle[site2]==0)){
		  dist = sqrt((double)((p-k)*(p-k)+(q-i)*(q-i)+(r-j)*(r-j)));
		  if(dist<=dist_merge){//if centrosomes are closer than dist_merge, then merge them.
		    if(label_particle[site2]==0)
		      label_particle[site2]=1; //include this particle into the same cluster if distance is less
		    else particle_dummy[site1]=0;
		  }
		};
	      };
	    };
	  };
	};
      };
    };
  };
  
  Npole_inst = 0;
  for(site=1;site<=N;site++){
    if(particle_dummy[site]==1){
      Npole_inst += 1;
      Npoles += 1;
    };
  };  
}

////Subroutine to obtain the histogram of final number of centrosomal clusters
void OutputPolesHist(int ensemble){
  char fn[100];
  FILE *fp;
  
  int i,j,k,s,ncnt;
  int x,y,z;

  sprintf(fn,"HistPole%dcnt%dNch%1.2fMergeDistFcs_surf%1.2f_Fcs_RetFib%1.2f_d_RFSpot%1.2f.dat",N_cent,N_ch,dist_merge,fcs_surf,fcs_RetFib,d_fibfromOrigin);
  if(ensemble==1){
    fp = fopen(fn,"w");
  }else fp = fopen(fn,"a");
  
  fprintf(fp,"%d\n",Npole_inst);
  
  fclose(fp);
}

/***********************************************/
/*                                             */
/*             Mitotic Spindle                 */
/*      Centerosome-Chromosome Patterning      */
/*                                             */
/***********************************************/

main( int argc, char** argv )
{
  char fn[20], fn2[40], fn3[40], fn4[40];
  FILE *fp, *fp2, *fp3, *fp4;
  
  int i, j, k,r, x,y,n, idum, seed, ensemble, mcstp,mc_max;
  double p;
  
  /************************** PARAMETER TABLE ****************************/ 
  Lx=32;                    // Linear system size along X-axis
  Ly = Lx;
  Lz = Lx;
  Lxy = Lx*Ly;
  mc_max = 10000;           // Run time (variable) 
  N = Lx*Ly*Lz;             // Total number of lattice points
  N_cent = 8;               // Number of centrosomes
  N_ch = 46;                // Number of chromosomes
              
  MAX_ENSEMBLE=1000;        // No of ensemble
    
  N_total = N_cent+N_ch;    // Total number of particles
  Temp = 0.05;
  fcs_csd= -1.5;            // amplitude of inter-centrosomal attraction
  fcs_csk=  0.0;            // amplitude of inter-centrosomal repulsion
  fcs_kt=-4.0;              // amplitude of centrosome-kinetochore attraction
  fcs_ch=10.0;              // amplitude of centrosome-chromosome repulsion
  fch_ch=4.0;               // steric repulsion between pairs of chromosomes
  fkt_kt=-0.0;              // long distance attraction between KTs
  r_cut = 2.0;              // cut-off distance for steric repulsion between pairs of chromosomes
  fch_wal = 2.0;            // steric repulsion between wall and chromosome
  fcs_wal= 0.0;             // steric repulsion between wall and centrosome
  fcs_surf = -10.0/100.0;   // attraction between centrosome and cell surface
  fcs_RetFib = -75.0/100.0; // attraction between centrosome and Retraction fiber (RF) zones

  //Center of the cell
  X0 = Lx/2; Y0 = Ly/2; Z0 = Lz/2;  //center of the spheroidal cell
  Ax = 15; Ay=15; Az=15;            //semi-axes of the spheroidal cell 
  Ax_ini=Ax; Ay_ini=Ay; Az_ini=Az;  //Defining initial volume for the distribution of centrosomes and chromosomes
  
  Lav = 18.0;                       //range of CS-CH and CS-CS interactions
  Lav_ast = Lav/4.0;                //range of CS and cortex (or RFs) interaction

  dist_merge = 1.5;                 //merging distance below which centrosomes are considered clustered

  label_RF=2;                       //label_RF = 1, 2, 3, and 4 corresponds to 1F-elongated, 2F-elongated (and also MF-elongated), 2D-rounded, and 2F-kite scenarios, respectively.
  NoOfFiber=2;                      //number of external fibers (for 2F-elongated case)---For 1F-elongated case, NoOfFiber=1;
  d_fibfromOrigin=5.5;              //Distance between the external fibers
  delang= 10.0*pi/180.0;            //corresponds to 1/2 times the value of d\theta (or d\phi) in the SI Appendix text

  /********* Memory allocation ***********/
  occ        = ivector (0,N); 
  particle      = ivector (0,N);
  neighb     = imatrix (0,N, 0,6); 
  label_particle = ivector (0,N);
  particle_dummy = ivector (0,N);
    
  Npoles=0;
  NearestNeighbAllocation();        //Define nearest neighbor allocation
  
  for(ensemble = 1; ensemble<=MAX_ENSEMBLE; ensemble ++){     //starting of ensemble loop 
    OccConfig();
    n =0 ;
    p = 0.0;
    t1 = Pw(p);    
    for(mcstp=1; mcstp<=mc_max; mcstp++){                     //starting of Monte Carlo loop
      Monte_Carlo_Update();
      //Visualization_povray(mcstp, ensemble);
    };
    IdentifyPoles();
    Final_configuration(ensemble,mc_max);
    OutputPolesHist(ensemble);
  }; 
}














