/* === This program performs a Molecular dynamics simulation of N beads, interacting through a van der Waals attraction/repulsion, V=4E( (s/r)^12 - (s/r)^6 ), within a cutoff distance. Dynamics is generated using the Leap-frog MD algorithm. PBC conditions taken into account === */
#include<iostream> 
#include<math.h>  
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>                                                                                   

#define DIM 3 // Number of dimensions (X,Y,Z)
#define PI 3.141592653589793   

class Catom  /* === ATOM CLASS === */
{
private:
  double r[DIM],  v_shifted_half[DIM] , f[DIM]; // positions, velocity (half time-step before), and force vectors
  double m, m_inv; // mass , 1/m will be used so it is good idea to compute it once  
  int ai; // atom index (in the array, atom must now which index it is)
  double s, e; // sigma and epsilon for the van der Waals potential
  double s6, e_s6_24; // auxiliary  optimization variables s^6 and e_s6_24; 
  
public:
  void initialize(double *r0, int ai0, double m0, double s0, double e0) ;
  void v_shifted_minus_half(double *v0, double dt) ;
  void computeforce( Catom  *atom , int Natoms , double *BOX, double cutoff2);
  void  move(double dt, double *BOX ) ; 
  double getr(int i) {return r[i];};
  double getv(int i) {return v_shifted_half[i];};  // note you could print any other "private" variable of atom in this way.
};

// Read initial positions (t=0), atom index, mass, sigma, and epsilon
void Catom::initialize(double *r0, int ai0, double m0, double s0, double e0) {
  for(int i=0;i<DIM;i++) r[i]=r0[i]; // r(t=0)
  ai=ai0;   m=m0;   m_inv=1.0/m;   s=s0;   e=e0;   s6=pow(s,6);   e_s6_24=e*s6*24; 
}

// compute velocity at time t=0 -dt/2  from v(t=0)
void Catom::v_shifted_minus_half(double *v0, double dt) {
  for(int i=0;i<DIM;i++) v_shifted_half[i]=v0[i]-m_inv*f[i]*dt*0.5 ; 
}

  // short range Van der Waals force. All against all atoms  (Catom array with Natoms members), in a BOX, within a cutoff (it reads cutoff2).
void Catom::computeforce( Catom  *atom , int Natoms , double *BOX, double cutoff2){
  // V=4e[ (s/rij)^12 - (s/rij)^6 ]
  // f= -4e [ 12s^12/|rij|^13 - 6s^6/|rij|^7  ] rij / |rij|
  // f= 24*e*s^6* [1/|rij|^8] [ 1 - 2s^6/|rij|^6 ] rij

  // define auxiliary variables
  double DR[DIM], DRaux[DIM] ; // ri - rj 
  double translation[DIM] ; //  translation vector to go through all periodic contiguous boxes
  double rij2,rij2aux;     // |rij|^2
  double rij2_inv , rij6_inv , rij8_inv ; // 1/|rij|^2, 1/|rij|^6 ,  1/|rij|^8, Auxiliary optimization variables
  int i,n,pi,pj,pk;  // counters

  for(i=0;i<DIM;i++)  f[i]=0.0 ;   // reset force 
  
  for(n=0 ; n< Natoms ; n++ ) // loop over all atoms
    if(n!=ai) { // avoid interaction with itself

      // compute DR vector and rij squared (Note that smallest distance could be through the periodic boundary, this is checked by computing distances through the 27 contiguous boxes and taking the smallest )
      rij2=1e6; // ridiculous large distance to be minimized
	  for(pi=-1;pi<=1;pi++)
	    for(pj=-1;pj<=1;pj++)
	      for(pk=-1;pk<=1;pk++){ // loop over 27 contigous boxes
		translation[0]=pi*BOX[0]; translation[1]=pj*BOX[1]; translation[2]=pk*BOX[2];
		rij2aux=0;
		for(i=0;i<DIM;i++) {
		  DRaux[i]=atom[n].getr(i)+translation[i] - r[i];
		  rij2aux+=DRaux[i]*DRaux[i];
		}

		if(rij2aux<rij2){ // take the minimum distance 
		  for(i=0;i<DIM;i++) DR[i]=DRaux[i];
		  rij2=rij2aux;
		}
	      }
       // continue only if distance is within cutoff
      if (rij2< cutoff2 ) {    
	// compute 1/|rij|^6 and 1/|rij|^8
	rij2_inv=1.0/rij2; 
	rij6_inv= rij2_inv * rij2_inv * rij2_inv ;
	rij8_inv= rij2_inv * rij2_inv * rij2_inv * rij2_inv ;
	// update force and acceleration
	for(i=0;i<DIM;i++) {
	  f[i]+= e_s6_24 * rij8_inv * ( 1.0 - 2*s6*rij6_inv ) * DR[i];
	}
      }
    }
}

// Move with Leap-frog taking PBC
void Catom:: move(double dt, double *BOX ){
  int i ;
  for(i=0;i<DIM;i++) {
    v_shifted_half[i]+=dt*m_inv*f[i]; // V(t+dt/2)=v(t-dt/2) + dt*F(t)/m. 
    r[i]+=dt*v_shifted_half[i]; // r(t+dt)=r(t)+dt*v(t+dt/2)

    if(r[i]>BOX[i]) r[i]-=BOX[i] ; // PBC right
    if(r[i]<0) r[i]+=BOX[i] ;  // PBC left
  }
}

using namespace std;
int main(){             
  // Parameter declaration
  int N=125; // Number of atoms
  Catom gas[N]; // an array of N atoms called gas
  double m=39.948 , s=0.340100  , e=0.978638 ; // constants for Argon m[amu] , sigma[nm] , epsilon[kJ/mol]
  double cutoff=1.0  , cutoff2=cutoff*cutoff ; //nm and nm^2,respectively
  double dt=0.001 ; //(ps) it yielded a reasonable speed distribution after 1e4 steps
  int nsteps=1e4 , out_freq=1e1 ;  // MD will run for nsteps printing positions every out_freq steps 
  double BOX[DIM]; for(int i=0;i<DIM;i++) BOX[i]=1.8242  ; // Simulation box. The chosen value gives liquid argon density for 125 atoms at 87 K 
  double r0[DIM], v0[DIM], speed0=0.233;  // r0=initial position, v0=initial velocity, speed0= initial speed (same for all atoms. 0.233 nm/ps corresponds to ~ T=87 K for Argon) 
  double theta , phi ; //auxiliary angles for initial velocities
  int n, t , i ;  // auxiliary counters
 
 // ==== initialize ===  
 for(n=0;n<N;n++) {
   r0[0]=(n%5+0.5)*BOX[0]/5  ;  r0[1]= (n%25/5+0.5)*BOX[1]/5;   r0[2]= (int(n/25)+0.5)*BOX[1]/5; // initial positions: ordered array 
   gas[n].initialize(r0,n,m,s,e) ;  // r(t=0), v(t=-dt/2), atom index, mass, sigma, epsilon
 }
 // ==== v(t=0-dt/2) ===
 for(n=0;n<N;n++) gas[n].computeforce(gas, N,BOX, cutoff2) ;
 for(n=0;n<N;n++) {
   theta=drand48()*PI ;  phi=drand48()*2*PI;  // pseudo-random angles thetha [0,180°] and phi [0,360°]
   v0[0]=speed0*sin(theta)*cos(phi);    v0[1]=speed0*sin(theta)*sin(phi);    v0[2]=speed0*cos(phi); // initial velocities (same magnitude, different direction)
   gas[n].v_shifted_minus_half(v0,dt);
 }

 // === DYNAMICS ===
 for (t=0 ; t<nsteps ; t++ ) {
     if( t%out_freq==0  || t==0 ) {    // print position in a pdb format every out_freq steps and at the begining 
       printf( "REMARK    GENERATED BY MD\nTITLE     box t=   %f\nREMARK    THIS IS A SIMULATION BOX\nCRYST1  %7.3f  %7.3f  %7.3f  90.00  90.00  90.00 P 1           1\nMODEL        %d\n", t*dt, BOX[0]*10, BOX[1]*10, BOX[2]*10 , int(t/out_freq) ) ; // pdb header per snapshot
       
       for(n=0 ; n<N ; n++ )  printf("%-6s%5d %4s %3s %1s        %8.3f%8.3f%8.3f%6.2f%6.2f\n", "ATOM ", n, " AR " , " AR", " " , gas[n].getr(0)*10.0 , gas[n].getr(1)*10.0  , gas[n].getr(2)*10.0  , 0.0 , 0.0) ; // coordinates in angstrom per snapshot
       
       printf("TER\nENDMDL\n"); // pdb footer per snapshot
     }  

     for(n=0;n<N;n++) gas[n].computeforce(gas, N,BOX, cutoff2) ; 
     for(n=0;n<N;n++) gas[n].move(dt, BOX) ; 
 } // end of loop over nsteps. END OF DYNAMICS!

 // final velocity distribution ( at tend-dt*0.5)
 double Vnorm, av_kin=0 ; // norm of velocity and  average kinetic energy
 for(n=0 ; n<N ; n++ ) {
   Vnorm= sqrt(gas[n].getv(0)*gas[n].getv(0)+gas[n].getv(1)*gas[n].getv(1)+gas[n].getv(2)*gas[n].getv(2));
   av_kin+=0.5*m*Vnorm*Vnorm ;
   cout<<"REMARK |V|= "<<Vnorm <<endl;
 }
 cout<<"REMARK <K>= "<<av_kin/N<<endl; 
  return 0;
}
