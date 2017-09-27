#include <iostream>
#include <cmath>

using namespace std;

int main(){

  static int N=10;
  double t=0.;
  double dt=0.000001;
  double x[N];
  double x0[N];
  double vx[N];
  double X[N];

  //Initial positions
  double pos=0.;
  for(int i=0;i<N;i++){
    x0[i]=pos;
    x[i]=pos;//in units of the interparticle spacing
    vx[i]=0;//in units of the interparticle spacing
    X[i]=0.;
    pos+=2*pow(10,-5);
  }

  //vx[0]=1000;

  double wp=0.05/dt;//3.98*pow(10,10);
  double cwp=cos(wp);
  double swp = sin(wp);

  /*
  x[0]=x[0]+vx[0]*swp*dt-X[0]*(1-cwp*dt);
  X[0]=x[0]-x0[0];
  cout << x[0]  << "   " << X[0] << endl;

  x[0]=x[0]+vx[0]*swp*dt-X[0]*(1-cwp*dt);
  X[0]=x[0]-x0[0];
  cout << x[0]  << "   " << X[0] << endl;

  x[0]=x[0]+vx[0]*swp*dt-X[0]*(1-cwp*dt);
  X[0]=x[0]-x0[0];
  cout << x[0]  << "   " << X[0] << endl;
  */

  
  //Next positions
  bool out=false;
  while(out==false){
    for(int i=0;i<N;i++){//estou so a por ate N-1 por causa do if que esta ali em baixo mas esta mal
      x[i]=x[i]+vx[i]*swp*dt-X[i]*(1-cwp*dt);
      vx[i]=vx[i]*cwp*dt-wp*X[i]*swp*dt;
      X[i]=x[i]-x0[i];

      cout << x[0] << "   " << x[1] << endl;

      if(i>0){
	if(x[i-1]>x[i]){
	  out=true;
	  cout << "crossing: " << i << endl;
	}
      }

    }
    t+=dt;
  }
  

  cout << "time: " << t << endl;



  return 0;

}
