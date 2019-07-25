#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"

using namespace std;

const int N=48;

const double Zeta  =0.1786178958448091;
const double Lambda=-0.2123418310626054;
const double Xi    =-0.06626458266981849;
const double uno_m2LAMBDA_S2=(1-2*Lambda)/2;
const double uno_m2_XIplusCHI=(1-2*(Xi+Zeta));

class Cuerpo;
class Colisionador;

//----------------------Clase Cuerpo---------------
class Cuerpo{

private:
	vector3D r, V, F;
	double m, q,R;

public:
  void Inicie(double x0, double y0,double Vx0,double Vy0, double q0, double m0, double R0);
  void BorreFuerza(void);
  void AgregueFuerza(vector3D F0);
  void Muevase_r(double dt, double Constante);
  void Muevase_V(double dt, double Constante);
  void Dibujese(void);
  double Getx(void){return r.x();}; //Función Inline
  double Gety(void){return r.y();}; //Función Inline
  
  double Getq(void){return q;}; //Función Inline
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0, double y0, double Vx0, double Vy0, double q0, double m0, double R0){
  r.cargue(x0, y0, 0); V.cargue(Vx0, Vy0, 0); q=q0; m=m0; R=R0;
}
void Cuerpo::BorreFuerza(void){
	F.cargue(0,0,0);
}

void Cuerpo::AgregueFuerza(vector3D F0){
	F+=F0;
}

void Cuerpo::Muevase_r(double dt, double Constante){
	r+=V*(Constante*dt);
}

void Cuerpo::Muevase_V(double dt, double Constante){
	V+=F*(Constante*dt/m);
}

void Cuerpo::Dibujese(void){
  cout<<", "<< Getx() <<"+"<<R<<"*cos(t),"<< Gety() <<"+"<<R<<"*sin(t)";
}
//-----------------------Clase Colisionador--------------

class Colisionador{
	private:
		double K, Kc, ga;
	public:
		void Inicie(double k0, double Kcentral, double Gamma);
		void CalculeTodasLasFuerzas(Cuerpo * Carga);
		void CalculeLaFuerzaEntre(Cuerpo & Carga1, Cuerpo & Carga2);
	
};

void Colisionador::Inicie(double k0, double Kcentral, double Gamma){
	K=k0;  Kc= Kcentral; ga=Gamma;
	}

void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Carga){
	int i,j;
	
	for(i=0;i<N;i++)  Carga[i].BorreFuerza();
	
	for(i=0;i<N;i++)
		for(j=0;j<i;j++)
			CalculeLaFuerzaEntre(Carga[i],Carga[j]);
	
	for(i=0;i<N;i++){
		Carga[i].AgregueFuerza(-Kc*Carga[i].r); //Fuerza central
		Carga[i].AgregueFuerza(-ga*Carga[i].V); //Fuerza Viscosa
		}
	
	}
void Colisionador::CalculeLaFuerzaEntre(Cuerpo & Cargai, Cuerpo & Cargaj){
	vector3D F1, rij=Cargaj.r-Cargai.r;
	double aux=K*Cargai.Getq()*Cargaj.Getq()*pow(norma2(rij),-1.5);
	F1=rij*aux;		Cargaj.AgregueFuerza(F1);		Cargai.AgregueFuerza((-1)*F1);
	}

//----------------------Funciones Globales---------------
/*
void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl; 
  cout<<"set output 'Carga.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange [-110:110]"<<endl;
  cout<<"set yrange [-110:110]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
}
void TermineCuadro(void){
    cout<<endl;
}
*/
//--------------------Programa Principal------------------
int main(void){
	
	Colisionador Newton;
	double k0=1000, Kcentral=1, Gamma=0.5;
	Newton.Inicie(k0, Kcentral, Gamma);
	
	Cuerpo Carga[N];
	
	//CARGAS.
	double x0, y0, Vx0, Vy0, q0 = 1, m0=1,  R0=3.0, Ran=40;  int i;
	
	//GENERADOR DE VELOCIDADES ALEATORIAS
	Crandom ran64(1);  double Vmag=1, theta;
	
	
	//GENERA CONDICIONES INICIALES
	for(i=0;i<N;i++){
		x0 = Ran*cos(i*(2.0*M_PI)/N);  y0 = Ran*sin(i*(2.0*M_PI)/N);
		theta=2*M_PI*ran64.r();  Vx0 = Vmag*cos(theta);  Vy0 = Vmag*sin(theta);
		Carga[i].Inicie(x0, y0, Vx0, Vy0, q0, m0, R0);
		}
	
	//InicieAnimacion();
	
	//---PARAMETROS-DE-TIEMPO-Y-ANIMACION---
	double T=2*M_PI*sqrt(m0/Kcentral), t, dt=0.001, tmax=20*T;
	double tdibujo; int ndibujos=100;
	
	for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
		/*
		if(tdibujo>tmax/ndibujos){
			InicieCuadro();
			for(i=0; i<N; i++) Carga[i].Dibujese();
			TermineCuadro();
			tdibujo=0;
		}
		*/
		
		//cout<< t << " " << Carga[0].Getx() << endl;
		
		for(i=0;i<N;i++) Carga[i].Muevase_r(dt, Zeta);
		Newton.CalculeTodasLasFuerzas(Carga);	for(i=0;i<N;i++) Carga[i].Muevase_V(dt, uno_m2LAMBDA_S2);
		for(i=0;i<N;i++) Carga[i].Muevase_r(dt, Xi);
		Newton.CalculeTodasLasFuerzas(Carga); for(i=0;i<N;i++) Carga[i].Muevase_V(dt, Lambda);
		for(i=0;i<N;i++) Carga[i].Muevase_r(dt, uno_m2_XIplusCHI);
		Newton.CalculeTodasLasFuerzas(Carga);	for(i=0;i<N;i++) Carga[i].Muevase_V(dt, Lambda);
		for(i=0;i<N;i++) Carga[i].Muevase_r(dt, Xi);
		Newton.CalculeTodasLasFuerzas(Carga);	for(i=0;i<N;i++) Carga[i].Muevase_V(dt, uno_m2LAMBDA_S2);
		for(i=0;i<N;i++) Carga[i].Muevase_r(dt, Zeta);
	}
	
	//DISTANCIA AL ORIGEN
	double Distancia;
	double Distancia_desp = 0;
	for(i=0;i<N;i++){
		Distancia = sqrt(Carga[i].Getx()*Carga[i].Getx() + Carga[i].Gety()*Carga[i].Gety());
		cout << Distancia << endl;
		
		if(Distancia > Distancia_desp){
			Distancia_desp = Distancia;
			}
		}
	cout << "El mayor valor es " << Distancia_desp << endl;
	
	return 0;
}

//Grafica con Gnuplot (1.pdf):
/*
gnuplot> set title 'POSICION DE UNA PARTICULA COMO FUNCION DEL TIEMPO CON FUERZA DISIPATIVA'
gnuplot> plot 'datos.dat' w l
gnuplot> set xlabel 't'
gnuplot> set ylabel 'x'
gnuplot> set grid
gnuplot> unset key
gnuplot> set terminal pdf
gnuplot> set output '2.pdf'
gnuplot> replot 
gnuplot> quit
*/

//Para generar el gif, escribir en consola:
/*
g++ Cristales_de_wigner_1.cpp
./a.out | gnuplot
*/
