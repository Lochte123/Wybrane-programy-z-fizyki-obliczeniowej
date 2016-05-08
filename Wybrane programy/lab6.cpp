#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <array>
using namespace std;
//autor: £ukasz Kokosza
//-std=c++11 do kompilacji
double h=1; double m=0.5;
double tolerancja=pow(10,-10);
double w=0.15;
#define x0 20
#define y0 20
double dx=0.8; double dy=0.8;

struct funkcje
{
       double zBiegunami(int x,int y, double B);
};

double funkcje::zBiegunami(int xi,int yi,double B)
{
       double phi;
       double x=xi*dx;
       double y=yi*dy;
       if(x-x0>tolerancja && y-y0>-tolerancja)
       {
             phi=atan((y-y0)/(x-x0));
       }
       else if(x-x0>tolerancja && y-y0<-tolerancja){
       phi=atan((y-y0)/(x-x0))+2*M_PI;
       }
       else if(x-x0<-tolerancja){
       phi=atan((y-y0)/(x-x0))+M_PI;
       }
       else if(fabs(x-x0)<tolerancja && y-y0>tolerancja){
       phi=M_PI/2;
       }
       else if(fabs(x-x0)<tolerancja && y-y0<-tolerancja){
       phi=3*M_PI/2;
       }
       else{
       phi=0;
       }
       return 0.5*m*w*w*((x-x0)*(x-x0)+(y-y0)*(y-y0))*(1+(2./7.)*cos(3*phi)); 
}

struct rownanie_Schrodingera
{
	void stanPodstawowy(int x,int y);
	void stanyWzbudzone(int x,int y);
};

void rownanie_Schrodingera::stanPodstawowy(int x,int y)
{
        funkcje f;
	//double dx,dy;
	double k,talk,E;
	k=0; talk=0; E=0;
	double norma=0;
	double **psi=new double*[x];
	double **Hpsi=new double*[x];
	double **V=new double*[51];
	
	for(int i=0; i<x; i++)
	{
		psi[i]=new double[y];
		Hpsi[i]=new double[y];
	}
	for(int i=0; i<51; i++)
	{
		V[i]=new double[51];
	}

	for(int i=0; i<51; i++)
	{
		for(int j=0; j<51; j++)
		{
			V[i][j]=f.zBiegunami(i,j,4);
		}
	}

	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			if(i==0 || j==0 || i==50 || j==50)
			{
				psi[i][j]=0;
				Hpsi[i][j]=0;
			}
			else
			{
				psi[i][j]=rand();
				Hpsi[i][j]=rand();
			}
		}
	}
	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			norma+=pow(fabs(psi[i][j]),2)*dx*dy;//normalizacja
		}
	}
	cout<<norma<<endl;
	norma=sqrt(norma);
	cout<<norma<<endl;
	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			psi[i][j]=psi[i][j]/norma;
		}
	}
	//Hpsi
	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			Hpsi[i][j]=0;
		}
	}
	for(int i=1; i<x-1; i++)
	{
		for(int j=1; j<y-1; j++)
		{
			Hpsi[i][j]=(-(h*h)/(2*m))*((psi[i+1][j]+psi[i-1][j]-2*psi[i][j])/pow(dx,2))+(-(h*h)/(2*m))*((psi[i][j+1]+psi[i][j-1]-2*psi[i][j])/pow(dy,2))+(V[i][j]*psi[i][j]);
		}
	}
	//energie
	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			E+=psi[i][j]*Hpsi[i][j]*dx*dy;
		}
	}
	cout<<E<<endl;

	double E2=0;
	talk=0;
	double dtalk=0.005;
	int iter=0;
	fstream plik1,plik0;
	plik1.open("zad1.dat",ios::out | ios::trunc);
	if(plik1.good()==true)
	{
	while(fabs(E-E2)>tolerancja)
	{
		talk+=dtalk;
		cout<<iter<<endl;
		norma=0;
		iter++;
		E2=E;
		for(int i=0; i<x; i++)
		{
			for(int j=0; j<y; j++)
			{
				psi[i][j]=psi[i][j]-(dtalk/h)*Hpsi[i][j];
			}
		}
		for(int i=0; i<x; i++)
		{
			for(int j=0; j<y; j++)
			{
				norma+=pow(fabs(psi[i][j]),2)*dx*dy;//normalizacja
			}
		}
		norma=sqrt(norma);
		for(int i=0; i<x; i++)
		{
			for(int j=0; j<y; j++)
			{
				psi[i][j]=psi[i][j]/norma;
			}
		}
		for(int i=1; i<x-1; i++)
		{
			for(int j=1; j<y-1; j++)
			{
				Hpsi[i][j]=(-(h*h)/(2*m))*((psi[i+1][j]+psi[i-1][j]-2*psi[i][j])/pow(dx,2))+(-(h*h)/(2*m))*((psi[i][j+1]+psi[i][j-1]-2*psi[i][j])/pow(dy,2))+(V[i][j]*psi[i][j]);
			}
		}
		E=0;
		for(int i=0; i<x; i++)
		{
			for(int j=0; j<y; j++)
			{
				E+=psi[i][j]*Hpsi[i][j]*dx*dy;
			}
		}
		plik1<<E<<" "<<iter<<endl;
	}
	}
	else{ cout<<"BRAK DOSTEPU DO PLIKU!"<<endl;}
plik1.close();
norma=0;
	for(int i=0; i<41; i++)
	{
		for(int j=0; j<41; j++)
		{
			V[i][j]=f.zBiegunami(i,j,4);
		}
	}

	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			if(i==0 || j==0 || i==50 || j==50)
			{
				psi[i][j]=0;
				Hpsi[i][j]=0;
			}
			else
			{
				psi[i][j]=rand();
				Hpsi[i][j]=rand();
			}
		}
	}
	
	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			norma+=pow(fabs(psi[i][j]),2)*dx*dy;//normalizacja
		}
	}
	cout<<norma<<endl;
	norma=sqrt(norma);
	cout<<norma<<endl;
	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			psi[i][j]=psi[i][j]/norma;
		}
	}
	//Hpsi
	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			Hpsi[i][j]=0;
		}
	}
	for(int i=1; i<x-1; i++)
	{
		for(int j=1; j<y-1; j++)
		{
			Hpsi[i][j]=(-(h*h)/(2*m))*((psi[i+1][j]+psi[i-1][j]-2*psi[i][j])/pow(dx,2))+(-(h*h)/(2*m))*((psi[i][j+1]+psi[i][j-1]-2*psi[i][j])/pow(dy,2))+(V[i][j]*psi[i][j]);
		}
	}
	//energie
	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			E+=psi[i][j]*Hpsi[i][j]*dx*dy;
		}
	}
	E2=0;
	talk=0;
	dtalk=0.05;
	iter=0;
	fstream plik2;
	plik2.open("zad1b.dat",ios::out | ios::trunc);
	if(plik2.good()==true)
	{
	while(fabs(E-E2)>tolerancja)
	{
		talk+=dtalk;
		cout<<iter<<endl;
		norma=0;
		E2=E;
		iter++;
		for(int i=0; i<x; i++)
		{
			for(int j=0; j<y; j++)
			{
				psi[i][j]=psi[i][j]-(dtalk/h)*Hpsi[i][j];
			}
		}
		for(int i=0; i<x; i++)
		{
			for(int j=0; j<y; j++)
			{
				norma+=pow(fabs(psi[i][j]),2)*dx*dy;//normalizacja
			}
		}
		norma=sqrt(norma);
		for(int i=0; i<x; i++)
		{
			for(int j=0; j<y; j++)
			{
				psi[i][j]=psi[i][j]/norma;
			}
		}
		for(int i=1; i<x-1; i++)
		{
			for(int j=1; j<y-1; j++)
			{
				Hpsi[i][j]=(-(h*h)/(2*m))*((psi[i+1][j]+psi[i-1][j]-2*psi[i][j])/pow(dx,2))+(-(h*h)/(2*m))*((psi[i][j+1]+psi[i][j-1]-2*psi[i][j])/pow(dy,2))+(V[i][j]*psi[i][j]);
			}
		}
		E=0;
		for(int i=0; i<x; i++)
		{
			for(int j=0; j<y; j++)
			{
				E+=psi[i][j]*Hpsi[i][j]*dx*dy;
			}
		}
		plik2<<E<<" "<<iter<<endl;
	}
	}
	else{ cout<<"BRAK DOSTEPU DO PLIKU!"<<endl;}
plik2.close();

norma=0;
	for(int i=0; i<41; i++)
	{
		for(int j=0; j<41; j++)
		{
			V[i][j]=f.zBiegunami(i,j,4);
		}
	}

	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			if(i==0 || j==0 || i==50 || j==50)
			{
				psi[i][j]=0;
				Hpsi[i][j]=0;
			}
			else
			{
				psi[i][j]=rand();
				Hpsi[i][j]=rand();
			}
		}
	}
	
	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			norma+=pow(fabs(psi[i][j]),2)*dx*dy;//normalizacja
		}
	}
	cout<<norma<<endl;
	norma=sqrt(norma);
	cout<<norma<<endl;
	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			psi[i][j]=psi[i][j]/norma;
		}
	}
	//Hpsi
	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			Hpsi[i][j]=0;
		}
	}
	for(int i=1; i<x-1; i++)
	{
		for(int j=1; j<y-1; j++)
		{
			Hpsi[i][j]=(-(h*h)/(2*m))*((psi[i+1][j]+psi[i-1][j]-2*psi[i][j])/pow(dx,2))+(-(h*h)/(2*m))*((psi[i][j+1]+psi[i][j-1]-2*psi[i][j])/pow(dy,2))+(V[i][j]*psi[i][j]);
		}
	}
	//energie
	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			E+=psi[i][j]*Hpsi[i][j]*dx*dy;
		}
	}
	E2=0;
	talk=0;
	dtalk=0.2;
	iter=0;
	fstream plik;
	plik.open("zad1c.dat",ios::out | ios::trunc);
	if(plik.good()==true)
	{
	while(fabs(E-E2)>tolerancja)
	{
		talk+=dtalk;
		cout<<iter<<endl;
		norma=0;
		E2=E;
		iter++;
		for(int i=0; i<x; i++)
		{
			for(int j=0; j<y; j++)
			{
				psi[i][j]=psi[i][j]-(dtalk/h)*Hpsi[i][j];
			}
		}
		for(int i=0; i<x; i++)
		{
			for(int j=0; j<y; j++)
			{
				norma+=pow(fabs(psi[i][j]),2)*dx*dy;//normalizacja
			}
		}
		norma=sqrt(norma);
		for(int i=0; i<x; i++)
		{
			for(int j=0; j<y; j++)
			{
				psi[i][j]=psi[i][j]/norma;
			}
		}
		for(int i=1; i<x-1; i++)
		{
			for(int j=1; j<y-1; j++)
			{
				Hpsi[i][j]=(-(h*h)/(2*m))*((psi[i+1][j]+psi[i-1][j]-2*psi[i][j])/pow(dx,2))+(-(h*h)/(2*m))*((psi[i][j+1]+psi[i][j-1]-2*psi[i][j])/pow(dy,2))+(V[i][j]*psi[i][j]);
			}
		}
		E=0;
		for(int i=0; i<x; i++)
		{
			for(int j=0; j<y; j++)
			{
				E+=psi[i][j]*Hpsi[i][j]*dx*dy;
			}
		}
		plik<<E<<" "<<iter<<endl;
	}
	}
	else{ cout<<"BRAK DOSTEPU DO PLIKU!"<<endl;}
plik.close();
}

void rownanie_Schrodingera::stanyWzbudzone(int x,int y)
{
        funkcje f;
	double k,talk,E;
	k=0; talk=0; E=0;
	double norma=0;

	double ***psi=new double**[x];
	double **Hpsi=new double*[x];
	double **V=new double*[51];

	typedef array<double,6> iloczyn;
	iloczyn skalarny;

	for(int i=0; i<x; i++)
	{
		psi[i]=new double*[y];
		Hpsi[i]=new double[y];
	}
	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			psi[i][j]=new double[6];//6stanow
		}
	}
	for(int i=0; i<51; i++)
	{
		V[i]=new double[51];
	}

	for(int i=0; i<51; i++)
	{
		for(int j=0; j<51; j++)
		{
			V[i][j]=f.zBiegunami(i,j,4);
		}
	}

	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			if(i==0 || j==0 || i==50 || j==50)
			{
				psi[i][j][0]=0;
				Hpsi[i][j]=0;
			}
			else
			{
				psi[i][j][0]=rand();
				Hpsi[i][j]=rand();
			}
		}
	}
	
	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			norma+=pow(fabs(psi[i][j][0]),2)*dx*dy;//normalizacja
		}
	}
	cout<<norma<<endl;
	norma=sqrt(norma);
	cout<<norma<<endl;
	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			psi[i][j][0]=psi[i][j][0]/norma;
		}
	}
	//Hpsi
	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			Hpsi[i][j]=0;
		}
	}
	for(int i=1; i<x-1; i++)
	{
		for(int j=1; j<y-1; j++)
		{
			Hpsi[i][j]=(-(h*h)/(2*m))*((psi[i+1][j][0]+psi[i-1][j][0]-2*psi[i][j][0])/pow(dx,2))+(-(h*h)/(2*m))*((psi[i][j+1][0]+psi[i][j-1][0]-2*psi[i][j][0])/pow(dy,2))+(V[i][j]*psi[i][j][0]);
		}
	}
	//energie
	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			E+=psi[i][j][0]*Hpsi[i][j]*dx*dy;
		}
	}
	cout<<E<<endl;
	double E2=0;
	talk=0;
	double dtalk=0.05;
	int iter=0;
	fstream plik3;
	plik3.open("zad2a.dat",ios::out | ios::trunc);
	if(plik3.good()==true)
	{
	while(fabs(E-E2)>tolerancja)
	{
		talk+=dtalk;
		cout<<iter<<endl;
		norma=0;
		iter++;
		for(int i=0; i<x; i++)
		{
			for(int j=0; j<y; j++)
			{
				psi[i][j][0]=psi[i][j][0]-(dtalk/h)*Hpsi[i][j];
			}
		}
		for(int i=0; i<x; i++)
		{
			for(int j=0; j<y; j++)
			{
				norma+=pow(fabs(psi[i][j][0]),2)*dx*dy;//normalizacja
			}
		}
		norma=sqrt(norma);
		for(int i=0; i<x; i++)
		{
			for(int j=0; j<y; j++)
			{
				psi[i][j][0]=psi[i][j][0]/norma;
			}
		}
		for(int i=1; i<x-1; i++)
		{
			for(int j=1; j<y-1; j++)
			{
				Hpsi[i][j]=(-(h*h)/(2*m))*((psi[i+1][j][0]+psi[i-1][j][0]-2*psi[i][j][0])/pow(dx,2))+(-(h*h)/(2*m))*((psi[i][j+1][0]+psi[i][j-1][0]-2*psi[i][j][0])/pow(dy,2))+(V[i][j]*psi[i][j][0]);
			}
		}
		E2=E;
		E=0;
		for(int i=0; i<x; i++)
		{
			for(int j=0; j<y; j++)
			{
				E+=psi[i][j][0]*Hpsi[i][j]*dx*dy;
			}
		}
	}
	}
	else{ cout<<"BRAK DOSTEPU DO PLIKU!"<<endl;}
	for(int i=0; i<x; i++)
	{
	     for(int j=0; j<y; j++)
	     {
	          plik3<<i<<" "<<j<<" "<<psi[i][j][0]<<endl;
	     }
	     plik3<<endl;
	}
plik3.close();

//
//
//

	fstream plik4,plik5,plik6,plik7,plik8;
	plik4.open("zad2b.dat",ios::out | ios::trunc);
	plik5.open("zad2c.dat",ios::out | ios::trunc);
	plik6.open("zad2d.dat",ios::out | ios::trunc);
	plik7.open("zad2e.dat",ios::out | ios::trunc);
	plik8.open("zad2f.dat",ios::out | ios::trunc);
	for(int i=0; i<51; i++)
	{
		for(int j=0; j<51; j++)
		{
			V[i][j]=f.zBiegunami(i,j,4);
		}
	}
for(int n=1; n<6; n++)
{
	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			if(i==0 || j==0 || i==50 || j==50)
			{
				psi[i][j][n]=0;
				Hpsi[i][j]=0;
			}
			else
			{
				psi[i][j][n]=rand();
				Hpsi[i][j]=rand();
			}
		}
	}
	norma=0;
	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			norma+=pow(fabs(psi[i][j][n]),2)*dx*dy;//normalizacja
		}
	}
	cout<<norma<<endl;
	norma=sqrt(norma);
	cout<<norma<<endl;
	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			psi[i][j][n]=psi[i][j][n]/norma;
		}
	}
	//Hpsi
	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			Hpsi[i][j]=0;
		}
	}
	for(int i=1; i<x-1; i++)
	{
		for(int j=1; j<y-1; j++)
		{
			Hpsi[i][j]=(-(h*h)/(2*m))*((psi[i+1][j][n]+psi[i-1][j][n]-2*psi[i][j][n])/pow(dx,2))+(-(h*h)/(2*m))*((psi[i][j+1][n]+psi[i][j-1][n]-2*psi[i][j][n])/pow(dy,2))+(V[i][j]*psi[i][j][n]);
		}
	}
	//energie
	for(int i=0; i<x; i++)
	{
		for(int j=0; j<y; j++)
		{
			E+=psi[i][j][n]*Hpsi[i][j]*dx*dy;
		}
	}
	cout<<E<<endl;

	double E2=0;
	talk=0;
	double suma=0;
	double dtalk=0.05;
	int iter=0;

	while(fabs(E-E2)>tolerancja)
	{
		talk+=dtalk;
		//cout<<iter<<endl;
		norma=0;
		iter++;
		for(int i=0; i<x; i++)
		{
			for(int j=0; j<y; j++)
			{
				psi[i][j][n]=psi[i][j][n]-(dtalk/h)*Hpsi[i][j];
			}
		}
		for(int k=0; k<n; k++)
		{
		suma=0;
		for(int i=0; i<x; i++)
		{
			for(int j=0; j<y; j++)
			{
				suma+=psi[i][j][n]*psi[i][j][k];
			}
		}
				skalarny[k]=suma*dx*dy;
				//cout<<skalarny[k]<<endl;
		}
		
                for(int k=0; k<n; k++){
		for(int i=1; i<x-1; i++)
		{
			for(int j=1; j<y-1; j++)
			{
				psi[i][j][n]=psi[i][j][n]-skalarny[k]*psi[i][j][k];
			}
		}
		}
		for(int i=0; i<x; i++)
		{
			for(int j=0; j<y; j++)
			{
				norma+=pow(fabs(psi[i][j][n]),2)*dx*dy;//normalizacja
			}
		}
		norma=sqrt(norma);
		for(int i=0; i<x; i++)
		{
			for(int j=0; j<y; j++)
			{
				psi[i][j][n]=psi[i][j][n]/norma;
			}
		}
		for(int i=1; i<x-1; i++)
		{
			for(int j=1; j<y-1; j++)
			{
				Hpsi[i][j]=(-(h*h)/(2*m))*((psi[i+1][j][n]+psi[i-1][j][n]-2*psi[i][j][n])/pow(dx,2))+(-(h*h)/(2*m))*((psi[i][j+1][n]+psi[i][j-1][n]-2*psi[i][j][n])/pow(dy,2))+(V[i][j]*psi[i][j][n]);
			}
		}
		E2=E;
		E=0;
		for(int i=0; i<x; i++)
		{
			for(int j=0; j<y; j++)
			{
				E+=psi[i][j][n]*Hpsi[i][j]*dx*dy;
			}
		}
	}
	cout<<"n="<<n<<"E="<<E<<endl;
	}
	for(int i=0; i<x; i++)
	{
	     for(int j=0; j<y; j++)
	     {
	          plik4<<i<<" "<<j<<" "<<psi[i][j][1]<<endl;
	     }
	     plik4<<endl;
	}
	for(int i=0; i<x; i++)
	{
	     for(int j=0; j<y; j++)
	     {
	          plik5<<i<<" "<<j<<" "<<psi[i][j][2]<<endl;
	     }
	     plik5<<endl;
	}
	for(int i=0; i<x; i++)
	{
	     for(int j=0; j<y; j++)
	     {
	          plik6<<i<<" "<<j<<" "<<psi[i][j][3]<<endl;
	     }
	     plik6<<endl;
	}
	for(int i=0; i<x; i++)
	{
	     for(int j=0; j<y; j++)
	     {
	          plik7<<i<<" "<<j<<" "<<psi[i][j][4]<<endl;
	     }
	     plik7<<endl;
	}
	for(int i=0; i<x; i++)
	{
	     for(int j=0; j<y; j++)
	     {
	          plik8<<i<<" "<<j<<" "<<psi[i][j][5]<<endl;
	     }
	     plik8<<endl;
	}
plik4.close();
plik5.close();
plik6.close();
plik7.close();
plik8.close();
}
int main()
{
	rownanie_Schrodingera r;
	r.stanPodstawowy(51,51);
	r.stanyWzbudzone(51,51);
	return 0;
}
