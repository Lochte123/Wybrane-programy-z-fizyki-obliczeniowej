#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
//autor: Łukasz Kokosza
struct zerowanie
{
	int licznik1,licznik2,t,spr;
	double h;
};

struct dyfuzja
{
	void izolacja();
	void nieszczelnosc();
};

void dyfuzja::izolacja()
{
	int imin,i1,i2,i3,imax,jmin,j1,j2,j3,j4,jmax;
	imin=0;
	i1=20;
	i2=60;
	i3=80;
	imax=90;
	jmin=0;
	j1=10;
	j2=20;
	j3=40;
	j4=50;
	jmax=60;
	double dx,dy,dt;  dx=1;  dy=dx;  dt=2.5;
	double p,c,k;	k=1; 	p=k;	c=p;
	double T0,T1,T2;
	double zewnatrz1=0,wewnatrz1=0,zewnatrz2=0,wewnatrz2=0;
	double tolerancja=pow(10,-5);
	T0=0; T1=20; T2=25;
	zerowanie z;
	z.licznik1=0;
	z.licznik2=0;
	z.spr=0;
	z.t=0;
	z.h=0.0;
	double **tab1=new double*[91];
	double **tab2=new double*[91];
	for(int i=imin; i<101; i++)
	{
		tab1[i]=new double[61];
		tab2[i]=new double[61];
	}
	for(int i=imin; i<91; i++)
	{
		for(int j=jmin; j<61; j++)
		{
			tab1[i][j]=T2;
			tab2[i][j]=T2;
			if(i>=i1 && i<=i2 && j==jmin) {tab2[i][j]=T0; tab1[i][j]=T0;}
			if(i>=i1 && i<=i2 && j==jmax) {tab2[i][j]=T0; tab1[i][j]=T0;}
			if(i==imax && j>=j3 && j<=jmax) {tab1[i][j]=T1; tab2[i][j]=T1;}
			if(i==imin && j>=j2 && j<=j4) {tab1[i][j]=T2; tab2[i][j]=T2;}
			if(i==imax && j>=jmin && j<=j3) {tab1[i][j]=T2; tab2[i][j]=T2;}
		}
	}
	fstream plik2;
	plik2.open("zad1inne.dat",ios::out | ios::trunc);
	if(plik2.good()==true)
	{
		while(z.licznik2<=6000)//zewnetrzna petla dla badania zbieznosci
		{
			z.licznik1=0;
			zewnatrz1=zewnatrz2;
			zewnatrz2=0;
			while(z.licznik1<=50)
			{
				z.licznik1++;
				//cout<<"iteracja nr.: "<<z.licznik1<<endl;
				for(int i=imin; i<=imax; i++)
				{
					for(int j=jmin; j<=jmax; j++)
					{
						if(i>=i1 && i<=i2 && j==jmin) {tab2[i][j]=T0; tab1[i][j]=T0;}
						if(i>=i1 && i<=i2 && j==jmax) {tab2[i][j]=T0; tab1[i][j]=T0;}
						if(i==imax && j>=j3 && j<=jmax) {tab1[i][j]=T1; tab1[i][j]=T1;}
						if(i==imin && j>=j2 && j<=j4) {tab1[i][j]=T2; tab2[i][j]=T2;}
						if(i==imax && j>=jmin && j<=j3) {tab1[i][j]=T2; tab2[i][j]=T2;}
						//KRAWĘDZIE
						if(i==imin && j>=jmin && j<j2) tab2[i][j]=(((k/dx)*tab2[i+1][j])+z.h*T0)/(z.h+(k/dx));//PIONOWA LEWA DOLNA
						if(i==imin && j>j4 && j<jmax) tab2[i][j]=(((k/dx)*tab2[i+1][j])+z.h*T0)/(z.h+(k/dx));//PIONOWA LEWA GORNA
						if(j==jmin && i>imin && i<i1) tab2[i][j]=(((k/dx)*tab2[i][j+1])+z.h*T0)/(z.h+(k/dx));//POZIOMA LEWA DOLNA
						if(j==jmin && i>i2 && i<imax) tab2[i][j]=(((k/dx)*tab2[i][j+1])+z.h*T0)/(z.h+(k/dx));//POZIOMA PRAWA DOLNA
						if(j==jmax && i>imin && i<i1)  tab2[i][j]=(((k/dx)*tab2[i][j-1])+z.h*T0)/(z.h+(k/dx));//POZIOMA LEWA GORNA
						if(j==jmax && i>i2 && i<imax)  tab2[i][j]=(((k/dx)*tab2[i][j-1])+z.h*T0)/(z.h+(k/dx));//POZIOMA PRAWA GORNA
						//KANTY
						if(i==imin && j==jmin) tab2[i][j]=(((k/(sqrt(2)*dx))*tab2[i+1][j+1])+z.h*T0)/(z.h+(k/sqrt(2)*dx));//1 KANT
						if(i==imin && j==jmax) tab2[i][j]=(((k/(sqrt(2)*dx))*tab2[i+1][j-1])+z.h*T0)/(z.h+(k/sqrt(2)*dx));//2 KANT
						if(i==imax && j==j3) tab2[i][j]=T1;
					}
						
				}
				zewnatrz1=zewnatrz2;
				zewnatrz2=0;
				for(int i=1; i<imax; i++)
				{
					for(int j=1; j<jmax; j++)
					{
						tab2[i][j]=(1/(1+((4*k*dt)/(2*p*c*dx*dx))))*(tab1[i][j]+((k*dt)/(2*p*c*dx*dx))*(tab1[i][j-1]+tab1[i][j+1]+tab1[i-1][j]+tab1[i+1][j]-4*tab1[i][j]+tab2[i][j-1]+tab2[i][j+1]+tab2[i-1][j]+tab2[i+1][j]));
						zewnatrz2+=fabs(tab2[i][j]);
					}
				}
				//z.licznik1++;
			}
			for(int i=imin; i<=imax; i++)
			{
				for(int j=jmin; j<=jmax; j++)	
				{
					tab1[i][j]=tab2[i][j];
				}
			}
			z.t+=dt;//zmiana kroku czasowego
			z.licznik2++;
		}
	}
	else
	{
		cout<<"BRAK DOSTEPU DO PLIKU!"<<endl;
	}
	for(int i=imin; i<91; i++)
	{
		for(int j=jmin; j<61; j++)
		{
			plik2<<i<<" "<<j<<" "<<tab2[i][j]<<endl;
		}
		plik2<<endl;
	}
	plik2<<endl;
	plik2.close();
	cout<<endl;
}

void dyfuzja::nieszczelnosc()
{
	int imin,i1,i2,i3,imax,jmin,j1,j2,j3,j4,jmax;
	imin=0;
	i1=20;
	i2=60;
	i3=80;
	imax=90;
	jmin=0;
	j1=10;
	j2=20;
	j3=40;
	j4=50;
	jmax=60;
	double dx,dy,dt;  dx=1;  dy=dx;  dt=2.5;
	double p,c,k;	k=1; 	p=k;	c=p;
	double T0,T1,T2;
	double zewnatrz1=0,wewnatrz1=0,zewnatrz2=0,wewnatrz2=0;
	double tolerancja=pow(10,-5);
	T0=0; T1=20; T2=25;
	zerowanie z;
	z.licznik1=0;
	z.licznik2=0;
	z.spr=0;
	z.t=0;
	z.h=0.1;
	double **tab1=new double*[91];
	double **tab2=new double*[91];
	for(int i=imin; i<101; i++)
	{
		tab1[i]=new double[61];
		tab2[i]=new double[61];
	}
	for(int i=imin; i<91; i++)
	{
		for(int j=jmin; j<61; j++)
		{
			tab1[i][j]=T2;
			tab2[i][j]=T2;
			if(i>=i1 && i<=i2 && j==jmin) {tab2[i][j]=T0; tab1[i][j]=T0;}
			if(i>=i1 && i<=i2 && j==jmax) {tab2[i][j]=T0; tab1[i][j]=T0;}
			if(i==imax && j>=j3 && j<=jmax) {tab1[i][j]=T1; tab2[i][j]=T1;}
			if(i==imin && j>=j2 && j<=j4) {tab1[i][j]=T2; tab2[i][j]=T2;}
			if(i==imax && j>=jmin && j<=j3) {tab1[i][j]=T2; tab2[i][j]=T2;}
		}
	}
	fstream plik1;
	plik1.open("zad2inne.dat",ios::out | ios::trunc);
	if(plik1.good()==true)
	{
		while(z.licznik2<=6000)//zewnetrzna petla dla badania zbieznosci
		{
			z.licznik1=0;
			zewnatrz1=zewnatrz2;
			zewnatrz2=0;
			while(z.licznik1<=50)
			{
				z.licznik1++;
				//cout<<"iteracja2 nr.: "<<z.licznik1<<endl;
				for(int i=imin; i<=imax; i++)
				{
					for(int j=jmin; j<=jmax; j++)
					{
						if(i>=i1 && i<=i2 && j==jmin) {tab2[i][j]=T0; tab1[i][j]=T0;}
						if(i>=i1 && i<=i2 && j==jmax) {tab2[i][j]=T0; tab1[i][j]=T0;}
						if(i==imax && j>=j3 && j<=jmax) {tab1[i][j]=T1; tab1[i][j]=T1;}
						if(i==imin && j>=j2 && j<=j4) {tab1[i][j]=T2; tab2[i][j]=T2;}
						if(i==imax && j>=jmin && j<=j3) {tab1[i][j]=T2; tab2[i][j]=T2;}
						//KRAWĘDZIE
						if(i==imin && j>=jmin && j<j2) tab2[i][j]=(((k/dx)*tab2[i+1][j])+z.h*T0)/(z.h+(k/dx));//PIONOWA LEWA DOLNA
						if(i==imin && j>j4 && j<jmax) tab2[i][j]=(((k/dx)*tab2[i+1][j])+z.h*T0)/(z.h+(k/dx));//PIONOWA LEWA GORNA
						if(j==jmin && i>imin && i<i1) tab2[i][j]=(((k/dx)*tab2[i][j+1])+z.h*T0)/(z.h+(k/dx));//POZIOMA LEWA DOLNA
						if(j==jmin && i>i2 && i<imax) tab2[i][j]=(((k/dx)*tab2[i][j+1])+z.h*T0)/(z.h+(k/dx));//POZIOMA PRAWA DOLNA
						if(j==jmax && i>imin && i<i1)  tab2[i][j]=(((k/dx)*tab2[i][j-1])+z.h*T0)/(z.h+(k/dx));//POZIOMA LEWA GORNA
						if(j==jmax && i>i2 && i<imax)  tab2[i][j]=(((k/dx)*tab2[i][j-1])+z.h*T0)/(z.h+(k/dx));//POZIOMA PRAWA GORNA
						//KANTY
						if(i==imin && j==jmin) tab2[i][j]=(((k/(sqrt(2)*dx))*tab2[i+1][j+1])+z.h*T0)/(z.h+(k/sqrt(2)*dx));//1 KANT
						if(i==imin && j==jmax) tab2[i][j]=(((k/(sqrt(2)*dx))*tab2[i+1][j-1])+z.h*T0)/(z.h+(k/sqrt(2)*dx));//2 KANT
						if(i==imax && j==j3) tab2[i][j]=T1;
					}
						
				}
				zewnatrz1=zewnatrz2;
				zewnatrz2=0;
				for(int i=1; i<imax; i++)
				{
					for(int j=1; j<jmax; j++)
					{
						tab2[i][j]=(1/(1+((4*k*dt)/(2*p*c*dx*dx))))*(tab1[i][j]+((k*dt)/(2*p*c*dx*dx))*(tab1[i][j-1]+tab1[i][j+1]+tab1[i-1][j]+tab1[i+1][j]-4*tab1[i][j]+tab2[i][j-1]+tab2[i][j+1]+tab2[i-1][j]+tab2[i+1][j]));
						zewnatrz2+=fabs(tab2[i][j]);
					}
				}
				//z.licznik1++;
			}
			for(int i=imin; i<=imax; i++)
			{
				for(int j=jmin; j<=jmax; j++)	
				{
					tab1[i][j]=tab2[i][j];
				}
			}
			z.t+=dt;//zmiana kroku czasowego
			z.licznik2++;
		}
	}
	else
	{
		cout<<"BRAK DOSTEPU DO PLIKU!"<<endl;
	}
	for(int i=imin; i<91; i++)
	{
		for(int j=jmin; j<61; j++)
		{
			plik1<<i<<" "<<j<<" "<<tab2[i][j]<<endl;
		}
		plik1<<endl;
	}
	plik1<<endl;
	plik1.close();
	cout<<endl;
}

int main()
{
	dyfuzja d;
	d.izolacja();
	d.nieszczelnosc();
	return 0;
}
