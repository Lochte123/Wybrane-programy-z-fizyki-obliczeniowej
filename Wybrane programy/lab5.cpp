#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
//autor: £ukasz Kokosza

struct licznik
{
	int licznik1;
};

struct struna
{
	double x;
	void schemat_Verleta();
	void granica_osrodkow_a();
	void granica_osrodkow_b();
	void granica_osrodkow_c();
	void drgania_tlumione();
	void drgania_wymuszone();
	void energia(double t1,double t2);
	double ro(double w);
	double sila_f(double q,double t);
	double sila_f2(double q,double t,double w);
};

double struna::ro(double w)
{
	double x0=0.75;
	double ro0=3;
	if(w<x0) return 1;
	if(w>x0) return ro0;
}

double struna::sila_f(double q, double t)
{
	double x0=0.5;
	double w=M_PI/2;
	return cos(w*t)*exp(-pow(10,5)*pow((q-x0),2));
}
double struna::sila_f2(double q, double t, double w)
{
	double x0=0.5;
	return cos(w*t)*exp(-pow(10,5)*pow((q-x0),2));
}

void struna::schemat_Verleta()
{
	licznik l;
	struna x;
	l.licznik1=0;
	double t=0.0;
	double dt=0.005;
	int n=101;
	double pi=M_PI;
	double *v,*u,*a,*uanalityczne,*vpom;
	v=new double[n];
	u=new double[n];
	a=new double[n];
	uanalityczne=new double[n];
	vpom=new double[n];
	double dx=0.01;
	x.x=0.0;
	for(int i=0; i<n; i++)
	{
		u[i]=sin(x.x*pi)-((0.5)*sin(2*pi*x.x));
		v[i]=0;
		a[i]=-(pi*pi)*cos(pi*t)*sin(pi*x.x)+(2*pi*pi)*cos(2*pi*t)*sin(2*pi*x.x);
		uanalityczne[i]=0;
		vpom[i]=0;		
		x.x+=dx;
	}
	u[0]=0;
	u[n-1]=0;
	fstream plik1;
	if(plik1.good()==true)
	{
	while(t<2)
	{
		x.x=0.0;
		l.licznik1++;
		cout<<l.licznik1<<endl;
		for(int i=0; i<n; i++)
		{
			uanalityczne[i]=cos(pi*t)*sin(pi*x.x)-(0.5)*cos(2*pi*t)*sin(2*pi*x.x);
			x.x+=dx;
			//cout<<t<<" "<<uanalityczne[i]<<endl;
		}
		for(int i=0; i<n; i++)
		{
			vpom[i]=v[i]+(dt/2)*a[i];
		}
		for(int i=0; i<n; i++)
		{
			u[i]=u[i]+dt*vpom[i];
		}
		for(int i=1; i<n-1; i++)
		{
			a[i]=(u[i+1]+u[i-1]-2*u[i])/(dx*dx);
		}
		for(int i=0; i<n; i++)
		{
			v[i]=vpom[i]+(dt/2)*a[i];
		}
		if(l.licznik1==50){
		plik1.open("zad1a.dat",ios::out | ios::trunc);
		for(int i=0; i<n; i++)
		{
			plik1<<i*dx<<" "<<u[i]<<" "<<uanalityczne[i]<<endl;
		}
		plik1.close();
		}
		if(l.licznik1==100){
		plik1.open("zad1b.dat",ios::out | ios::trunc);
		for(int i=0; i<n; i++)
		{
			plik1<<i*dx<<" "<<u[i]<<" "<<uanalityczne[i]<<endl;
		}
		plik1.close();
		}
		if(l.licznik1==150){
		plik1.open("zad1c.dat",ios::out | ios::trunc);
		for(int i=0; i<n; i++)
		{
			plik1<<i*dx<<" "<<u[i]<<" "<<uanalityczne[i]<<endl;
		}		
		plik1.close();
		}
		if(l.licznik1==200){
		plik1.open("zad1d.dat",ios::out | ios::trunc);
		for(int i=0; i<n; i++)
		{
			plik1<<i*dx<<" "<<u[i]<<" "<<uanalityczne[i]<<endl;
		}
		plik1.close();
		}
		if(l.licznik1==250){
		plik1.open("zad1e.dat",ios::out | ios::trunc);
		for(int i=0; i<n; i++)
		{
			plik1<<i*dx<<" "<<u[i]<<" "<<uanalityczne[i]<<endl;
		}
		plik1.close();
		}
		if(l.licznik1==300){
		plik1.open("zad1f.dat",ios::out | ios::trunc);
		for(int i=0; i<n; i++)
		{
			plik1<<i*dx<<" "<<u[i]<<" "<<uanalityczne[i]<<endl;
		}
		plik1.close();
		}
		if(l.licznik1==350){
		plik1.open("zad1g.dat",ios::out | ios::trunc);
		for(int i=0; i<n; i++)
		{
			plik1<<i*dx<<" "<<u[i]<<" "<<uanalityczne[i]<<endl;
		}
		plik1.close();
		}
		if(l.licznik1==400){
		plik1.open("zad1h.dat",ios::out | ios::trunc);
		for(int i=0; i<n; i++)
		{
			plik1<<i*dx<<" "<<u[i]<<" "<<uanalityczne[i]<<endl;
		}
		plik1.close();
		}
		t+=dt;
	}
	}
}

void struna::granica_osrodkow_a()
{
//a
	licznik l;
	struna x;
	l.licznik1=0;
	double t=0.0;
	double dt=0.005;
	int n=101;
	double pi=M_PI;
	double *v,*u,*a,*uanalityczne,*vpom;
	v=new double[n];
	u=new double[n];
	a=new double[n];
	uanalityczne=new double[n];
	vpom=new double[n];
	double dx=0.01;
	x.x=0.0;
	for(int i=0; i<n; i++)
	{
		u[i]=exp(-100*pow(x.x-0.5,2));
		v[i]=0;
		a[i]=-(pi*pi)*cos(pi*t)*sin(pi*x.x)+(2*pi*pi)*cos(2*pi*t)*sin(2*pi*x.x);
		uanalityczne[i]=0;
		vpom[i]=0;		
		x.x+=dx;
	}
	u[0]=0;
	u[n-1]=0;
	fstream plik2;
	plik2.open("zad2a.dat",ios::out | ios::trunc);
	if(plik2.good()==true)
	{
	while(t<4)
	{
		x.x=0.0;
		l.licznik1++;
		cout<<l.licznik1<<endl;
		for(int i=0; i<n; i++)
		{
			uanalityczne[i]=cos(pi*t)*sin(pi*x.x)-(0.5)*cos(2*pi*t)*sin(2*pi*x.x);
			x.x+=dx;
			//cout<<t<<" "<<uanalityczne[i]<<endl;
		}
		for(int i=0; i<n; i++)
		{
			vpom[i]=v[i]+(dt/2)*a[i];
		}
		for(int i=0; i<n; i++)
		{
			u[i]=u[i]+dt*vpom[i];
		}
		for(int i=1; i<n-1; i++)
		{
			a[i]=(u[i+1]+u[i-1]-2*u[i])/(dx*dx);
		}
		for(int i=0; i<n; i++)
		{
			v[i]=vpom[i]+(dt/2)*a[i];
		}
		for(int i=0; i<n; i++)
		{
			plik2<<t<<" "<<i*dx<<" "<<u[i]<<endl;
		}
		plik2<<endl;		

		t+=dt;
	}
	}
	plik2.close();

}

void struna::granica_osrodkow_b()
{
	licznik l;
	struna x;
	l.licznik1=0;
	double t=0.0;
	double dt=0.005;
	int n=101;
	double pi=M_PI;
	double *v,*u,*a,*uanalityczne,*vpom;
	v=new double[n];
	u=new double[n];
	a=new double[n];
	uanalityczne=new double[n];
	vpom=new double[n];
	double dx=0.01;
	x.x=0.0;
	for(int i=0; i<n; i++)
	{
		u[i]=exp(-100*pow(x.x-0.5,2));
		v[i]=0;
		a[i]=-(pi*pi)*cos(pi*t)*sin(pi*x.x)+(2*pi*pi)*cos(2*pi*t)*sin(2*pi*x.x);
		uanalityczne[i]=0;
		vpom[i]=0;		
		x.x+=dx;
	}
	u[0]=0;
	u[n-1]=0;
	fstream plik3;
	plik3.open("zad2b.dat",ios::out | ios::trunc);
	if(plik3.good()==true)
	{
	while(t<4)
	{
		x.x=0.0;
		l.licznik1++;
		cout<<l.licznik1<<endl;
		for(int i=0; i<n; i++)
		{
			uanalityczne[i]=cos(pi*t)*sin(pi*x.x)-(0.5)*cos(2*pi*t)*sin(2*pi*x.x);
			x.x+=dx;
			//cout<<t<<" "<<uanalityczne[i]<<endl;
		}
		for(int i=0; i<n; i++)
		{
			vpom[i]=v[i]+(dt/2)*a[i];
		}
		for(int i=0; i<n; i++)
		{
			u[i]=u[i]+dt*vpom[i];
		}
		u[0]=u[1];
		u[n-1]=u[n-2];	
		for(int i=1; i<n-1; i++)
		{
			a[i]=(u[i+1]+u[i-1]-2*u[i])/(dx*dx);
		}
		for(int i=0; i<n; i++)
		{
			v[i]=vpom[i]+(dt/2)*a[i];
		}
		for(int i=0; i<n; i++)
		{
			plik3<<t<<" "<<u[0]<<" "<<u[n-1]<<endl;
		}	

		t+=dt;
	}
	}
	plik3.close();

}

void struna::granica_osrodkow_c()
{
	licznik l;
	struna x;
	l.licznik1=0;
	double t=0.0;
	double dt=0.005;
	int n=101;
	double pi=M_PI;
	double *v,*u,*a,*uanalityczne,*vpom;
	v=new double[n];
	u=new double[n];
	a=new double[n];
	uanalityczne=new double[n];
	vpom=new double[n];
	double dx=0.01;
	x.x=0.0;
	for(int i=0; i<n; i++)
	{
		u[i]=exp(-100*pow(x.x-0.5,2));
		v[i]=0;
		a[i]=-(pi*pi)*cos(pi*t)*sin(pi*x.x)+(2*pi*pi)*cos(2*pi*t)*sin(2*pi*x.x);
		uanalityczne[i]=0;
		vpom[i]=0;		
		x.x+=dx;
	}
	u[0]=0;
	u[n-1]=0;
	fstream plik4;
	plik4.open("zad2c.dat",ios::out | ios::trunc);
	if(plik4.good()==true)
	{
	while(t<4)
	{
		x.x=0.0;
		l.licznik1++;
		cout<<l.licznik1<<endl;
		for(int i=0; i<n; i++)
		{
			uanalityczne[i]=cos(pi*t)*sin(pi*x.x)-(0.5)*cos(2*pi*t)*sin(2*pi*x.x);
			x.x+=dx;
			//cout<<t<<" "<<uanalityczne[i]<<endl;
		}
		for(int i=0; i<n; i++)
		{
			vpom[i]=v[i]+(dt/2)*a[i];
		}
		for(int i=0; i<n; i++)
		{
			u[i]=u[i]+dt*vpom[i];
		}
		for(int i=1; i<n-1; i++)
		{
			a[i]=(u[i+1]+u[i-1]-2*u[i])/(dx*dx)*(1/x.ro(i*dx));
		}
		for(int i=0; i<n; i++)
		{
			v[i]=vpom[i]+(dt/2)*a[i];
		}
		for(int i=0; i<n; i++)
		{
			plik4<<t<<" "<<i*dx<<" "<<u[i]<<endl;
		}
		plik4<<endl;		

		t+=dt;
	}
	}
	plik4.close();

}

void struna::drgania_tlumione()
{/*
	double B=0.2;
	licznik l;
	struna x;
	l.licznik1=0;
	double t=0.0;
	double dt=0.005;
	int n=101;
	double pi=M_PI;
	double *v,*u,*a,*uanalityczne,*vpom,*upom;
	v=new double[n];
	u=new double[n];
	a=new double[n];
	uanalityczne=new double[n];
	vpom=new double[n];
	upom=new double[n];
	double dx=0.01;
	x.x=0.0;
	for(int i=0; i<n; i++)
	{
		u[i]=exp(-100*pow(x.x-0.5,2));
		v[i]=0;
		a[i]=-(pi*pi)*cos(pi*t)*sin(pi*x.x)+(2*pi*pi)*cos(2*pi*t)*sin(2*pi*x.x);
		uanalityczne[i]=0;
		vpom[i]=0;		
		x.x+=dx;
	}
	u[0]=0;
	u[n-1]=0;
	fstream plik5;
	plik5.open("zad3a.dat",ios::out | ios::trunc);
	if(plik5.good()==true)
	{
	while(t<4)
	{
		x.x=0.0;
		l.licznik1++;
		cout<<l.licznik1<<endl;
		for(int i=0; i<n; i++)
		{
			uanalityczne[i]=cos(pi*t)*sin(pi*x.x)-(0.5)*cos(2*pi*t)*sin(2*pi*x.x);
			x.x+=dx;
			//cout<<t<<" "<<uanalityczne[i]<<endl;
		}
		for(int i=0; i<n; i++)
		{
			vpom[i]=v[i]+(dt/2)*a[i];
		}
		for(int i=0; i<n; i++)
		{
			upom[i]=u[i];
			u[i]=u[i]+dt*vpom[i];
		}
		for(int i=1; i<n-1; i++)
		{
			a[i]=a[i]=(u[i+1]+u[i-1]-2*u[i])/(dx*dx)-((2*B)*(u[i]-upom[i])/dt);
		}
		for(int i=0; i<n; i++)
		{
			v[i]=vpom[i]+(dt/2)*a[i];
		}
		for(int i=0; i<n; i++)
		{
			plik5<<t<<" "<<i*dx<<" "<<u[i]<<endl;
		}
		plik5<<endl;		

		t+=dt;
	}
	}
	plik5.close();*/

/*
	double B=1;
	licznik l;
	struna x;
	l.licznik1=0;
	double t=0.0;
	double dt=0.005;
	int n=101;
	double pi=M_PI;
	double *v,*u,*a,*uanalityczne,*vpom,*upom;
	v=new double[n];
	u=new double[n];
	a=new double[n];
	uanalityczne=new double[n];
	vpom=new double[n];
	upom=new double[n];
	double dx=0.01;
	x.x=0.0;
	for(int i=0; i<n; i++)
	{
		u[i]=exp(-100*pow(x.x-0.5,2));
		v[i]=0;
		a[i]=-(pi*pi)*cos(pi*t)*sin(pi*x.x)+(2*pi*pi)*cos(2*pi*t)*sin(2*pi*x.x);
		uanalityczne[i]=0;
		vpom[i]=0;		
		x.x+=dx;
	}
	u[0]=0;
	u[n-1]=0;
	fstream plik5;
	plik5.open("zad3b.dat",ios::out | ios::trunc);
	if(plik5.good()==true)
	{
	while(t<4)
	{
		x.x=0.0;
		l.licznik1++;
		cout<<l.licznik1<<endl;
		for(int i=0; i<n; i++)
		{
			uanalityczne[i]=cos(pi*t)*sin(pi*x.x)-(0.5)*cos(2*pi*t)*sin(2*pi*x.x);
			x.x+=dx;
			//cout<<t<<" "<<uanalityczne[i]<<endl;
		}
		for(int i=0; i<n; i++)
		{
			vpom[i]=v[i]+(dt/2)*a[i];
		}
		for(int i=0; i<n; i++)
		{
			upom[i]=u[i];
			u[i]=u[i]+dt*vpom[i];
		}
		for(int i=1; i<n-1; i++)
		{
			a[i]=(u[i+1]+u[i-1]-2*u[i])/(dx*dx)-((2*B)*(u[i]-upom[i])/dt);
		}
		for(int i=0; i<n; i++)
		{
			v[i]=vpom[i]+(dt/2)*a[i];
		}
		for(int i=0; i<n; i++)
		{
			plik5<<t<<" "<<i*dx<<" "<<u[i]<<endl;
		}
		plik5<<endl;		

		t+=dt;
	}
	}
	plik5.close();*/




	double B=3;
	licznik l;
	struna x;
	l.licznik1=0;
	double t=0.0;
	double dt=0.005;
	int n=101;
	double pi=M_PI;
	double *v,*u,*a,*uanalityczne,*vpom,*upom;
	v=new double[n];
	u=new double[n];
	a=new double[n];
	uanalityczne=new double[n];
	vpom=new double[n];
	upom=new double[n];
	double dx=0.01;
	x.x=0.0;
	for(int i=0; i<n; i++)
	{
		u[i]=exp(-100*pow(x.x-0.5,2));
		v[i]=0;
		a[i]=-(pi*pi)*cos(pi*t)*sin(pi*x.x)+(2*pi*pi)*cos(2*pi*t)*sin(2*pi*x.x);
		uanalityczne[i]=0;
		vpom[i]=0;		
		x.x+=dx;
	}
	u[0]=0;
	u[n-1]=0;
	fstream plik5;
	plik5.open("zad3c.dat",ios::out | ios::trunc);
	if(plik5.good()==true)
	{
	while(t<4)
	{
		x.x=0.0;
		l.licznik1++;
		cout<<l.licznik1<<endl;
		for(int i=0; i<n; i++)
		{
			uanalityczne[i]=cos(pi*t)*sin(pi*x.x)-(0.5)*cos(2*pi*t)*sin(2*pi*x.x);
			x.x+=dx;
			//cout<<t<<" "<<uanalityczne[i]<<endl;
		}
		for(int i=0; i<n; i++)
		{
			vpom[i]=v[i]+(dt/2)*a[i];
		}
		for(int i=0; i<n; i++)
		{
			upom[i]=u[i];
			u[i]=u[i]+dt*vpom[i];
		}
		for(int i=1; i<n-1; i++)
		{
			a[i]=a[i]=(u[i+1]+u[i-1]-2*u[i])/(dx*dx)-((2*B)*(u[i]-upom[i])/dt);
		}
		for(int i=0; i<n; i++)
		{
			v[i]=vpom[i]+(dt/2)*a[i];
		}
		for(int i=0; i<n; i++)
		{
			plik5<<t<<" "<<i*dx<<" "<<u[i]<<endl;
		}
		plik5<<endl;		

		t+=dt;
	}
	}
	plik5.close();
}

void struna::drgania_wymuszone()
{
		double B=1;
	licznik l;
	struna x;
	l.licznik1=0;
	double t=0.0;
	double dt=0.005;
	int n=101;
	double pi=M_PI;
	double *v,*u,*a,*uanalityczne,*vpom,*upom,*af;
	v=new double[n];
	af=new double[n];
	u=new double[n];
	a=new double[n];
	uanalityczne=new double[n];
	vpom=new double[n];
	upom=new double[n];
	double dx=0.01;
	x.x=0.0;
	for(int i=0; i<n; i++)
	{
		u[i]=0;
		v[i]=0;
		a[i]=0;	
		x.x+=dx;
	}
	u[0]=0;
	u[n-1]=0;
	fstream plik6;
	plik6.open("zad4.dat",ios::out | ios::trunc);
	if(plik6.good()==true)
	{
	while(t<=10)
	{
		x.x=0.0;
		l.licznik1++;
		cout<<l.licznik1<<endl;
		for(int i=0; i<n; i++)
		{
			uanalityczne[i]=cos(pi*t)*sin(pi*x.x)-(0.5)*cos(2*pi*t)*sin(2*pi*x.x);
			x.x+=dx;
			//cout<<t<<" "<<uanalityczne[i]<<endl;
		}
		for(int i=0; i<n; i++)
		{
			vpom[i]=v[i]+(dt/2)*a[i];
		}
		for(int i=0; i<n; i++)
		{
			upom[i]=u[i];
			u[i]=u[i]+dt*vpom[i];
		}
		for(int i=1; i<n-1; i++)
		{
			a[i]=a[i]=(u[i+1]+u[i-1]-2*u[i])/(dx*dx)-((2*B)*(u[i]-upom[i])/dt)+sila_f(i*dx,t);
		}
		for(int i=0; i<n; i++)
		{
			v[i]=vpom[i]+(dt/2)*a[i];
		}
		for(int i=0; i<n; i++)
		{
			plik6<<t<<" "<<i*dx<<" "<<u[i]<<endl;
		}
		plik6<<endl;		

		t+=dt;
	}
	}
	plik6.close();
}

void struna::energia(double t1,double t2)
{
	double dw=0.5;
	double du=0;
	double w=0;
	double B=1;
	licznik l;
	struna x;
	l.licznik1=0;
	double t=0.0;
	double v2=0;
	double u2=0;
	double E;
	double Es=0;
	double dt=0.005;
	int n=101;
	double pi=M_PI;
	double *v,*u,*a,*uanalityczne,*vpom,*upom,*af;
	v=new double[n];
	af=new double[n];
	u=new double[n];
	a=new double[n];
	uanalityczne=new double[n];
	vpom=new double[n];
	upom=new double[n];
	double dx=0.01;
	x.x=0.0;
	for(int i=0; i<n; i++)
	{
		u[i]=0;
		v[i]=0;
		a[i]=0;	
		x.x+=dx;
	}
	u[0]=0;
	u[n-1]=0;
	fstream plik7;
	plik7.open("zad5.dat",ios::out | ios::trunc);
	if(plik7.good()==true)
	{
	while(w<=10*M_PI){
	Es=0;
	t=0;
	while(t<=20)
	{
		x.x=0.0;
		//l.licznik1++;
		//cout<<l.licznik1<<endl;
		for(int i=0; i<n; i++)
		{
			uanalityczne[i]=cos(pi*t)*sin(pi*x.x)-(0.5)*cos(2*pi*t)*sin(2*pi*x.x);
			x.x+=dx;
			//cout<<t<<" "<<uanalityczne[i]<<endl;
		}
		for(int i=0; i<n; i++)
		{
			vpom[i]=v[i]+(dt/2)*a[i];
		}
		for(int i=0; i<n; i++)
		{
			upom[i]=u[i];
			u[i]=u[i]+dt*vpom[i];
		}
		for(int i=1; i<n-1; i++)
		{
			a[i]=a[i]=(u[i+1]+u[i-1]-2*u[i])/(dx*dx)-((2*B)*(u[i]-upom[i])/dt)+sila_f2(i*dx,t,w);
		}
		v2=0;
		du=0;
		for(int i=0; i<n; i++)
		{
			v[i]=vpom[i]+(dt/2)*a[i];
			v2+=pow(v[i],2);
			if(i<n-1) du+=pow((u[i+1]-u[i])/dx,2);
		}

		E=(0.5*v2*dx)+(0.5*du*dx);	
		if(t>t1) {Es+=E;}
		
		t+=dt;
	}
	Es=(1/(t2-t1))*Es*dt;
	plik7<<w<<" "<<Es<<endl;
	w+=dw;
}
	}
	plik7.close();
}

int main()
{
	struna s;
	s.schemat_Verleta();
	s.granica_osrodkow_a();
	s.granica_osrodkow_b();
	s.granica_osrodkow_c();
	s.drgania_tlumione();
	s.drgania_wymuszone();
	s.energia(16,20);
	return 0;
}
