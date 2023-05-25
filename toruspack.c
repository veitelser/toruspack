#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define DMAX 30


double (**x)[DMAX],(**xA)[DMAX],(**xR)[DMAX],(*xB)[DMAX];
double **weight,**conerr;

double **x1,*w1;

double *val;
int *rank;

int dim,num;

int monostep;
double *errseq;

double size,hsize,beta,gam;

char cmdfile[20],errfile[20],runfile[20],solfile[20],weightfile[20];

	
void allocvars()
	{
	int i,m;
	
	x=malloc(num*sizeof(double*[DMAX]));
	xA=malloc(num*sizeof(double*[DMAX]));
	xR=malloc(num*sizeof(double*[DMAX]));
	xB=malloc(num*sizeof(double[DMAX]));
	
	weight=malloc(num*sizeof(double*));
	conerr=malloc(num*sizeof(double*));
	
	x1=malloc(dim*sizeof(double*));
	w1=malloc((num-1)*sizeof(double));
	
	val=malloc((num-1)*sizeof(double));
	rank=malloc((num-1)*sizeof(int));
	
	if(monostep)
		errseq=malloc((monostep+1)*sizeof(double));
	
	for(i=0;i<num;++i)
		{
		x[i]=malloc(num*sizeof(double[DMAX]));
		xA[i]=malloc(num*sizeof(double[DMAX]));
		xR[i]=malloc(num*sizeof(double[DMAX]));
		
		weight[i]=malloc(num*sizeof(double));
		conerr[i]=malloc(num*sizeof(double));
		}
		
	for(m=0;m<dim;++m)
		x1[m]=malloc((num-1)*sizeof(double));
	}
		
	
void center(double u[DMAX])
	{
	int m;
	
	for(m=0;m<dim;++m)
		u[m]-=size*rint(u[m]/size);
	}
	
	
double norm2(int len,double u[DMAX])
	{
	double nu;
	int l;
	
	nu=0.;
	for(l=0;l<len;++l)
		nu+=u[l]*u[l];
		
	return nu;
	}
	
	
void normalize(int len,double nv,double u[DMAX],double v[DMAX])
	{
	int l;
	double rescale;
	
	rescale=sqrt(nv/norm2(len,u));
	
	for(l=0;l<len;++l)
		v[l]=rescale*u[l];
	}
	

void satnormalize(int len,double nv0,double u0[DMAX],double v0[DMAX])
	{
	double u1[DMAX],v1[DMAX],nu1;
	int l,k,sat[DMAX];
	
	normalize(len,nv0,u0,v0);
	
	k=0;
	nu1=nv0;
	
	for(l=0;l<len;++l)
		{
		if(fabs(v0[l])>hsize)
			{
			sat[l]=1;
			
			v0[l]=v0[l]<0. ? -hsize : hsize;
			
			nu1-=hsize*hsize;
			}
		else
			{
			sat[l]=0;
			
			u1[k++]=v0[l];
			}
		}
		
	if(k==len)
		return;
		
	satnormalize(k,nu1,u1,v1);
	
	k=0;
	for(l=0;l<len;++l)
		if(!sat[l])
			v0[l]=v1[k++];
	}
	
	
void projA()
	{
	double u[DMAX],v[DMAX],nu,dx;
	int i,j,m;
	
	for(i=0;i<num-1;++i)
	for(j=i+1;j<num;++j)
		{
		for(m=0;m<dim;++m)
			u[m]=x[i][j][m]-x[j][i][m];
			
		center(u);
		
		nu=norm2(dim,u);
		
		if(nu>1.)
			{
			for(m=0;m<dim;++m)
				{
				xA[i][j][m]=x[i][j][m];
				xA[j][i][m]=x[j][i][m];
				}
			}
		else
			{
			satnormalize(dim,1.,u,v);
					
			for(m=0;m<dim;++m)
				{
				dx=.5*(v[m]-u[m]);
				
				xA[i][j][m]=x[i][j][m]+dx;
				xA[j][i][m]=x[j][i][m]-dx;
				}
			}
		}
	}
	
	
void reflA()
	{
	int i,j,m;
	
	for(i=0;i<num;++i)
	for(j=0;j<num;++j)
		{
		if(i==j)
			continue;
		
		for(m=0;m<dim;++m)
			xR[i][j][m]=2.*xA[i][j][m]-x[i][j][m];
			
		center(xR[i][j]);
		}
	}
	
	
int valinc(const void *a,const void *b)
	{
    double arg1 = val[*(const int*)a];
    double arg2 = val[*(const int*)b];
 
    return (arg1 > arg2) - (arg1 < arg2);
	}
	
	
double centroid1(double *u,double *w)
	{
	int k,r;
	double c,cmin,d,dmin;
	
	c=0.;
	for(k=0;k<num-1;++k)
		{
		if(u[k] > 0.)
			u[k]-=size;
			
		c+=w[k]*u[k];
		
		rank[k]=k;
		val[k]=u[k];
		}
		
	d=0.;
	for(k=0;k<num-1;++k)
		d+=w[k]*(u[k]-c)*(u[k]-c);
		
	qsort(rank,num-1,sizeof(int),valinc);	
	
	dmin=d;
	cmin=c;
	for(r=0;r<num-2;++r)
		{
		k=rank[r];
		
		d+=(w[k]-w[k]*w[k])*size*size+2.*w[k]*size*(u[k]-c);
		c+=w[k]*size;
		
		if(d<dmin)
			{
			dmin=d;
			cmin=c;
			}
		}
		
	if(cmin>hsize)
		cmin-=size;
		
	return cmin;
	}
	
	
void projB()
	{
	int i,j,k,m;
	double wsum;

	for(i=0;i<num;++i)
		{
		k=0;
		wsum=0.;
		for(j=0;j<num;++j)
			{
			if(i==j)
				continue;
				
			for(m=0;m<dim;++m)
				x1[m][k]=xR[i][j][m];
				
			w1[k]=(i<j) ? weight[i][j] : weight[j][i];
			
			wsum+=w1[k];
			
			++k;
			}
			
		for(k=0;k<num-1;++k)
			w1[k]/=wsum;
			
		for(m=0;m<dim;++m)
			xB[i][m]=centroid1(x1[m],w1);
		}
	}
	
	
double RRR()
	{
	double u[DMAX],u2,err,averr;
	int i,j,m;
	
	for(i=0;i<num-1;++i)
	for(j=i+1;j<num;++j)
		conerr[i][j]=0.;
		
	projA();
	reflA();
	projB();
	
	err=0.;
	for(i=0;i<num;++i)
	for(j=0;j<num;++j)
		{
		if(i==j)
			continue;
			
		for(m=0;m<dim;++m)
			u[m]=xB[i][m]-xA[i][j][m];
			
		center(u);
		
		for(m=0;m<dim;++m)
			{
			x[i][j][m]+=beta*u[m];
			
			u2=u[m]*u[m];
			
			err+=u2;
			
			if(i<j)
				conerr[i][j]+=u2;
			else
				conerr[j][i]+=u2;
			}
			
		center(x[i][j]);
		}
		
	averr=0.;
	for(i=0;i<num-1;++i)
	for(j=i+1;j<num;++j)
		averr+=conerr[i][j];
		
	averr/=.5*num*(num-1.);
	
	for(i=0;i<num-1;++i)
	for(j=i+1;j<num;++j)
		weight[i][j]+=gam*(conerr[i][j]/averr-weight[i][j]);
				
	return sqrt(err/num);
	}
	
	
double urand()
	{
	return ((double)rand())/RAND_MAX;
	}
	
	
void randstart()
	{
	int i,j,m;
	
	for(i=0;i<num;++i)
		{
		for(m=0;m<dim;++m)
			xB[i][m]=size*urand()-hsize;
		
		for(j=0;j<num;++j)
			{
			if(i==j)
				continue;
			
			for(m=0;m<dim;++m)
				x[i][j][m]=xB[i][m];
			}
		}		
	}
	
	
void printtopweight()
	{
	int i,j,itop,jtop;
	double top;
	FILE *fp;
	
	top=0.;
	
	for(i=0;i<num-1;++i)
	for(j=i+1;j<num;++j)
		if(weight[i][j]>top)
			{
			top=weight[i][j];
			
			itop=i;
			jtop=j;
			}
	
	fp=fopen(weightfile,"a");
	fprintf(fp,"%10d%10d%15.1lf\n",itop,jtop,top);
	fclose(fp);
	}
	
	
int solve(int itermax,int errstride,double errstop)
	{
	int iter,i,j,m,monostart;
	double err,errmax;
	FILE *fp;
	
	for(i=0;i<num-1;++i)
	for(j=i+1;j<num;++j)
		weight[i][j]=1.;
		
	randstart();
	
	if(monostep)
		{
		for(m=0;m<=monostep;++m)
			errseq[m]=0.;
	
		errmax=0.;
		}
	
	for(iter=1;iter<=itermax;++iter)
		{
		err=RRR();
		
		if(monostep)
			{
			if(err>errmax)
				{
				errmax=err;
				monostart=iter;
				}
		
			for(m=0;m<monostep;++m)
				errseq[m]=errseq[m+1];
			
			errseq[monostep]=err;
		
			if(iter >= monostart+monostep)
				if(err>errseq[0])
					return 0;
			}

		if(err<errstop || iter%errstride==0)
			{	
			fp=fopen(errfile,"a");
			fprintf(fp,"%10d%20.10lf\n",iter,err);
			fclose(fp);
			
			printtopweight();
			}
		
		if(err<errstop)
			return iter;
		}
		
	return 0;
	}
	
	
void printsol()
	{
	FILE *fp;
	int i,m;
	double u[DMAX];
	
	fp=fopen(solfile,"a");
	
	for(i=0;i<num;++i)
		{
		for(m=0;m<dim;++m)
			u[m]=xB[i][m]-xB[0][m]-hsize;
			
		center(u);
		
		for(m=0;m<dim;++m)
			fprintf(fp,"%17.14lf",u[m]+hsize);
		
		fprintf(fp,"\n");
		}
	
	fprintf(fp,"\n");
	fclose(fp);
	}
	
	
int main(int argc,char* argv[])
	{
	char *id;
	int a,iter,itermax,errstride,trials,t,success;
	double diter,errstop,iterave,iter2ave,itererr;
	FILE *fp;
	
	if(argc==12)
		{
		dim=atoi(argv[1]);
		num=atoi(argv[2]);
		size=atof(argv[3]);
		beta=atof(argv[4]);
		gam=atof(argv[5]);
		monostep=atoi(argv[6]);
		itermax=atoi(argv[7]);
		errstride=atoi(argv[8]);
		errstop=atof(argv[9]);
		trials=atoi(argv[10]);
		id=argv[11];
		}
	else
		{
		printf("expected 11 arguments: dim num size beta gam monostep itermax errstride errstop trials id\n");
		return 1;
		}
	
	hsize=size/2.;
	
	allocvars();
	
	sprintf(cmdfile,"%s.cmd",id);
	fp=fopen(cmdfile,"w");
	for(a=0;a<argc;++a)
    	fprintf(fp,"%s ",argv[a]);
    fprintf(fp,"\n");
	fclose(fp);
	
	sprintf(runfile,"%s.run",id);
	fp=fopen(runfile,"w");
	fclose(fp);
	
	sprintf(errfile,"%s.err",id);
	fp=fopen(errfile,"w");
	fclose(fp);
	
	sprintf(solfile,"%s.sol",id);
	fp=fopen(solfile,"w");
	fclose(fp);
	
	sprintf(weightfile,"%s.wgt",id);
	fp=fopen(weightfile,"w");
	fclose(fp);
	
	srand(time(NULL));
	
	success=0;
	iterave=0.;
	iter2ave=0.;
	
	for(t=1;t<=trials;++t)
		{
		iter=solve(itermax,errstride,errstop);
		
		fp=fopen(runfile,"a");
		fprintf(fp,"%4d%12d\n",t,iter);
		fclose(fp);
		
		if(iter)
			{
			++success;
			printsol();
			
			diter=(double)iter;
			iterave+=diter;
			iter2ave+=diter*diter;
			}
		}
	
	iterave/=(double)success;
	iter2ave/=(double)success;
	itererr=sqrt((iter2ave-iterave*iterave)/(double)success);
	
	fp=fopen(runfile,"a");
	fprintf(fp,"\n%d / %d solutions%20.2lf +/-%13.2lf iterations/solution\n",success,trials,iterave,itererr);
	fclose(fp);
		
	return 0;
	}
	