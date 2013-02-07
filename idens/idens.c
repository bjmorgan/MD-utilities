/* This program calcultates ionic densities as a function of a given coordinate direction (x,y,z)
from the positions' file from Paul Madden's simulation code (cartesian coordinates). */
/* One can put limit conditions for the 2 directions that are not considered for the density 
to look at the density in a particular region of the cell. */

/* To execute the program, you just have to launch it with an input file giving
the following information (program < inputfile)
first line : name of the positions's file
second line : number of configurations in the positions' file and number of sectors in time
third line : direction along which the ionic densities will be calculated (x, y or z)
following lines : length of the box in the 3 cartesian directions
following line : number of different species
following lines : number of ions for each species 
following line : number of boxes into which the simulation cell will be divided 
following line : normalization by the bulk density (T or F) 
following line : limits (T or F), if T then complete the last lines, if F, it is the last line
following line : position of the center x / space / y / space / z (xc, yc, zc) ]
last line : maximum distance between the center and each ion for the calculation of the density in bohrs (Rmax). */


/* The code is written without considering the periodic boundary conditions, so be careful or modify the code ! 
   Rmax represents a distance in the plane perpendicular to the direction specified in the third line of the input file. */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#define L 1000
#define pi 3.14159265


char nom[50],dirxyz[50],nom2[50],test[50],testlim[50];
int nions,nconfigs,nsectors,nbox,*nspecies,diffspecies;
double Lbox[4],**pos,**boxes,**boxeslim,*bulk,Rmax,center[4];


void read();
void write();
int find(int j);
void normalize(int s);


char **cmatrix(int nl, int nc);
int **imatrix(int nl, int nc);
double *dvector(int n);
int *ivector(int n);
double **dmatrix(int nl, int nc);

int main ( void )
{
    printf("Hello\n");
    read();
    write();
    printf("The output files are ionic_density.out and bloc_density.out if you put limits. \n");    
    return 0;
} // end of main()

/*Reads the input file and allocates the vectors and matrices*/
void read()
{
    int i;    
    scanf("%s",nom);                    /* name of the positions' file */
    scanf("%d %d",&nconfigs,&nsectors);            /* number of configurations in the positions's file */ /* and number of sectors in time */
    nconfigs/=nsectors;
    scanf("%s",dirxyz);                        /* direction of divisions */
    scanf("%lf",&Lbox[1]);                       /* length of the simulation box */
    scanf("%lf",&Lbox[2]);                        /* length of the simulation box */
    scanf("%lf",&Lbox[3]);                         /* length of the simulation box */
    scanf("%d",&diffspecies);                  /* number of different species */    
    nspecies=ivector(diffspecies+1);
    for(i=1;i<=diffspecies;i++) scanf("%d",&nspecies[i]);     /* number of ions for each species */
    scanf("%d",&nbox);                         /* number of boxes in the chosen direction */
    scanf("%s",test);                    /* T or F for the normalization by the bulk density */
    scanf("%s",testlim);                    /* T or F for limits */    
    if(strcmp(testlim,"T")==0)
    {
        scanf("%lf %lf %lf",&center[1],&center[2],&center[3]);    
        scanf("%lf",&Rmax);      
        boxeslim=dmatrix(diffspecies+1,nbox+1);    
    }    
    nions=0;
    for(i=1;i<=diffspecies;i++) nions+=nspecies[i];    /* Calculation of the total number of ions */
    nspecies[0]=nions;    
    boxes=dmatrix(diffspecies+1,nbox+1);            /* Boxes to fill in with the number of ions in each box */
    pos=dmatrix(4,nions+1);    
    printf("Reading : OK !\n");
} // end of read()


/* Calculates the number densities and output them */
void write()
{
    FILE *in,*out;
    char ligne[L];
    int i,j,k,s,dir,dir1,dir2;
    double R,dx,dy;    
    /* Transformation of the direction into a given number */
    if (strcmp(dirxyz,"x")==0) dir=1;
    if (strcmp(dirxyz,"y")==0) dir=2;
    if (strcmp(dirxyz,"z")==0) dir=3;    
    /* Transformation of the direction of limits into numbers if there are limits */
    if(strcmp(testlim,"T")==0)
    {
        if(dir==1) {dir1=2; dir2=3;}
        if(dir==2) {dir1=1; dir2=3;}
        if(dir==3) {dir1=1; dir2=2;}
    }    
    /* Opening of the positions'file only one time */
    in=fopen(nom,"r");
     
    /* Filling of boxes */
    for(s=1; s<=nsectors; s++)
    {
        for(i=1; i<=nconfigs; i++)
        {
            if ( i%100 == 0 ) printf("Config %d over %d\n",(s-1)*nconfigs+i,nconfigs*nsectors);
            for(j=1; j<=nions; j++) 
            {
                fgets(ligne,L,in);
                sscanf(ligne,"%lf %lf %lf",&pos[1][j],&pos[2][j],&pos[3][j]);
                k=1;
                while( pos[dir][j] > (k*Lbox[dir]/nbox) ) k++;    /* Find the box in which the ion is. */
                boxes[0][k]++; 
                boxes[find(j)][k]++;        /* One more ion in the box where the ion is. */
                if(strcmp(testlim,"T")==0)
                {    
                    dx=pos[dir1][j]-center[dir1];
                    dy=pos[dir2][j]-center[dir2];
                    R=sqrt(dx*dx+dy*dy);
                    if(R<=Rmax)          
                    {
                        while( pos[dir][j] > (k*Lbox[dir]/nbox) ) k++;
                        boxeslim[0][k]++; boxeslim[find(j)][k]++;
                    }        
                }
            }
        }
    /* Average over the number of configurations, and obtention of the ionic density dividing by the volume of each box */
        for(k=1;k<=nbox;k++) 
        {
              for(j=0;j<=diffspecies;j++)          
              {          
                  boxes[j][k]/=(Lbox[1]*Lbox[2]*Lbox[3]*nconfigs/nbox);
                  if(strcmp(testlim,"T")==0) boxeslim[j][k]/=(Lbox[dir]*Rmax*Rmax*pi*nconfigs/nbox);
              }
        }
    /* Filling of the output files */
        sprintf(nom2,"ionic_density.out-%d",s);
        out=fopen(nom2,"w");
        for(k=1;k<=nbox;k++)
        {
            for(j=0;j<=diffspecies;j++)
            {         
                if ( j==0 ) fprintf(out,"%e    %e",k*Lbox[dir]/nbox,boxes[j][k]);         
                if ( (j!=0) && (j!=diffspecies) ) fprintf(out,"    %e",boxes[j][k]);         
                if ( j==diffspecies) fprintf(out,"    %e\n",boxes[j][k]);
            }
        }  
        fclose(out);    
        if(strcmp(testlim,"T")==0)
        {
            sprintf(nom2,"bloc_density.out-%d",s);
            out=fopen(nom2,"w");
            for(k=1;k<=nbox;k++)
            {
                for(j=0;j<=diffspecies;j++)
                {
                    if (j==0) fprintf(out,"%e    %e", k*Lbox[dir] / nbox, boxeslim[j][k]);
                    if ((j!=0)&&(j!=diffspecies)) fprintf(out,"    %e",boxeslim[j][k]);
                    if (j==diffspecies) fprintf(out,"    %e\n",boxeslim[j][k]);
                }
            }
            fclose(out);
        }
    /* Normalization */
        if(strcmp(test,"T")==0) 
        {
            normalize(s); 
            printf("Normalization of the ionic densities are written in norm-idens-s.outi.\n");
        }    
    }
    fclose(in);
} // end of write()

void normalize(int s)
{
    FILE *out;
    int i, j, k, dir;  
    /* Transformation of the direction into a given number */
    if (strcmp(dirxyz,"x")==0) dir=1;
    if (strcmp(dirxyz,"y")==0) dir=2;
    if (strcmp(dirxyz,"z")==0) dir=3;   
    bulk = dvector(diffspecies+1);    
    /* Bulk density = average density */
    for( j=0; j<=diffspecies; j++) 
    {
        bulk[j]=nspecies[j]/(Lbox[1]*Lbox[2]*Lbox[3]);  
    }
    /* Normalization by the bulk density */
    for(k=1;k<=nbox;k++) 
    {
        for(j=0;j<=diffspecies;j++) { boxes[j][k]/=bulk[j] }; 
    }
    
    /* Filling of the output files */
    for(j=0;j<=diffspecies;j++)
    {
        sprintf(nom2,"norm-idens-%d.out%d",s,j);
        out=fopen(nom2,"w");
        for(k=1;k<=nbox;k++) fprintf(out,"%lf    %lf\n",k*Lbox[dir]*0.52918/nbox,boxes[j][k]);
        fclose(out);
    }
}

/* Find the species ion j belongs to */
int find(int j)
{
    int i, sum;  
    i=1; 
    sum=nspecies[i];
    while(j>sum) 
    {
        i++; 
        sum+=nspecies[i];
    }   
    return i;
} // find

/********************************************************************************************/
/*Allocation dynamique*/

char **cmatrix(int nl, int nc)
{
    int i;
    char **m;
    m=(char **) malloc(nl*sizeof(char*));
    if (m)
    { 
        m[0]=(char *) malloc(nl*nc*sizeof(char));
        if (m[0]==NULL) return NULL;
        for (i=1;i<nl;i++) 
        {
            m[i]=m[i-1]+nc;
        }
    }
    return m;
}


int **imatrix(int nl, int nc)
{
    int i;
    int **m;
    m=(int **) malloc(nl*sizeof(int*));
    if (m) 
    {
        m[0]=(int *) malloc(nl*nc*sizeof(int));
        if (m[0]==NULL) return NULL;
        for (i=1;i<nl;i++) 
        {
            m[i]=m[i-1]+nc;
        }
    }
    return m;
}

double *dvector(int n)
{
    return (double *) malloc(n*sizeof(double));
}

int *ivector(int n)
{
    return (int *) malloc(n*sizeof(int));
}

double **dmatrix(int nl, int nc)
{
    int i;
    double **m;
    m=(double **) malloc(nl*sizeof(double*));
    if (m)
    { 
        m[0]=(double *) malloc(nl*nc*sizeof(double));
        if (m[0]==NULL) return NULL;
        for (i=1;i<nl;i++)
        {
            m[i]=m[i-1]+nc;
        }
    }
    return m;
}
