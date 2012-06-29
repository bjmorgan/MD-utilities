/*  This file is part of xbs.
   
 *  xbs is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.

 *  xbs is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.

 *  You should have received a copy of the GNU General Public License
 *  along with xbs.  If not, see <http://www.gnu.org/licenses/>. */

/* ----- various useful routines ------- */
void rx(char *msg)
{
  printf("Error: %s\n", msg);
  exit(1);
}

void cross(float a[3], float b[3], float c[3])
{
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];
  return;
}

float sp(float a[3], float b[3])
{
  float sp1;
  sp1=a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
  return sp1;
}

void vscal(float a[3], float ca, float v[3])
{
  v[0]=a[0]*ca;
  v[1]=a[1]*ca;
  v[2]=a[2]*ca;
  return;
}

void vsum(float a[3], float b[3], float ca, float cb, float v[3])
{
  v[0]=ca*a[0]+cb*b[0];
  v[1]=ca*a[1]+cb*b[1];
  v[2]=ca*a[2]+cb*b[2];
  return;
}


int parse_args (str,w)
char str[],w[8][41];
{
  int i, num, reading, quoted;
  char *p;

  p=str;
  reading=quoted=0;
  num=-1;
  while (*p != 0) {
    if (reading) {
      if ((quoted&&(*p=='\'')) ||
          ((!quoted) && (*p==' ')) ) { i++; w[num][i]=0; reading=0; }
      else { i++; w[num][i]=*p; }
    }
    else {
      if (*p != ' ') {
        num++;
        if (*p == '\'') { quoted=1; i=-1; }
        else { i=0; w[num][i]=*p; }
        reading=1;
      }
    }
    p++;
  }
  if (reading) { i++; w[num][i]=0; }
  num++;
  return num;
}

void strip (str1,str)
char str[],str1[];
{
  int l,i,i1,i2;
  l=strlen(str);

  i1=0;
  for (i=0; i<l; i++) 
    if ((str[i]!=' ') && (str[i]!='\n')) { i1=i; break; }
  i2=0;
  for (i=l-1; i>=0; i--) 
    if ((str[i]!=' ') && (str[i]!='\n')) { i2=i+1; break; }
  for (i=i1;i<i2;i++) str1[i-i1]=str[i];
  str1[i2-i1]=0;
/*  printf (" l=%d i1=%d i2=%d <%s> <%s>\n", l, i1, i2, str, str1);*/
}

int abbrev (str,ab,nchar)
char str[],ab[];
int nchar;
{
  int i,nc;
  if (strlen(str) > strlen(ab)) return 0;
  nc=strlen(str);
  if (nc<nchar) nc=nchar;
  for (i=0;i<nc;i++) if (str[i] != ab[i]) return 0;
  return 1;
}

/*  set extension on a file identifier, output=fid1 */
/*  force=1 forces change even if fid already has an extension */
/*  force=0 does not change the extension if there already is one */
void strext (fid1, fid, ext, force)
char fid[],ext[],fid1[];
int force;
{
  int i,l;
  char *p,*q;

  strcpy (fid1, fid);
  l=strlen(fid1);
  p=fid1;
  for (i=0;i<l;i++) 
    if (fid1[i]=='/') p=fid1+i;

  if (!force) {
    q=strchr(p,'.');
    if (q && (q!=fid1+strlen(fid1)-1)) return;
  }
  if (!strchr(p,'.')) strcat (fid1,".");
  q=strchr(p,'.'); 
  if (strlen(ext)>0) q++;
  *q = 0;
  strcat(fid1,ext);

}


/* ----- match ------- */
int match (str, pat)
char str[], pat[];
{
  char *p,*s;
  p=pat;
  s=str; 

  while (*p != 0) {

    if (*p == '*') {           /* found wildcard '*' in pattern */
      p++;
      while (*p == '*') p++;
      if (*p == 0) return 1;   /* trailing '*' matches all */
      for (;;) {               /* find match to char after '*' */
        if (*s == 0) return 0;
        if ((*s == *p) || (*p == '+')) 
          if (match(s+1,p+1)) return 1;   /* ok if rest matches */
        s++;
      }
      return 0;                /* tried all cases but none worked */
    }
    
    else {                     /* no wildcard -- char must match */   
      if (*s == 0) return 0;
      if ((*p != *s) && (*p != '+')) return 0;
      s++;
    }
    p++;
  }

  if (*s != 0) return 0;       /* pattern but not string exhausted */
  return 1;
}

/* ----- get_extent ----- */
void get_extent(xx1,xx2,yy1,yy2,zz1,zz2)
float *xx1,*xx2,*yy1,*yy2,*zz1,*zz2;
{
  float big,x1,x2,y1,y2,z1,z2;
  int i;

  big=1000000; if(nbas==0) big=0;
  x1=y1=z1=big;
  x2=y2=z2=-big;
  for (i=0;i<nbas;i++) {
    if (p[i][0]<x1) x1=p[i][0];
    if (p[i][0]>x2) x2=p[i][0];
    if (p[i][1]<y1) y1=p[i][1];
    if (p[i][1]>y2) y2=p[i][1];
    if (p[i][2]<z1) z1=p[i][2];
    if (p[i][2]>z2) z2=p[i][2];
  }
  *xx1=x1+center[0];  *xx2=x2+center[0];
  *yy1=y1+center[1];  *yy2=y2+center[1];
  *zz1=z1+center[2];  *zz2=z2+center[2];
}

/* ----- atompos: position and radius on paper for an atom --- */
void atompos( float fac, float p[3], float rad,
             float zp[2], float *zr)
{
  float y[3],q[3],v1[3],v2[3];
  float xxx,za1,za2,zb1,zb2,a,b;

  if (pmode==1) {
    zp[0]=fac*p[0];
    zp[1]=fac*p[1];
    *zr=fac*rad;
    *zr=MAXRAD;
    if (dist0-p[2]>0) *zr=fac*rad*dist0/(dist0-p[2]);
    if (*zr > MAXRAD) *zr=MAXRAD;
    return;
  }

  vscal(p, 1.0, q);
  q[2] = q[2]-dist;
  vscal(p, 1.0, y);
  xxx = -sp(y,q)/sp(q,q);
  vsum(y, q, 1.0, xxx, y);
  if(sp(y,y)<=1e-3) { y[0]=1.0; y[1]=0.0; y[2]=0.0; }

  a = -rad*rad / sp(q,q);
  b = rad*sqrt((1.0+a) / sp(y,y));
  vsum(q, y, a, b, v1);
  vsum(q, y, a, -b, v2);
  vsum(p, v1, 1.0, 1.0, v1);
  vsum(p, v2, 1.0, 1.0, v2);
  za1 = fac*v1[0]*dist / (dist-v1[2]);
  za2 = fac*v1[1]*dist / (dist-v1[2]);
  zb1 = fac*v2[0]*dist / (dist-v2[2]);
  zb2 = fac*v2[1]*dist / (dist-v2[2]);
  zp[0] = 0.5*(za1+zb1);
  zp[1] = 0.5*(za2+zb2);
  *zr = (zb1-za1)*(zb1-za1) + (zb2-za2)*(zb2-za2);
  *zr = 0.5 * sqrt(*zr);
}

/* ----- readclusterdata ---- */
int readclusterdata(char infile[])
{
  FILE *fp;
  char str[257],token[81];
  char xxx [81];
  char *p;
  int l,i,nn;

  if(( fp = fopen( infile, "r" ) ) == NULL) return 0;
  
  fgets(str, 257, fp); // read first line of infile
  while(!feof(fp)) {
    l=strlen(str);
    str[l+1]='\0';
    strcpy(token,"SNOT");
    sscanf(str,"%s",token);   
    if ( p = strstr(str,"frame") ) { // read in one frame of atom data from FILE.mv // p = first character of "frame"
      if (nframe * nbas > FBMAX) rx("increase internal dimension FBMAX");
      if (nframe > NFRMAX) rx("increase internal dimension NFRMAX");
      p = p + 6;
      sprintf(frstr[nframe], "%-80s", p); // write frame number to frstr
      for (i=0; i<strlen(frstr[nframe]); i++) {
        if (frstr[nframe][i] == '\n') frstr[nframe][i]='\0';
      }
      for (i=0;i<nbas;i++) {
        nn = nframe * nbas + i;
        fgets(str, 257, fp);
        flags[nn] = -1;
        sscanf(str,"%f %f %f %i", &frame[0][nn], &frame[1][nn], &frame[2][nn], &flags[nn]);
      }
      nframe++;
      if ( nframe%50 == 0) {    
        sprintf (xxx, ": frame %d         ", nframe);
        showline (win, 10, 8, "Reading ", inmv, xxx);
        XFlush (dpy);
      }
      // fgets(str, 257, fp);
    } else {
      readclusterline(str,0);
    } 
    fgets(str, 257, fp);
  }
  sprintf(curf,"%s ", inf);
  return 1;
}

/* ----- readclusterline ------ */
  int readclusterline(char str[], int helpme)
  {
    char token[81], cname[81];
    float rval, bval, gval, gray;
    int l,n,i,nn, flag;

    l=strlen(str); if(l<1) return 0;
    strcpy(token, "SNOT");
    sscanf(str,"%s",token);   
    if(!strcmp(token,"SNOT")) return 0;  /* empty line */
    if(token[0]=='*') return 0;          /* comment -- no error */

    if (helpme) {
      if (abbrev(token,"spec",4))
        sprintf(gmsg,"Usage: spec label radius color  - define species");
      else if (abbrev(token,"atom",4))
        sprintf(gmsg,"Usage: atom label x y z (flag) (px py pz)- place atom at (x,y,z), optional flag and polarisation (x, y, z)");
      else if (abbrev(token,"bonds",5))
        sprintf(gmsg,"Usage: bonds pat1 pat2 min max rad color  - select bonds");
      else if (abbrev(token,"light",5))
        sprintf(gmsg,"Usage: light vx vy vz  - light along vector in bw mode)");
      else if (abbrev(token,"inc",3))
        sprintf(gmsg,"Usage: inc degrees  - angle increment for rotation");
      else if (abbrev(token,"dist",5))
        sprintf(gmsg,"Usage: dist d  - set distance for perspective");
      else if (abbrev(token,"frm",3))
        sprintf(gmsg,"Usage: frm n  - goto frame n");
      else if (abbrev(token,"step",4))
        sprintf(gmsg,"Usage: step n  - set step for frames");
      else if (abbrev(token,"gramp",5))
        sprintf(gmsg,"Usage: gramp slope [middle]  - set gray ramp in bw mode)");
      else if (abbrev(token,"scale",5))
        sprintf(gmsg,"Usage: scale x  - set overall scale factor");
      else if (abbrev(token,"rfac",4))
        sprintf(gmsg,"Usage: rfac x  - scale all sphere radii by x");
      else if (abbrev(token,"bfac",4))
        sprintf(gmsg,"Usage: bfac x  - scale all bond radii by x");
      else if (abbrev(token,"pos",3))
        sprintf(gmsg,"Usage: pos px py  - set position on page");
      else if (abbrev(token,"dpos",4))
        sprintf(gmsg,"Usage: dpos x  - set increment for position");
      else
        sprintf(gmsg,"No help available on %s", token);
      return 0;
    }

    if(!strcmp(token,"spec")) {
      sscanf(str,"%*s %s %f %n", spec[nspec].lab, &spec[nspec].rad, &n);
      strip (spec[nspec].cname, str+n); 
      nspec++; 
      if (nspec>NSPMAX) rx("increase internal dimension NSPMAX");
      return 2;
    }

    if(!strcmp(token,"atom")) {
      atom[nbas].pol[0] = 0;
      atom[nbas].pol[1] = 0;
      atom[nbas].pol[2] = 0;
      atom[nbas].flag = 1;
      sscanf(str,"%*s %s %f %f %f %i %f %f %f", atom[nbas].lab, 
             &atom[nbas].pos[0], &atom[nbas].pos[1], &atom[nbas].pos[2], &atom[nbas].flag,
             &atom[nbas].pol[0], &atom[nbas].pol[1], &atom[nbas].pol[2]); 
      nbas++; if (nbas>NAMAX) rx("increase internal dimension NAMAX");
      return 2;
    }

    if(!strcmp(token,"bonds")) {
      sscanf(str,"%*s %s %s %f %f %f %n", bonds[nbonds].lab1, 
        bonds[nbonds].lab2, &bonds[nbonds].min, &bonds[nbonds].max,
             &bonds[nbonds].rad, &n);
      strip (bonds[nbonds].cname, str+n); 
      nbonds++; 
      if (nbonds>NBTMAX) rx("increase internal dimension NBTMAX");
      return 2;
    }

    if(!strcmp(token,"tmat")) {
      sscanf(str, "%*s %f %f %f %f %f %f %f %f %f",
             &tmat[0][0], &tmat[0][1], &tmat[0][2], 
             &tmat[1][0], &tmat[1][1], &tmat[1][2], 
             &tmat[2][0], &tmat[2][1], &tmat[2][2]);
      return 1;
    }

    if(!strcmp(token,"dist") || !strcmp(token,"d")) {
      sscanf(str, "%*s %f", &dist0);
      return 1;
    }

    if(!strcmp(token,"inc")) {sscanf(str, "%*s %f", &dalfa); return 3;}

    if(!strcmp(token,"frm")) {
      sscanf(str, "%*s %d", &nn); 
      if (nn>nframe || nn<1) {
        sprintf (emsg, "No frame %d available", nn);
        return 0;
      }
      if (iframe==nn-1) return 0;
      iframe=nn-1;
      return 2;
    }

    if(!strcmp(token,"light")) {
      light[0]=light[1]=light[2]=0.0;
      sscanf(str,"%*s %f %f %f", &light[0], &light[1], &light[2]);
      if (light[0]*light[0]+light[1]*light[1]+light[2]*light[2]<.01) {
        sprintf (gmsg, "Use standard coloring");
        gmode=G_STD;;
        return 2;
      }
      if (color) {
        sprintf(emsg, "light: only works in b/w mode"); 
        return 0;
      }
      gmode=G_LIGHT;
      return 2;
    }

    if(!strcmp(token,"step"))  {sscanf(str,"%*s %d", &fstep); return 1;}
    if(!strcmp(token,"scale")) {
      sscanf(str,"%*s %f", &scale); 
      scale=scale*igs;
      return 1;
    }
    if(!strcmp(token,"rfac"))  {sscanf(str,"%*s %f", &radfac); return 1;}
    if(!strcmp(token,"bfac"))  {sscanf(str,"%*s %f", &bndfac); return 1;}
    if(!strcmp(token,"amp"))   {sscanf(str,"%*s %f", &amp); return 2;}
    if(!strcmp(token,"pos")) 
      {sscanf(str,"%*s %f %f", &taux, &tauy); return 1;}
    if(!strcmp(token,"dpos"))   {
      sscanf(str,"%*s %f", &dtaux); dtauy=dtaux; chginfo=1; return 0;}

    if(!strcmp(token,"gramp"))   {
      gslope=gz0=0; 
      sscanf(str, "%*s %f %f", &gslope, &gz0); 
      if (gslope*gslope<0.1) {
        sprintf (gmsg, "Use standard coloring");
        gmode=G_STD;
        return 2;
      }
      if (color) {
        sprintf(emsg, "gramp: only works in b/w mode"); 
        return 0;
      }
      gmode=G_RAMP; 
      return 2; 
    }

    if(!strcmp(token,"switches"))   {
      sscanf(str, "%*s %d %d %d %d %d %d %d %d %d",
             &usepixmap,&numbers,&grayvalues,&bline,&wire,
             &withbonds,&recenter,&pmode,&shadow);
      return 2;
    }

    sprintf(emsg,"Undefined command: %s", token); 
    if (startup) printf ("Cannot understand line: %s\n", str);
    return 0;
  }    

/* ----- writeclusterdata ---- */
void writeclusterdata(char outfile[], int svstep, int svrgb)
{
  FILE *fp;
  char nm[81];
  int i,nn,n,nfrm;
  time_t  ltime;
  char    timestr[41];

  if((fp = fopen (outfile,"w")) == NULL) {
    sprintf(emsg, "Cannot open file %s\n", outfile); return;}
  
  time(&ltime);              
  strcpy (timestr,  ctime(&ltime));
  timestr[24]=0;
  fprintf (fp, "* Saved %s from %s\n\n", timestr, inf);

  for (i=0;i<nbas;i++) 
    fprintf(fp, "atom %6s  %10.3f %10.3f %10.3f %2d\n",
            atom[i].lab, atom[i].pos[0], atom[i].pos[1], atom[i].pos[2], atom[i].flag);

  fprintf(fp, "\n");
  for (i=0;i<nspec;i++) 
    if (svrgb || reverse) 
      fprintf(fp, "spec %6s %10.3f   %.2f %.2f %.2f\n",
              spec[i].lab, spec[i].rad, spec[i].r, spec[i].g, spec[i].b);
    else
      fprintf(fp, "spec %6s %10.3f   %s\n",
              spec[i].lab, spec[i].rad, spec[i].cname);
  
  fprintf(fp, "\n");
  for(i=0;i<nbonds;i++) 
    if (svrgb || reverse) 
      fprintf(fp, "bonds %5s %5s %8.3f %8.3f %8.3f   %.2f %.2f %.2f\n",
              bonds[i].lab1, bonds[i].lab2, bonds[i].min, bonds[i].max,
              bonds[i].rad, bonds[i].r, bonds[i].g, bonds[i].b);
    else
      fprintf(fp, "bonds %5s %5s %8.3f %8.3f %8.3f   %s\n",
              bonds[i].lab1, bonds[i].lab2, bonds[i].min, bonds[i].max,
              bonds[i].rad, bonds[i].cname);
  
  fprintf(fp,"\ntmat %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",
          tmat[0][0], tmat[0][1], tmat[0][2], tmat[1][0], tmat[1][1],
          tmat[1][2], tmat[2][0], tmat[2][1], tmat[2][2]);

  fprintf(fp, "dist  %8.3f\n", dist0);
  fprintf(fp, "inc   %8.3f\n", dalfa);
  fprintf(fp, "scale %8.3f\nrfac %.2f\nbfac %.2f\n", scale,radfac,bndfac);
  fprintf(fp, "pos %8.3f %8.3f\n", taux, tauy);
  if (gmode==G_RAMP) 
    fprintf(fp, "gramp %8.3f %8.3f\n", gslope, gz0);
  if (gmode==G_LIGHT) 
    fprintf(fp, "light %8.3f %8.3f %8.3f\n", light[0], light[1], light[2]);
    
  fprintf(fp, "switches %d %d %d %d %d %d %d %d %d\n",
          usepixmap,numbers,grayvalues,bline,wire,
          withbonds,recenter,pmode,shadow);

  fclose(fp);
  sprintf(gmsg, "Saved in %s", outfile);

  if (nframe>1) {
    strext (nm, outfile, "mv", 1);
    if((fp = fopen (nm,"w")) == NULL) {
      sprintf(emsg,"Cannot open file %s\n", nm); return;}
    nfrm=1;
    for (i=1;i<nframe;i=i+svstep) {
      nfrm++;
      fprintf(fp,"frame %s\n", frstr[i]);
      for(n=0;n<nbas;n++) {
        nn=i*nbas+n;
        fprintf (fp, "%.3f %.3f %.3f %.3i", frame[0][nn], frame[1][nn], frame[2][nn], flags[nn]);
      }
      fprintf (fp, "\n\n");
    }
    fclose(fp);
    sprintf(gmsg, "Saved %d frames in %s and %s", nfrm, outfile, nm);
    
  }
  
  return;
}

/* ----- parse_all_colors ------- */
void parse_all_colors ()
{
  int i;
  
  for (i=0;i<nspec;i++) 
    parse_color (spec[i].cname, &spec[i].r, &spec[i].g, &spec[i].b,
                  &spec[i].gray);
                  
  for (i=0;i<nbonds;i++) 
    parse_color (bonds[i].cname, &bonds[i].r, &bonds[i].g, &bonds[i].b,
                  &bonds[i].gray);
}                  

/* ----- set_auto_colors ------ */
void set_auto_colors()
{
  int k, i, up1;
  char id[21];
  char *p;
  
  for (k=0; k<nspec; k++) {
    /* separate off the identifier and change to upper case */
    up1=isupper(spec[k].lab[0]);
    p = &spec[k].lab[0];
    i=0;
    while (*p != 0) {
      id[i]=toupper(*p); p++; i++;
      if (!isalpha(*p)) break;
      if (up1 && isupper(*p)) break;
    }
    id[i]=0;
    
    /* Set default colors here. Use rgb values or color name */
    
    strcpy(spec[k].cname, ".68 .85 .90"); 

    if (!strcmp(id,"H"))   strcpy(spec[k].cname, "1.0 0.2 0.2"); 
    if (!strcmp(id,"C"))   strcpy(spec[k].cname, "0.65 0.7 0.7");
/*    if (!strcmp(id,"O"))   strcpy(spec[k].cname, "0.2 0.2 1.0"); */
    if (!strcmp(id,"O"))   strcpy(spec[k].cname, "blue");
    if (!strcmp(id,"N"))   strcpy(spec[k].cname, "0.8 0.0 1.0");
    if (!strcmp(id,"P"))   strcpy(spec[k].cname, "0.0 0.8 0.0");
    if (!strcmp(id,"CL"))  strcpy(spec[k].cname, "0.9 0.9 0.7");
    if (!strcmp(id,"TI"))  strcpy(spec[k].cname, "0.2 1.0 1.0");
  }

  for (k=0; k<nbonds; k++) {
    strcpy(bonds[k].cname, "0.8 0.8 0.8");
  }
}

/* ----- ball_list ----- */
int ball_list(struct ballstr ball[], int jpr)
{
  int i,j,k,m;
  float top,bot,sp;
  
  for(i=0;i<nbas;i++) {
    k=-1;
    for(j=0;j<nspec;j++) if(!strcmp(spec[j].lab, atom[i].lab)) k=j;
    if(k==-1) {
      if (!startup) {
        sprintf(emsg,"Undefined species %s ", atom[i].lab); 
      } else 
        printf("Undefined species %s\n", atom[i].lab);
      continue;
    }
    for(m=0;m<3;m++) ball[i].pos[m]=atom[i].pos[m]+amp*atom[i].pol[m];
    ball[i].flag = atom[i].flag;
    ball[i].rad=spec[k].rad;
    ball[i].gray=spec[k].gray;
    ball[i].r=spec[k].r;
    ball[i].g=spec[k].g;
    ball[i].b=spec[k].b;
    strcpy(ball[i].lab,spec[k].lab);
    ball[i].special=0;
    ball[i].col=spec[k].col;
    if(spec[k].gray<-0.1) {ball[i].gray=1.0; ball[i].special=1;}
  }
  
  if (gmode==G_LIGHT) {
    top=-1000.0;
    bot= 1000.0;
    for (i=0;i<nbas;i++) {
      if (! ball[i].special) {
        sp = ball[i].pos[0]*light[0] + ball[i].pos[1]*light[1] + ball[i].pos[2]*light[2];
        if (sp>top) top=sp;
        if (sp<bot) bot=sp;
      }
    }
    for (i=0;i<nbas;i++) {
      if (! ball[i].special) {
        sp=ball[i].pos[0]*light[0]+ball[i].pos[1]*light[1]
          +ball[i].pos[2]*light[2];
        ball[i].gray=(sp-bot)/(top-bot);
      }
    }
  }

  return nbas;
}

/* ----- stick_list ------ */
int stick_list(struct ballstr ball[], struct stickstr stick[])
{
  int i,j,k,l,m,nbond,kb;
  float dis,dd;
  
  i=-1;
  for(k=0;k<nbas;k++)
    for(l=k+1;l<nbas;l++) {
      kb=-1;
      for(j=0;j<nbonds;j++) {
        if(match(ball[k].lab,bonds[j].lab1) &&
           match(ball[l].lab,bonds[j].lab2)) kb=j;
        if(match(ball[l].lab,bonds[j].lab1) &&
           match(ball[k].lab,bonds[j].lab2)) kb=j;
/*        if( (!strcmp(bonds[j].lab1,ball[k].lab)) &&
            (!strcmp(bonds[j].lab2,ball[l].lab)) ) kb=j;
        if( (!strcmp(bonds[j].lab1,ball[l].lab)) &&
            (!strcmp(bonds[j].lab2,ball[k].lab)) ) kb=j; */
      }
      if (kb>-1) {
        dis=0.0;
        for (m=0;m<3;m++) { 
          dd = ball[k].pos[m] - ball[l].pos[m];
          dis=dis+dd*dd; 
        }
        dis=alat*sqrt(dis);
        if( (dis>=bonds[kb].min) && (dis<=bonds[kb].max) ) {
          i++; if (i>NBMAX) rx("increase internal dimension NBMAX");
          stick[i].start = k;
          stick[i].end = l;
          stick[i].rad = bonds[kb].rad;
          stick[i].gray = bonds[kb].gray;
          stick[i].col = bonds[kb].col;
        }
      }
    }

  nbond=i+1;
  return nbond;
}

/* ----- duplicate_atoms ---- */
int duplicate_atoms( sh, helpme )
float sh[3][6];
int helpme;
{
  int k,l,iv,nbas1,ndup,fr,nn,nn1;
  float cx,cy,cz;

  if (helpme) {
    sprintf (gmsg,"Usage: dup vx vy vz  - duplicate shifted by vector");
    return 0;
  }

  ndup=0;
  for (iv=0;iv<6;iv++) {
    cx=sh[0][iv]; cy=sh[1][iv]; cz=sh[2][iv];
    if (cx*cx+cy*cy+cz*cz>0.001) ndup++;
  }

  if (ndup==0) {
    sprintf (emsg, "Cannot dup for (0,0,0)");
    return 0;
  }
  
  nbas1=nbas*(ndup+1);
  if (nframe*nbas1>FBMAX) {
    sprintf(emsg,"Cannot dup, internal dimension FBMAX too small");
    return 0;
  }
  if (nbas*(ndup+1)>NAMAX) {
    sprintf(emsg,"Cannot dup, internal dimension NAMAX too small");
    return 0;
  }
  
  for (iv=0;iv<ndup;iv++) {
    for (k=0;k<nbas;k++) {
      l=k+nbas*(iv+1);
      strcpy (atom[l].lab, atom[k].lab);
      atom[l].pos[0]=atom[k].pos[0]+sh[0][iv];
      atom[l].pos[1]=atom[k].pos[1]+sh[1][iv];
      atom[l].pos[2]=atom[k].pos[2]+sh[2][iv];
      atom[l].pol[0]=atom[k].pol[0];
      atom[l].pol[1]=atom[k].pol[1];
      atom[l].pol[2]=atom[k].pol[2];
			atom[l].flag = atom[k].flag;
    }
  }
  
  for (fr=nframe-1; fr>=0; fr-- ) {
    for ( k=0; k<nbas; k++ ) {
      nn = fr * nbas + k;
      nn1 = fr * nbas1 + k;
      frame[0][nn1]=frame[0][nn];
      frame[1][nn1]=frame[1][nn];
      frame[2][nn1]=frame[2][nn];
      flags[nn1]=flags[nn];
      for (iv=0;iv<ndup;iv++) {
        nn1=nn1+nbas;
        frame[0][nn1]=frame[0][nn]+sh[0][iv];
        frame[1][nn1]=frame[1][nn]+sh[1][iv];
        frame[2][nn1]=frame[2][nn]+sh[2][iv];
      }
    }
  }
  
  sprintf (gmsg, "Increased from %d to %d atoms", nbas, nbas1);
  nbas=nbas1;
  return 1;
}
   
/* ----- cut_atoms ---- */
int cut_atoms(cut, cut1, cut2, helpme)
float cut[3], cut1, cut2;
int helpme;
{
  int i,j,nbas1,fr,nn,nn1;
  float c2,s,fuzz;
  int ip[NAMAX];
  
  if (helpme) {
    sprintf(gmsg,"Usage: cut vx vy vz a b  - cut along vector at a and b");
    return 0;
  }

  fuzz=0.001;
  c2=cut[0]*cut[0]+cut[1]*cut[1]+cut[2]*cut[2];
  if (c2<0.01) {
    sprintf(emsg,"cut: invalid vector (%.2f,%.2f,%.2f)", 
            cut[0],cut[1],cut[2],cut1,cut2);
    return 0;
  }

  nbas1=0;
  for (i=0; i<nbas; i++) {
    s=atom[i].pos[0]*cut[0]+atom[i].pos[1]*cut[1]+atom[i].pos[2]*cut[2];
    if (c2*cut1-fuzz<s && s<c2*cut2+fuzz) {
      ip[nbas1]=i;
      nbas1=nbas1+1;
    }
  }
  if (nbas1==nbas) {
    sprintf (gmsg, "No atoms were cut");
    return 0;
  }

  for (i=0;i<nbas1;i++) {
    j=ip[i];
    atom[i].pos[0]=atom[j].pos[0];
    atom[i].pos[1]=atom[j].pos[1];
    atom[i].pos[2]=atom[j].pos[2];
    atom[i].pol[0]=atom[j].pol[0];
    atom[i].pol[1]=atom[j].pol[1];
    atom[i].pol[2]=atom[j].pol[2];
    atom[i].flag = atom[j].flag;
    strcpy(atom[i].lab,atom[j].lab);
  }   

  for (fr=0; fr<nframe; fr++) {
    for (i=0;i<nbas1;i++) {
      j=ip[i];
      nn  = fr*nbas+j;
      nn1 = fr*nbas1+i;
      frame[0][nn1]=frame[0][nn];
      frame[1][nn1]=frame[1][nn];
      frame[2][nn1]=frame[2][nn];
	    flags[nn1]=flags[nn];
    }
  }
  
  sprintf(gmsg, "Reduced from %d to %d atoms", nbas, nbas1); 
  nbas=nbas1;
  return 1;
}
  
/* ----- select bonds to plot ----- */
int selectbonds(int natom, struct ballstr ball[], float blen,
          float rad, float gray, struct stickstr stick[])
{
  int i,k,l,m,nbond;
  float dis,dd;
                 
  i=-1;
  for(k=0;k<natom;k++)
    for(l=k+1;l<natom;l++) {
      dis=0.0;
      for (m=0;m<3;m++) { 
        dd = ball[k].pos[m] - ball[l].pos[m];
        dis=dis+dd*dd; 
      }
      dis=sqrt(dis);
      if(dis<=blen) {
        i++;
        stick[i].start = k;
        stick[i].end = l;
        stick[i].rad = rad;
        stick[i].gray = gray;
      }
      
    }
  printf("bond   start  end    radius    gray\n");
  nbond=i+1;
  for(i=0;i<nbond;i++) 
    printf("%3d    %3d   %3d    %7.3f  %6.2f\n", 
           i+1, stick[i].start+1, stick[i].end+1, 
           stick[i].rad,stick[i].gray);
  return nbond;
}


/* ----- rotmat ------ */
void rotmat(int ixyz, float alfa)
{
  int i,j,k;
  float rot[3][3], w[3][3];

  switch (ixyz) {
  case 3:
    rot[0][0]=cos(alfa);  rot[0][1]=-sin(alfa); rot[0][2]=0.0;
    rot[1][0]=sin(alfa);  rot[1][1]=cos(alfa);  rot[1][2]=0.0;
    rot[2][0]=0.0;        rot[2][1]=0.0;        rot[2][2]=1.0;
    break;
  case 1:
    rot[0][0]=cos(alfa);  rot[0][1]=0.0;        rot[0][2]=sin(alfa);
    rot[1][0]=0.0;        rot[1][1]=1.0;        rot[1][2]=0.0;
    rot[2][0]=-sin(alfa); rot[2][1]=0.0;        rot[2][2]=cos(alfa);
    break;
  case 2:
    rot[0][0]=1.0;        rot[0][1]=0.0;        rot[0][2]=0.0;
    rot[1][0]=0.0;        rot[1][1]=cos(alfa);  rot[1][2]=-sin(alfa);
    rot[2][0]=0.0;        rot[2][1]=sin(alfa);  rot[2][2]=cos(alfa);
    break;
  default:
    return;
  }

  for(i=0;i<3;i++) for(j=0;j<3;j++) {
    w[i][j]=0.0;
    for(k=0;k<3;k++) w[i][j]=w[i][j]+rot[i][k]*tmat[k][j];
  }
   
  for(i=0;i<3;i++) for(j=0;j<3;j++) tmat[i][j]=w[i][j];
  return;
}
/* ----- eumat ------ */
void eumat(float alfa, float beta, float gama)
{
  int i,j,k;
  float gam[3][3], bet[3][3], alf[3][3], w[3][3];
/*  printf("Euler angles  %7.1f %7.1f %7.1f\n", alfa,beta,gama);*/

  gam[0][0]=1.0;        gam[0][1]=0.0;        gam[0][2]=0.0;
  gam[1][0]=0.0;        gam[1][1]=cos(gama);  gam[1][2]=-sin(gama);
  gam[2][0]=0.0;        gam[2][1]=sin(gama);  gam[2][2]=1.0;

  bet[0][0]=cos(beta);  bet[0][1]=0.0;        bet[0][2]=sin(beta);
  bet[1][0]=0.0;        bet[1][1]=1.0;        bet[1][2]=0.0;
  bet[2][0]=-sin(beta); bet[2][1]=0.0;        bet[2][2]=cos(beta);

  alf[0][0]=cos(alfa);  alf[0][1]=-sin(alfa); alf[0][2]=0.0;
  alf[1][0]=sin(alfa);  alf[1][1]=cos(alfa);  alf[1][2]=0.0;
  alf[2][0]=0.0;        alf[2][1]=0.0;        alf[2][2]=1.0;

  for(i=0;i<3;i++) for(j=0;j<3;j++) {
    w[i][j]=0.0;
    for(k=0;k<3;k++) w[i][j]=w[i][j]+bet[i][k]*alf[k][j];
  }
   
  for(i=0;i<3;i++) for(j=0;j<3;j++) {
    tmat[i][j]=0.0;
    for (k=0;k<3;k++) tmat[i][j]=tmat[i][j]+gam[i][k]*w[k][j];
  } 
  return;
}

/* ----- dbond ------- */
void dbond(float gray, float m1[6], float m2[6])
{
  float dfac=0.8;
  float rfac=2.0;

  float r,ax,ay,bx,by,a1,b1,dx,dy,alf,d1,bb,x;

  m2[4]=dfac*m2[4]+(1-dfac)*m1[4];
  m2[5]=dfac*m2[5]+(1-dfac)*m1[5];
  
  ax=m2[0];
  ay=m2[1];
  bx=m2[2];
  by=m2[3];

  b1=sqrt(bx*bx+by*by);
  a1=sqrt(ax*ax+ay*ay);
  r=rfac*b1;
  m2[0] = -r*ax/a1;
  m2[1] = -r*ay/a1;
  m2[2] = r*bx/b1;
  m2[3] = r*by/b1;
  
  dx=m2[4]-m1[4];
  dy=m2[5]-m1[5];

  if(m2[0]*dx<0 || m2[1]*dy<0) { 
    m1[0]=-m1[0]; m1[1]=-m1[1]; m1[2]=-m1[2]; m1[3]=-m1[3];
    m2[0]=-m2[0]; m2[1]=-m2[1]; m2[2]=-m2[2]; m2[3]=-m2[3];
  }

  d1=sqrt(dx*dx+dy*dy);
  bb=sqrt(m1[2]*m1[2]+m1[3]*m1[3]);
  if(r-bb<d1) {
    x=bb*d1/(r-bb);
    alf=asin(bb/x)*57.3;
    if(hardcopy)
      hardcopy_xdbond(gray, m1, m2, alf);
    else
      printf("PSWxdbond.. not changed yet\n");
/*      PSWxdbond (gray, m1, m2, alf); */
  }
  else {
    if(hardcopy)
      hardcopy_ydbond (gray, m1, m2, alf);
    else
      printf("PSWydbond.. not changed yet\n");
/*      PSWydbond (gray, m1, m2, alf); */
  }
}


/* ----- getframe ------ */
void getframe(struct ballstr ball[], int fnum)
{
  register int m,n,nn;
  float sum;

  if (recenter) {
    for (m=0;m<3;m++) {
      sum=0.0;
      for (n=0;n<nbas;n++) {
        nn=fnum*nbas+n;
        sum=sum+frame[m][nn];
      }
      center[m]=sum/nbas;
    }
  }
  
  for(n=0;n<nbas;n++) {
    nn=fnum*nbas+n;
    ball[n].pos[0] = frame[0][nn]-center[0];
    ball[n].pos[1] = frame[1][nn]-center[1];
    ball[n].pos[2] = frame[2][nn]-center[2];
    if ( flags[nn] == -1 ) {
      ball[n].flag = atom[n].flag;
    } else {
      ball[n].flag = flags[nn];
    }
  }
}


/* ----- putframe ------ */
void putframe(struct ballstr ball[], int fnum)
{
  int   n,nn;
  for(n=0;n<nbas;n++) {
    nn=fnum*nbas+n;
    frame[0][nn] = ball[n].pos[0];
    frame[1][nn] = ball[n].pos[1];
    frame[2][nn] = ball[n].pos[2];
    flags[nn] = ball[n].flag;
  }
}


/* ----- prframes  ------ */
void prframes()
{
  printf ("Number of frames: %d\n", nframe );
  printf ("frame %-5d <%s>\n", 1, frstr[0]);
  if (nframe>1) printf ("frame %-5d <%s>\n", 2, frstr[1]);
  if (nframe>2) printf ("frame %-5d <%s>\n", nframe, frstr[nframe-1]);

  return;
}

/* ----- draw_axes ----- */
void draw_axes ()
{
  float e0[3],e1[3],e2[3],z0[2],z1[2],z2[2];
  float fac,r0,r1,r2,tx,ty,zz[3];
  int m,i,j,i0;

  tx  =  (70       -midx)/PSFAC;
  ty  = -(igh-120  -midy)/PSFAC;
  fac = 30;

  for (m=0;m<3;m++) {
    e0[m]=tmat[m][0];
    e1[m]=tmat[m][1];
    e2[m]=tmat[m][2];
  }

  atompos(fac, e0, 1.0, z0, &r0);
  atompos(fac, e1, 1.0, z1, &r1);
  atompos(fac, e2, 1.0, z2, &r2);

  /* sort vectors (clumsily) back to front */
  zz[0]=e0[2];
  zz[1]=e1[2];
  zz[2]=e2[2];

  for (i=0;i<3;i++) {
    i0=0;
    if (zz[1]<zz[i0]) i0=1;
    if (zz[2]<zz[i0]) i0=2;
    if (i0==0) DrawArrow (tx, ty, z0[0]+tx, z0[1]+ty, r0, "100");
    if (i0==1) DrawArrow (tx, ty, z1[0]+tx, z1[1]+ty, r1, "010");
    if (i0==2) DrawArrow (tx, ty, z2[0]+tx, z2[1]+ty, r2, "001");
    zz[i0]=1e10;
  }
  
}

/* ----- bs_transform ------ */
void bs_transform(int natom, struct ballstr ball[])
{
  register int m,n;
  
  for(n=0;n<natom;n++) {
    for(m=0;m<3;m++)
      p[n][m]=tmat[m][0]*ball[n].pos[0]
        +tmat[m][1]*ball[n].pos[1]
          +tmat[m][2]*ball[n].pos[2];
  }
}

/* ----- bs_kernel ------ */
void bs_kernel(int natom, struct ballstr ball[], int nbond,
               struct stickstr stick[])
{
  int   flag[NAMAX], ip[NAMAX], j,k,kk,m,n,ibot;
  float br, xx,bx,by,rk,rkk,th1,th2,cth1,cth2,sth1,sth2;
  float w,ww,bb,aa,crit1,crit2,fac,beta,gray;
  char label[81]; 
  int   ib,note;
  float zp[NAMAX][2], zr[NAMAX];
  float q1[3], q2[3], b[3], d[3], bot, big;
  float fudgefac = 0.6,bmidx,bmidy,dd;
  float m1[6],m2[6];
  int nbx, ibx, jbx, kbx[100], abx[100], pbx[100], fbx[100]; /* faster bond search */
  dist=dist0;
  if (pmode==0 || pmode==1) dist=10000.0;
  d[0]=0.0; d[1]=0.0; d[2]=dist;
  fac=scale;

/* ------- sort atoms back to front ----- */
  for(k=0;k<natom;k++) flag[k]=0;
  for(n=0;n<natom;n++) {
    bot=1.0e10; ibot=0;
    for(k=0;k<natom;k++) 
      if(p[k][2]<bot && !flag[k]) { bot=p[k][2]; ibot=k; }
    ip[n]=ibot; flag[ibot]=1;
  }
  
/* ------- make list of sphere centers and radii ---- */
  big=1000000; if(nbas==0) big=0;
  xbot=ybot=big;
  xtop=ytop=-big;
  for(k=0;k<natom;k++) {
    atompos(fac, p[k], ball[k].rad, zp[k], &zr[k]);
    zr[k] = radfac * zr[k];
    if (zp[k][0]-zr[k] < xbot) xbot = zp[k][0] - zr[k];
    if (zp[k][0]+zr[k] > xtop) xtop = zp[k][0] + zr[k];
    if (zp[k][1]-zr[k] < ybot) ybot = zp[k][1] - zr[k];
    if (zp[k][1]+zr[k] > ytop) ytop = zp[k][1] + zr[k];
  }
/*  printf ("bounds x %.3f %.3f  y %.3f %.3f\n", xbot,xtop,ybot,ytop); */
  
/* ------- start loop over atoms; plot ball first ----- */ 
  for (n=0;n<natom;n++) {
    k=ip[n];
    if (ball[k].flag == 0) continue; // do not draw this atom
    rk=ball[k].rad;

    if(!ball[k].special) {
      beta=exp(gslope*(p[k][2]-gz0)*gslope) ;
      if (grayvalues) {
        if (gmode==G_RAMP) 
          gray=beta*ball[k].gray +(1-beta)*GRAY0;
        else
          gray=ball[k].gray;
      }
      else gray=1.0;
      if(hardcopy) 
        hardcopy_ball(gray,ball[k].r,ball[k].g,ball[k].b,
                      zp[k][0]+taux,zp[k][1]+tauy,zr[k]);
      else 
        DrawBall(gray, ball[k].col, zp[k][0]+taux, zp[k][1]+tauy, zr[k]); 

      if(numbers || coords) {
        if (numbers == 1) sprintf(label, "%d", k+1); 
        if (numbers == 2) sprintf(label, "%s", ball[k].lab);
        if (coords) 
          sprintf(label, "(%.2f, %.2f, %.2f)", 
                  ball[k].pos[0]+center[0], ball[k].pos[1]+center[1], ball[k].pos[2]+center[2]);
        if(hardcopy)
          hardcopy_label(zp[k][0]+taux,zp[k][1]+tauy, label); 
        else
        LabelBG(zp[k][0]+taux,zp[k][1]+tauy-2, 0.0, 1.0, label);  
      }        
    }

/*  ------ make list of bonds to this atom ----- */
    if (!withbonds) continue;
    nbx=0; 
    for (j=0;j<nbond;j++) {
      // do not draw bonds if either terminating atom is not displayed
      if ( ball[stick[j].start].flag == 0 || ball[stick[j].end].flag == 0) continue;
      if(k == stick[j].start) { 
        kbx[nbx]=j;
        abx[nbx]=stick[j].end;
        nbx++; 
      }
      else if(k == stick[j].end) { 
        kbx[nbx]=j; 
        abx[nbx]=stick[j].start;
        nbx++; 
      }
    }
    if (nbx==0) continue;

    for (m=0;m<nbx;m++) fbx[m]=0;   /* sort mini-list */
    for (m=0;m<nbx;m++) {
      bot=1.0e10; ibot=0;
      for (j=0;j<nbx;j++) 
        if (p[abx[j]][2]<bot && !fbx[j]) { bot=p[abx[j]][2] ; ibot=j; }
      pbx[m]=ibot; fbx[ibot]=1;
    }

/*  ------ inner loop over bonds ----- */
    for (ibx=0;ibx<nbx;ibx++) {
      jbx=pbx[ibx];
      kk=abx[jbx];
      ib=kbx[jbx];
      if (ib<0) printf("this cannot happen\n");
      rkk = ball[kk].rad;  
      
      if (ib>=0) {
        br = bndfac*stick[ib].rad;
        bx = zp[kk][0]-zp[k][0];
        by = zp[kk][1]-zp[k][1];
        xx = sqrt(bx*bx+by*by);
        if ( xx*xx < 0.0001 ) continue;
        bx=bx/xx;
        by=by/xx;
        vsum(d, p[k],  1.0, -1.0, q1);
        vsum(d, p[kk], 1.0, -1.0, q2);
        vsum( p[kk], p[k], 1.0, -1.0, b);
        cth1 =  sp(q1,b) / sqrt(sp(q1,q1)*sp(b,b));
        cth2 = -sp(q2,b) / sqrt(sp(q2,q2)*sp(b,b));
        th1=acos(cth1);
        th2=acos(cth2);
        crit1 = asin(br/rk) * fudgefac;
        if (crit1<0.0) crit1=0.0;
        crit2 = asin(br/rkk) * fudgefac;
        if (crit2<0.0) crit2=0.0;
        note=0;
        if(th2-0.5*PI>crit2 && k<kk) note=1;
        if(th1-0.5*PI<crit1 && k>kk) note=2;
/*        if(th2-0.5*PI>crit2 && n<nn) note=1;
        if(th1-0.5*PI<crit1 && n>nn) note=2; */

/* ------- plot a stick ------ */
        if (note==1 || note==2) {
          w = sqrt(rk*rk - br*br);
          sth1 = sqrt(1.0-cth1*cth1);
          ww = w*sth1*zr[k]/rk;
          bb = br*zr[k]/rk;
          aa = br*cth1*zr[k]/rk;
          m1[0]=bx*aa;   m1[1]=by*aa;
          m1[2]=-by*bb;  m1[3]=bx*bb;
          m1[4]=zp[k][0]+bx*ww+taux;
          m1[5]=zp[k][1]+by*ww+tauy; 
          w = sqrt(rkk*rkk-br*br);
          sth2 = sqrt(1.0-cth2*cth2);
          ww=w*sth2*zr[kk]/rkk;
          bb=br*zr[kk]/rkk;
          aa=br*cth2*zr[kk]/rkk;
          m2[0]=bx*aa;   m2[1]=by*aa;
          m2[2]=-by*bb;  m2[3]=bx*bb;
          m2[4]=zp[kk][0]-bx*ww+taux; 
          m2[5]=zp[kk][1]-by*ww+tauy; 

          beta=exp(gslope*(0.5*(p[k][2]+p[kk][2])-gz0)*gslope);
          if (grayvalues) {
            if (gmode==G_STD) 
              gray = stick[ib].gray;
            else if (gmode==G_LIGHT)    
              gray=0.5*(ball[k].gray+ball[kk].gray);  
            else if(gmode==G_RAMP) 
              gray=beta*stick[ib].gray+(1-beta)*GRAY0;
          }
          else gray=1.0;
          
          if (grayvalues && bline && gray>0.7) gray=0.7;
          if (!grayvalues && bline) gray=0.0;
          if (grayvalues && wire) gray=0.0;
          if (bline) gray=0.0;    /* overrides... black if lines */
          if (ball[k].special) 
            dbond(1.0, m2, m1);
          else if(ball[kk].special) 
            dbond(1.0, m1, m2); 
          else {
            if(hardcopy)
              hardcopy_stick(gray, m1, m2);
            else
/*              printf ("draw stick %d  col0=%d\n", ib, stick[ib].col); */
              DrawStick(gray, stick[ib].col, m1, m2);   
          }

/*  next part writes bond lengths onto the sticks */
          if (bondnums) {
            bmidx=0.5*(zp[k][0]+zp[kk][0])+taux;
            bmidy=0.5*(zp[k][1]+zp[kk][1])+tauy;
            dd=0;
            for (m=0;m<3;m++)
              dd=dd+pow(ball[k].pos[m]-ball[kk].pos[m],2);
            dd=sqrt(dd);
            sprintf(label, "%.2f", dd);
            if(hardcopy)
              hardcopy_label (bmidx, bmidy, label);
            else
            LabelBG (bmidx, bmidy, 0.0, 1.0, label); 
          }

        }
      } /* if (ib!=0) */
    }   /* end loop over nn */
  }     /* end loop over n */

  if (showaxes) draw_axes();

}
