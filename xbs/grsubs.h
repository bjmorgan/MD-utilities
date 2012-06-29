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

/*  procedures for graphics */

#define MAXRGB   65535

#define HIGHGRAY 65535     /* top of gray scale */
#define LOWGRAY      0     /* bottom of gray scale */


/* ------- parse_color ------- */
int parse_color (str, rval, gval, bval, gray)
char str[];
float *rval, *gval, *bval, *gray;
{
  int isname,i,j,j0,rc;
  XColor color,exact;
  float r,g,b;
  char cname[81];

  rc=1;
  strcpy (cname, str);
  j0=0;
  isname=0;
  for (j=0;j<strlen(cname);j++) {
    if (isalpha(cname[j])) isname=1;
    if (cname[j] != ' ') j0=j;
  }
  
  if (isname) {
    if (j0>0) cname[j0+1]=0;
    if (j0>0) str[j0+1]=0;
    if (XLookupColor (dpy, cmap, cname, &exact, &color)) {
      r=exact.red   / (float) MAXRGB;
      g=exact.green / (float) MAXRGB;
      b=exact.blue  / (float) MAXRGB;
    } else {
/*      printf("parse_color: Cannot look up color <%s>\n", cname); */
      r=g=b=1.0;
      rc=0;
    }
  }
  else {
    r=b=g=999.9;
    sscanf (cname, "%f %f %f", &r, &g, &b);
    if (r>999.0) r=0.0;
    if (g>999.0) g=r;
    if (b>999.0) b=r;
  }

  if (reverse) { r=1.0-r; g=1.0-g; b=1.0-b; }
  *rval=r;
  *gval=g;
  *bval=b;
  *gray=(r+g+b)/3.0;
/*  printf ("parse_color: rgb %.3f %.3f %.3f  gray %.3f  <%s>\n", 
          *rval, *gval, *bval, *gray, cname);  */
  return rc;
}


/* ------- GetColorGC ------- */
int GetColorGC (gpixel)
unsigned long gpixel;
{
  int i,j,j0,icol;

  j0=-1;
  for (j=0;j<ncol;j++) if (gpx[j]==gpixel) j0=j;
  
  if (j0 >=0) {
/*    printf("GetColorGC: already have pixel %d in slot %d\n",gpixel,j0);  */
    return j0;
  }

  if (ncol>=NCOL) {
    sprintf(emsg, "no more room for colors.. increase dimension NCOL"); 
    return 0;
  }
  
  /* find a slot and make a graphics context */
  icol=ncol;
  for (j=0;j<ncol;j++)
    if (gpx[j] == -1) icol=j;
 
/*  if (icol<ncol) 
    printf ("Make GC in old slot %d (%d) for pixel %d\n", icol,ncol,gpixel);
  else
    printf ("Make GC in new slot %d (%d) for pixel %d\n", icol,ncol,gpixel); */

  if (icol==ncol) ncol++;
  graygc[icol] = XCreateGC (dpy, win, 0, 0);
  XSetForeground (dpy, graygc[icol], gpixel); 
  gpx[icol]=gpixel;
  return icol;
  
}

/* ------- FreeColorGC: free GC if not needed ------- */
FreeColorGC (icol)
int icol;
{
  int i,num;

  num=0;
  for (i=0;i<nspec;i++) 
    if (spec[i].col == icol) num++;
  for (i=0;i<nbonds;i++) 
    if (bonds[i].col == icol) num++;

  if (num==0) {
    XFreeGC (dpy, graygc[icol]);
    gpx[icol]=-1;
  }
  
}

/* ------- SetColors ------- */
void SetColors ()
{
  int i,j,j0;
  int colptr;
  unsigned long gpixel,bground;
  XColor color,col1;
  
  bground = WhitePixel (dpy, screen);
  
  ncol=0;
  for (i=0;i<nspec+nbonds;i++) {
    if (i<nspec) {
      col1.red   = spec[i].r * MAXRGB;
      col1.green = spec[i].g * MAXRGB;
      col1.blue  = spec[i].b * MAXRGB;
    }
    else {
      col1.red   = bonds[i-nspec].r * MAXRGB;
      col1.green = bonds[i-nspec].g * MAXRGB;
      col1.blue  = bonds[i-nspec].b * MAXRGB;
    }

    color=col1;
    if (XAllocColor(dpy, cmap, &color)) {
      gpixel = color.pixel;
    }
    else {
      printf("Cannot allocate %6d %6d %6d, use background instead\n", 
             col1.red, col1.blue, col1.green);
      gpixel=bground;
    }
    
    colptr = GetColorGC(gpixel);
    if (i < nspec) 
      spec[i].col=colptr;
    else 
      bonds[i-nspec].col=colptr;
  }
}


/* ------- NewSpecColor ------- */ 
int NewSpecColor (pat, cname, helpme)
char pat[],cname[];
int helpme;
{
  int i,icol,nmatch;
  unsigned long gpixel;
  float f;
  char list[151];
  XColor col1,col2;

  if (helpme) {
    sprintf (gmsg, "Usage: color pattern [color]  - query or set color");
    return 0;
  }
  
  if (strlen(pat)==0) {
    sprintf (emsg, "color: no pattern specified");
    return 0;
  }

  nmatch=0;
  strcpy (list, "color");
  for (i=0;i<nspec;i++) {
    if (match(spec[i].lab, pat)) {
      nmatch++;
      strcat (list, " ");
      strcat (list, spec[i].lab);
      if (strlen(cname)==0) {
        sprintf (gmsg, "Species %s has color <%s> rgb %.2f %.2f %.2f", 
                 spec[i].lab, spec[i].cname, spec[i].r, spec[i].g, spec[i].b); 
        return 0;
      }
      
      if (! parse_color(cname, &spec[i].r, &spec[i].g, 
                        &spec[i].b, &spec[i].gray)) {
        sprintf (emsg, "Cannot identify color <%s>", cname);
        return 0;
      }

      strcpy (spec[i].cname, cname);
      if (color) {
        col1.red   = spec[i].r*MAXRGB;
        col1.green = spec[i].g*MAXRGB;
        col1.blue  = spec[i].b*MAXRGB;
        col2=col1;
        if (XAllocColor(dpy, cmap, &col1)) 
          gpixel = col1.pixel;
        else {
          sprintf(emsg, "Cannot allocate color <%s>", cname);
          return 0;
        }
        
        icol=spec[i].col;
        spec[i].col=-1;
        FreeColorGC (icol);
        icol=GetColorGC (gpixel);
        spec[i].col=icol;
      }
    }
  }

  if (nmatch==0 && (strlen(cname)==0)) {
    sprintf (emsg, "No species matches \'%s\'", pat);
    return 0;
  }

  f=MAXRGB;
  sprintf (gmsg, "%s <%s> rgb %.2f %.2f %.2f shown as %.2f %.2f %.2f",
           list, cname, col2.red/f,col2.green/f,col2.blue/f,
           col1.red/f,col1.green/f,col1.blue/f);  
  if (nmatch==0) return 0;
  return 1;
}
        
/* ------- SetSmoothGrays --- --- */
void SetSmoothGrays ()
{
  int i,allocated[NCOL];
  int lastgnum;
  unsigned long gpixel;
  XColor color;
  int g1,g2,gnum;
  float dg;

  ngray=21;
  g1=HIGHGRAY;
  g2=LOWGRAY;
  dg=(g2-g1)/(ngray-1.0);

  lastgnum=0;
  for (i=0;i<ngray;i++) {
    gnum = g1 + i*dg; 
    color.red=color.blue=color.green=gnum;
    allocated[i]=XAllocColor(dpy, cmap, &color);
    if (allocated[i]) {
      lastgnum=gnum;
    }
    else {
      printf("%7d: cannot allocate color.. use %d\n", gnum, lastgnum);
      color.red=color.blue=color.green=lastgnum;
      XAllocColor(dpy, cmap, &color);
    }
    graygc[i] = XCreateGC (dpy, win, 0, 0);
    gpixel = color.pixel;
    XSetForeground (dpy, graygc[i], gpixel); 
  }
}

/* ------- SetStippled4x4 ------- */
void SetStippled4x4 ()
{
  int i,j,n,depth;
  int gnum1,gnum2;
  unsigned long fg,bg,lg,dg,g1,g2;
  XColor color1,color2;
  Pixmap tile;
  int tw=4, th=4;
  static char t0[]  = {0x00, 0x00, 0x00, 0x00};
  static char t2[]  = {0x01, 0x00, 0x04, 0x00}; 
  static char t4[]  = {0x01, 0x04, 0x01, 0x04};
  static char t6[]  = {0x05, 0x02, 0x05, 0x08};
  static char t8[]  = {0x05, 0x0a, 0x05, 0x0a};
  static char t10[] = {0xfa, 0xfd, 0xfa, 0xf7};
  static char t12[] = {0xfe, 0xfb, 0xfe, 0xfb};
  static char t14[] = {0xfe, 0xff, 0xfb, 0xff};
  
  depth  = XDefaultDepth (dpy, screen);
  fg = BlackPixel (dpy, screen);
  bg = WhitePixel (dpy, screen);

  gnum1 = 0.33 * HIGHGRAY + 0.67 * LOWGRAY;
  color1.red=color1.blue=color1.green=gnum1;
  if (XAllocColor(dpy, cmap, &color1)) dg = color1.pixel;
  else printf("Could not allocate gray value %d\n", gnum1);

  gnum2 = 0.67 * HIGHGRAY + 0.33 * LOWGRAY;
  color2.red=color2.blue=color2.green=gnum2;
  if (XAllocColor(dpy, cmap, &color2)) lg = color2.pixel;
  else printf("Could not allocate gray value %d\n", gnum2);
/*  printf("Gray values: %d %d\n", color1.red, color2.red); */
  
  n=0;
  for (j=0; j<3; j++) {
    if (j==2) { g1=fg; g2=dg; }
    if (j==1) { g1=dg; g2=lg; }
    if (j==0) { g1=lg; g2=bg; }
    for (i=0; i<8; i++) {
      graygc[n] = XCreateGC (dpy, win, 0, 0);
      if (i==0) tile=XCreatePixmapFromBitmapData(dpy,win,t0,tw,th,g1,g2,depth);
      if (i==1) tile=XCreatePixmapFromBitmapData(dpy,win,t2,tw,th,g1,g2,depth);
      if (i==2) tile=XCreatePixmapFromBitmapData(dpy,win,t4,tw,th,g1,g2,depth);
      if (i==3) tile=XCreatePixmapFromBitmapData(dpy,win,t6,tw,th,g1,g2,depth);
      if (i==4) tile=XCreatePixmapFromBitmapData(dpy,win,t8,tw,th,g1,g2,depth);
      if (i==5) tile=XCreatePixmapFromBitmapData(dpy,win,t10,tw,th,g1,g2,depth);
      if (i==6) tile=XCreatePixmapFromBitmapData(dpy,win,t12,tw,th,g1,g2,depth);
      if (i==7) tile=XCreatePixmapFromBitmapData(dpy,win,t14,tw,th,g1,g2,depth);
      XSetTile (dpy, graygc[n], tile);
      XSetFillStyle (dpy, graygc[n], FillTiled);
      n++;
    }
  }
  ngray=n;
}

/* ------- SetStipple4x6    ------- */
void SetStippled4x6 ()
{
  int i,depth;
  unsigned long fg,bg,gr;
  XColor color;
  Pixmap tile;
  int tw=4, th=6;
  static char t0[]  = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
  static char t6[]  = {0x01, 0x04, 0x01, 0x04, 0x01, 0x04};
  static char t12[] = {0x05, 0x0a, 0x05, 0x0a, 0x05, 0x0a};
  static char t18[] = {0xfe, 0xfb, 0xfe, 0xfb, 0xfe, 0xfb};
  static char t24[] = {0xff, 0xff, 0xff, 0xff, 0xff, 0xff};
  
  depth  = XDefaultDepth (dpy, screen);
  fg = BlackPixel (dpy, screen);
  bg = WhitePixel (dpy, screen);

  if (XAllocColor(dpy, cmap, &color)) { 
    gr = color.pixel;
    printf(" gray is %d %d %d\n", color.red, color.blue, color.green);
  }
  else printf("Could not allocate gray value\n");

  ngray=9;
  for (i=0;i<ngray;i++) {
    graygc[i] = XCreateGC (dpy, win, 0, 0);

    if (i==8) {
      tile=XCreatePixmapFromBitmapData (dpy,win,t0,tw,th,gr,bg,depth);
      XSetTile (dpy, graygc[i], tile);
    }
    if (i==7) {
      tile=XCreatePixmapFromBitmapData (dpy,win,t6,tw,th,gr,bg,depth);
      XSetTile (dpy, graygc[i], tile);
    }
    if (i==6) {
      tile=XCreatePixmapFromBitmapData (dpy,win,t12,tw,th,gr,bg,depth);
      XSetTile (dpy, graygc[i], tile);
    }
    if (i==5) {
      tile=XCreatePixmapFromBitmapData (dpy,win,t18,tw,th,gr,bg,depth);
      XSetTile (dpy, graygc[i], tile);
    }

    if (i==4) {
      tile=XCreatePixmapFromBitmapData (dpy,win,t0,tw,th,fg,gr,depth);
      XSetTile (dpy, graygc[i], tile);
    }
    if (i==3) {
      tile=XCreatePixmapFromBitmapData (dpy,win,t6,tw,th,fg,gr,depth);
      XSetTile (dpy, graygc[i], tile);
    }
    if (i==2) {
      tile=XCreatePixmapFromBitmapData (dpy,win,t12,tw,th,fg,gr,depth);
      XSetTile (dpy, graygc[i], tile);
    }
    if (i==1) {
      tile=XCreatePixmapFromBitmapData (dpy,win,t18,tw,th,fg,gr,depth);
      XSetTile (dpy, graygc[i], tile);
    }
    if (i==0) {
      tile=XCreatePixmapFromBitmapData (dpy,win,t24,tw,th,fg,gr,depth);
      XSetTile (dpy, graygc[i], tile);
    }

    XSetFillStyle (dpy, graygc[i], FillTiled);
  }
}

/* ------ ChooseColor ------ */
GC ChooseColor (gray, col0)
float gray;
int col0;
{
  int igray;

  if (! grayvalues) return gcbg;
  if (color) {
    return graygc[col0];
  } 
  else {
    igray=(1.0-gray)*ngray-0.49;
    if (igray>ngray-1) igray=ngray-1;
    if (igray<0) igray=0;
    return graygc[igray];
  }
}

/* ------ showline ------ */
void showline(drawable, x, y, s1, s2, s3)
Drawable      drawable;
int           x,y;
char          s1[],s2[],s3[];
{
  char sx[1001];
  strcpy (sx, s1);
  strcat (sx, s2);
  strcat (sx, s3);
  XFillRectangle(dpy, drawable, gcbg, x, igh-y-9, x+6*strlen(sx), 12); 
  XDrawString (dpy, drawable, labelgc, x, igh-y, sx, strlen(sx)); 

}     

/* ------ clearline ------ */
void clearline(drawable, x, y)
Drawable      drawable;
int           x,y;
{
  XFillRectangle(dpy, drawable, gcbg, x, igh-y-11, igw-x-20, 14); 
}     

/* ------ DrawArrow ------ */
void DrawArrow(x0, y0, x1, y1, rad, str)
float x0,y0,x1,y1,rad;
char str[];
{
  int xx1,yy1,xx0,yy0,xx,yy,dd;

  xx1=midx+PSFAC*x1;
  yy1=midy-PSFAC*y1;
  xx0=midx+PSFAC*x0;
  yy0=midy-PSFAC*y0;
  xx=0.9*xx0+0.1*xx1;
  yy=0.9*yy0+0.1*yy1;
  XDrawLine (dpy, drw, shadowgc, xx, yy, xx1, yy1);
  XDrawLine (dpy, drw, gc, xx0, yy0, xx1, yy1);

  xx=xx1+0.3*(xx1-xx0);
  yy=yy1+0.3*(yy1-yy0)+4;
  dd=strlen(str)*6;
  XDrawLine   (dpy, drw, labbggc, xx+2-dd/2, yy-4, xx+dd/2-2, yy-4);
  XDrawString (dpy, drw, labelgc, xx-dd/2, yy, str, strlen(str)); 

  
}


/* ------ DrawBall ------ */
void DrawBall(gray, col0, x, y, rad)
float gray, x, y, rad;
int   col0;
{
  int xx,yy,rr,icol;
  GC gcfill;
  
  rr=PSFAC*rad;
  xx=midx+PSFAC*x;
  yy=midy-PSFAC*y;
  if (shadow) 
    XDrawArc (dpy, drw, shadowgc, xx-rr, yy-rr, 2*rr, 2*rr, 0, 360*64);
  if (! wire) {
    gcfill=ChooseColor(gray,col0);
    XFillArc (dpy, drw, gcfill, xx-rr, yy-rr, 2*rr, 2*rr, 0, 360*64);
  }
  XDrawArc (dpy, drw, gc, xx-rr, yy-rr, 2*rr, 2*rr, 0, 360*64);
}

/* ------ DrawStick ------ */
void DrawStick(gray, col0, m1, m2)
float gray, m1[6], m2[6];
int col0;
{
  int x1,y1,x2,y2,i,igray;
  float x,y;
  XPoint pp[NPOINTS*2+1];
  GC gcfill;

  if (bline) {
    x1=midx+PSFAC*m1[4];
    y1=midy-PSFAC*m1[5];
    x2=midx+PSFAC*m2[4];
    y2=midy-PSFAC*m2[5];
    if (shadow) 
      XDrawLine (dpy, drw, shadowgc, x1, y1, x2, y2);
    XDrawLine (dpy, drw, gc, x1, y1, x2, y2);
    return;
  }

  for (i=0;i<NPOINTS;i++) {
    x=m1[0]*arc[i][0]+m1[2]*arc[i][1]+m1[4];
    y=m1[1]*arc[i][0]+m1[3]*arc[i][1]+m1[5];
    pp[i].x=midx+PSFAC*x;
    pp[i].y=midy-PSFAC*y;
  }
  for (i=0;i<NPOINTS;i++) {
    x=-m2[0]*arc[i][0]+m2[2]*arc[i][1]+m2[4];
    y=-m2[1]*arc[i][0]+m2[3]*arc[i][1]+m2[5];
    pp[2*NPOINTS-i-1].x=midx+PSFAC*x;
    pp[2*NPOINTS-i-1].y=midy-PSFAC*y;
  }
  pp[2*NPOINTS]=pp[0];

  if (shadow) 
    XDrawLines (dpy, drw, shadowgc, &pp[0], 2*NPOINTS+1, CoordModeOrigin);
  if (! wire) {
    gcfill=ChooseColor(gray,col0);
    XFillPolygon (dpy, drw, gcfill, &pp[0], 2*NPOINTS+1, 
                  Nonconvex, CoordModeOrigin); 
  }
  XDrawLines (dpy, drw, gc, &pp[0], 2*NPOINTS+1, CoordModeOrigin);

}

/* ------ LabelBG ------ */
void LabelBG (x, y, g1, g2, str)
float x,y,g1,g2;
char str[];
{
  int xx,yy,dd;
  
  xx=midx+PSFAC*x;
  yy=midy-PSFAC*y;
/*  dd=strlen(str)*6+5; */
/*  XFillRectangle(dpy, drw, gcbg, xx-2, yy-10, dd, 13); */
  dd=strlen(str)*6;

  XDrawLine (dpy, drw, labbggc, xx+2-dd/2, yy-4, xx+dd/2-2, yy-4);


  XDrawString (dpy, drw, labelgc, xx-dd/2, yy, str, strlen(str)); 

}



