/* ------- subroutines to set ball and stick mode ---------  */
#include <stdlib.h>

void HCballFull()
{
  if (color) 
    fprintf(outfp,
            "/BALL  { 0 360 arc"
            "    gsave setrgbcolor fill grestore stroke } bind def\n");
  else
    fprintf(outfp,
            "/BALL  { 0 360 arc"
            "    gsave setgray fill grestore stroke } bind def\n");
}

void HCballWire()
{
  fprintf(outfp,
          "/BALL  { 0 360 arc pop stroke } bind def\n");
}

void HCballShadow(float gray, float wid)
{
  fprintf(outfp,
          "/BALL  { 0 360 arc\n "
          "   gsave %.2f setgray %.2f setlinewidth stroke grestore\n"
          "   gsave setgray fill grestore\n"
          "   stroke } bind def\n", gray, wid);
}

void HCstickFull()
{
  fprintf(outfp,
          "/STICK { CP gsave setgray fill grestore stroke } bind def\n");
}

void HCstickShadow(float gray, float wid)
{
  fprintf(outfp,
          "/STICK { CP\n"
          "gsave %.2f setgray %.2f setlinewidth stroke grestore\n"
          "gsave setgray fill grestore\n"
          "stroke } bind def\n", gray, wid);
}

void HCstickWire()
{
  fprintf(outfp,
          "/STICK { CPW pop stroke } bind def\n");
}

void HCstickLine()
{
  fprintf(outfp,
          "/STICK { CP0 gsave setgray stroke grestore newpath} bind def\n");

}

void HCstickLineShadow(float gray, float wid)
{
  fprintf(outfp,
          "/STICK { CP0 \n"
          "gsave %.2f setgray %.2f setlinewidth stroke grestore\n"
          "gsave setgray stroke grestore \n"
          "newpath} bind def\n", gray, wid);
}

void HCdbond()
{
  fprintf(outfp, "/XDB { } bind def\n");
}

void HCdbondShadow(float gray, float wid)
{
  fprintf(outfp,
          "/XDB { gsave %.2f setgray %.2f setlinewidth\n"
          "stroke grestore } bind def\n", gray, wid);
}

/* ----- hardcopy_init ------ */
void hardcopy_init(char outname[])
{
  if((outfp = fopen (outname,"w")) == NULL) {
    printf("Cannot open output file\n"); exit (1);}

  fprintf(outfp,
          "%%!PS-Adobe-3.0\n"
          "%%%%Creator: bs\n"
          "%%%%Pages: 1\n"
          "%%%%EndComments\n"
          "%%%%BeginPageSetup\n"
          "%%%%Page: 1 1\n"
          " 15 15 translate\n"
          "%%%%PageOrientation: Landscape\n"
          " 90 rotate 0 -580 translate\n"
          "/CP { %% Full cylinder path. Args: [CM1] [CM2]\n"
          "  exch matrix currentmatrix exch concat 0 0 1 90 270 arc\n"
          "  setmatrix matrix currentmatrix exch concat 0 0 1 -90 90 arc\n"
          "  closepath setmatrix } bind def\n"
          "/CP0 {\n"
          "  exch matrix currentmatrix exch concat 0 0 moveto\n"
          "  setmatrix matrix currentmatrix exch concat\n"
          "  0 0 lineto setmatrix } bind def\n"
          "/CPW {\n"
          "exch matrix currentmatrix exch concat 0 0 1 90 630 arc\n"
          "setmatrix matrix currentmatrix exch concat\n"
          "0 0 1 -90 450 arc closepath setmatrix } bind def\n"
          "/LAB { gsave 3 copy /Courier 9 selectfont\n"
          "1 setgray 10 setlinewidth 1 setlinecap 3 add exch 2 add exch\n"
          "moveto stringwidth exch 4 sub exch rlineto stroke\n"
          "0 setgray moveto show grestore} bind def\n");

/*          "/LAB { gsave 1 setgray 2 copy 1 sub exch 1 sub exch 20 8 rectfill\n"
          "   /Courier 9 selectfont 0 setgray\n"
          "   moveto show grestore} bind def\n"); */


  if (wire) {
    HCballWire();
    HCstickWire();
  }
  else {
    if (shadow) {
      HCballShadow(1.0, 4.0);
      if (bline) HCstickLineShadow(1.0, 4.0);
      else       HCstickShadow(1.0, 4.0);
    }
    else {
      HCballFull();
      if (bline) HCstickLine();
      else       HCstickFull();
    }
  }
  if (shadow) HCdbondShadow(1.0, 4.0);
  else        HCdbond();
  fprintf(outfp, "%%%%EndPageSetup\n");
}


/* ----- hardcopy_redefine ------ */
void hardcopy_redefine()
{
  if (wire) {
    HCballWire();
    HCstickWire();
  }
  else {
    if (shadow) {
      HCballShadow(1.0, 4.0);
      if (bline) HCstickLineShadow(1.0, 4.0);
      else       HCstickShadow(1.0, 4.0);
    }
    else {
      HCballFull();
      if (bline) HCstickLine();
      else       HCstickFull();
    }
  }
  if (shadow) HCdbondShadow(1.0, 4.0);
  else        HCdbond();
}

/* ------ hardcopy_ball, hardcopy_stick ----- */
void hardcopy_ball(float gray, float r, float g, float b, float x, float y, float rad)
{
  if (color)
    fprintf(outfp, "%.2f %.2f %.2f %8.1f%8.1f%8.1f  BALL\n", r, g, b,
            x+midx/PSFAC, y+600-midy/PSFAC, rad);
  else
    fprintf(outfp, "%6.2f%8.1f%8.1f%8.1f  BALL\n", gray, 
            x+midx/PSFAC, y+600-midy/PSFAC, rad);
}


void hardcopy_label(float x, float y, char *str)
{
  int sh;
  sh = 2.2*strlen(str);
  fprintf(outfp, "(%s) %7.2f %7.2f LAB\n", str, 
          x+midx/PSFAC-sh,y+600-midy/PSFAC-2);
}

void hardcopy_stick(float gray, float m1[6], float m2[6])
{
  fprintf(outfp,
          "%6.2f [%7.1f%7.1f%7.1f%7.1f%7.1f%7.1f ]\n",
          gray, 
          m1[0], m1[1], m1[2], m1[3], 
          m1[4]+midx/PSFAC, m1[5]+600-midy/PSFAC);
  fprintf(outfp,
          "       [%7.1f%7.1f%7.1f%7.1f%7.1f%7.1f ]  STICK\n",
          m2[0], m2[1], m2[2], m2[3], 
          m2[4]+midx/PSFAC, m2[5]+600-midy/PSFAC);
}


/* ----- hardcopy_close ------ */
void hardcopy_close()
{
  fprintf(outfp,
"showpage\n"
"%%%%EOF\n");
 fclose(outfp);
}

void hardcopy_xdbond(float gray, float m1[6], float m2[6], float alf)
{
  fprintf(outfp," matrix currentmatrix");
  fprintf(outfp," [%7.1f%7.1f%7.1f%7.1f%7.1f%7.1f ] concat\n",
          m1[0], m1[1], m1[2], m1[3], m1[4], m1[5]);
  fprintf(outfp," 0 0 1 %7.2f %7.2f arc setmatrix\n", 
          90+alf, 270-alf);
  fprintf(outfp," matrix currentmatrix");
  fprintf(outfp," [%7.1f%7.1f%7.1f%7.1f%7.1f%7.1f ] concat\n",
          m2[0], m2[1], m2[2], m2[3], m2[4], m2[5]);
  fprintf(outfp," 0 0 1 %7.2f %7.2f arc closepath setmatrix XDB \n", 
          -90-alf, 90+alf);
  fprintf(outfp," gsave %6.2f setgray fill grestore stroke\n",
          gray);
}

void hardcopy_ydbond(float gray, float m1[6], float m2[6], float alf)
{
  fprintf(outfp," matrix currentmatrix");
  fprintf(outfp," [%7.1f%7.1f%7.1f%7.1f%7.1f%7.1f ] concat\n",
          m2[0], m2[1], m2[2], m2[3], m2[4], m2[5]);
  fprintf(outfp," 0 0 1 0 360 arc setmatrix XDB\n");
  fprintf(outfp," gsave %6.2f setgray fill grestore stroke\n",
          gray);
}



