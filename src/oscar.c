#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <tom.h>
#include <zeit.h>
#include <fftw.h>
#include <rfftw.h>
#include <mpi.h>


MPI_Status status;
MPI_Request request, request1, request2, request3, request4;

int myrank, mysize, tag = 99;

main (int argc, char *argv[])
{
  struct em_file inputdata1;
  struct em_file inputdata2;
  struct em_file inputdata3;
  struct em_file inputdata4;
  struct em_file outputdata;

  fftw_real *Vol_tmpl_sort, *Volume, *e3, *PointCorr, *sqconv;
  fftw_complex *C3, *PointVolume, *PointSq;
  rfftwnd_plan p3, pi3, r3, ri3;
  fftw_real scale;

  struct tm *zeit;
  struct tm start;
  char name[200];
  int Rx_max, Ry_max, Rz_max;
  int Rx_min, Ry_min, Rz_min;
  int Vx_min, Vy_min, Vz_min;
  int Vx_max, Vy_max, Vz_max;
  float Phi, Psi, Theta, winkel_lauf;
  float *Rot_tmpl, *Vol_tmpl;
  int i, j, k, tmpx, tmpy, tmpz,lauf_pe, ksub;
  int ijk;
  int lauf, n;
  float max, eps;
  time_t lt;
  float Ctmp, Ctmpim, Dtmp, Dtmpim;
  int dim_fft;
  int sub[3],range[3],range_sub[3],subc[3],offset[3],dimarray[3];
  int FullVolume_dims[3];
  int nr[3];
  int area[3];
  
/* MPI Variablen */
  int winkel_max, winkel_min;
  int winkel_max_pe, winkel_min_pe;
  int winkel_step_pe;
  int Phi_max, Psi_max, Theta_max;
  int Phi_min, Psi_min, Theta_min;
  int Phi_step, Psi_step, Theta_step;
  int Theta_winkel_start, Psi_winkel_start, Phi_winkel_start;
  int Theta_winkel_nr, Psi_winkel_nr, Phi_winkel_nr;
  int Theta_winkel_end, Psi_winkel_end, Phi_winkel_end;
  int Theta_steps, Psi_steps, Phi_steps;
  float Theta_winkel_rest_nr, Psi_winkel_rest_nr, Phi_winkel_rest_nr;
  int in_max;
  float rms_wedge, tempccf;
  float *Ergebnis, *conv;
  float cycles;
  int cycle;
  
/* MPI Variablen Ende*/

  if (argc < 15)
    {
      printf ("\n\n");
      printf (" 'OSCAR' is an Optimized SCanning AlgoRithm for \n");
      printf (" local correlation.\n");
      printf (" All files in EM-V-int4 format !!!\n\n");
      printf (" Input: Volume to be searched, Template mask for local \n ");
      printf ("   correlation, pointspread function and angular search \n");
      printf ("   range. \n");
      printf (" Output: locally normalized X-Correlation Function Out.ccf.norm, \n");
      printf ("   non-normalized X-Correlation Function Out.ccf, and Out.ang \n");
      printf ("   with the corresponding angles.\n\n");
      printf (" usage: oscar Volume Template Out ...\n");
      printf ("         ... Phi_min Phi_max Phi_step Psi_min Psi_max Psi_step The_min The_max The_step\n");
      printf ("    ... Poinspread-function mask-file dim_of_fft\n\n");
      printf (" with Message Passing Interface (MPI)\n");
      printf (" the total number of angles must be modulo\n");
      printf (" of used processors!\n\n");
      printf (" Linux:   	1.'lamboot' to start MPI\n");
      printf ("		2.'mpirun -np 2 oscar Volume Templ Out 30 180 30 30 180 30 30 180 30 Poinspread-function mask-file 256'\n\n");
      printf (" In this version asymmetric masks can be used ! \n");
      printf (" last revision  ,  11.11.03, Friedrich Foerster");
      printf (" \n\n");
      exit (1);
    }

  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &mysize);
  MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
  /* Dimensionen auslesen */
  // Dimension of fft
  dim_fft = atoi (argv[15]);
  nr[0]=1;
  nr[1]=1;
  nr[2]=1;
  area[0]=dim_fft;
  area[1]=dim_fft;
  area[2]=dim_fft;
  read_em_header(argv[1], &inputdata1); /* Searchvolume */
  read_em (argv[2], &inputdata2); /* Template */
  FullVolume_dims[0]=inputdata1.dims[0];
  FullVolume_dims[1]=inputdata1.dims[1];
  FullVolume_dims[2]=inputdata1.dims[2];
  Rx_min = 1;
  Ry_min = 1;
  Rz_min = 1;
  Rx_max = (inputdata2.dims[0]);
  Ry_max = (inputdata2.dims[1]);
  Rz_max = (inputdata2.dims[2]);
  Vx_min = 1;
  Vy_min = 1;
  Vz_min = 1;
  Vx_max = dim_fft;
  Vy_max = dim_fft;
  Vz_max = dim_fft;
  p3 = rfftw3d_create_plan (Vx_max, Vy_max, Vz_max, FFTW_REAL_TO_COMPLEX,
			    FFTW_MEASURE | FFTW_IN_PLACE);	/*FFTW_ESTIMATE FFTW_MEASURE */
  pi3 = rfftw3d_create_plan (Vx_max, Vy_max, Vz_max, FFTW_COMPLEX_TO_REAL,
			     FFTW_MEASURE | FFTW_IN_PLACE);
  r3 = rfftw3d_create_plan (Rx_max, Rx_max, Rx_max, FFTW_REAL_TO_COMPLEX,
  		    FFTW_MEASURE | FFTW_IN_PLACE);	/*FFTW_ESTIMATE FFTW_MEASURE */
  ri3 = rfftw3d_create_plan (Rx_max, Rx_max, Rx_max, FFTW_COMPLEX_TO_REAL,
  		     FFTW_MEASURE | FFTW_IN_PLACE);
  if (myrank == 0)
    {
      printf("Plans for FFTW created \n");fflush(stdout);
    }
  Volume = (fftw_real *) calloc (Vx_max * Vx_max * 2 * (Vx_max / 2 + 1),sizeof (fftw_real) );
  Rot_tmpl = (float *) malloc (sizeof (float) * Rx_max * Ry_max * Rz_max);
  Vol_tmpl = (float *) malloc (sizeof (float) * Vx_max * Vy_max * Vz_max);
  conv = (float *) malloc (sizeof (float) * Vx_max * Vy_max * Vz_max);
  sqconv = (fftw_real *) calloc(Vz_max * Vy_max * 2 * (Vx_max / 2 + 1), sizeof (fftw_real));
  if (!
      (inputdata1.floatdata =
       (float *) malloc (sizeof (float) * Vx_max * Vy_max * Vz_max)))
    {
      printf ("Memory allocation  failure in inputdata1.floatdata!!!");
      fflush (stdout);
      exit (1);
    }
  if (!
      (outputdata.floatdata =
       (float *) malloc (sizeof (float) * Vx_max * Vy_max * Vz_max)))
    {
      printf ("Memory allocation  failure in outputdata.floatdata!!!");
      fflush (stdout);
      exit (1);
    }
  
  if (!
      (Vol_tmpl_sort = (fftw_real *) calloc (Vz_max*Vy_max*2*(Vx_max / 2 + 1), sizeof (fftw_real) )))
    {
      printf ("Memory allocation  failure in Volume_tmpl_sort!!!");
      printf ("Nx = %i, Ny = %i, Nz = %i, bytes = %i \n",2 *(Vx_max / 2 + 1),Vy_max, Vz_max, sizeof (fftw_real));
      fflush (stdout);
      exit (1);
    }
  Ergebnis = (float *) calloc (Vz_max * Vy_max * Vx_max, sizeof (float));
   /* Winkelraum */
  Phi_min = atof (argv[4]);
  Phi_max = atof (argv[5]);
  Phi_step = atof (argv[6]);
  Psi_min = atof (argv[7]);
  Psi_max = atof (argv[8]);
  Psi_step = atof (argv[9]);
  Theta_min = atof (argv[10]);
  Theta_max = atof (argv[11]);
  Theta_step = atof (argv[12]);
  /* Pointspread Function*/
  read_em (argv[13], &inputdata3);
  /* mask function */
  read_em (argv[14], &inputdata4);
  Phi_steps = (Phi_max - Phi_min) / Phi_step + 1;
  Psi_steps = (Psi_max - Psi_min) / Psi_step + 1;
  Theta_steps = (Theta_max - Theta_min) / Theta_step + 1;
  winkel_max = Phi_steps * Psi_steps * Theta_steps;
  winkel_min = 0;
  range[0]=dim_fft-1;
  range[1]=dim_fft-1;
  range[2]=dim_fft-1;
  range_sub[0]=range[0]-Rx_max;
  range_sub[1]=range[1]-Rx_max;
  range_sub[2]=range[2]-Rx_max;
  sub[0]=1;
  sub[1]=1;
  sub[2]=1;
  cycles=(int)(FullVolume_dims[2]/(dim_fft-Rx_max)+0.5);
  cycles=(int)(FullVolume_dims[1]/(dim_fft-Rx_max)+0.5)*cycles;
  cycles=(int)(FullVolume_dims[0]/(dim_fft-Rx_max)+0.5)*cycles;
  cycle=0;
  if (myrank == 0)
    {
      printf ("\n oscar starts to run ... ");tack (&start);fflush (stdout);
      /* prepare Output */
      strcpy (name, argv[3]);
      strcat (name, ".ccf");
      printf ("\nCreate outputfile: %s ... \n", name);fflush(stdout);
      create_em (name, FullVolume_dims);
      strcpy (name, argv[3]);
      strcat (name, ".ang");
      printf ("Create outputfile: %s ... \n", name);fflush(stdout);
      create_em (name, FullVolume_dims);
      strcpy (name, argv[3]);
      strcat (name, ".ccf.norm");
      printf ("Create outputfile: %s ... \n", name);fflush(stdout);
      create_em (name, FullVolume_dims);
    }
  for (sub[2]=1; sub[2] < FullVolume_dims[2]-Rz_max;sub[2]=sub[2]+dim_fft-Rz_max)
    {		  
    if (myrank == 0)
	{
	tack (&start);
	printf ("%f%%..", (float) (cycle / cycles * 100));
	fflush (stdout);
	}

      for (sub[1]=1; sub[1] < FullVolume_dims[1]-Ry_max;sub[1]=sub[1]+dim_fft-Ry_max)
        {
	  for (sub[0]=1; sub[0] < FullVolume_dims[0]-Rx_max;sub[0]=sub[0]+dim_fft-Rx_max)
	    {
	      cycle=cycle+1;
	      subc[0]=sub[0];
	      subc[1]=sub[1];
	      subc[2]=sub[2]; 
	      if (sub[2] + range[2] > FullVolume_dims[2]) subc[2]=FullVolume_dims[2]-range[2];  /* we are at the corner ?!*/
	      if (sub[1] + range[1] > FullVolume_dims[1]) subc[1]=FullVolume_dims[1]-range[1];  /* we are at the corner ?!*/
	      if (sub[0] + range[0] > FullVolume_dims[0]) subc[0]=FullVolume_dims[0]-range[0];  /* we are at the corner ?!*/
	      read_em_subregion (argv[1], &inputdata1,subc,range);
	      read_em_subregion (argv[1], &outputdata,subc,range);
	      /* Umsortieren der Daten */
	      lauf = 0;
	      for (k = 0; k < Vz_max; k++)
		{
		  for (j = 0; j < Vy_max; j++)
		    {
		      for (i = 0; i < Vx_max; i++)
			{
			  /* square - needed for normalization */
			  sqconv[i + 2 * (Vx_max / 2 + 1) * (j + Vy_max * k)] = inputdata1.floatdata[lauf]*inputdata1.floatdata[lauf];
			  Volume[i + 2 * (Vx_max / 2 + 1) * (j + Vy_max * k)] = inputdata1.floatdata[lauf];
			  inputdata1.floatdata[lauf] = -1.0; /* kleine Zahl wg Max-Op , hier kommen die CCFs rein*/
			  outputdata.floatdata[lauf] = -1.0; /* hier kommen die Winkel rein*/
			  lauf++;
			}
		    }
		}
	      rfftwnd_one_real_to_complex (p3, &Volume[0], NULL); /* einmalige fft von Suchvolumen */
	      rfftwnd_one_real_to_complex (p3, &sqconv[0], NULL); /* FFT of square*/
	      winkel_step_pe = (int) winkel_max / mysize;
	      winkel_min_pe = myrank * winkel_step_pe;
	      winkel_max_pe = winkel_min_pe + winkel_step_pe;
	      Theta_winkel_nr = (int) winkel_min_pe / (Psi_steps * Phi_steps);
	      Theta_winkel_rest_nr = winkel_min_pe - Theta_winkel_nr * (Psi_steps * Phi_steps);
	      Psi_winkel_nr = (int) Theta_winkel_rest_nr / (Phi_steps);
	      Psi_winkel_rest_nr = Theta_winkel_rest_nr - Psi_winkel_nr * (Phi_steps);
	      Phi_winkel_nr = (int) Psi_winkel_rest_nr;
	      Theta = Theta_winkel_nr * Theta_step + Theta_min;
	      Phi = Phi_winkel_nr * Phi_step + Phi_min - Phi_step;
	      Psi = Psi_winkel_nr * Psi_step + Psi_min;
	      eps = 0.001;
	      n = 0;
	      //Friedrich -> Zaehlung der voxels 
	      n = countvoxel(inputdata4.dims[0], inputdata4.floatdata, eps);
	      eps = 0.001;
	      for (winkel_lauf = winkel_min_pe; winkel_lauf < winkel_max_pe;winkel_lauf++)
		{
		  if (Phi < Phi_max)
		    Phi = Phi + Phi_step;
		  else
		    {
		      Phi = Phi_min;
		      Psi = Psi + Psi_step;
		    }
		  if (Psi > Psi_max)
		    {
		      Psi = Psi_min;
		      Theta = Theta + Theta_step;
		    }
		  tom_rotate3d (&Rot_tmpl[0], &inputdata2.floatdata[0], Phi, Psi, Theta, Rx_max, Ry_max, Rz_max);
		  /*calculate Ref variance */
		  rms_wedge = energizer (Rx_min, Rx_max, n, &Rot_tmpl[0], &inputdata3.floatdata[0], &inputdata4.floatdata[0], r3, ri3); 
		  pastes (&Rot_tmpl[0], &Vol_tmpl[0], 1, 1, 1, Rx_max, Ry_max, Rz_max, Vx_max);
		  scale = 1.0 / ((double)Vx_max * (double)Vy_max * (double)Vz_max * ((double) rms_wedge) );
		  //printf("hippo1: scale = %.10f \n",scale);
		  sort4fftw(&Vol_tmpl_sort[0],&Vol_tmpl[0],Vx_max, Vy_max, Vz_max);
		  rfftwnd_one_real_to_complex (p3, &Vol_tmpl_sort[0], NULL);
		  PointVolume = (fftw_complex *) & Volume[0];
		  C3 = (fftw_complex *) & Vol_tmpl_sort[0];
		  /* Correlation */
		  correl(&PointVolume[0], &C3[0], Vx_max, Vy_max, Vz_max, scale);
		  /* back to real space */
		  rfftwnd_one_complex_to_real (pi3, &C3[0], NULL);
		  PointCorr = (fftw_real *) & C3[0];
		  /* Umsortieren der Daten */
		  sortback4fftw( &PointCorr[0], &Ergebnis[0], Vx_max, Vy_max, Vz_max);
		  // crossen 
		  cross(&Ergebnis[0], Vx_max);
		  /* 3rd: divide */
		  lauf = 0;
		  for (k = 0 ; k < Vz_max  ; k++)
		    {
		      for (j = 0; j < Vy_max; j++)
			{
			  for (i = 0; i < Vx_max; i++)
			    {
			      if (inputdata1.floatdata[lauf] < Ergebnis[lauf] )
				{
				  inputdata1.floatdata[lauf] = Ergebnis[lauf];
				  outputdata.floatdata[lauf] = (int) winkel_lauf;
				}
			      lauf++;
			    }
			}
		    }
		}				/* Ende winkel_lauf */
	      //FF
	      MPI_Barrier (MPI_COMM_WORLD);
	      /* Ergebnisse einsammeln (myrank 0)*/
	      if (myrank == 0)
		{
		  for (lauf_pe = 1; lauf_pe < mysize; lauf_pe++)
		    {
		      MPI_Recv (&Ergebnis[0], Vx_max * Vy_max * Vz_max, MPI_FLOAT, lauf_pe,
				99, MPI_COMM_WORLD, &status);
		      MPI_Recv (&conv[0], Vx_max * Vy_max * Vz_max, MPI_FLOAT,
				lauf_pe, 98, MPI_COMM_WORLD, &status);
		      /* use conv as temporary memory for angles  */
		      for (lauf = 0; lauf < Vx_max * Vy_max * Vz_max; lauf++)
			{
			  if (inputdata1.floatdata[lauf] < Ergebnis[lauf])
			    {
			      inputdata1.floatdata[lauf] = Ergebnis[lauf];
			      outputdata.floatdata[lauf] = conv[lauf];
			    }
			}
		    }
		  /*Ergebnisse eingesammelt */
		  
		}
	      // myrank > 0: Ergebnisse senden
	      else
		{
		  MPI_Send (inputdata1.floatdata, Vx_max * Vy_max * Vz_max, MPI_FLOAT, 0,
			    99, MPI_COMM_WORLD);
		  MPI_Send (outputdata.floatdata, Vx_max * Vy_max * Vz_max, MPI_FLOAT, 0,
			    98, MPI_COMM_WORLD);
		}
	      MPI_Barrier (MPI_COMM_WORLD);
	      // nicht normalisiertes Volumen und Winkel rausschreiben
	      subc[0]=subc[0]+Rx_max/2;
	      subc[1]=subc[1]+Rx_max/2;
	      subc[2]=subc[2]+Rx_max/2;
	      if (myrank==0)
		{
		  offset[0]=Rx_max/2;
		  offset[1]=Rx_max/2;
		  offset[2]=Rx_max/2;
		  dimarray[0]=dim_fft;
		  dimarray[1]=dim_fft;
		  dimarray[2]=dim_fft;
		  strcpy (name, argv[3]);
		  strcat (name, ".ccf");
		  write_em_subsubregion (name, &inputdata1,subc,range_sub,offset,dimarray); 
		  strcpy (name, argv[3]);
		  strcat (name, ".ang");
		  write_em_subsubregion (name, &outputdata,subc,range_sub,offset,dimarray);
		  /* ------------------- normalization - here only PE 0 ---------- */
		  pastes (&inputdata4.floatdata[0], &Vol_tmpl[0], 1, 1, 1, Rx_max, Ry_max, Rz_max, Vx_max); /* paste mask into zero volume*/
		  /* 1st local mean */
		  sort4fftw(&Vol_tmpl_sort[0], &Vol_tmpl[0], Vx_max, Vy_max, Vz_max);
		  rfftwnd_one_real_to_complex (p3, &Vol_tmpl_sort[0], NULL);
		  C3 = (fftw_complex *) & Vol_tmpl_sort[0];
		  /* Convolution of volume and mask */
		  scale = 1.0 / ((double)Vx_max * (double)Vy_max * (double)Vz_max );
		  convolve( &PointVolume[0], &C3[0], Vx_max, Vy_max, Vz_max, scale);
		  rfftwnd_one_complex_to_real (pi3, &C3[0], NULL);
		  PointCorr = (fftw_real *) & C3[0];
		  /* Umsortieren der Daten */
		  sortback4fftw( &PointCorr[0], &conv[0], Vx_max, Vy_max, Vz_max);
		  /* 2nd : convolution of square and resorting*/
		  pastes (&inputdata4.floatdata[0], &Vol_tmpl[0], 1, 1, 1, Rx_max, Ry_max, Rz_max, Vx_max); /* paste mask into zero volume*/
		  sort4fftw( &Vol_tmpl_sort[0], &Vol_tmpl[0], Vx_max, Vy_max, Vz_max);
		  rfftwnd_one_real_to_complex (p3, &Vol_tmpl_sort[0], NULL);
		  C3 = (fftw_complex *) & Vol_tmpl_sort[0];
		  PointSq = (fftw_complex *) & sqconv[0];// set pointer to FFT of square
		  convolve( &PointSq[0], &C3[0], Vx_max, Vy_max, Vz_max, scale);
		  rfftwnd_one_complex_to_real (pi3, &C3[0], NULL);
		  PointCorr = (fftw_real *) &C3[0];
		  //FF
		  lauf = 0;
		  for (k = 0; k < Vz_max; k++)
		    {
		      for (j = 0; j < Vy_max; j++)
			{
			  for (i = 0; i < Vx_max; i++)
			    {
			      conv[lauf] = sqrt(PointCorr[i + 2 * (Vx_max / 2 + 1) * (j + Vy_max * k)] - conv[lauf]*conv[lauf]/((float) n) ) ;/*local variance*/
			      lauf++;
			    }
			}
		    }
		  cross(&conv[0], Vx_max);
		  /* perform division */
		  for (lauf = 0; k < Vz_max*Vy_max*Vz_max; lauf++)
		    {
		      if (conv[lauf] > eps)
			{
			  inputdata1[lauf].floatdata = inputdata1[lauf].floatdata/conv[lauf];
			}
		      else
			{
			  inputdata1[lauf].floatdata = 0;
			}
		    }
		  strcpy (name, argv[3]);
		  strcat (name, ".ccf.norm");
		  write_em_subsubregion (name, &inputdata1,subc,range_sub,offset,dimarray);
		}
	      MPI_Barrier (MPI_COMM_WORLD);
	    }
	} /* these are the new brackets from the subregion_read , SN */
    }
  free(Ergebnis);
  free(inputdata1.floatdata);  
  free(inputdata2.floatdata);
  free(inputdata3.floatdata);
  free(inputdata4.floatdata);
  rfftwnd_destroy_plan(p3);
  rfftwnd_destroy_plan(pi3);
  rfftwnd_destroy_plan(r3);
  rfftwnd_destroy_plan(ri3);
  free(Volume);
  free(sqconv);
  free(conv);
  free(Rot_tmpl);
  free(Vol_tmpl_sort);
  free(outputdata.floatdata);
  if (myrank==0)
    {
      printf ("oscar finished. ");
      tack (&start); fflush(stdout);
    }
  MPI_Finalize();

  /* end main */
}

/* Variablen:

p3 und pi3: 	 Plan fuer die FFT bzw. IFFT (Volume)
r3 und ri3: 	 Plan fuer die FFT bzw. IFFT (Template Vol)
Volume:		 Suchvolumen (sortiert fuer FFT)
PointVolume:	 komplexer Pointer auf Suchvolumen fft(Volume)
sqconv:          Quadrat des Volumens (sortiert fuer FFT)
PointSq:         komplexer Pointer auf fft(sqconv)
conv:            Volumen fuer convolution Vol - mask / Temp fuer Winkel
Rot_tmpl:	 Rotiertes Tmpl fuer Suche
Vol_tmpl:	 Volumen mit rot. Template
Vol_tmpl_sort:   Volumen mit rot. Template sortiert fuer FFT
Vx_max...:	 Dimensionen des Suchvolumens
Rx_max...:	 Dimensionen des Templatevolumens
PointCorr:	 realer Pointer auf Correlationsergebnis
C3:		 komplexer Pointer auf fft(Vol_tmpl_sort)
Ergebnis:        max CCF, spaeter local mean
Ergebnis_winkel: Winkel Index 
i,j,k,lauf:	 Laufvariablen
*/













