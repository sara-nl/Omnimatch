void cross(float *volinout, int Vx_max)
     /* procedure exchanges quadrants 
	only for cubic volumes of even dimensions
	FF 01/13/02
     */
{
  int lauf, i, j, k, tmpx, tmpy, tmpz;
  float tempval;

  lauf = 0;
  for (i = 0; i < Vx_max/2; i++)
    {
      for (j = 0 ; j < Vx_max; j++)
	{
	  for (k = 0; k < Vx_max; k++)
	    {
	      if(i < (Vx_max/2) ) tmpx=i+Vx_max/2;
	      else tmpx = i - Vx_max/2;
	      if(j<Vx_max/2) tmpy=j+Vx_max/2;
	      else tmpy=j-Vx_max/2;
	      if(k<Vx_max/2) tmpz=k+Vx_max/2;
	      else tmpz=k-Vx_max/2;
	      //printf ("tmpz: %d ... \n", tmpz);fflush(stdout); 
	      tempval = volinout[tmpx*Vx_max*Vx_max+tmpy*Vx_max+tmpz];
	      volinout[tmpx*Vx_max*Vx_max+tmpy*Vx_max+tmpz] = volinout[lauf];
	      volinout[lauf] = tempval;
	      lauf++;
	    }
	}
    }
}
