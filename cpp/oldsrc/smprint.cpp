void timeseries::smprint(std::string nf,int c,bool usegmt) const
{
  int i,s=size();

  double* ave = new double[6];
  average(ave,0);

  switch (usegmt) {
    /* DO NOT USE GMT FORMAT; DATES ARE PRINTED AS DECIMAL YEARS */
    case (0):
      if (nf=="") {/* print to stream */
        switch (c) {
          case (0):
            for (i=0;i<s;i++) {
              printf("\n%10.5f %8.5f %8.6f %1i %1i",time[i],x.vals[i],x.sigmas[i],x.flags[i],soltype[i]);
              printf(" %8.5f %8.6f %1i",y.vals[i],y.sigmas[i],y.flags[i]);
              printf(" %8.5f %8.6f %1i",z.vals[i],z.sigmas[i],z.flags[i]);
            }
            break;
          case (1):
            for (i=0;i<s;i++) printf("\n%10.5f %8.5f %8.6f %1i %1i",time[i],x.vals[i],x.sigmas[i],x.flags[i],soltype[i]);
            break;
          case (2):
            for (i=0;i<s;i++) printf("\n%10.5f %8.5f %8.6f %1i %1i",time[i],y.vals[i],y.sigmas[i],y.flags[i],soltype[i]);
            break;
          case (3):
            for (i=0;i<s;i++) printf("\n%10.5f %8.5f %8.6f %1i %1i",time[i],z.vals[i],z.sigmas[i],z.flags[i],soltype[i]);
            break;
        }
      } else {/* print to file */
        FILE* pFile;
        pFile = fopen(nf.c_str(),"w");
	        if (!pFile) {
		        std::cerr<<"\n*** Could not write time-series file!";
	        } else {
                switch (c) {
                  case (0):
                    for (i=0;i<s;i++) {
                      printf("\n%10.5f %8.5f %8.6f %1i %1i",time[i],x.vals[i],x.sigmas[i],x.flags[i],soltype[i]);
                      printf(" %8.5f %8.6f %1i",y.vals[i],y.sigmas[i],y.flags[i]);
                      printf(" %8.5f %8.6f %1i",z.vals[i],z.sigmas[i],z.flags[i]);
                    }
                    break;
                  case (1):
                    for (i=0;i<s;i++) fprintf(pFile,"\n%10.5f %8.5f %8.6f %1i %1i",time[i],x.vals[i],x.sigmas[i],x.flags[i],soltype[i]);
                    break;
                  case (2):
                    for (i=0;i<s;i++) fprintf(pFile,"\n%10.5f %8.5f %8.6f %1i %1i",time[i],y.vals[i],y.sigmas[i],y.flags[i],soltype[i]);
                    break;
                  case (3):
                    for (i=0;i<s;i++) fprintf(pFile,"\n%10.5f %8.5f %8.6f %1i %1i",time[i],z.vals[i],z.sigmas[i],z.flags[i],soltype[i]);
                    break;
              }
	        	fclose(pFile);
		        printf("\n$$ t-series file written: %20s",nf.c_str());
	        }
      }
      break;
    /* USE GMT FORMAT; DATES ARE PRINTED AS DD-MM-YYYY */
    case (1):
      int iy,im,id;
      std::string sm;
      if (nf=="") {/* print to stream */
        switch (c) {
          case (0):
            for (i=0;i<s;i++) {
              fy2cal(time[i],iy,im,id);
              sm = m2m(im);
              printf("\n%02i-%3s-%4i %8.5f %8.6f %1i %1i",id,sm.c_str(),iy,x.vals[i],x.sigmas[i],x.flags[i],soltype[i]);
              printf(" %8.5f %8.6f %1i",y.vals[i],y.sigmas[i],y.flags[i]);
              printf(" %8.5f %8.6f %1i",z.vals[i],z.sigmas[i],z.flags[i]);
            }
            break;
          case (1):
            for (i=0;i<s;i++) {
              fy2cal(time[i],iy,im,id);
              sm = m2m(im);
              printf("\n%02i-%3s-%4i %8.5f %8.6f %1i %1i",id,sm.c_str(),iy,x.vals[i],x.sigmas[i],x.flags[i],soltype[i]);
            }
            break;
          case (2):
            for (i=0;i<s;i++) {
              fy2cal(time[i],iy,im,id);
              sm = m2m(im);
              printf("\n%02i-%3s-%4i %8.5f %8.6f %1i %1i",id,sm.c_str(),iy,y.vals[i],y.sigmas[i],y.flags[i],soltype[i]);
            }
            break;
          case (3):
            for (i=0;i<s;i++) {
              fy2cal(time[i],iy,im,id);
              sm = m2m(im);
              printf("\n%02i-%3s-%4i %8.5f %8.6f %1i %1i",id,sm.c_str(),iy,z.vals[i],z.sigmas[i],z.flags[i],soltype[i]);
            }
            break;
        }
      } else {/* print to file */
        FILE* pFile;
        pFile = fopen(nf.c_str(),"w");
	if (!pFile) {
	        std::cerr<<"\n*** Could not write time-series file!";
	} else {
                switch (c) {
          case (0):
            for (i=0;i<s;i++) {
              fy2cal(time[i],iy,im,id);
              sm = m2m(im);
              fprintf(pFile,"\n%02i-%3s-%4i %8.5f %8.6f %1i",id,sm.c_str(),iy,x.vals[i],x.sigmas[i],x.flags[i]);
              fprintf(pFile," %8.5f %8.6f %1i",y.vals[i],y.sigmas[i],y.flags[i]);
              fprintf(pFile," %8.5f %8.6f %1i",z.vals[i],z.sigmas[i],z.flags[i]);
            }
            break;
          case (1):
            for (i=0;i<s;i++) {
              fy2cal(time[i],iy,im,id);
              sm = m2m(im);
              double limit = 3.0 * ave[1];
              if (fabs(ave[0]-x.vals[i])<limit) {
                fprintf(pFile,"\n%02i-%3s-%4i %8.5f %8.6f %1i %1i",id,sm.c_str(),iy,x.vals[i],x.sigmas[i],x.flags[i],soltype[i]);
              }
            }
            break;
          case (2):
            for (i=0;i<s;i++) {
              fy2cal(time[i],iy,im,id);
              sm = m2m(im);
              double limit = 3.0 * ave[3];
              if (fabs(ave[2]-y.vals[i])<limit) {
                fprintf(pFile,"\n%02i-%3s-%4i %8.5f %8.6f %1i %1i",id,sm.c_str(),iy,y.vals[i],y.sigmas[i],y.flags[i],soltype[i]);
              }
            }
            break;
          case (3):
            for (i=0;i<s;i++) {
              fy2cal(time[i],iy,im,id);
              sm = m2m(im);
              double limit = 3.0 * ave[5];
              if (fabs(ave[4]-z.vals[i])<limit) {
              fprintf(pFile,"\n%02i-%3s-%4i %8.5f %8.6f %1i %1i",id,sm.c_str(),iy,z.vals[i],z.sigmas[i],z.flags[i],soltype[i]);
              }
            }
            break;
              }
	        	fclose(pFile);
		        printf("\n$$ t-series file written: %20s",nf.c_str());
	        }
      }
      break;
  }

  return;
}
