c
c   OLEM compiled for size 648x648x400 using 64 processors
c
c  : make_param,v 1.2 1998/10/16 21:21:25 dwd Exp $
c
      integer idim,jdim,kdim,idim1,jdim1,kdim1,trcdim,vardim,
     *  nsize,ijdimt2,fltdim,dscdim,varoutdim,kdimt2,jtotal2,
     *  jtotal,nprocs,jtotal1,idim2,jdim2,kdim2,ijk2dim,ij2dim,
     *  itotal,itotal1,itotal2,nprocx,nprocy
      parameter(nprocx=18,nprocy=18,itotal=648,jtotal=648)
      parameter(idim=itotal/nprocx,jdim=jtotal/nprocy,kdim=400,trcdim=1,
     *     vardim=129,fltdim=4000,dscdim=3,nprocs=nprocx*nprocy,
     *     idim1=idim+1,jdim1=jdim+1,kdim1=kdim+1,nsize=2000,
     *     varoutdim=20,ijdimt2=idim*jdim*2,idim2=idim+2,jdim2=jdim+2,
     *     kdim2=kdim+2,ijk2dim=idim2*jdim2*kdim2,jtotal1=jtotal+1,
     *     ij2dim=idim2*jdim2,kdimt2=2*kdim,jtotal2=jtotal+2,
     *     itotal1=itotal+1, itotal2=itotal+2)
