c
c   OLEM compiled for size 1250x960x20 using 16 processors
c
c  : make_param,v 1.2 1998/10/16 21:21:25 dwd Exp $
c
      integer idim,jdim,kdim,idim1,jdim1,kdim1,trcdim,vardim,
     *  nsize,ijdimt2,fltdim,dscdim,varoutdim,kdimt2,jtotal2,
     *  jtotal,nprocs,jtotal1,idim2,jdim2,kdim2,ijk2dim,ij2dim,
     *  idimd,jdimd,kdimd,idimd1,jdimd1,kdimd1,idimd2,jdimd2,
     *  kdimd2,jtotald,jtotald1
      parameter(nprocs=24,jtotal=960,jtotald=1920)
      parameter(idim=1250,jdim=jtotal/nprocs,kdim=20,trcdim=1,
     *     vardim=127,fltdim=4000,dscdim=3,
     *     idim1=idim+1,jdim1=jdim+1,kdim1=kdim+1,nsize=2000,
     *     varoutdim=20,ijdimt2=idim*jdim*2,idim2=idim+2,jdim2=jdim+2,
     *     kdim2=kdim+2,ijk2dim=idim2*jdim2*kdim2,jtotal1=jtotal+1,
     *     jdimd=2*jdim,idimd=2*idim,kdimd=2*kdim,
     *     jdimd1=2*jdim+1,idimd1=2*idim+1,kdimd1=2*kdim+1,
     *     idimd2=idimd+2,jdimd2=jdimd+2,kdimd2=kdimd+2,
     *     ij2dim=idim2*jdim2,kdimt2=2*kdim,jtotal2=jtotal+2,
     *     jtotald1=jtotald+1)
