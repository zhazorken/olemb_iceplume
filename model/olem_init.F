      subroutine olem_init(ier)
c
#include "param.inc"
#include "olemp.inc"
#include "olemtke.inc"
#include "communicate.inc"
c
      character value*80
      integer unts(3),ier,i,j,k,len,lnblk
      real xtemp
      logical flgs(3),mkstst
c
c initialize sgs component
c
      ier=0
      sgscom1='$Id: olem_init.F,v 1.3 1998/08/04 18:12:50 dwd Exp $'
      sgscom2='TKE based diffusion coefficients'
      sgscom3='Laplacian diffusion, cyclic B.C. version'
c
c parse sgs flags   sgsfstr
c
c values are
c    flag for using sub-grid model, usesgs
c    flag for using tke, usedeardorff
c         (if usesgs .and. .not.usedeardorff,  then metais)
c    flag for turning off advection, advtke
c
      call parse(sgsfstr,value)
      call mks(value,xtemp,flgs,unts)
      if(mkstst(flgs,unts,.false.,.false.,.true.,0,0,1))goto 2000
      sgsflgs(1)=int(xtemp)
      usesgs=sgsflgs(1).eq.1
c
      call parse(sgsfstr,value)
      call mks(value,xtemp,flgs,unts)
      if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 2000
      sgsflgs(2)=int(xtemp)
      usedeardorff=sgsflgs(2).eq.1
c
      call parse(sgsfstr,value)
      call mks(value,xtemp,flgs,unts)
      if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 2000
      sgsflgs(3)=int(xtemp)
      advtke=sgsflgs(3).eq.1
c
c parse sgs values  sgsvstr
c 
c values are 
c    momentum background, kmb
c    heat background, khb
c    error diffusion, nu0
c 
      call parse(sgsvstr,value)
      call mks(value,kmb,flgs,unts)
      if(mkstst(flgs,unts,.true.,.false.,.true.,2,0,-1))goto 2000
      sgsvals(1)=kmb
c
      call parse(sgsvstr,value)
      call mks(value,khb,flgs,unts)
      if(mkstst(flgs,unts,.true.,.false.,.true.,2,0,-1))goto 2000
      sgsvals(2)=khb
c
      call parse(sgsvstr,value)
      call mks(value,nu0,flgs,unts)
      if(mkstst(flgs,unts,.true.,.false.,.true.,2,0,-1))goto 2000
      sgsvals(3)=nu0
c
c interpret sgs values
c
      if(iope) then
        len=lnblk(sgscom1,120)
        write(lunit,1) sgscom1(:len)
 1      format('SGS Comments:'/1x,a)
        len=lnblk(sgscom2,120)
        write(lunit,5)sgscom2(:len)
 5      format(1x,a)
        len=lnblk(sgscom3,120)
        write(lunit,5)sgscom3(:len)
c
        write(lunit,4)
 4      format(/'SGS flags:')
        if(usesgs)then
           write(lunit,3)1,'Use Sub-grid parameterization'
        else
           write(lunit,3)1,'Use kmb and khb for constant diffusion'
        endif
        if(usedeardorff)then
           write(lunit,3)2,'Use Deardorff parameterization'
 3         format(4x,i2,'. ',a)
        else
           write(lunit,3)2,'Use Metais paramterization'
        endif
        if(advtke) then
           write(lunit,3)3,'Advect TKE'
        else
           write(lunit,3)3,'Do not advect TKE'
        endif
c
        write(lunit,2) kmb,khb,nu0
 2      format(/'SGS values:'/
     *5x,'1.    background momentum diff, kmb = ',1pg12.3,' (m2/sec)'/
     *5x,'2. background heat/scalar diff, khb = ',g12.3,' (m2/sec)'/
     *5x,'3.                maximum diff, nu0 = ',g12.3,' (m2/sec)'/)
      endif
c
      if(.not.isrestart)then
c
c compute initial diffusion coefficients
c
         do 100 k=0,kdim1
            do 101 j=0,jdim1
               do 102 i=0,idim1
                  kmom(i,j,k) = kmb
                  kheat(i,j,k) = khb
 102           continue
 101        continue
 100     continue
      endif
c
      return
 2000 write(*,999)value(:30)
 999  format('Units error on SGS card for entry: ',a)
      ier=1
      return
      end
