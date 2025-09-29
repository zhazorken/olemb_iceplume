      subroutine flt_read(ffile, ier)
c
c  $Id: flt_read.F,v 1.6 1998/10/16 21:21:21 dwd Exp $
c
c read the float configuration file
c
#include "param.inc"
#include "olemp.inc"
#include "float.inc"
#include "moddef.inc"
c
	  character ffile*80
      character str*80,value*80,pattern*80,name*80
      character curline*80
      logical flgs(3),mkstst
      integer unts(3),line,ier,num,strtflt,stopflt,type,i,j,k
      real xtemp,xsize,ysize,zsize,xloc,yloc,zloc
      real delta,deltax,deltay,deltaz
c
      open(funit,file=ffile,status='old')
      nflt=0
      ndesc=0
      ier=0
      fltcomment=' '
      line=1
c
 20   read(funit,30,end=100)str
 30   format(a)
      if(str(1:1).eq.'#')goto 20
      curline=str
c
      if(line.eq.1)then
c
c first line read, comment
c comment = character string describing floats
c
         fltcomment=str
         line=2
      else if(line.eq.2)then
c
c second line read, fltdt
c fltdt = time steps between output
c
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.true.,0,0,1))goto 1000
c
         fltdt=int(xtemp/delt+0.5)
         line=3
c
      else
c
c additional pairs of lines
c
c first line 
c pattern, xsize, ysize, zsize, no. floats
c pattern = pattern of floats (point, xline, yline, zline
c                              xyplane, xzplane, yzplane, cube)
c xsize   = size of pattern in meters
c ysize
c zsize
c no. floats= number of floats in EACH dimension
c
         call parse(str,pattern)
c
         call parse(str,value)
         call mks(value,xsize,flgs,unts)
         if(mkstst(flgs,unts,.true.,.false.,.false.,1,0,0))goto 1000
c
         call parse(str,value)
         call mks(value,ysize,flgs,unts)
         if(mkstst(flgs,unts,.true.,.false.,.false.,1,0,0))goto 1000
c
         call parse(str,value)
         call mks(value,zsize,flgs,unts)
         if(mkstst(flgs,unts,.true.,.false.,.false.,1,0,0))goto 1000
c
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.false.,0,0,0))goto 1000
         num=int(xtemp)
c
         read(funit,30,end=1000)str
         curline=str
c
c second line
c xloc, yloc, zloc, start, end, type, name
c xloc = centroid position in meters
c yloc
c zloc
c start = deployment time of floats
c end   = retrieval time of floats
c type  = type of float (neutral, drogue)
c name  = character string designating float name
c
         call parse(str,value)
         call mks(value,xloc,flgs,unts)
         if(mkstst(flgs,unts,.true.,.false.,.false.,1,0,0))goto 1000
c
         call parse(str,value)
         call mks(value,yloc,flgs,unts)
         if(mkstst(flgs,unts,.true.,.false.,.false.,1,0,0))goto 1000
c
         call parse(str,value)
         call mks(value,zloc,flgs,unts)
         if(mkstst(flgs,unts,.true.,.false.,.false.,1,0,0))goto 1000
c
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.true.,0,0,1))goto 1000
         strtflt=int(xtemp/delt+0.5)
c
         call parse(str,value)
         call mks(value,xtemp,flgs,unts)
         if(mkstst(flgs,unts,.false.,.false.,.true.,0,0,1))goto 1000
         stopflt=int(xtemp/delt+0.5)
c
         call parse(str,value)
         if(value.eq."neutral")then
            type = MDFLTXYZ
         else
            type = MDFLTXY
         endif
c
         call parse(str,name)
         ndesc=ndesc+1
         if(name.ne.' ') then
           fltroot(ndesc)=name
         else
           fltroot(ndesc)="flt1"
         endif
c
c build float structures
c
         fltstrt(ndesc)=strtflt
         fltstop(ndesc)=stopflt
         flttype(ndesc)=type
         fltname(ndesc)=name
         if(pattern.eq."point")then
c
c create point
c
            nflt=nflt+1
            if(nflt.gt.fltdim) goto 2000
            fltx(ndesc,nflt)   =xloc
            flty(ndesc,nflt)   =yloc
            fltz(ndesc,nflt)   =zloc
            fltdescsize(ndesc)=num
         else if(pattern.eq.'xline')then
c
c create xline pattern
c
            delta=xsize/float(num-1)
            xloc = xloc - xsize*0.5
            do 200 i=1,num
               nflt=nflt+1
               if(nflt.gt.fltdim) goto 2000
               fltx(ndesc,nflt)=xloc+delta*(i-1)
               flty(ndesc,nflt)=yloc
               fltz(ndesc,nflt)=zloc
c
 200        continue
            fltdescsize(ndesc)=num
         else if(pattern.eq.'yline')then
c
c create yline pattern
c
            delta=ysize/float(num-1)
            yloc = yloc - ysize*0.5
            do 210 i=1,num
               nflt=nflt+1
               if(nflt.gt.fltdim) goto 2000
               fltx(ndesc,nflt)=xloc
               flty(ndesc,nflt)=yloc+delta*(i-1)
               fltz(ndesc,nflt)=zloc
c
 210        continue
            fltdescsize(ndesc)=num
         else if(pattern.eq.'zline')then
c
c create zline pattern
c
            delta=zsize/float(num-1)
            zloc = zloc - zsize*0.5
            do 220 i=1,num
               nflt=nflt+1
               if(nflt.gt.fltdim) goto 2000
               fltx(ndesc,nflt)=xloc
               flty(ndesc,nflt)=yloc
               fltz(ndesc,nflt)=zloc+delta*(i-1)
c
 220        continue
            fltdescsize(ndesc)=num
         else if(pattern.eq.'xyplane')then
c
c create xyplane pattern
c
            deltax=xsize/float(num-1)
            deltay=ysize/float(num-1)
            xloc = xloc - xsize*0.5
            yloc = yloc - ysize*0.5
            do 300 j=1,num
               do 301 i=1,num
                  nflt=nflt+1
                  if(nflt.gt.fltdim) goto 2000
                  fltx(ndesc,nflt)=xloc+deltax*(i-1)
                  flty(ndesc,nflt)=yloc+deltay*(j-1)
                  fltz(ndesc,nflt)=zloc
c
 301           continue
 300        continue
            fltdescsize(ndesc)=num*num
         else if(pattern.eq.'xzplane')then
c
c create xzplane pattern
c
            deltax=xsize/float(num-1)
            deltaz=zsize/float(num-1)
            xloc = xloc - xsize*0.5
            zloc = zloc - zsize*0.5
            do 310 k=1,num
               do 311 i=1,num
                  nflt=nflt+1
                  if(nflt.gt.fltdim) goto 2000
                  fltx(ndesc,nflt)=xloc+deltax*(i-1)
                  flty(ndesc,nflt)=yloc
                  fltz(ndesc,nflt)=zloc+deltaz*(k-1)
c
 311           continue
 310        continue
            fltdescsize(ndesc)=num*num
         else if(pattern.eq.'yzplane')then
c
c create yzplane pattern
c
            deltay=ysize/float(num-1)
            deltaz=zsize/float(num-1)
            yloc = yloc - ysize*0.5
            zloc = zloc - zsize*0.5
            do 320 k=1,num
               do 321 j=1,num
                  nflt=nflt+1
                  if(nflt.gt.fltdim) goto 2000
                  fltx(ndesc,nflt)=xloc
                  flty(ndesc,nflt)=yloc+deltay*(j-1)
                  fltz(ndesc,nflt)=zloc+deltaz*(k-1)
c
 321           continue
 320        continue
            fltdescsize(ndesc)=num*num
         else if(pattern.eq.'cube')then
c
c create cube pattern
c
            deltax=xsize/float(num-1)
            deltay=ysize/float(num-1)
            deltaz=zsize/float(num-1)
            xloc = xloc - xsize*0.5
            yloc = yloc - ysize*0.5
            zloc = zloc - zsize*0.5
            do 400 k=1,num
               do 401 j=1,num
                  do 402 i=1,num
                     nflt=nflt+1
                     if(nflt.gt.fltdim) goto 2000
                     fltx(ndesc,nflt)=xloc+deltax*(i-1)
                     flty(ndesc,nflt)=yloc+deltay*(j-1)
                     fltz(ndesc,nflt)=zloc+deltaz*(k-1)
c
 402              continue
 401           continue
 400        continue
            fltdescsize(ndesc)=num*num*num
         endif
      endif
      goto 20
 2000 write(*,997)fltdim
 997  format('Exceeded maximum number of floats = ',i4)
      ier=2
      goto 100
 1000 write(*,999)curline
 999  format('Format error in float file'/'line: ',a)
      ier=1
 100  return
      end
