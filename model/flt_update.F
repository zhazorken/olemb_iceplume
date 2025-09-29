      subroutine flt_update(n)
      integer n
c
c     $Id: flt_update.F,v 1.3 1998/07/24 23:15:28 dwd Exp $
c
#include "param.inc"
#include "olemp.inc"
#include "olemf.inc"
#include "float.inc"
#include "cyclic.inc"
#include "moddef.inc"
c
      integer l,i,j,k,nn
      real tmpxu,tmpxw,tmpyu,tmpyv,tmpzu,tmpzw,xi,yi,zi
      real ui,u1,u2,u3,u4,u5,u6,u7,u8
      real vi,v1,v2,v3,v4,v5,v6,v7,v8
      real wi,w1,w2,w3,w4,w5,w6,w7,w8
      real xsize,ysize,zsize
c
      xsize=dx*idim
      ysize=dy*jdim
      zsize=dz*kdim
c
      do 10 l=1,ndesc
         if(n.ge.fltstrt(l) .and. n.le.fltstop(l))then
            do 11 nn=1,fltdescsize(l)

            tmpxu = fltx(l,nn)/dx+1.0
            tmpxw = fltx(l,nn)/dx+0.5
            tmpyu = flty(l,nn)/dy+0.5
            tmpyv = flty(l,nn)/dy+1.0
            tmpzu = kdim-fltz(l,nn)/dz+0.5
            tmpzw = kdim-fltz(l,nn)/dz+1.0
c
c u velocity
c
            i = ifix(tmpxu)
            j = ifix(tmpyu)
            k = ifix(tmpzu)
            xi = tmpxu - i
            yi = tmpyu - j
            zi = tmpzu - k
c
#if defined(STOKES)
           if(k.le.0)then
               u1=u(i,j,1)
               u2=u(i+1,j,1)
               u5=u(i,j+1,1)
               u6=u(i+1,j+1,1)
               ui=u1+(u2-u1)*xi+(u5-u1)*yi+(u6-u5+u1-u2)*xi*yi
            else if(k.ge.kdim)then
               u1=um(i,j,kdim)+us(kdim1)
               u2=um(i+1,j,kdim)+us(kdim1)
               u5=um(i,j+1,kdim)+us(kdim1)
               u6=um(i+1,j+1,kdim)+us(kdim1)
               ui=u1+(u2-u1)*xi+(u5-u1)*yi+(u6-u5+u1-u2)*xi*yi
            else
               u1=u(i,j,k)
               u2=u(i+1,j,k)
               u3=u(i+1,j,k+1)
               u4=u(i,j,k+1)
               u5=u(i,j+1,k)
               u6=u(i+1,j+1,k)
               u7=u(i+1,j+1,k+1)
               u8=u(i,j+1,k+1)
               ui=u1+(u4-u1)*zi+(u2-u1)*xi+(u1-u2+u3-u4)*zi*xi+
     *              (u5-u1)*yi+(u6-u5+u1-u2)*xi*yi+(u8-u5+u1-u4)*zi*yi+
     *              (u7-u8+u5-u6+u4-u3+u2-u1)*xi*yi*zi
            endif

#else
            if(k.le.0)then
               u1=u(i,j,1)
               u2=u(i+1,j,1)
               u5=u(i,j+1,1)
               u6=u(i+1,j+1,1)
               ui=u1+(u2-u1)*xi+(u5-u1)*yi+(u6-u5+u1-u2)*xi*yi
            else if(k.ge.kdim)then
               u1=u(i,j,kdim)
               u2=u(i+1,j,kdim)
               u5=u(i,j+1,kdim)
               u6=u(i+1,j+1,kdim)
               ui=u1+(u2-u1)*xi+(u5-u1)*yi+(u6-u5+u1-u2)*xi*yi
            else
               u1=u(i,j,k)
               u2=u(i+1,j,k)
               u3=u(i+1,j,k+1)
               u4=u(i,j,k+1)
               u5=u(i,j+1,k)
               u6=u(i+1,j+1,k)
               u7=u(i+1,j+1,k+1)
               u8=u(i,j+1,k+1)
               ui=u1+(u4-u1)*zi+(u2-u1)*xi+(u1-u2+u3-u4)*zi*xi+
     *              (u5-u1)*yi+(u6-u5+u1-u2)*xi*yi+(u8-u5+u1-u4)*zi*yi+
     *              (u7-u8+u5-u6+u4-u3+u2-u1)*xi*yi*zi
            endif
#endif
c
c v velocity
c
#if defined(STOKES)
            i = ifix(tmpxw)
            j = ifix(tmpyv)
            k = ifix(tmpzu)
            xi = tmpxw - i
            yi = tmpyv - j
            zi = tmpzu - k
c
            if(k.le.0)then
               v1=v(i,j,1)
               v2=v(i+1,j,1)
               v5=v(i,j+1,1)
               v6=v(i+1,j+1,1)
               vi=v1+(v2-v1)*xi+(v5-v1)*yi+(v6-v5+v1-v2)*xi*yi
            else if(k.ge.kdim)then
               v1=vm(i,j,kdim)+vs(kdim1)
               v2=vm(i+1,j,kdim)+vs(kdim1)
               v5=vm(i,j+1,kdim)+vs(kdim1)
               v6=vm(i+1,j+1,kdim)+vs(kdim1)
               vi=v1+(v2-v1)*xi+(v5-v1)*yi+(v6-v5+v1-v2)*xi*yi
            else
               v1=v(i,j,k)
               v2=v(i+1,j,k)
               v3=v(i+1,j,k+1)
               v4=v(i,j,k+1)
               v5=v(i,j+1,k)
               v6=v(i+1,j+1,k)
               v7=v(i+1,j+1,k+1)
               v8=v(i,j+1,k+1)
               vi=v1+(v4-v1)*zi+(v2-v1)*xi+(v1-v2+v3-v4)*zi*xi+
     *              (v5-v1)*yi+(v6-v5+v1-v2)*xi*yi+(v8-v5+v1-v4)*zi*yi+
     *              (v7-v8+v5-v6+v4-v3+v2-v1)*xi*yi*zi
            endif
#else
            i = ifix(tmpxw)
            j = ifix(tmpyv)
            k = ifix(tmpzu)
            xi = tmpxw - i
            yi = tmpyv - j
            zi = tmpzu - k
c
            if(k.le.0)then
               v1=v(i,j,1)
               v2=v(i+1,j,1)
               v5=v(i,j+1,1)
               v6=v(i+1,j+1,1)
               vi=v1+(v2-v1)*xi+(v5-v1)*yi+(v6-v5+v1-v2)*xi*yi
            else if(k.ge.kdim)then
               v1=v(i,j,kdim)
               v2=v(i+1,j,kdim)
               v5=v(i,j+1,kdim)
               v6=v(i+1,j+1,kdim)
               vi=v1+(v2-v1)*xi+(v5-v1)*yi+(v6-v5+v1-v2)*xi*yi
            else
               v1=v(i,j,k)
               v2=v(i+1,j,k)
               v3=v(i+1,j,k+1)
               v4=v(i,j,k+1)
               v5=v(i,j+1,k)
               v6=v(i+1,j+1,k)
               v7=v(i+1,j+1,k+1)
               v8=v(i,j+1,k+1)
               vi=v1+(v4-v1)*zi+(v2-v1)*xi+(v1-v2+v3-v4)*zi*xi+
     *              (v5-v1)*yi+(v6-v5+v1-v2)*xi*yi+(v8-v5+v1-v4)*zi*yi+
     *              (v7-v8+v5-v6+v4-v3+v2-v1)*xi*yi*zi
            endif
#endif
            if(flttype(l).eq. MDFLTXYZ)then
c
c w velocity, float is type neutral
c
               i = ifix(tmpxw)
               j = ifix(tmpyu)
               k = ifix(tmpzw)
               xi = tmpxw - i
               yi = tmpyv - j
               zi = tmpzw - k
c     
               w1=w(i,j,k)
               w2=w(i+1,j,k)
               w3=w(i+1,j,k+1)
               w4=w(i,j,k+1)
               w5=w(i,j+1,k)
               w6=w(i+1,j+1,k)
               w7=w(i+1,j+1,k+1)
               w8=w(i,j+1,k+1)
c
               wi=w1+(w4-w1)*zi+(w2-w1)*zi+(w1-w2+w3-w4)*zi*xi+
     *              (w5-w1)*yi+(w6-w5+w1-w2)*xi*yi+(w8-w5+w1-w4)*zi*yi+
     *              (w7-w8+w5-w6+w4-w3+w2-w1)*xi*yi*zi
            else
c
c float is not neutral
c
               wi=0.0
            endif
c
c vertical coordinate is depth (NOT Z!!!)
c
            fltx(l,nn)=fltx(l,nn)+ui*delt
            flty(l,nn)=flty(l,nn)+vi*delt
            fltz(l,nn)=fltz(l,nn)-wi*delt
c
c test for wrap and bounds
c
            if(fltx(l,nn).lt.0) fltx(l,nn) = xsize+fltx(l,nn)
            if(fltx(l,nn).gt.xsize) fltx(l,nn) = fltx(l,nn)-xsize
            if(flty(l,nn).lt.0) flty(l,nn) = ysize+flty(l,nn)
            if(flty(l,nn).gt.ysize) flty(l,nn) = flty(l,nn)-ysize
            if(fltz(l,nn).lt.0) fltz(l,nn) = 0.01*dz
            if(fltz(l,nn).gt.zsize) fltz(l,nn) = zsize-0.01*dz
 11         continue
         endif
 10   continue
      return
      end
