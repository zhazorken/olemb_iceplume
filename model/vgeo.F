#if defined(FRONT)
c
c apply thermal wind
c
         ml = kdim-20./dz
         do k=1,kdim
           if(k.gt.ml) then
             coef = 1.0
           else
             coef = exp((k-ml)/5.)
           endif

           do j=1,jdim
             do i=1,idim
               ii = iblkstart(my_pe+1)+i
               jj = jblkstart(my_pe+1)+j
               if(ii.lt.itotal/2) then
                 vz = coef*grav*
     *               0.5/cosh((real(ii-itotal/3.))/300.)**2
               else
                 vz = -grav*coef*
     *               0.5/cosh((real(ii-2.*itotal/3))/300.)**2
c     *               0.5/cosh((real(ii-2.*itotal/3))/80.)**2
               endif
               vz = 2.e-4*vz/(rhob(k)*tomegaz)
               vm(i,j,k) = vm(i,j,k-1)+vz*dz
               v(i,j,k) = vm(i,j,k)
               vmm(i,j,k) = vm(i,j,k)
             enddo
           enddo
         enddo 
#endif
