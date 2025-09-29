#if defined(FRONT)
c
c         if(iope) then
           thinmin = th(1,1,kdim)-1.0
           ml = kdim-20./dz
           do k=0,kdim1
c             if(k.gt.ml) then
c               xof = (500.*(k - ml)/(kdim-ml))/dx
c               if(iope) write(*,*)"xof = ",xof,k
c             endif
             do j=1,jdim
               do i=1,idim
                
c
                 ii = iblkstart(my_pe+1)+i
                 jj = jblkstart(my_pe+1)+j
c
c double front warm in the middle
c
               if(k.gt.ml) then
                 coef = 1.0
               else
                 coef = exp((k-ml)/5.)
               endif
               if(ii.lt.itotal/2) then
                 phi(i,j,k) = th(i,j,k)+
     *               0.1*tanh((real(ii-itotal/3.))/80.)*coef
c     *               0.04*tanh((float(i-idim/3)+xof)/100.)
c     *               0.2*tanh(float(i-idim/4)/100.)
               else
                 phi(i,j,k) = th(i,j,k)-
     *            0.1*tanh((real(ii-2.*itotal/3.))/80.)*coef
c     *             0.04*tanh((float(i-idim/2-(idim/2-idim/3))-xof)/100.)
c     *             0.2*tanh(float(i-idim/2-idim/4)/100.)
               endif
c               if(i+iblkstart(my_pe+1).gt.509.and.
c     *            i+iblkstart(my_pe+1).lt.514
c     *             .and.k.ge.kdim-2) then
c                 tke(i,j,k) = 1.0
c               else if(i.gt.330.and.i.lt.340
c     *             .and.k.gt.23) then
c                 tke(i,j,k) = 1.0
c               else
c                 tke(i,j,k) = 0.0
c               endif
               enddo
             enddo
           enddo
c
c turn off front
c
           do k=0,kdim1
           do j=1,jdim
           do i=1,idim
             if(phi(i,j,k).gt.thinmin) then
               th(i,j,k) = phi(i,j,k)
             endif
           enddo
           enddo
           enddo
c         call imag(tke)
#endif
