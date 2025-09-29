C Use projections of U and V onto the legendre polynomilas
C Coefficients and the order of the approximation
C should be read in for now they are hard coded (boo!)
C
      subroutine initializeUV(UA,VA)
#include 'param.inc'
#include 'UVCoef.inc'
        real, dimension(kdim),intent(inout) :: UA, VA
        real :: z
        integer :: i,l,m
        real, dimension(lmax,lmax) :: LegPolyCoef
C Set the Legendre Polynomial Coefficients
        LegPolyCoef = 0.0
C P_0(z) = 1
        LegPolyCoef(1,1) = +1.0
C P_2(z) = z
        LegPolyCoef(2,2) = +1.0
        do l = 1,lmax-2
         do m = 1,l+2
C Generate the Legendre Polynomial Coefficients for order >1
C (n+1) * P_{n+1}(z) = (2*n+1) * z * P_{n}(z) - n * P_{n-1}(z)
C n = l+1
         if(m.gt.1) then
          LegPolyCoef(l+2,m)= (2*l+1)*LegPolyCoef(l+1,m-1)
	 endif
         LegPolyCoef(l+2,m)= LegPolyCoef(l+2,m)-l*LegPolyCoef(l,m)
         LegPolyCoef(l+2,m)= LegPolyCoef(l+2,m)/(l+1)
         enddo
        enddo
C Fill in U and V
       do i=1,kdim
C Calculate z for the legendre polynomial z in [-1,1]
         z = -1+2*(i-.5)/(kdim+.5)
C Sum over legendre ploynomials
        do l=1,lmax
         do m=1,lmax
          UA(i) = UA(i)+UCoef(l)*LegPolyCoef(l,m)*(z**(m-1))
          VA(i) = VA(i)+VCoef(l)*LegPolyCoef(l,m)*(z**(m-1))
         enddo
        enddo
       enddo
      end subroutine
