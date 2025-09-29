          deptho = 0.5*dz
          cd = 0.16/(log(deptho/0.10))**2
          uw(i,j,olev+1) = randn(i,j)*cd*um(i,j,1)*
     *       abs(um(i,j,1))
          vw(i,j,olev+1) = randn(i,j)*cd*(vm(i,j,1))*
     *       abs(vm(i,j,1))
          endif
          uistar(i,j) = sqrt(abs(uw(i,j,2))+abs(vw(i,j,2)))

