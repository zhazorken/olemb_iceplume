        integer lmax
        parameter(lmax=6)
        real  UCoef(lmax),VCoef(lmax)
C projection of U velocity from 2011-Nov-24 upper 150 m
        data Ucoef/+0.090021,-0.324402,+0.048989,+0.143456,
     *             +0.028982,-0.119650/
C projection of V velocity from 2011-Nov-24 upper 150 m
        data Vcoef/-0.0278478, +0.1446929, +0.0751749, +0.0310151,
     *              +0.0487805, -0.0327094/
