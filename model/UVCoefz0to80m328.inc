C velocity data from /media/mhoecker/8982053a-3b0f-494e-84a1-98cdce5e67d9/Dynamo/Observations/netCDF/ADCP/adcp150_filled_with_140.nc
        integer parameter :: lmax=5
        real dimension(lmax) :: UCoef,VCoef
C projection of U velocity from upper 80 m
        parameter(Ucoef=(/0.18596   -0.22003   -0.13806 -0.0050221    0.15151/))
C projection of V velocity from upper 80 m
        parameter(Vcoef=(/-0.069816    0.21766   -0.05713 -0.0060056   -0.11861/))
