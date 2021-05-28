function dfunc(x)
    implicit none
     double precision dfunc, x
    dfunc = -8.0*x/sqrt(900.0-x**2)-8.0*x/sqrt(400.0-x**2)-(4*x**3-2600.0*x)/sqrt(360000.0-1300.0*x**2+x**4)
end function
