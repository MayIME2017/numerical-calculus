function func(x)
    implicit none
     double precision func, x
    func = 8.0*sqrt(900.0-x**2)+8.0*sqrt(400.0-x**2)-sqrt(360000.0-1300.0*x**2+x**4)
end function
