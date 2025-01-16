real function init_p(x, y)
    real x, y
    !init_p = x+3*y
    init_p = x**2+y**2+2*x*y+x+y
    !init_p = x**3+5*y**3+x**2*y+y**2*x+x*y+x+y
    !init_p = x**4+y**4+x**2+y**2
end function

subroutine calc_gradP_exact(x, y, gradP)
    real x, y
    real gradP(2)

    !gradP(1) = 1
    !gradP(2) = 3

    !gradP(1) = 1+2*x+2*y
    !gradP(2) = 1+2*x+2*y

    gradP(1) = 1+3*x**2+y+y**2+2*x*Y
    gradP(2) = 1+x+x**2+2*x*y+15*y**2


end subroutine

real function calc_laplP_exact(x,y)
    real x,y
    calc_laplP_exact = 4
    !calc_laplP_exact = 8*(x+4*y)
    !calc_laplP_exact = 12*(x**2+Y**2)+4
end function


subroutine init_V(x,y,V)
    real x, y
    real V(2)

    !V(1) = 1
    !V(2) = 1

    !V(1) = y
    !V(2) = -3*x

    !V(1) = y**2-4*x*y+y
    !V(2) = x**2+x*y+x

    V(1) = -3*y**3-x*y-y
    V(2) = x**3-x*y+4*x

end subroutine

real function calc_divV_exact(x,y)
    real x,y
    !calc_divV_exact = 4
    calc_divV_exact = 2+4*x+4*y
    !calc_divV_exact = 2 + x + 4*x**2+y+4*x*y+16*y**2
end function

real function calc_rotV_exact(x,y)
    real x, y

    !calc_rotV_exact = -4
    !calc_rotV_exact = 6*x-y
    calc_rotV_exact = 3*x**2+9*y**2+x-y+5
end function



! Линейная интерполяция: выч знач в (.) между 2 изв знач
! d1, d2 - расстояния от (.) интерп до 2 изв (.)
! p1, p2 - знач в 2 изв (.)
real function RLinearInterp(d1,d2,p1,p2)
    real :: d1, d2, p1, p2
    RLinearInterp = (d1*p2 + d2*p1)/(d1+d2)
End Function