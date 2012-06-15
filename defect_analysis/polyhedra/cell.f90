module cell

    double precision, dimension(3,3) :: h
    double precision, dimension(3) :: boxlen, cboxlen, halfboxlen, halfcboxlen, cpplane

contains

    pure double precision function relr(r,i)
        integer, intent(in) :: i
        double precision, dimension(3), intent(in) :: r
        double precision, dimension(3) :: temp_relr
    
        temp_relr(1) = r(1) - r(2)*h(2,1)/h(2,2) - r(3)*h(3,1)/h(3,3)
        temp_relr(2) = r(2) - r(1)*h(1,2)/h(1,1) - r(3)*h(3,2)/h(3,3)
        temp_relr(3) = r(3) - r(1)*h(1,2)/h(1,1) - r(2)*h(2,3)/h(2,2)

        relr = temp_relr(i)

    end function relr

end module cell
