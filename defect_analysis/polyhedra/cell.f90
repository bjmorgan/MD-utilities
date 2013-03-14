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

    pure function r_as_minimum_image( r )
        double precision, dimension(3) :: temp_r
        double precision, dimension(3), intent(in) :: r
        double precision, dimension(3) :: r_as_minimum_image
        integer :: k

        temp_r = r
        
        do k=1, 3
            temp_r(:) = temp_r(:) - cboxlen(k) * h(k,:) * int( relr( temp_r(:), k ) / halfcboxlen(k) )
        end do

        r_as_minimum_image = temp_r
    end function r_as_minimum_image

    pure function dr( r1, r2 )
        double precision, dimension(3), intent(in) :: r1, r2
        double precision, dimension(3) :: dr

        dr = r1 - r2

    end function dr

end module cell
