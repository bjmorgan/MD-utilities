module faces

    use atoms

    type face
        double precision, dimension(3) :: normal
        type(atom), dimension(3) :: vertex
    contains
        procedure :: define_normal
    end type face

    contains
  
    subroutine define_normal(this, point_outside)
        implicit none
        type point
            double precision, dimension(3) :: r
        end type point
        class(face) :: this
        type(point), dimension(4) :: p
        double precision, dimension(3), intent(in) :: point_outside
        double precision :: invfac
        double precision, dimension(3) :: normal

        p(1)%r = this%vertex(1)%r
        p(2)%r = this%vertex(2)%r
        p(3)%r = this%vertex(3)%r
        p(4)%r = point_outside

        normal(1) =   p(2)%r(2)*p(3)%r(3) - p(2)%r(3)*p(3)%r(2) & ! this can be generalised into a function f(2,3)
                    - p(1)%r(2)*p(3)%r(3) + p(1)%r(3)*p(3)%r(2) &
                    + p(1)%r(2)*p(2)%r(3) - p(1)%r(3)*p(2)%r(2)

        normal(2) = - p(2)%r(1)*p(3)%r(3) + p(2)%r(3)*p(3)%r(1) & ! = f(3,1)
                    + p(1)%r(1)*p(3)%r(3) - p(1)%r(3)*p(3)%r(1) &
                    - p(1)%r(1)*p(2)%r(3) + p(1)%r(3)*p(2)%r(1)

        normal(3) =   p(2)%r(1)*p(3)%r(2) - p(2)%r(2)*p(3)%r(1) & ! = f(1,2)
                    - p(1)%r(1)*p(3)%r(2) + p(1)%r(2)*p(3)%r(1) &
                    + p(1)%r(1)*p(2)%r(2) - p(1)%r(2)*p(2)%r(1)

        invfac = -sign(1.0, normal(1) * (p(4)%r(1)-p(1)%r(1)) &
                          + normal(2) * (p(4)%r(2)-p(1)%r(2)) &
                          + normal(3) * (p(4)%r(3)-p(1)%r(3)))

        this%normal = normal * invfac
    end subroutine define_normal

end module faces
