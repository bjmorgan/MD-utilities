module class_Tetrahedron

    implicit none

    integer, parameter :: num_v = 4

    type surface_normal
        double precision, dimension(3) :: vec
        double precision, dimension(3) :: inplane
    end type surface_normal

    type Tetrahedron
        integer, dimension( num_v ) :: atom_id
        double precision, dimension( num_v, 3 ) :: vertex
        double precision, dimension( 3 ) :: centre
        logical, dimension( 3)  :: shift
        type (surface_normal), dimension( num_v ) :: snorm
        logical :: occupied
        integer :: occnum
        integer :: id
    contains
        procedure :: find_snorms
        procedure :: point_inside
        procedure :: occupied_by
end type Tetrahedron

contains

subroutine find_snorms( this )

    implicit none
    class(Tetrahedron) :: this
    integer :: i
    double precision, dimension(4,3) :: p
    double precision :: invfac

    do i=1, 4
        p(1,:) = this%vertex(mod(i,4)+1,:)
        p(2,:) = this%vertex(mod(i+1,4)+1,:)
        p(3,:) = this%vertex(mod(i+2,4)+1,:)
        p(4,:) = this%vertex(mod(i+3,4)+1,:)
  
        this%snorm(i)%inplane(:) = p(1,:)

        this%snorm(i)%vec(1) =   p(2,2)*p(3,3) - p(2,3)*p(3,2) & ! this can be generalised into a function f(2,3)
                               - p(1,2)*p(3,3) + p(1,3)*p(3,2) &
                               + p(1,2)*p(2,3) - p(1,3)*p(2,2)

        this%snorm(i)%vec(2) = - p(2,1)*p(3,3) + p(2,3)*p(3,1) & ! = f(3,1)
                               + p(1,1)*p(3,3) - p(1,3)*p(3,1) &
                               - p(1,1)*p(2,3) + p(1,3)*p(2,1)

        this%snorm(i)%vec(3) =   p(2,1)*p(3,2) - p(2,2)*p(3,1) & ! = f(1,2)
                               - p(1,1)*p(3,2) + p(1,2)*p(3,1) &
                               + p(1,1)*p(2,2) - p(1,2)*p(2,1)

        invfac = -sign(1.0,this%snorm(i)%vec(1) * (p(4,1)-p(1,1)) &
                         + this%snorm(i)%vec(2) * (p(4,2)-p(1,2)) &
                         + this%snorm(i)%vec(3) * (p(4,3)-p(1,3)))
    
        this%snorm(i)%vec(:) = this%snorm(i)%vec(:) * invfac
    enddo

end subroutine find_snorms

logical function point_inside( this, r )

    implicit none
    class(Tetrahedron) :: this
    integer :: k, dotsum
    double precision, dimension(3), intent(in) :: r
    double precision, dimension(3) :: testp
    double precision :: dotprod

    dotsum = 0
    do k=1, 4 ! loop over each face of the tetrahedron
        testp(:) = r(:) - this%snorm(k)%inplane(:) 
        dotprod = testp(1) * this%snorm(k)%vec(1) &
                + testp(2) * this%snorm(k)%vec(2) &
                + testp(3) * this%snorm(k)%vec(3)
        dotsum = dotsum - sign(1.0,dotprod)
    enddo
    if (dotsum == 4) then
        point_inside = .true.
    else
        point_inside = .false.
    endif

end function point_inside

subroutine occupied_by( this, id )

    implicit none
    class(Tetrahedron) :: this
    integer :: id
    this%occupied = .true.
    if (this%occnum /= 0) then
        write(6,*) "warning! double occupation of tetrahedron ", this%id
        write(6,*) "by ions ", this%occnum, id
    endif
    this%occnum = id
end subroutine occupied_by

end module class_Tetrahedron

module class_Octahedron

    integer, parameter :: num_v = 6

    type octahedron
        integer, dimension( num_v ) :: atom_id
        double precision, dimension( num_v , 3 ) :: vertex
        double precision, dimension(3) :: centre
        logical, dimension(3) :: shift
        logical :: occupied
        integer :: occnum
    end type octahedron
end module class_Octahedron
