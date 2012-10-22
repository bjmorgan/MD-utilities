module polyhedra
    
    use atoms
    use cell

    implicit none
    integer, parameter :: string_length = 11
 
    type :: polyhedron
        double precision, dimension(3) :: centre
        type(atom), dimension(:), allocatable :: vertex
        integer, dimension(:), allocatable :: vertex_ids
        integer :: num_vert
        character(len=string_length) :: string
    contains
        procedure :: set_vertex
        procedure :: set_vertices
        procedure :: alloc_vertices
        procedure :: enforce_pbc
        procedure, private :: set_centre
    end type polyhedron

    contains

    subroutine alloc_vertices( this )
        class(polyhedron) :: this
        allocate( this%vertex(this%num_vert) )
        allocate( this%vertex_ids(this%num_vert) )
    end subroutine alloc_vertices

    subroutine set_vertex( this, v, this_atom )
        class(polyhedron) :: this
        integer :: v
        type(atom) :: this_atom

        this%vertex(v) = this_atom
        this%vertex_ids(v) = this_atom%id
    end subroutine set_vertex

    subroutine set_vertices( this, atoms )
        class(polyhedron) :: this
        integer :: i
        type(atom), dimension(:) :: atoms
        do i=1, this%num_vert
            this%vertex(i) = atoms(i)
            this%vertex_ids(i) = atoms(i)%id
        end do
        call this%set_centre
    end subroutine set_vertices

    subroutine set_centre( this )
        class(polyhedron) :: this
        integer :: l
        forall (l=1:3)
            this%centre(l) = sum( this%vertex(:)%r(l) ) / this%num_vert
        end forall
    end subroutine set_centre

    subroutine enforce_pbc( this )
        class(polyhedron) :: this
        type point
            double precision, dimension(3) :: r
        end type point
        integer :: i, k
        type(point), dimension(this%num_vert) :: p
        double precision :: tetspread
        double precision, dimension(this%num_vert) :: dr

        forall (k=1:this%num_vert) p(k)%r = this%vertex(k)%r
        do i=1, 3
            forall (k=1:this%num_vert) dr(k) = relr(p(k)%r, i)
            tetspread = maxval(dr)-minval(dr)
            if (tetspread > halfboxlen(i)) then
                do k=1, this%num_vert
                    if (relr(p(k)%r, i) < halfboxlen(i)) then
                        p(k)%r = p(k)%r + h(i,:) * cboxlen
                    end if
                end do
            end if
        end do
        forall (k=1:this%num_vert) this%vertex(k)%r = p(k)%r
        call this%set_centre ! update centre point
    end subroutine enforce_pbc

end module polyhedra
