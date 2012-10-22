module tetrahedra
    
    use polyhedra

    implicit none
    
    integer :: ntet = 0
    integer, parameter :: ntet_vertices = 4

    type, extends(polyhedron) :: tetrahedron
        integer :: orientation
    contains
        procedure :: init => init_tet
    end type tetrahedron

    contains

    subroutine init_tet( this )
        class(tetrahedron) :: this
        this%num_vert = ntet_vertices
        call this%alloc_vertices
    end subroutine init_tet

end module tetrahedra
