module octahedra
    
    use polyhedra

    implicit none
    
    integer :: noct = 0
    integer, parameter :: noct_vertices = 6
    integer, parameter :: noct_faces = 8
    integer :: occoct
    character(len=string_length), parameter :: oct_string = 'octahedron'

    type, extends(polyhedron) :: octahedron
    contains
        procedure :: init => init_oct
    end type octahedron

    contains

    subroutine init_oct( this )
        class(octahedron) :: this
        this%num_vert = 6
        call this%alloc_vertices
    end subroutine init_oct

    logical function oct_exists( vertex_ids, octs )
        integer, dimension(6), intent(in) :: vertex_ids
        class(octahedron), dimension(:), intent(in) :: octs
        integer i
        do i=1, size(octs)
            if (count(sort(vertex_ids) .eq. sort(octs(i)%vertex_ids)) == 6) then
                oct_exists = .true.
                return
            end if
        end do
        oct_exists = .false.
    end function oct_exists

    function sort(x)
        implicit none

        ! returns a sorted array using an insertion sort

        integer, dimension(:), intent(in) :: x
        integer, dimension(size(x)) :: sort
        integer :: lindex, i, j

        sort = x
        do i=2, size(x)
!    If the Ith element is out of order with the preceeding
!    element search for its position in the portion of the
!    array that is already ordered.
            if ( sort(i) < sort(i-1) ) then
                lindex = 1
                do j=i-2,1,-1
                    if(sort(j) < sort(i)) then
                        lindex = j+1
                        exit 
                    endif
                enddo
                sort(lindex:i) = cshift(sort(lindex:i),-1)
            endif
        enddo
    end function sort 

end module octahedra
