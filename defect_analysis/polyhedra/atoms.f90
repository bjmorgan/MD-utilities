module atoms

    implicit none

    type atom
        double precision :: r(3)  ! atom position
        logical, dimension(:), allocatable :: neigh !neighbour list
        integer :: nneigh ! number of neighbours
        integer :: id
        integer, dimension(:), allocatable :: neighbour_id
    contains
        procedure :: set_neighbour_ids
    end type atom

contains

    subroutine set_neighbour_ids( this )
        class(atom) :: this
        integer, allocatable, dimension(:) :: id_list 
        integer :: i

        allocate( this%neighbour_id(this%nneigh) )
        allocate( id_list(size(this%neigh)) )

        forall (i=1:size(this%neigh)) id_list(i) = i

        this%neighbour_id = 0
        this%neighbour_id = pack(id_list, this%neigh)

    end subroutine set_neighbour_ids        

end module atoms
