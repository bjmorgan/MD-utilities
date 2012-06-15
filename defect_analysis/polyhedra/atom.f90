module atom

    type atom
        double precision :: r(3)  ! atom position
        logical, dimension(:), allocatable :: neigh !neighbour list
        integer :: nneigh ! number of neighbours
        integer :: id
    end type atom

end module atom
