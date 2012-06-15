program find_polyhedra

use atoms
use cell
use polyhedra

implicit none

type (atom), dimension(:), allocatable :: part
type (octahedron), dimension(:), allocatable :: octa
type (tetrahedron), dimension(:), allocatable :: tetra
character(len=20) :: posfile, inptfile
integer :: natoms, i, j, k, l, m, n, v, fout, dotprod, temp
integer :: fin
integer, dimension(:), allocatable :: v_mask, v_index
integer, dimension(6) :: v_list
double precision :: rcut, rcutsq
double precision :: rij(3), rijsq
logical, dimension(:), allocatable :: pair_list

interface

    function diagonal( square_matrix )
        double precision, dimension(:,:), intent(in) :: square_matrix
        double precision, dimension( size(square_matrix, 1) ) :: diagonal
    end function diagonal

    pure function ortho_pbc( r, boxlen )
        double precision, dimension(3), intent(in) :: r, boxlen
        double precision, dimension(3) :: ortho_pbc
    end function ortho_pbc

    subroutine read_input( inptfile, posfile, natoms, cboxlen, rcut, h, cpplane )
        double precision, dimension(3), intent(out) :: cboxlen, cpplane
        character(len=20), intent(out) :: posfile
        character(len=*), intent(in) :: inptfile
        integer, intent(out) :: natoms
        double precision, intent(out) :: rcut
        double precision, intent(out) :: h(3,3)
    end subroutine read_input

end interface

inptfile = "getocts.inpt"
call read_input( inptfile, posfile, natoms, cboxlen, rcut, h, cpplane )

rcutsq = rcut * rcut
boxlen = cboxlen * diagonal(h)
halfboxlen = boxlen/2
halfcboxlen = boxlen/2 !warning. Should this be cboxlen/2?

allocate (part(natoms))
allocate (octa(natoms))
allocate (tetra(natoms*2))
allocate (pair_list(natoms), v_mask(natoms), v_index(natoms))

open(file=posfile, status='old', form='formatted', newunit=fin)
do i=1, natoms
    allocate (part(i)%neigh(natoms))
    read(fin,*) part(i)%r(1:3)
end do

forall (i=1:natoms) 
    part(i)%neigh(:) = .false.
    part(i)%neigh(i) = .true.
    part(i)%id = i
    ! map to an orthorhombic cell, assuming rhombohedral input
    part(i)%r = ortho_pbc(part(i)%r, boxlen) 
    v_index(i) = i
end forall

! create neighbour list
! assuming the input file is in lab coordinates
do i=1, natoms
    do j=i+1, natoms
        rij = part(i)%r - part(j)%r
        ! minimum image convention
        do k=1, 3
            rij(:) = rij(:)-cboxlen(k) * h(k,:) * int(relr(rij(:), k) / halfcboxlen(k))
        end do
        rijsq = sum(rij*rij)
        if (rijsq < rcutsq) then
            part(i)%neigh(j) = .true. 
            part(j)%neigh(i) = .true.    
        end if
    end do
end do


forall (i=1:natoms) part(i)%nneigh = count(part(i)%neigh(:))

if (any(part%nneigh /= 13)) then
    write(*) "Not closed-packed, or rcut incorrect"
    stop
end if

do i=1, natoms
    call part(i)%set_neighbour_ids
end do

! find octahedra
do i=1, natoms-1
    do j=i+1, natoms
        pair_list = ( part(i)%neigh .and. part(j)%neigh )
        if (count(pair_list) == 4) then
            pair_list(i) = .true.
            pair_list(j) = .true.
            v_list = pack(v_index, pair_list)
            if (.not.oct_exists(v_list, octa(1:noct))) then
                noct = noct + 1
                associate( oct => octa(noct) )
                    call oct%init
                    call oct%set_vertices( pack(part, pair_list) )
                end associate
            end if
        end if
    end do
end do

!find tetrahedra
!tetrahedra are defined by sets of four atoms, where any triplet are neighbours of the fourth
do i=1, natoms-3
    do j=i+1, natoms-2
        if ( .not.part(j)%neigh(i) )cycle
        do k=j+1, natoms-1
            if ( .not. ( part(k)%neigh(j) .and. part(k)%neigh(i) ) ) cycle
            do l=k+1, natoms
                if ( .not. ( part(l)%neigh(k) .and. part(l)%neigh(j) .and. part(l)%neigh(i) ) ) cycle
                pair_list = .false.
                pair_list(i) = .true.
                pair_list(j) = .true.
                pair_list(k) = .true.
                pair_list(l) = .true.
                ntet = ntet + 1
                associate( tet => tetra(ntet) )
                    call tet%init
                    call tet%set_vertices( pack(part, pair_list) )
                end associate
            end do
        end do
    end do
end do

write(6,*) noct,"octahedra found"
write(6,*) ntet,"tetrahedra found"

do j=1, ntet
    call tetra(j)%enforce_pbc
end do

do i=1, noct
    call octa(i)%enforce_pbc
end do

do i=1, ntet
    tetra(i)%orientation = 0
    do j=1, 4
        temp = -int( sign( 1.0,sum( (tetra(i)%vertex(j)%r-tetra(i)%centre)*cpplane ) ) )
        tetra(i)%orientation = tetra(i)%orientation + temp
    end do
    tetra(i)%orientation = sign(1,tetra(i)%orientation)
end do

call write_output(tetra, octa)

contains

subroutine write_output( tetra, octa )

    use polyhedra

    implicit none

    type(tetrahedron), dimension(:), intent(in) :: tetra
    type(octahedron), dimension(:), intent(in) :: octa
    integer :: ntet_up = 0, ntet_down = 0
    integer :: fcent, ftet1, ftet2, foct
    integer i

    open(file='centres.out', newunit=fcent)
    open(file='tet1.list', newunit=ftet1)
    open(file='tet2.list', newunit=ftet2)
    open(file='oct.list', newunit=foct)

    do i=1, size(tetra)    
        write(fcent,*) tetra(i)%centre
        if (tetra(i)%orientation == 1) then
            ntet_up = ntet_up + 1
            write(ftet1,*) tetra(i)%vertex%id
        else
            ntet_down = ntet_down + 1
            write(ftet2,*) tetra(i)%vertex%id
        end if
    end do

    do i=1, size(octa)
        write(fcent,*) octa(i)%centre
        write(foct,*) octa(i)%vertex%id
    end do

    close(fcent)
    close(ftet1)
    close(ftet2)
    close(foct)

    write(6,*) ntet_up,'tet1',ntet_down,'tet2',noct,'oct'

end subroutine write_output
    
end program find_polyhedra

subroutine read_input( inptfile, posfile, natoms, cboxlen, rcut, h, cpplane )
    implicit none
    double precision, dimension(3), intent(out) :: cboxlen, cpplane
    character(len=20), intent(out) :: posfile
    character(len=*), intent(in) :: inptfile
    integer, intent(out) :: natoms
    double precision, intent(out) :: rcut
    double precision, intent(out) :: h(3,3)
    integer :: fin

    open(file=inptfile, status='old', form='formatted', newunit=fin)
        read(fin,*) posfile
        read(fin,*) natoms
        read(fin,*) cboxlen(1:3)
        read(fin,*) rcut
        read(fin,*) h(1:3,1:3)
        read(fin,*) cpplane(1:3)
    close(fin)
end subroutine read_input
 
function diagonal( square_matrix )
    implicit none
    double precision, dimension(:,:), intent(in) :: square_matrix
    double precision, dimension(size(square_matrix)) :: temp_array
    double precision, dimension(size(square_matrix, 1)) :: diagonal
    temp_array = pack(square_matrix, .true.)
    diagonal = temp_array(1::size(diagonal)+1)
end function diagonal

pure function ortho_pbc( r, boxlen )
    implicit none
    double precision, dimension(3), intent(in) :: r
    double precision, dimension(3), intent(in) :: boxlen
    double precision, dimension(3) :: ortho_pbc
    integer :: i
    forall(i=1:3) ortho_pbc(i) = r(i) - (int(r(i)/boxlen(i)) * boxlen(i)) ! map to orthorhombic cell
end function ortho_pbc
