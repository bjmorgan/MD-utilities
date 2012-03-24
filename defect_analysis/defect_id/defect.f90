module class_Atom
    
    implicit none

    type Atom
        double precision :: r(3) = 0.0 ! atom position
        double precision, dimension(:), allocatable :: neigh !neighbour list
        integer :: nneigh ! number of neighbours
        integer :: intet = 0
        integer :: inoct = 0
        integer :: id
    contains
        procedure :: site_id
    end type Atom

contains

integer function site_id( this ) 

    class(Atom), intent(in) :: this

    if ( this%intet > 0 ) then ! atom is occupying a tetrahedron
        site_id = this%intet
    else if ( this%inoct > 0 ) then ! atom is occupying an octahedron
        site_id = -this%inoct
    else ! atom is unassigned.
        site_id = 0
    endif
  
end function site_id

end module class_Atom

program defect_new

use class_atom
use class_tetrahedron
use class_octahedron
use class_file
 
implicit none

type species
    type (atom), dimension(:), allocatable :: particle
end type species  

interface 
    double precision function relr(r,h,i)
        implicit none
        integer :: i
        double precision, dimension(3) :: r
        double precision, dimension(3,3) :: h
    end function relr

    double precision function site_spread(r)
        implicit none
        double precision, dimension(:) :: r
    end function site_spread
end interface

type( species ), dimension(:), target, allocatable :: spec
type( tetrahedron ), dimension(:), target, allocatable :: tetra
type( octahedron ), dimension(:), allocatable :: octa
character( len=20 ) :: posfile_name, inptfile_name, tetfile_name, octfile_name
double precision :: boxlen(3), cboxlen(3), halfboxlen(3), halfcboxlen(3), h(3,3)
double precision :: shiftvec(3), thisp(3)

type(tetrahedron), pointer :: p_tetra
type(atom), pointer :: p_atom

integer :: ntet, noct, i, j, k, l, m, nv
integer :: nspec, tetspec, occtet
integer, allocatable, dimension(:) :: nsp
integer :: tet_id, oct_id
double precision, dimension(3) :: testp
double precision, dimension(4,3) :: p
double precision, dimension(4) :: r
double precision :: dotprod
double precision :: rcut = 10.0
integer :: nconfigs, nstep
integer :: natomsout ! number of mobile atoms ?
integer, allocatable, dimension(:) :: tetlist
character( len=20 ) :: fmtout, fmtout2
logical :: oct_analysis = .false.
integer :: io_status

type( iofile ) :: inptfile, posfile, tetfile, octfile

!double precision :: rij(3), rijsq 
!integer :: dotprod, indshift

inptfile_name = "defect_new.inpt"

call inptfile%init( inptfile_name )
call inptfile%fopen( "read"  )

read( inptfile%unit, * ) posfile_name
read( inptfile%unit, * ) nconfigs
read( inptfile%unit, * ) nspec
allocate (nsp(nspec))
allocate (spec(nspec))
do i=1, nspec
    read( inptfile%unit, * ) nsp(i)
    allocate (spec(i)%particle(nsp(i)))
    do j=1, nsp(i)
        spec(i)%particle(j)%id = j
    enddo
enddo
read( inptfile%unit, * ) tetspec
read( inptfile%unit, * ) ntet
read( inptfile%unit, * ) cboxlen(:)
read( inptfile%unit, * ) h(:,:)
read( inptfile%unit, *, iostat=io_status ) oct_analysis

if ( io_status < 0 ) then ! end of file reached. Assume old input file
    write(6,*) "Reached the end of the input file"
    write(6,*) "Assuming old input file, and only analysing tetrahedra"
else if ( io_status > 0 ) then ! read error
    write(6,*) "Error reading input file"
    stop
endif

if (oct_analysis) then
    read( inptfile%unit, * ) octfile_name
    read( inptfile%unit, * ) noct
endif

close( inptfile%unit )

!generate cell lengths for orthorhombic cell from boxlen(:)
boxlen(1) = cboxlen(1)*h(1,1)
boxlen(2) = cboxlen(2)*h(2,2)
boxlen(3) = cboxlen(3)*h(3,3)

halfcboxlen = cboxlen/2
halfboxlen = boxlen/2

allocate (tetra(ntet))
allocate (octa(noct))

natomsout = 0
do i=1, nspec
    if (i /= tetspec) then
        natomsout = natomsout+nsp(i)
    endif
enddo

allocate(tetlist(natomsout))

tetra(:)%occupied = .false.

tetfile_name = "fort.21"
call tetfile%init( tetfile_name )
call tetfile%fopen( "read" )
do i=1, ntet
    read( tetfile%unit, * ) tetra(i)%atom_id(1:4) ! read tetrahedral vertex ions
    tetra(i)%id = i
enddo
close( tetfile%unit )

if ( oct_analysis ) then
  call octfile%init( octfile_name )
  call octfile%fopen( "read" )
do i=1, noct
    read( octfile%unit, * ) octa(i)%atom_id(1:6) ! read octahedral vertex ions
enddo
close( octfile%unit )
endif

call posfile%init( posfile_name )
call posfile%fopen( "read" )

open(42, file="ntetocc.dat")

do nstep=1, nconfigs

tetra(:)%occnum = 0

    do i=1, nspec ! read ion positions
        do j=1, nsp(i)
            read( posfile%unit, * ) spec(i)%particle(j)%r(1:3)
        enddo
    enddo

    occtet = 0

    do tet_id = 1, ntet
        p_tetra => tetra(tet_id)
        do nv = 1, 4
            p(nv,:) = spec(tetspec)%particle(p_tetra%atom_id(nv))%r(:)
        enddo

        do i=1, 3
            do nv = 1, 4
                r(nv) = relr(p(nv,:), h, i)
            enddo

            if ( site_spread(r) > halfcboxlen(i) ) then
                do nv = 1, 4
                    if (relr(p(nv,:), h, i) < halfboxlen(i)) then
                        p(nv,:) = p(nv,:) + h(i,:) * cboxlen(:)
                    endif
                enddo
            endif
        enddo
 
        do nv = 1, 4
            p_tetra%vertex(nv,:) = p(nv,:) ! ions belong to multiple tetrahedra which may have different periodic boundary requirements
        enddo

        call p_tetra%find_snorms()

        do i=1, nspec
            if (i /= tetspec) then
                do j=1, nsp(i) 
                    p_atom => spec(i)%particle(j)
                    do l=0, 7 ! loop over periodic images
                        shiftvec(:) = mod(l,2) * cboxlen(1) * h(1,:) &
                                    + mod(int(l/2),2) * cboxlen(2) * h(2,:) &
                                    + mod(int(l/4),2) * cboxlen(3) * h(3,:)
                        if ( p_tetra%point_inside( p_atom%r(:) + shiftvec(:) ) ) then
                            call p_tetra%occupied_by( p_atom%id ) 
                            p_atom%intet = tet_id
                        endif
                    enddo
                    nullify(p_atom)
                enddo
            endif
        enddo

        nullify( p_tetra )
    enddo

    write(6,*) "step ", nstep, ": number of occupied tetrahedra: ",occtet 
    write(42,*) nstep, occtet

    if (oct_analysis) then ! perform analysis of labelled octahedra
        do oct_id = 1, noct
           do nv = 1, 6
               p(nv,:) = spec( tetspec )%particle( octa( oct_id )%atom_id( nv ))%r(:)
            enddo

            do i=1, 3
               do nv = 1, 6
                   r(nv) = relr(p(nv,:), h, i)
               enddo
    
               if ( site_spread(r) > halfcboxlen(i) ) then
                   do nv = 1, 6
                       if (relr(p(nv,:), h, i) < halfboxlen(i)) then
                           p(nv,:) = p(nv,:) + h(i,:) * cboxlen(:)
                       endif
                   enddo
                endif
            enddo
 
            do nv = 1, 4
                octa(oct_id)%vertex(nv,:) = p(nv,:) ! ions belong to multiple octahedra which may have different periodic boundary requirements
            enddo
        enddo 
    endif

    natomsout = 0
    tetlist(:) = 0
    do i=1, nspec
        if (i /= tetspec) then
            do j=1, nsp(i)
                natomsout = natomsout + 1
                tetlist(natomsout) = spec(i)%particle(j)%site_id() 
            enddo
        endif
    enddo

    write(fmtout, '(A4,I4,A7)') "(I5,",natomsout,"(I5,X))" ! internal write to define output formatting
    write(fmtout2, '(A4,I4,A7)') "(I5,",ntet,"(I5,X))" ! internal write to define output formatting
    write(40,fmtout) nstep, tetlist(:)
    write(41,fmtout2) nstep, tetra(:)%occnum
    write(401) nstep, tetlist(:)
    write(411) nstep, tetra(:)%occnum

enddo !ends loop over nconfig steps

close(42)

stop

end program defect_new

double precision function relr(r,h,i)

    implicit none
    integer :: i
    double precision, dimension(3) :: r
    double precision, dimension(3,3) :: h
    double precision, dimension(3) :: temp_relr

    temp_relr(1) = r(1) - r(2)*h(2,1)/h(2,2) - r(3)*h(3,1)/h(3,3)
    temp_relr(2) = r(2) - r(1)*h(1,2)/h(1,1) - r(3)*h(3,2)/h(3,3)
    temp_relr(3) = r(3) - r(1)*h(1,2)/h(1,1) - r(2)*h(2,3)/h(2,2)

    relr = temp_relr(i)

end function relr

double precision function site_spread( r )

    implicit none
    double precision, dimension(:) :: r

    site_spread = maxval(r) - minval(r)

end function site_spread
