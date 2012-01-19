module atom_def
    type atom
        double precision :: r(3)  ! atom position
        double precision, dimension(:), allocatable :: neigh !neighbour list
        integer :: nneigh ! number of neighbours
        integer :: intet
end type atom
end module atom_def

module tet_def
    type surface_normal
        double precision, dimension(3) :: vec
        double precision, dimension(3) :: inplane
    end type surface_normal

    type tetrahedron
        integer, dimension(4) :: atom
        double precision, dimension(4,3) :: vertex
        double precision, dimension(3) :: centre
        logical, dimension(3) :: shift
        type (surface_normal), dimension(4) :: snorm
        logical :: occupied
        integer :: occnum
    end type tetrahedron
end module tet_def

program defect_new

use atom_def
use tet_def
 
implicit none

type species
    type (atom), dimension(:), allocatable :: part
end type species  

double precision :: relr

type (species), dimension(:), allocatable :: spec
type (tetrahedron), dimension(:), allocatable :: tetra
character(len=20) :: posfile
character(len=20) :: inptfile
double precision :: boxlen(3), cboxlen(3), halfboxlen(3), halfcboxlen(3), h(3,3)
double precision :: shiftvec(3), thisp(3)

integer :: ntet, i, j, k, l, m
integer :: nspec, tetspec, occtet
integer, allocatable, dimension(:) :: nsp
integer :: thistet, dotsum
double precision, dimension(3) :: p1, p2, p3, p4, testp, tetspread(3)
double precision :: r1, r2, r3, r4
double precision :: invfac, dotprod
double precision :: rcut = 10.0
integer :: nconfigs, nstep
integer :: natomsout ! number of mobile atoms ?
integer, allocatable, dimension(:) :: tetlist
character(len=20) :: fmtout, fmtout2

!double precision :: rij(3), rijsq 
!integer :: dotprod, indshift

inptfile = "defect_new.inpt"

open(11, file=inptfile, status="old")
read(11,*) posfile
read(11,*) nconfigs
read(11,*) nspec
allocate (nsp(nspec))
allocate (spec(nspec))
do i=1, nspec
    read(11,*) nsp(i)
    allocate (spec(i)%part(nsp(i)))
enddo
read(11,*) tetspec
read(11,*) ntet
read(11,*) cboxlen(:)
read(11,*) h(:,:)
close(11)

!generate cell lengths for orthorhombic cell from boxlen(:)
boxlen(1) = cboxlen(1)*h(1,1)
boxlen(2) = cboxlen(2)*h(2,2)
boxlen(3) = cboxlen(3)*h(3,3)

halfcboxlen = cboxlen/2
halfboxlen = boxlen/2

allocate (tetra(ntet))

natomsout = 0
do i=1, nspec
    if (i /= tetspec) then
        natomsout = natomsout+nsp(i)
    endif
enddo

allocate(tetlist(natomsout))

tetra(:)%occupied = .false.

open(21, file="fort.21", status="old")
do i=1, ntet
    read(21,*)  tetra(i)%atom(1:4)
enddo
close(21)

open(10, file=posfile, status="old")
open(42, file="ntetocc.dat")

do nstep=1, nconfigs

tetra(:)%occnum = 0

    do i=1, nspec
        do j=1, nsp(i)
            spec(i)%part(j)%intet = 0
        enddo
    enddo

    do i=1, nspec
        do j=1, nsp(i)
            read(10,*) spec(i)%part(j)%r(1:3)
        enddo
    enddo

    occtet = 0

    do thistet = 1, ntet
        p1(:) = spec(tetspec)%part(tetra(thistet)%atom(1))%r(:)
        p2(:) = spec(tetspec)%part(tetra(thistet)%atom(2))%r(:)
        p3(:) = spec(tetspec)%part(tetra(thistet)%atom(3))%r(:)
        p4(:) = spec(tetspec)%part(tetra(thistet)%atom(4))%r(:)

        tetspread = 0.0
        do i=1, 3
            r1 = relr(p1(:), h, i)
            r2 = relr(p2(:), h, i)
            r3 = relr(p3(:), h, i)
            r4 = relr(p4(:), h, i)

            tetspread(i) = (max(r1 ,r2 ,r3 ,r4)-min(r1 ,r2 ,r3 ,r4))
            if (tetspread(i) > halfcboxlen(i)) then
                if (relr(p1(:), h, i) < halfboxlen(i)) then
                    p1(:) = p1(:) + h(i,:) * cboxlen(:)
                endif
        
                if (relr(p2(:), h, i) < halfboxlen(i)) then
                    p2(:) = p2(:) + h(i,:) * cboxlen(:)
                endif
        
                if (relr(p3(:), h, i) < halfboxlen(i)) then
                    p3(:) = p3(:) + h(i,:) * cboxlen(:)
                endif
        
                if (relr(p4(:), h, i) < halfboxlen(i)) then
                    p4(:) = p4(:) + h(i,:) * cboxlen(:)
                endif
            endif
        enddo
 
        tetra(thistet)%vertex(1,:) = p1(:) ! ions belong to multiple tetrahedra which may have different
        tetra(thistet)%vertex(2,:) = p2(:) ! periodic boundary requirements
        tetra(thistet)%vertex(3,:) = p3(:)
        tetra(thistet)%vertex(4,:) = p4(:)

    enddo

    do thistet = 1, ntet
        do i=1, 4

            p1(:) = tetra(thistet)%vertex(mod(i,4)+1,:)
            p2(:) = tetra(thistet)%vertex(mod(i+1,4)+1,:)
            p3(:) = tetra(thistet)%vertex(mod(i+2,4)+1,:)
            p4(:) = tetra(thistet)%vertex(mod(i+3,4)+1,:)
  
            tetra(thistet)%snorm(i)%inplane(:) = p1(:)

            tetra(thistet)%snorm(i)%vec(1) =   p2(2)*p3(3) - p2(3)*p3(2) & ! this can be generalised into a function f(2,3)
                                             - p1(2)*p3(3) + p1(3)*p3(2) &
                                             + p1(2)*p2(3) - p1(3)*p2(2)

            tetra(thistet)%snorm(i)%vec(2) = - p2(1)*p3(3) + p2(3)*p3(1) & ! = f(3,1)
                                             + p1(1)*p3(3) - p1(3)*p3(1) &
                                             - p1(1)*p2(3) + p1(3)*p2(1)

            tetra(thistet)%snorm(i)%vec(3) =   p2(1)*p3(2) - p2(2)*p3(1) & ! = f(1,2)
                                             - p1(1)*p3(2) + p1(2)*p3(1) &
                                             + p1(1)*p2(2) - p1(2)*p2(1)

            invfac = -sign(1.0,tetra(thistet)%snorm(i)%vec(1) * (p4(1)-p1(1)) &
                             + tetra(thistet)%snorm(i)%vec(2) * (p4(2)-p1(2)) &
                             + tetra(thistet)%snorm(i)%vec(3) * (p4(3)-p1(3)))
    
            tetra(thistet)%snorm(i)%vec(:) = tetra(thistet)%snorm(i)%vec(:) * invfac
        enddo

        do i=1, nspec
            if (i /= tetspec) then
                do j=1, nsp(i) 
                    do l=0, 7 ! loop over periodic images
                        shiftvec(:) = mod(l,2) * cboxlen(1) * h(1,:) &
                                    + mod(int(l/2),2) * cboxlen(2) * h(2,:) &
                                    + mod(int(l/4),2) * cboxlen(3) * h(3,:)
                        dotsum = 0
                        do k=1, 4 ! loop over each face of the tetrahedron
                            testp(:) = spec(i)%part(j)%r(:) + shiftvec(:) - tetra(thistet)%snorm(k)%inplane(:) 
                            dotprod = testp(1) * tetra(thistet)%snorm(k)%vec(1) &
                                    + testp(2) * tetra(thistet)%snorm(k)%vec(2) &
                                    + testp(3) * tetra(thistet)%snorm(k)%vec(3)
                            dotsum = dotsum - sign(1.0,dotprod)
                        enddo
                        if (dotsum == 4) then
                            tetra(thistet)%occupied = .true.
                            occtet = occtet + 1
                            if (tetra(thistet)%occnum /= 0) then
                                write(6,*) "warning! double occupation of tetrahedron ",thistet
                                write(6,*) "by ions ", tetra(thistet)%occnum, j
                            endif
                            tetra(thistet)%occnum = j
                            spec(i)%part(j)%intet = thistet
                        endif
                    enddo
                enddo
            endif
        enddo

    enddo

    write(6,*) "step ", nstep, ": number of occupied tetrahedra: ",occtet 
    write(42,*) nstep, occtet

    natomsout = 0
    tetlist(:) = 0
    do i=1, nspec
        if (i /= tetspec) then
            do j=1, nsp(i)
                natomsout = natomsout + 1
                tetlist(natomsout) = spec(i)%part(j)%intet
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

