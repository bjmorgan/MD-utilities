module atom_def
    type atom
        double precision :: r(3)  ! atom position
        double precision, dimension(:), allocatable :: neigh !neighbour list
        integer :: nneigh ! number of neighbours
    end type atom
end module atom_def

module tet_def
    type tetrahedron
        integer, dimension(4) :: atom
        double precision, dimension(4,3) :: vertex
        double precision, dimension(3) :: centre
        logical, dimension(3) :: shift
        integer :: orientation 
    end type tetrahedron
end module tet_def

program gettet

! Finds tetrahedra in a close-packed system.
! Generates neighour lists for all atoms of a given species
! These are given as vector neighbour lists e.g. v(a)=[a b c ... ]
! Which includes the central ion in the list.
! Tetrahedra are given by those sets of four neighbour lists
! where |v(1).v(2).v(3).v(4)| = 4

use atom_def
use tet_def

implicit none

type (atom), dimension(:), allocatable :: part
type (tetrahedron), dimension(:), allocatable :: tetra
character(len=20) :: posfile
character(len=20) :: inptfile
double precision :: rcut, rcutsq
double precision :: boxlen(3), cboxlen(3), halfboxlen(3), halfcboxlen(3), h(3,3)
double precision, dimension(3) :: p1, p2, p3, p4, tetspread
double precision :: r1, r2, r3, r4
integer :: natoms, ntet, i, j, k, l, m, thistet
double precision :: rij(3), rijsq 
integer :: dotprod, tetcnt
double precision :: rijsqav
integer :: nrijsq
logical :: clpackflag
integer :: temp
double precision, dimension(3) :: tempv
double precision, dimension(3) :: cpplane

double precision :: relr

ntet = 0
inptfile = "gettet.inpt"

open(11, file=inptfile, status="old")
read(11,*) posfile
read(11,*) natoms
read(11,*) cboxlen(1:3)
read(11,*) rcut
read(11,*) h(1:3,1:3)
read(11,*) cpplane(1:3)

rcutsq = rcut*rcut

!generate cell lengths for orthorhombic cell from boxlen(:)
boxlen(1) = cboxlen(1)*h(1,1) 
boxlen(2) = cboxlen(2)*h(2,2) 
boxlen(3) = cboxlen(3)*h(3,3) 

halfboxlen = boxlen/2
halfcboxlen = boxlen/2

allocate ( part( natoms ) )
allocate ( tetra( natoms*2 ) )

do i=1, natoms
    allocate (part(i)%neigh(natoms))
    part(i)%neigh(:) = 0
    part(i)%neigh(i) = 1
enddo

open(10, file=posfile, status="old")
do i=1, natoms
    read(10,*) part(i)%r(1:3)
enddo

! create neighbur list
! assuming the input file is in cartesian coordinates (which makes implementing this a pain)
! it would probably be easier to transform to fractional coordinates?

rijsqav = 0.0
nrijsq = 0

do i=1, natoms-1
    do j=i+1, natoms

        rij = part(i)%r - part(j)%r

        rij(:)=rij(:)-(cboxlen(:) * h(1,:) * int((relr(rij(:), h, 1) / halfcboxlen(1)))) ! Minimum image convention
        rij(:)=rij(:)-(cboxlen(:) * h(2,:) * int((relr(rij(:), h, 2) / halfcboxlen(2))))
        rij(:)=rij(:)-(cboxlen(:) * h(3,:) * int((relr(rij(:), h, 3) / halfcboxlen(3))))

        rijsq = (rij(1)*rij(1) + &
                 rij(2)*rij(2) + &
                 rij(3)*rij(3))

        !write(6,*) rijsq

        !rijsq = sum( rij(:) * rij(:) )

        !write(6,*) rijsq

        if (rijsq < rcutsq) then 
            part(i)%neigh(j) = 1    ! construct neighbour list
            part(j)%neigh(i) = 1    
            
            rijsqav = rijsqav + rijsq
            nrijsq = nrijsq + 1
        endif

    enddo
enddo

write(6,*) "average (rij < rcut) =", sqrt(rijsqav/nrijsq)

clpackflag = .false.
do i=1, natoms
    part(i)%nneigh = sum(part(i)%neigh(:))
    if (part(i)%nneigh /= 13) then
       clpackflag = .true.
    endif
enddo

if (clpackflag) then
    do i=1, natoms
        write(50,*) i, part(i)%nneigh
    enddo
    write(6,*) "rcut not set correctly for identifying all close-packed neighbours"
    write(6,*) "rcut = ", rcut
    write(6,*) "average number of neighbours = ", sum(part(:)%nneigh)/natoms
    stop
endif

do i=1, natoms-3
    write(6,*) i
    do j=i+1, natoms-2
        do k=j+1, natoms-1
            do l=k+1, natoms
                dotprod = part(j)%neigh(i) * part(k)%neigh(i) * part(l)%neigh(i) &
                        + part(i)%neigh(j) * part(k)%neigh(j) * part(l)%neigh(j) &
                        + part(i)%neigh(k) * part(j)%neigh(k) * part(l)%neigh(k) &
                        + part(i)%neigh(l) * part(j)%neigh(l) * part(k)%neigh(l) 
                if (dotprod == 4) then
                    ntet = ntet + 1
 
                    tetra(ntet)%atom(1) = i
                    tetra(ntet)%atom(2) = j
                    tetra(ntet)%atom(3) = k
                    tetra(ntet)%atom(4) = l

                    tetra(ntet)%vertex(1,:) = part(i)%r
                    tetra(ntet)%vertex(2,:) = part(j)%r
                    tetra(ntet)%vertex(3,:) = part(k)%r
                    tetra(ntet)%vertex(4,:) = part(l)%r
                
                endif
            enddo
        enddo
    enddo
enddo

! enforce periodic boundaries for tetrahedra. 
! ensures tetrahedra are all "real space"

do thistet=1, ntet
    p1(:) = tetra(thistet)%vertex(1,:) 
    p2(:) = tetra(thistet)%vertex(2,:) 
    p3(:) = tetra(thistet)%vertex(3,:) 
    p4(:) = tetra(thistet)%vertex(4,:) 

    if (thistet == 1) then
        write(6,*) p1
        write(6,*) p2
        write(6,*) p3
        write(6,*) p4
    endif

    tetspread = 0.0
    do i=1, 3
        r1 = relr(p1(:), h, i)
        r2 = relr(p2(:), h, i)
        r3 = relr(p3(:), h, i)
        r4 = relr(p4(:), h, i)

        tetspread = (max(r1, r2, r3, r4)-min(r1, r2, r3, r4))
        if (tetspread(i) > halfboxlen(i)) then
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

    tetra(thistet)%vertex(1,:) = p1(:)
    tetra(thistet)%vertex(2,:) = p2(:)
    tetra(thistet)%vertex(3,:) = p3(:)
    tetra(thistet)%vertex(4,:) = p4(:)

    if (thistet == 1) then
        write(6,*) p1
        write(6,*) p2
        write(6,*) p3
        write(6,*) p4
    endif
enddo

! calculate tetrahedra centres. 
write(6,*) ntet,"tetrahedra found"

do i=1, ntet
    tetra(i)%shift = .false.
    do l=1, 3
        tetra(i)%centre(l) = (sum(tetra(i)%vertex(:,l)))/4
    enddo
   
    if (i == 1) then
        write(6,*) tetra(i)%centre
    endif

    tetra(i)%orientation = 0
    do j=1, 4
        temp = int(sign(1.0,sum((tetra(i)%vertex(j,:)-tetra(i)%centre(:))*cpplane(:))))
        tetra(i)%orientation = tetra(i)%orientation - temp
    enddo

    tetra(i)%orientation = sign(1,tetra(i)%orientation)
   
enddo

do i=1, ntet
    do j=1, 4 ! swap vertices so that tetra(i)%atom(1) is the tetrahedron vertex (normal to close packed planes)
        if (int(sign(1.0,sum((tetra(i)%vertex(j,:)-tetra(i)%centre(:))*cpplane(:)))) == tetra(i)%orientation) then

            temp = tetra(i)%atom(1)
            tetra(i)%atom(1) = tetra(i)%atom(j)
            tetra(i)%atom(j) = temp
 
            tempv(:) = tetra(i)%vertex(1,:)
            tetra(i)%vertex(1,:) = tetra(i)%vertex(j,:)
            tetra(i)%vertex(j,:) = tempv(:)

        endif
    enddo
enddo

tetcnt = 0
do i=1, ntet    
    if (tetra(i)%orientation == 1) then
        tetcnt = tetcnt + 1
        write(20,*) tetra(i)%centre(:)
        write(21,*) tetra(i)%atom(:)
    endif
enddo
write(6,*) tetcnt

tetcnt = 0
do i=1, ntet
    if (tetra(i)%orientation == -1) then
        tetcnt = tetcnt + 1
        write(20,*) tetra(i)%centre(:)
        write(21,*) tetra(i)%atom(:)
    endif
enddo
write(6,*) tetcnt

end program gettet

double precision function relr(r,h,i)

integer :: i
double precision, dimension(3) :: r
double precision, dimension(3,3) :: h
double precision, dimension(3) :: temp_relr

temp_relr(1) = r(1) - r(2)*h(2,1)/h(2,2) - r(3)*h(3,1)/h(3,3)
temp_relr(2) = r(2) - r(1)*h(1,2)/h(1,1) - r(3)*h(3,2)/h(3,3)
temp_relr(3) = r(3) - r(1)*h(1,2)/h(1,1) - r(2)*h(2,3)/h(2,2)

relr = temp_relr(i)

end function relr
