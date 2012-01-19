program defectrdf

type centype
    double precision :: r(3)
    character(len=4) :: ctype
end type centype

type(centype), dimension(:), allocatable :: centre

character(len=20) :: infile = "defectrdf.inpt"
character(len=20) :: deffile = "defectpos.dat"
character(len=20) :: vacfile = "vaccentpos.dat"
character(len=30) :: outfile 

double precision :: boxlen(3), hboxlen(3)
double precision :: dr(3)
double precision :: rdfmax
double precision :: xnorm, rdfdr, rdfr
double precision :: cellvol
double precision :: xint

integer :: nframes
integer :: ncen
character(len=4) :: spec1, spec2 
integer :: i, j, k
integer :: nbins
integer :: norm1, norm2
integer :: axis
integer, allocatable, dimension(:) :: rdfout
double precision, allocatable, dimension(:) :: rdffunc
integer :: ibin

open(10, file=infile, status="old")
read(10,*) nframes
read(10,*) ndefects
read(10,*) nvacs
read(10,*) boxlen(1), boxlen(2), boxlen(3)
read(10,*) spec1
read(10,*) spec2
read(10,*) rdfmax
read(10,*) nbins
close(10)

ncen = ndefects + nvacs
allocate(centre(ncen))
allocate(rdfout(0:nbins))
allocate(rdffunc(0:nbins))

cellvol = boxlen(1) * boxlen(2) * boxlen(3)

outfile = "def_"//trim(spec1)//"_"//trim(spec2)//"_rdf.dat"

open(11, file=deffile, status="old")
open(12, file=vacfile, status="old")

norm1 = 0
norm2 = 0
hboxlen = boxlen/2
rdfout = 0

do i=1, nframes
    do j=1, ndefects
        read(11,*) centre(j)%ctype, centre(j)%r(:)
    enddo
    do j=ndefects+1, ncen
        read(12,*) centre(j)%ctype, centre(j)%r(:)
    enddo

    do j=1, ncen
        if (centre(j)%ctype == spec1) then
            norm1 = norm1 + 1
            do k=1, ncen
                if (centre(k)%ctype == spec2) then
                    norm2 = norm2 + 1
                    dr = centre(j)%r - centre(k)%r            
                    ! apply periodic boundary conditions
                    do axis = 1,3
                        if (dr(axis) > hboxlen(axis)) then
                            dr(axis) = dr(axis) - boxlen(axis)
                        elseif (dr(axis) < -hboxlen(axis)) then
                            dr(axis) = dr(axis) + boxlen(axis)
                        endif
                    enddo
                    dist = sqrt(sum(dr*dr))
                    if (dist <= rdfmax) then
                        ibin = int(nbins*dist/rdfmax)
                        rdfout(ibin) = rdfout(ibin) + 1
                    endif
                endif
            enddo
        endif
    enddo
enddo

rdfdr = rdfmax/float(nbins)

xint = 0.0d0
do i=1, nbins
    rdfr= float(i)*rdfdr
    xnorm = cellvol/(float(norm2) * 12.56637 * rdfr*rdfr*rdfdr)
    rdffunc(i) = float(rdfout(i))*xnorm
    xint = xint + float(rdfout(i))*xnorm*rdfdr
enddo

rdffunc = rdffunc / xint

open(15, file=outfile)
do i=1, nbins
    rdfr= float(i)*rdfdr
    write(15,*) rdfr, rdffunc(i)
enddo
close(15)
end
