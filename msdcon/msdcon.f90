program msdcalc

! November 20, 2011
! Rewritten to use dynamic memory and user-defined data structures

use fileio

implicit none 

integer, parameter :: xyz=3
integer, parameter :: x=1, y=2, z=3

double precision, allocatable, dimension(:,:,:) :: msdcoll
double precision, allocatable, dimension(:,:,:,:) :: ds 
integer, allocatable, dimension( : ) :: normtot

type ion_type
    double precision, dimension( xyz ) :: disp
    double precision, allocatable, dimension( :,: ) :: dispstore
    double precision, allocatable, dimension( :,: ) :: msd
    double precision, allocatable, dimension( :,: ) :: dispabs
    integer, allocatable, dimension( : ) :: norm    
end type ion_type

type species_type
    double precision :: suma
    double precision, dimension( xyz ) :: sum, sumabs
    double precision :: z
    integer :: num
    
    double precision, allocatable, dimension(:) :: msd_av
    double precision, allocatable, dimension(:,:) :: msd, msdabs, tot 
    
    type( ion_type ), allocatable, dimension(:) :: ion
end type species_type

type(species_type), allocatable, target, dimension(:) :: species
type(species_type), pointer :: p_spec, p_spec1, p_spec2
type(ion_type), pointer :: p_ion

integer :: nmsdlength, nmsdcalltime
integer :: mcorrtime = 1
integer :: mcorrmax

integer :: num
integer :: nstep = 0
integer :: nrun
double precision :: dtime
integer :: ncation1, ncation2, ncation3, nanion
integer :: nspecies
double precision :: rn, time, work, wnernst

character( len = 200 ) :: filein
logical :: overflow = .false.
logical :: restart = .false.
logical :: endrun = .false.
logical :: readfrominpt_log = .false.

integer :: i, j, k, l, i1, n, ns, ns2, i2, ispec, nt, ipt, ip, ip1, ip2, n1, n2
integer :: ipoint
double precision :: dum1, dum2
double precision, parameter :: time_au_to_ps = 2.418884d-5

real :: bin
integer :: nbin
integer :: err

external :: check_allocation

! iofiles
type(iofile) :: inptfile, dispfile, rstfile, nernstfile, workfile
type(iofile), allocatable, dimension(:) :: msdabsfile, msdfile

character(len=20) filename

!Set up io
call file_init( inptfile, "msd.inpt" )
call file_init( dispfile, "disp.out" )
call file_init( rstfile, "msdrst.dat" )
call file_init( workfile, "work.dat" )
call file_init( nernstfile, "nernst.dat" )

call file_open( inptfile, "read" )
call file_open( dispfile, "read" )
call file_open( workfile, "write" )
call file_open( nernstfile, "write" )

!==================================================
! Get parameters
!==================================================

open( inptfile%unit, file='msd.inpt' )
read( inptfile%unit, '(a)' ) filein
read( inptfile%unit, * ) num
read( inptfile%unit, * ) nmsdlength
read( inptfile%unit, * ) nspecies
allocate( species(nspecies), stat = err ); call check_allocation( err, "species(:)")  
do i=1, nspecies
    read( inptfile%unit,* ) species(i)%num
    read( inptfile%unit,* ) species(i)%z
enddo
read( inptfile%unit, * ) restart

!==================================================
! Allocate remaining arrays
!==================================================

allocate( normtot( 0:nmsdlength ), stat = err); call check_allocation( err, "normtot(:)" )
  
do i=1, nspecies
    p_spec => species(i)
    allocate( p_spec%ion( p_spec%num ), stat = err ); call check_allocation( err, "species(:)%ion(:)" )
    allocate( p_spec%msd_av( 0:nmsdlength ), stat = err ); call check_allocation( err, "species(:)%msd_av(:)" )
    allocate( p_spec%msd( 0:nmsdlength, xyz ), stat = err ); call check_allocation( err, "species(:)%msd(:,:)" )
    allocate( p_spec%msdabs( 0:nmsdlength, xyz ), stat = err ); call check_allocation( err, "species(:)%msdabs(:,:)" )
    allocate( p_spec%tot( 0:nmsdlength, xyz ), stat = err ); call check_allocation( err, "species(:)%tot(:,:)" )
    do j=1, p_spec%num
        p_ion => p_spec%ion(j)
        allocate( p_ion%dispstore( nmsdlength, xyz ), stat = err ); call check_allocation( err, "species(:)%ion(:)%dispstore(:)")
        allocate( p_ion%msd( 0:nmsdlength, xyz ), stat = err ); call check_allocation( err, "species(:)%ion(:)%msd(:)")   
        allocate( p_ion%dispabs( 0:nmsdlength, xyz ), stat = err ); call check_allocation( err, "species(:)%ion(:)%dispabs(:)")
        allocate( p_ion%norm( 0:nmsdlength ), stat = err ); call check_allocation( err, "species(:)%ion(:)%norm(:)")
        nullify( p_ion )
    enddo
    nullify( p_spec )
enddo

allocate( ds( 0:nmsdlength, nspecies, nspecies, xyz ), stat = err ); call check_allocation( err, "ds(:,:,:,:)" ) 
allocate( msdcoll( 0:nmsdlength, nspecies, nspecies ), stat = err ); call check_allocation( err, "msdcoll(:,:,:)")

allocate( msdabsfile(nspecies), stat = err ) ; call check_allocation( err, "msdabsfile(:)" )
allocate( msdfile(nspecies), stat = err ) ; call check_allocation( err, "msdfile(:)" )
  
!==================================================
! Zero arrays
!==================================================

ds = 0.0d0
ds = 0.0d0
do i=1, nspecies
    do j=1, species(i)%num
        p_ion => species(i)%ion(j)
        p_ion%msd = 0.0d0     
        p_ion%dispabs = 0.0d0
        p_ion%norm = 0
        nullify( p_ion )
    enddo
enddo

!================================================
! Open output files
!================================================


do ns=1, nspecies

    call file_init( msdabsfile(ns), 'msdabs'//char(ns+48)//'.dat' )
    call file_init( msdfile(ns), 'msd'//char(ns+48)//'.dat' )
    call file_open( msdabsfile(ns), "write" )
    call file_open( msdfile(ns), "write")

    do ns2=ns,nspecies
        filename='msdcollect'//char(ns+48)//char(ns2+48)//'.dat'
        open (70+(ns*nspecies)+ns2,file=filename)
    enddo
enddo



if (restart) then

    call file_open( rstfile, "read" )

    do i=0,nmsdlength
        do j=1, nspecies
            do n=1, species(j)%num
                p_ion => species(j)%ion(n)
                do l=x, z
                    read( rstfile%unit, * ) p_ion%msd( i, l )
                enddo
                read( rstfile%unit, * ) p_ion%norm( i ) 
                nullify(p_ion)
            enddo
        enddo
        do i1=1, nspecies
            do i2=1, nspecies
                do l=x, z
                    read( rstfile%unit, * ) ds(i,i1,i2,l)
                enddo
            enddo
        enddo
        read( rstfile%unit, * ) normtot(i)
    enddo
    close( rstfile%unit )

endif

read( inptfile%unit, * ) readfrominpt_log

if (readfrominpt_log) then
    read (inptfile%unit,*) nmsdcalltime
    read (inptfile%unit,*) dtime
    read (inptfile%unit,*) nrun

    read( dispfile%unit, * ) nbin, bin, nbin
    do i=1,num
        read( dispfile%unit, * ) nbin
    enddo
else
    read ( dispfile%unit, * ) nmsdcalltime, dtime, nrun
    do i=1,num
        read( dispfile%unit, * ) dum1
    enddo
endif

close( inptfile%unit )

nstep = 0
mcorrtime = 1

do while ( endrun == .false. )

    nstep = nstep + nmsdcalltime

    if (nstep.ge.nrun) endrun = .true.

    if (overflow) then
        mcorrmax = nmsdlength
    else
        mcorrmax = mcorrtime
    endif
 
    print *,'entering msdcalc',nstep
       
    do i=1,nspecies
        p_spec => species(i)
        do j=1, p_spec%num
            p_ion => p_spec%ion(j)
            p_ion%dispstore(mcorrtime,:) = 0.0d0
            read( dispfile%unit, * ) ( p_ion%disp(l), l=x, z )
            
            do k=1, mcorrmax
                p_ion%dispstore(k,:) = p_ion%dispstore(k,:) + p_ion%disp(:)
            enddo
            
            nullify( p_ion )
        enddo
        nullify( p_spec )
    enddo

     
    do j=1, mcorrmax
        
        if ( j <= mcorrtime ) then ! this should be moved to a separate function, and computed numerically
            nt = mcorrtime - j
        else
            nt = mcorrtime - j + nmsdlength
        endif
        
        normtot(nt) = normtot(nt) + 1
        
        do i=1, nspecies
            p_spec => species(i)
            p_spec%tot(j,:) = 0.0d0
            do ip=1, p_spec%num
                p_ion => p_spec%ion(ip)
              
                p_spec%tot(j,:) = p_spec%tot(j,:) + p_ion%dispstore(j,:)
                p_ion%msd(nt,:) = p_ion%msd(nt,:) + p_ion%dispstore(j,:)**2
                p_ion%dispabs(nt,:) = p_ion%dispabs(nt,:) + sign(p_ion%dispstore(j,:)**2, p_ion%dispstore(j,:))
                p_ion%norm(nt) = p_ion%norm(nt) + 1
                
                nullify( p_ion )
            enddo           
            nullify( p_spec )
        enddo

        do i1=1,nspecies
            do i2=i1,nspecies
                ds(nt,i1,i2,:) = ds(nt,i1,i2,:) + species(i1)%tot(j,:) * species(i2)%tot(j,:)
            enddo
        enddo

    enddo

    !========================================================================
    ! Update array counters
    !========================================================================
    if (mod(float(mcorrtime),float(nmsdlength)).eq.0) then
        overflow=.true.
    endif

    mcorrtime = int( mod( float(mcorrtime), float(nmsdlength) ) )
    mcorrtime = mcorrtime + 1

enddo

close( dispfile%unit )

!=========================================================
! Write out restart file
!=========================================================

do i=1, nspecies
    p_spec => species(i)
    p_spec%msd_av(:) = 0.0d0
    nullify( p_spec )
enddo
    
msdcoll = 0.0d0
  
!   Writing out values for restart file

call file_open( rstfile, "write" )

do i=0,nmsdlength
    do j=1, nspecies
        do n=1,species(j)%num
            p_ion => species(j)%ion(n)
            write( rstfile%unit, * ) p_ion%msd(i,x)
            write( rstfile%unit, * ) p_ion%msd(i,y)
            write( rstfile%unit, * ) p_ion%msd(i,z)
            write( rstfile%unit, * ) p_ion%norm(i)
            nullify( p_ion )
        enddo
    enddo
    do i1=1,nspecies
        do i2=1,nspecies
            write( rstfile%unit, * ) ds(i,i1,i2,x),"xds"
            write( rstfile%unit, * ) ds(i,i1,i2,y),"yds"
            write( rstfile%unit, * ) ds(i,i1,i2,z),"zds"
        enddo
    enddo

    write(rstfile%unit,*) normtot(i)
enddo

close(rstfile%unit)

!======================================================================
! Average over components and number of molecules
!======================================================================

do i=1, nspecies
    do k=1, species(i)%num
        p_ion => species(i)%ion(k)
        do j=0, nmsdlength
            if ( p_ion%norm(j) > 0) then
                p_ion%msd(j,:) = p_ion%msd(j,:) / float( p_ion%norm(j) )
                p_ion%dispabs(j,:) = p_ion%dispabs(j,:) / float( p_ion%norm(j))
            endif
        enddo
        nullify( p_ion )
    enddo
enddo

do j=0,nmsdlength
    do n1=1,nspecies
        do n2=1,nspecies
            if (normtot(j) > 0) then
                rn=sqrt(float(species(n1)%num*species(n2)%num))*float(normtot(j))
                ds(j,n1,n2,:) = ds(j,n1,n2,:)/rn
                msdcoll(j,n1,n2) = sum(ds(j,n1,n2,:))
            endif
        enddo
    enddo
enddo

do i=0,nmsdlength
    do j=1,nspecies
        p_spec => species(j)
        p_spec%suma = 0.0d0
        p_spec%sum = 0.0d0
        p_spec%sumabs = 0.0d0
        do k=1, p_spec%num
            p_ion => species(j)%ion(k)
            p_spec%suma = p_spec%suma + sum(p_ion%msd(i,:))
            p_spec%sum(:) = p_spec%sum(:) + p_ion%msd(i,:)
            p_spec%sumabs(:) = p_spec%sumabs(:) + p_ion%dispabs(i,:)
            nullify(p_ion)
        enddo
        nullify( p_spec )
    enddo

    do ip=1,nspecies
        p_spec => species(ip)
        if (p_spec%num > 0) then
            p_spec%msd_av(i) = p_spec%suma / float( p_spec%num )
            p_spec%msd(i,:) = p_spec%sum(:) / float( p_spec%num )
            p_spec%msdabs(i,:) = p_spec%sumabs(:) / float( p_spec%num )
        endif
        nullify( p_spec )
    enddo
enddo

dum1=0.0
dum2=0.0

do ip=1,nspecies
    write ( msdfile(ip)%unit, * ) dum1,dum1,dum1,dum1,dum1
enddo

write( workfile%unit,* )dum1,dum2
write( nernstfile%unit,* )dum1,dum2

do ip1=1,nspecies
    if (species(ip1)%num > 0) then
        do ip2 = ip1, nspecies
            write(70+(ip1*nspecies)+ip2,*) dum1,dum2
        enddo
    endif
enddo

do i=0,nmsdlength-1
    time=(dble(i)+1)*dble(nmsdcalltime)*dtime*2.418d-5
    do ip=1, nspecies
        p_spec => species(ip)
        write( msdfile(ip)%unit, '(5(e12.5,1x))') time, p_spec%msd_av(i), p_spec%msd(i,x), p_spec%msd(i,y), p_spec%msd(i,z)
        write( msdabsfile(ip)%unit, '(4(e12.5,1x))') time, p_spec%msdabs( i, x ), p_spec%msdabs( i, y ), p_spec%msdabs( i, z )
        nullify( p_spec )
    enddo

    work=0.0d0
    wnernst=0.0d0

    do ip1=1,nspecies
        p_spec1 => species(ip1)
        if (species(ip1)%num > 0) then
            do ip2=ip1,nspecies
                p_spec2 => species(ip2)
                write(70+(ip1*nspecies)+ip2,*) time, &
                               msdcoll(i,ip1,ip2)
                msdcoll(i,ip1,ip2)=msdcoll(i,ip1,ip2)* &
                                   p_spec1%z*p_spec2%z* &
                sqrt(float(p_spec1%num*p_spec2%num))/float(num)
                if (ip1.ne.ip2) then 
                    work=work+msdcoll(i,ip1,ip2)
                endif
                work=work+msdcoll(i,ip1,ip2)
                nullify( p_spec2 )
            enddo
        endif

        wnernst=wnernst+species(ip1)%msd_av(i)*p_spec1%z*p_spec1%z*float(p_spec1%num)/float(num)
        
        nullify( p_spec1 )

    enddo
    write( workfile%unit,* ) time,work
    write( nernstfile%unit,* ) time,wnernst
enddo
  
write (6,*)
write (6,*) '*** Mean squared displacements written out. ***'
write (6,*)

stop
end

subroutine check_allocation( err, string )
    implicit none
    integer, intent(in) :: err
    character(len=*) :: string
    
    if ( err /= 0 ) then
        print *,"Error allocating memory to "//string
        stop
    endif
end subroutine check_allocation