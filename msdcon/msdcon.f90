program msdcalc

use fileio
use ion_class, only : ion_type
use species_class, only : species_type

implicit none 
integer, parameter :: long = selected_real_kind(9, 99)
integer, parameter :: xyz=3
integer, parameter :: x=1, y=2, z=3

real(long), allocatable, dimension(:,:,:) :: msdcoll
real(long), allocatable, dimension(:,:,:,:) :: ds 
integer, allocatable, dimension( : ) :: normtot

type(species_type), allocatable, target, dimension(:) :: species

integer :: nmsdlength, nmsdcalltime
integer :: mcorrtime = 1
integer :: mcorrmax

integer :: num = 0
integer :: nstep = 0
integer :: nrun
integer :: nspecies
real(long) :: dtime, rn, time, work, wnernst

character( len = 200 ) :: filein
logical :: overflow = .false.
logical :: restart = .false.
logical :: endrun = .false.
logical :: readfrominpt_log = .false.

integer :: i, j, k, l, i1, n, i2, ispec, nt, ipt, ip, ip1, ip2, n1, n2
integer :: ipoint, nbin, err
real(long) :: dum1, dum2
real(long), parameter :: time_au_to_ps = 2.418884d-5
real :: bin

integer :: get_nt
external :: check_allocation, get_nt

! iofiles
type(iofile) :: inptfile
type(iofile) :: dispfile, rstfile, nernstfile, workfile
type(iofile), allocatable, dimension(:) :: msdabsfile, msdfile

character(len=20) filename

!Set up io
call inptfile%init( "msd.inpt" )
call rstfile%init( "msdrst.dat" )
call workfile%init( "work.dat" )
call nernstfile%init( "nernst.dat" )

call inptfile%fopen('r')
call workfile%fopen('w')
call nernstfile%fopen('w')

!==================================================
! Get parameters
!==================================================

read( inptfile%unit, '(a)' ) filein 
call dispfile%init( filein ) 
call dispfile%fopen('r')

read( inptfile%unit, * ) nmsdlength
read( inptfile%unit, * ) nspecies
allocate( species(nspecies), stat = err )
call check_allocation( err, "species(:)")  
do i=1, nspecies
    read( inptfile%unit,* ) species(i)%num
    if ( species(i)%num > 0 ) then 
        num = num + species(i)%num
    else ! if this is negative it is likely a misplaced charge
        print *,"Error in input"
        print *
        print *,"Example input file:"
        print *
        call example_input
        stop
    endif
    read( inptfile%unit,* ) species(i)%z
enddo
read( inptfile%unit, * ) restart

!==================================================
! Allocate remaining arrays
!==================================================

allocate( normtot( 0:nmsdlength ), stat = err); call check_allocation( err, "normtot(:)" )
  
do i=1, nspecies
    call species(i)%init( nmsdlength )
enddo

allocate( ds( 0:nmsdlength, nspecies, nspecies, xyz ), stat = err ); call check_allocation( err, "ds(:,:,:,:)" ) 
allocate( msdcoll( 0:nmsdlength, nspecies, nspecies ), stat = err ); call check_allocation( err, "msdcoll(:,:,:)")
allocate( msdabsfile(nspecies), stat = err ) ; call check_allocation( err, "msdabsfile(:)" )
allocate( msdfile(nspecies), stat = err ) ; call check_allocation( err, "msdfile(:)" )
  
ds = 0.0d0

!================================================
! Open output files
!================================================


do ip1=1, nspecies

    call msdabsfile(ip1)%init( 'msdabs'//char(ip1+48)//'.dat' )
    call msdabsfile(ip1)%fopen( 'w' )

    call msdfile(ip1)%init( 'msd'//char(ip1+48)//'.dat' )
    call msdfile(ip1)%fopen( 'w' )

    do ip2=ip1,nspecies
        filename='msdcollect'//char(ip1+48)//char(ip2+48)//'.dat'
        open (70+(ip1*nspecies)+ip2,file=filename)
    enddo
enddo

if (restart) then

   write(6,*) "Reading restart data"
   call rstfile%fopen( 'r' )

    do i=0,nmsdlength
        do j=1, nspecies
            do n=1, species(j)%num
                associate( p_ion => species(j)%ion(n) )
                    do l=x, z
                        read( rstfile%unit, * ) p_ion%msd( i, l )
                    enddo
                    read( rstfile%unit, * ) p_ion%norm( i ) 
                end associate
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
    do i=1, num
        read( dispfile%unit, * ) nbin
    enddo
else
    read ( dispfile%unit, * ) nmsdcalltime, dtime, nrun
    do i=1, num
        read( dispfile%unit, * ) dum1
    enddo
endif

close( inptfile%unit )

nstep = 0
mcorrtime = 1

do while ( endrun .eqv. .false. )

    nstep = nstep + nmsdcalltime

    if (nstep.ge.nrun) endrun = .true.

    if (overflow) then
        mcorrmax = nmsdlength
    else
        mcorrmax = mcorrtime
    endif
 
    print *,'entering msdcalc',nstep
       
    do i=1,nspecies
        associate( p_spec => species(i) )
        do j=1, p_spec%num
            associate( p_ion => p_spec%ion(j) )
                p_ion%dispstore(mcorrtime,:) = 0.0d0
                read( dispfile%unit, * ) ( p_ion%disp(l), l=x, z )
            
                do k=1, mcorrmax
                    p_ion%dispstore(k,:) = p_ion%dispstore(k,:) + p_ion%disp(:)
                enddo
            
            end associate
        enddo
        end associate
    enddo
    
    do j=1, mcorrmax
        
        nt = get_nt(mcorrtime, j, nmsdlength)
        normtot(nt) = normtot(nt) + 1
        
        do i=1, nspecies
            associate( p_spec => species(i) )
                p_spec%tot(j,:) = 0.0d0
                do ip=1, p_spec%num
                    associate( p_ion => p_spec%ion(ip) )              
                        p_spec%tot(j,:) = p_spec%tot(j,:) + p_ion%dispstore(j,:)
                        p_ion%msd(nt,:) = p_ion%msd(nt,:) + p_ion%dispstore(j,:)**2
                        p_ion%dispabs(nt,:) = p_ion%dispabs(nt,:) + sign(p_ion%dispstore(j,:)**2, p_ion%dispstore(j,:))
                        p_ion%norm(nt) = p_ion%norm(nt) + 1
                    end associate
                enddo           
            end associate
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
    if( mod( float(mcorrtime), float(nmsdlength) ) == 0 ) overflow = .true.

    mcorrtime = int( mod( float(mcorrtime), float(nmsdlength) ) ) + 1

enddo

close( dispfile%unit )

!=========================================================
! Write out restart file
!=========================================================

forall( i=1:nspecies ) species(i)%msd_av = 0.0d0    
msdcoll = 0.0d0
  
!   Writing out values for restart file

call rstfile%fopen('w')

do i=0,nmsdlength
    do j=1, nspecies
        do n=1,species(j)%num
            associate( p_ion => species(j)%ion(n) )
                write( rstfile%unit, * ) p_ion%msd(i,x)
                write( rstfile%unit, * ) p_ion%msd(i,y)
                write( rstfile%unit, * ) p_ion%msd(i,z)
                write( rstfile%unit, * ) p_ion%norm(i)
            end associate
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
        associate( p_ion => species(i)%ion(k) )
            do j=0, nmsdlength
                if ( p_ion%norm(j) > 0) then
                    p_ion%msd(j,:) = p_ion%msd(j,:) / float( p_ion%norm(j) )
                    p_ion%dispabs(j,:) = p_ion%dispabs(j,:) / float( p_ion%norm(j))
                endif
            enddo
        end associate
    enddo
enddo

do j=0, nmsdlength
    do n1=1, nspecies
        do n2=1, nspecies
            if (normtot(j) > 0) then
                rn = sqrt( float( species( n1 )%num*species( n2 )%num) ) * float( normtot( j ) )
                ds(j,n1,n2,:) = ds(j,n1,n2,:) / rn
                msdcoll(j,n1,n2) = sum(ds(j,n1,n2,:))
            endif
        enddo
    enddo
enddo

do i=0,nmsdlength
    do j=1,nspecies
        associate(p_spec => species(j))
        p_spec%suma = 0.0d0
        p_spec%sum = 0.0d0
        p_spec%sumabs = 0.0d0
        do k=1, p_spec%num
            associate( p_ion => species(j)%ion(k) )
                p_spec%suma = p_spec%suma + sum(p_ion%msd(i,:))
                p_spec%sum(:) = p_spec%sum(:) + p_ion%msd(i,:)
                p_spec%sumabs(:) = p_spec%sumabs(:) + p_ion%dispabs(i,:)
            end associate
        enddo
        end associate
    enddo

    do ip=1,nspecies
        associate(p_spec => species(ip))
            if (p_spec%num > 0) then
                p_spec%msd_av(i) = p_spec%suma / float( p_spec%num )
                p_spec%msd(i,:) = p_spec%sum(:) / float( p_spec%num )
                p_spec%msdabs(i,:) = p_spec%sumabs(:) / float( p_spec%num )
            endif
        end associate
    enddo
enddo

dum1=0.0
dum2=0.0

do ip=1,nspecies
    write ( msdfile(ip)%unit, '(5(e12.5,1x))') dum1, dum1, dum1, dum1, dum1
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
        associate(p_spec => species(ip))
            write( msdfile(ip)%unit, '(5(e12.5,1x))') time, p_spec%msd_av(i), p_spec%msd(i,x), p_spec%msd(i,y), p_spec%msd(i,z)
            write( msdabsfile(ip)%unit, '(4(e12.5,1x))') time, p_spec%msdabs( i, x ), p_spec%msdabs( i, y ), p_spec%msdabs( i, z )
        end associate
    enddo

    work=0.0d0
    wnernst=0.0d0

    do ip1=1,nspecies
        associate( p_spec1 => species(ip1) )
            if (species(ip1)%num > 0) then
                do ip2=ip1,nspecies
                    associate( p_spec2 => species(ip2) )
                        write(70+(ip1*nspecies)+ip2,*) time, msdcoll(i,ip1,ip2)
                        msdcoll(i,ip1,ip2) = msdcoll(i,ip1,ip2) * p_spec1%z*p_spec2%z &
                                          * sqrt( float(p_spec1%num*p_spec2%num) )/float(num)
                        if (ip1.ne.ip2) work=work+msdcoll(i,ip1,ip2)
                        work=work+msdcoll(i,ip1,ip2)
                    end associate
                enddo
            endif
            wnernst=wnernst+species(ip1)%msd_av(i)*p_spec1%z*p_spec1%z*float(p_spec1%num)/float(num)
        end associate
    enddo
    write( workfile%unit,* ) time,work
    write( nernstfile%unit,* ) time,wnernst
enddo
  
write (6,*)
write (6,*) '*** Mean squared displacements written out. ***'
write (6,*)

stop
end

subroutine example_input
    implicit none
    
    write(6,*) "disp.out <filein>"
    write(6,*) "150      <nmsdlength>"
    write(6,*) "2        <nspecies>"
    write(6,*) "448      <number of species 1>"
    write(6,*) "-1.0     <charge of species 1>"
    write(6,*) "448      <number of species 1>"
    write(6,*) "1.0      <charge of species 1>"
    write(6,*) ".false.  <restart?>"
    write(6,*) ".false.  <read from input?> if .true. then:"
    write(6,*) "?        <nmsdcalltime>"
    write(6,*) "?        <dtime>"
    write(6,*) "?        <nrun>"
end subroutine example_input

subroutine check_allocation( err, string )
    use iso_fortran_env, only : error_unit
    implicit none
    integer, intent(in) :: err
    character(len=*) :: string
    
    if ( err /= 0 ) then
        write( error_unit, * ) "Error allocating memory to "//string
        stop
    endif
end subroutine check_allocation

function get_nt(mcorrtime, j, nmsdlength)
    integer, intent(in) :: mcorrtime, j, nmsdlength
    integer :: get_nt

    if ( j <= mcorrtime ) then 
        get_nt = mcorrtime - j
    else
        get_nt = mcorrtime - j + nmsdlength
    endif   
end function get_nt 
