program framedisp

! Reads in disp.out, and outputs the total displacement for each species, per frame.

    use fileio

    implicit none
    
    character(len=200) :: filein

    type pos_type
        double precision, dimension(3) :: dr = 0.0d0
        double precision, dimension(3) :: r 
    end type pos_type

    type ion_type
        type(pos_type), allocatable, dimension(:) :: time
        double precision, dimension(3) :: oldr 
        double precision, dimension(3) :: total_disp
        logical :: mobile
        integer :: static_count ! the number of frames since an ion has been mobile
        logical :: monitor
    end type ion_type

    type species_type
        type(ion_type), allocatable, dimension(:) :: ion
        integer :: number_of_ions     ! number of ions
        logical :: count_displacement ! Are we considering this species for possible diffusion events?
        integer :: nmobile ! used to count the number of mobile ions at each frame
        real :: z ! species charge
    end type species_type

    type(species_type), target, allocatable, dimension(:) :: species
    type(species_type), pointer :: p_species
    type(ion_type), pointer :: p_ion
    type(pos_type), pointer :: p_time

    type(iofile) :: inptfile, dispfile, posfile, outfile, outfile_long, dispdistfile, monitorfile, chgdispfile

    double precision :: msd_tot, msd_mobile_tot
    double precision :: dummy
    double precision :: max_disp
    double precision :: this_disp_sq
    double precision :: disp_cut ! displacement cutoff for identifying "mobile" ions
    double precision, dimension(3) :: chg_disp, chg_disp_mobile, chg_disp_tot, chg_disp_mobile_tot, this_disp, this_disp_mobile

    integer :: nspecies, nsteps, istep, nsp
    integer :: nindex, mindex
    integer :: i, ios, err, nion
    integer :: i_dummy
    integer :: disp_length
    integer :: nbins, this_bin
    integer, allocatable, dimension(:) :: disp_bin
    
    integer :: ntotal_ions = 0
    integer :: ncounted_ions = 0

    double precision, external :: length_sq

    integer, parameter :: z = 3

    ! set up input/output files

    call file_init( inptfile, "frame_disp.inpt" )
    call file_init( dispfile, "disp.out" )
    call file_init( outfile, "frame_disp.dat" )
    call file_init( outfile_long, "frame_disp.out" )
    call file_init( dispdistfile, "disp_dist.dat" )
    call file_init( posfile, "poscart.out" )
    call file_init( monitorfile, "ion_disp.dat" )
    call file_init( chgdispfile, "chg_disp.dat" )

    call file_open( inptfile, "read" )
    call file_open( dispfile, "read" )
    call file_open( outfile, "write" )
    call file_open( outfile_long, "write" )
    call file_open( dispdistfile, "write" )
    call file_open( posfile, "read" )
    call file_open( monitorfile, "write" )
    call file_open( chgdispfile, "write" )

    ! read calculation parameters from parameters file, inptfile

    read(inptfile%iounit,*) nsteps
    read(inptfile%iounit,*) disp_length
    read(inptfile%iounit,*) max_disp
    read(inptfile%iounit,*) disp_cut
    read(inptfile%iounit,*) nbins
    allocate( disp_bin( nbins ), stat=err )
    if ( err /= 0 ) then
        print *,"Cannot allocate memory for disp_bins(:)" ; stop
    endif
    read(inptfile%iounit,*) nspecies
    allocate( species( nspecies ), stat=err ) 
    if ( err /= 0 ) then
        print *,"Cannot allocate memory for species(:)" ; stop
    endif

    do nsp=1, nspecies ! read number_of_ions and count_displacement for each species
        
        p_species => species( nsp ) ! p_species points to species( nsp )
        
        read ( inptfile%iounit, * ) p_species%number_of_ions ! read number of ions for each species
        ntotal_ions = ntotal_ions + p_species%number_of_ions ! add to the total number of ions in the simulation
        allocate( p_species%ion( p_species%number_of_ions ), stat=err )
        if ( err /= 0 ) then
            print *,"Cannot allocate memory for species ", nsp ; stop
        endif
        
        do nion = 1, p_species%number_of_ions
            allocate( p_species%ion(nion)%time(disp_length), stat=err )
            if (err /= 0) then
                print *,"Cannot allocate memory for ion ", nion, " in species ", nsp ; stop
            endif
        enddo
        
        read( inptfile%iounit, * ) p_species%count_displacement ! consider this species when summing the total displacement? ( .true. || .false. )
        if ( p_species%count_displacement ) ncounted_ions = ncounted_ions + p_species%number_of_ions ! add to the total number of contributing ions
        read( inptfile%iounit, * ) p_species%z ! read ion charge for this species
        
        nullify (p_species)
        
    enddo

    close( inptfile%iounit ) ! end parameter readin from inptfile
    
    disp_bin = 0
    chg_disp_tot = 0.0d0
    chg_disp_mobile_tot = 0.0d0
    
    do nsp = 1, nspecies
        species( nsp )%ion(:)%mobile = .false.
        species( nsp )%ion(:)%static_count = 0
    enddo
    
    print *,"Analysing", nsteps, "steps"
    do i=1, ntotal_ions+1 
        read( dispfile%iounit, * ) i_dummy ! read past header information in <disp.out>
    enddo
    
    do istep = 1, nsteps

        nindex = mod(istep-1, disp_length) + 1 ! index for time(:) arrays [displacements and positions]

        do nsp = 1, nspecies
            p_species => species( nsp ) ! p_species points to species(nsp)
            
            p_species%nmobile = 0
            
            do nion = 1, p_species%number_of_ions               
                p_ion => p_species%ion( nion )
                p_time => p_ion%time( nindex )
                
                p_ion%mobile = .false.
                p_ion%oldr = p_time%r
                read( dispfile%iounit, * ) p_time%dr ! read displacement for each ion
                read( posfile%iounit,* ) p_time%r   ! read positions for each ion
                
                nullify( p_ion )
                nullify( p_time )
                
            enddo
            
            nullify( p_species )
            
        enddo

        if (istep < disp_length) cycle
        
        msd_tot = 0.0d0
        msd_mobile_tot = 0.0d0
        chg_disp = 0.0d0
        chg_disp_mobile = 0.0d0

        do nsp = 1, nspecies
            
            p_species => species( nsp )
            this_disp = 0.0d0
            this_disp_mobile = 0.0d0
            
            do nion=1, p_species%number_of_ions
                p_ion => p_species%ion( nion )
                
                p_ion%total_disp(:) = 0.0d0
                do i=1, disp_length
                    p_ion%total_disp(:) = p_ion%total_disp(:) + p_ion%time(i)%dr(:)
                enddo

                this_disp(:) = this_disp(:) + p_ion%total_disp(:)
                this_disp_sq = length_sq( p_ion%total_disp(:) )
 
                chg_disp_tot(:) = chg_disp_tot(:) + p_ion%time(nindex)%dr(:) * p_species%z
                
                if ( sqrt( this_disp_sq ) > max_disp ) then
                    print *, "Warning, displacement exceeds max_disp: ", sqrt( this_disp_sq ), " > ", max_disp
                    stop
                else if ( sqrt( this_disp_sq ) > disp_cut ) then ! ion is mobile
                    p_ion%mobile = .true.
                    p_ion%static_count = 0
                    p_species%nmobile = p_species%nmobile + 1
                    this_disp_mobile(:) = this_disp_mobile(:) + p_ion%total_disp(:)

                    if (  p_ion%static_count >= disp_length ) then
  ! at least (disp_length) steps since this ion counted as being mobile, so we need to count all (disp_length) steps
                        chg_disp_mobile_tot(:) = chg_disp_mobile_tot(:) + p_ion%total_disp(:) * p_species%z
                    else    
                        chg_disp_mobile_tot(:) = chg_disp_mobile_tot(:) + p_ion%time(nindex)%dr(:) * p_species%z
                    endif
                else
                    p_ion%static_count = p_ion%static_count + 1    
                endif

                if ( p_species%count_displacement == .true. ) then                                    
                    msd_tot = msd_tot + this_disp_sq
                    this_bin = int( ( sqrt( this_disp_sq )/max_disp )*nbins ) + 1
                    disp_bin( this_bin ) = disp_bin( this_bin ) + 1
                    if ( p_ion%mobile ) then
                        msd_mobile_tot = msd_mobile_tot + this_disp_sq
                    endif
                endif
                
                nullify( p_ion )
                
                chg_disp(:) = chg_disp(:) + this_disp(:) * p_species%z
                chg_disp_mobile(:) = chg_disp_mobile(:) + this_disp_mobile(:) * p_species%z
                
            enddo
            
            nullify( p_species )
            
        enddo
        
        ! frame output
        write(outfile_long%iounit,'(a5,i5,a9,i3,a8)',advance="no") "Step ", istep, " Configs ", disp_length+1, " Mobile "

        do nsp = 1, nspecies
            
            p_species => species( nsp )
            
            write(outfile_long%iounit,'(i5)',advance="no") p_species%nmobile
            do nion=1, species(nsp)%number_of_ions
                
                p_ion => p_species%ion( nion )
                
                if (p_ion%monitor == .true.) then
                    write(monitorfile%iounit,*) length_sq(p_ion%total_disp)
                endif
                
                nullify( p_ion )
                
            enddo
            
            nullify( p_species )
            
        enddo
 
        write(outfile_long%iounit,*)
 
        do i=nindex, nindex-disp_length+1, -1
            mindex=i-disp_length*int((i-disp_length)/disp_length)
            do nsp=1, nspecies
                do nion=1, species(nsp)%number_of_ions
                    
                    p_ion => species(nsp)%ion(nion)
                    
                    if (p_ion%mobile == .true.) then
                         write(outfile_long%iounit,'(3(f9.5,1x),1x,i5)') p_ion%time(mindex)%r(1:3), nion
                    endif
                    
                    nullify( p_ion )
                    
                enddo
            enddo    
        enddo

        do nsp=1, nspecies
            do nion=1, species(nsp)%number_of_ions
                
                p_ion => species(nsp)%ion(nion)
                
                if (p_ion%mobile == .true.) then
                    write(outfile_long%iounit,'(3(f9.5,1x),1x,i5)') p_ion%oldr(1:3), nion
                endif
                
                nullify( p_ion )
                
            enddo
        enddo
        
        write(outfile%iounit, '(i5,1x,2(e10.3,1x))' ) istep, msd_tot/float(ncounted_ions), msd_mobile_tot/float(ncounted_ions) 

        write( chgdispfile%iounit, '(i5,1x,5(e10.3,1x))' ) istep, length_sq( chg_disp ), length_sq( chg_disp_mobile ), length_sq( chg_disp_tot ), length_sq( chg_disp_mobile_tot ), sign( chg_disp_mobile_tot(z)**2, chg_disp_mobile_tot(z) )
    enddo
    
    do this_bin=1, nbins
        write(dispdistfile%iounit,*) max_disp*this_bin/nbins, disp_bin(this_bin)/float(ncounted_ions*(nsteps-disp_length))/max_disp*nbins
    enddo

end program framedisp

double precision function length_sq( r )

    implicit none
    double precision, dimension(3), intent(in) :: r
    
    length_sq = sum(r*r)
    
end function length_sq
