module species_class

    use ion_class, only : ion_type

    implicit none

    integer, parameter :: long = selected_real_kind(9, 99)
    integer, parameter :: xyz = 3

    type species_type
        real(long) :: suma
        real(long), dimension( xyz ) :: sum, sumabs
        real(long) :: z
        integer :: num
    
        real(long), allocatable, dimension(:) :: msd_av
        real(long), allocatable, dimension(:,:) :: msd, msdabs, tot 
    
        type( ion_type ), allocatable, dimension(:) :: ion
		
	contains
		procedure :: init => species_init
		
    end type species_type

contains
	
	subroutine species_init(this, nmsdlength)
	implicit none
	class(species_type) :: this
	integer, intent(in) :: nmsdlength	
	integer :: err
	integer :: xyz = 3
	integer :: j
	
    allocate( this%ion( this%num ) )
    allocate( this%msd_av( 0:nmsdlength ) )
    allocate( this%msd( 0:nmsdlength, xyz ) )
    allocate( this%msdabs( 0:nmsdlength, xyz ) )
    allocate( this%tot( 0:nmsdlength, xyz ) )
	do j=1, this%num
		call this%ion(j)%init( nmsdlength )
	enddo
	
	end subroutine species_init
end module species_class