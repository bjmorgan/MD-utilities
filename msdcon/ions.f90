module ion_class

    implicit none

    integer, parameter :: long = selected_real_kind(9, 99)
    integer, parameter :: xyz = 3

    type ion_type
        real(long), dimension( xyz ) :: disp
        real(long), allocatable, dimension( :,: ) :: dispstore
        real(long), allocatable, dimension( :,: ) :: msd
        real(long), allocatable, dimension( :,: ) :: dispabs
        integer, allocatable, dimension( : ) :: norm
    contains
        procedure :: init => ion_init
    end type ion_type

contains

    subroutine ion_init(this, nmsdlength)
        implicit none
        class(ion_type) :: this
        integer, intent(in) :: nmsdlength
        integer :: xyz = 3
        integer :: err

        allocate( this%dispstore( nmsdlength, xyz ) )
        allocate( this%msd( 0:nmsdlength, xyz ) )
        allocate( this%dispabs( 0:nmsdlength, xyz ) )
        allocate( this%norm( 0:nmsdlength ) )
		
		this%msd = 0.0d0
		this%dispabs = 0.0d0
		this%norm = 0.0d0
		
    end subroutine ion_init

end module ion_class
