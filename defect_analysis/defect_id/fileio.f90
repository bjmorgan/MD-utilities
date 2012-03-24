module class_File

    implicit none

    type, public :: iofile
        character(len=20) :: filename
        integer :: unit
    contains
        procedure :: fopen => file_open
        procedure :: init => file_init
    end type iofile

contains

subroutine file_open( this, act )     
    
    class(iofile), intent(in) :: this
    character(len=*), intent(in) :: act
    integer :: ios

    open(unit=this%unit, file=this%filename, iostat=ios, action=act)
    if ( ios /= 0 ) then
        print *, "Error opening file " // this%filename 
        stop 
    endif
    
end subroutine file_open
    
subroutine file_init( this, filename )
    
    class(iofile), intent(out) :: this
    character(len=*), intent(in) :: filename
    integer, save :: unit = 10 

    this%filename = filename
    this%unit = unit
    
    unit = unit + 1
    
end subroutine file_init

end module class_File
