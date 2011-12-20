module fileio

    implicit none

    type iofile
        character(len=20) :: filename
        integer :: unit
    end type iofile

contains

subroutine file_open( the_file, act )     
    
    type(iofile), intent(in) :: the_file
    character(len=*), intent(in) :: act
    integer :: ios

    open(unit=the_file%unit, file=the_file%filename, iostat=ios, action=act)
    if ( ios /= 0 ) then
        print *, "Error opening file " // the_file%filename 
        stop 
    endif
    
end subroutine file_open
    
subroutine file_init( the_file, filename )
    
    type(iofile), intent(out) :: the_file
    character(len=*), intent(in) :: filename
    integer, save :: unit = 10 

    the_file%filename = filename
    the_file%unit = unit
    
    unit = unit + 1
    
end subroutine file_init

end module fileio
