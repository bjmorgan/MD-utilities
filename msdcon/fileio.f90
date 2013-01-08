module fileio

    implicit none

    type iofile
        character(len=20) :: filename
        integer :: unit
    contains
        procedure :: init => file_init
        procedure :: fopen => file_open
    end type iofile

contains

subroutine file_open( this, act )     
    
    class(iofile), intent(inout) :: this
    character, intent(in) :: act
    character(len=5) :: action_str
    integer :: ios
    
    select case (act)
        case ('r')
            action_str = 'read'
        case ('w')
            action_str = 'write'
        case default
            ! warning. action_str will be undefined!
    end select

    open(newunit=this%unit, file=this%filename, iostat=ios, action=action_str)
    if ( ios /= 0 ) then
        print *, "Error opening file " // this%filename 
        stop 
    endif
    
end subroutine file_open
    
subroutine file_init( this, filename )
    
    class(iofile) :: this
    character(len=*), intent(in) :: filename

    this%filename = filename
    
end subroutine file_init

end module fileio
