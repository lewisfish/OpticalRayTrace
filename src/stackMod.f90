module stackMod
    implicit none
    
    type :: pointtype
        real :: x, y, z
    end type pointtype
    
    type :: stack
        type(pointtype), allocatable :: data(:) ! stack
        integer                  :: size = 0
        contains
            procedure :: pop   => pop_fn
            procedure :: push  => push_sub
            procedure :: peek  => peek_fn
            procedure :: empty => empty_fn
            procedure :: zero  => zero_sub
    end type stack

    integer, parameter :: block_size = 4

    contains
    

    subroutine zero_sub(this)
        
        implicit none
        
        class(stack) :: this

        type(pointtype) :: tmp

        do while(.not. this%empty())
            tmp = this%pop()
        end do

    end subroutine zero_sub

    type(pointtype) function pop_fn(this)
    ! pop top enrty off stack
        implicit none

        class(stack) :: this

        if(this%size == 0 .or. .not. allocated(this%data))then
            !if nothing in stack send back garbage data
            pop_fn = pointtype(-999d0, -999.d0, -999.d0)
            return
        end if
        pop_fn = this%data(this%size)
        this%size = this%size - 1

    end function pop_fn


    type(pointtype) function peek_fn(this)

        implicit none

        class(stack) :: this

        if(this%size == 0 .or. .not. allocated(this%data))then
            peek_fn = pointtype(-999d0, -999.d0, -999.d0)
            return
        end if
        peek_fn = this%data(this%size)

    end function peek_fn

    logical function empty_fn(this)

        implicit none

        class(stack) :: this

        empty_fn = (this%size == 0 .or. .not. allocated(this%data))

    end function empty_fn

    subroutine push_sub(this, pt)
    ! add pt to stack
        implicit none

        class(stack) :: this

        type(pointtype), intent(IN)  :: pt
        type(pointtype), allocatable :: tmp(:)

        if(.not. allocated(this%data))then
            ! Allocate space if not yet done
            allocate(this%data(block_size))
        elseif(this%size == size(this%data))then
            ! Grow the allocated space
            allocate(tmp(size(this%data)+block_size))
            tmp(1:this%size) = this%data
            call move_alloc(tmp,this%data)
        end if

        ! Store the data in the stack
        this%size = this%size + 1
        this%data(this%size) = pt
    end subroutine push_sub


end module stackMod