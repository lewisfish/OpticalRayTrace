module stackMod

    use vector_class, only : vector

    implicit none
    
    type :: stack
        type(vector), allocatable :: data(:) ! stack
        integer                   :: size = 0
        contains
            procedure :: pop   => pop_fn
            procedure :: push  => push_sub
            procedure :: peek  => peek_fn
            procedure :: empty => empty_fn
            procedure :: zero  => zero_sub
            procedure :: write => write_sub
            procedure :: write_empty => write_empty_sub
    end type stack

    integer, parameter :: block_size = 4

    contains

    subroutine write_empty_sub(this, u)

        implicit none
        
        class(stack) :: this
        integer, intent(IN) :: u

        call this%zero()
        write(u,*)" "
        write(u,*)" "
        write(u,*)" "

    end subroutine write_empty_sub

    subroutine write_sub(this, u)

        implicit none
        
        class(stack) :: this
        integer, intent(IN) :: u

        do while(.not. this%empty())
            write(u,"(3(F10.7,1x))")this%pop()
        end do
        write(u,*)" "
        write(u,*)" "
        write(u,*)" "

    end subroutine write_sub

    subroutine zero_sub(this)
        
        implicit none
        
        class(stack) :: this

        type(vector) :: tmp

        do while(.not. this%empty())
            tmp = this%pop()
        end do

    end subroutine zero_sub

    type(vector) function pop_fn(this)
    ! pop top enrty off stack
        implicit none

        class(stack) :: this

        if(this%size == 0 .or. .not. allocated(this%data))then
            !if nothing in stack send back garbage data
            pop_fn = vector(-999d0, -999.d0, -999.d0)
            return
        end if
        pop_fn = this%data(this%size)
        this%size = this%size - 1

    end function pop_fn


    type(vector) function peek_fn(this)

        implicit none

        class(stack) :: this

        if(this%size == 0 .or. .not. allocated(this%data))then
            peek_fn = vector(-999d0, -999.d0, -999.d0)
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

        type(vector), intent(IN)  :: pt
        type(vector), allocatable :: tmp(:)

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