module random

    implicit none

    private
    public  :: ran2, ranu, rang

    contains


        real function ran2(seed)

            implicit none

            integer, intent(INOUT) :: seed

            call random_number(ran2)

        end function ran2


        real function ranu(a, b, iseed)

            implicit none

            real,    intent(IN)    :: a, b
            integer, intent(INOUT) :: iseed

            ranu = a + ran2(iseed) * (b - a)

        end function ranu

        subroutine rang(x, y, avg, sigma, iseed)

            implicit none

            real,    intent(IN) :: avg, sigma
            integer, intent(INOUT) :: iseed
            real, intent(OUT):: x,y
            
            real :: s, tmp

            s = 1.

            do while(s >= 1.)
                x = ranu(-1., 1., iseed)
                y = ranu(-1., 1., iseed)
                s = y**2 + x**2
            end do

            tmp = x*sqrt(-2.*log(s)/s)
            x = avg + sigma*tmp

            tmp = y*sqrt(-2.*log(s)/s)
            y = avg + sigma*tmp

        end subroutine rang

end module random