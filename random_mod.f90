module random

    implicit none

    private
    public  :: ran2, ranu, rang, init_rng

    contains

        subroutine init_rng(iseed, id)

            implicit none

            integer, intent(OUT) :: iseed
            integer, intent(IN)  :: id

            iseed = 2137834897 + id
            iseed = -abs(iseed)

        end subroutine init_rng

        real function ran2(seed)

            implicit none

            integer,      parameter     :: k4b=selected_int_kind(9)
            integer(k4b), intent(INOUT) :: seed

            integer(k4b), parameter :: ia=16807, im=2147483647, iq=127773, ir=2836
            real,         save      :: am
            integer(k4b), save      :: ix=-1, iy=-1,k

            if(seed <= 0 .or. iy < 0)then
                am = nearest(1.0, -1.0)/im
                iy = ior(ieor(888889999,abs(seed)),1)
                ix = ieor(777755555,abs(seed))
                seed=abs(seed)+1
            end if

            ix = ieor(ix,ishft(ix,13))
            ix = ieor(ix,ishft(ix,-17))
            ix = ieor(ix,ishft(ix,5))
            k = iy/iq
            iy = ia*(iy-k*iq)-ir*k
            if(iy < 0)iy = iy + im
            ran2 = am*ior(iand(im,ieor(ix,iy)),1)

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