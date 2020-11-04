module lensMod

    implicit none

    contains

    subroutine plano_convex(pos, dir, D1, Lt, R, u, iseed)
        use vector_class

        implicit none

        type(vector), intent(INOUT) :: pos, dir
        integer,      intent(INOUT) :: iseed, u
        real,         intent(IN)    :: D1, Lt, R

        type(vector) :: centre, flatNormal, curvedNormal
        real :: n1, n2, n3, d, t
        logical :: flag

        n1 = 1.0  !air
        n2 = 1.51 !https://www.filmetrics.com/refractive-index-database/Schott+N-BK7
        n3 = n1

        ! move to flat surface
        d = (D1 - pos%z) / dir%z
        pos = pos + dir * d

        ! centre of lens defining sphere
        centre = vector(0., 0., (D1 + Lt) - R)
        !check for refraction on flat surface
        flatNormal = vector(0., 0., -1.)
        flag = .false.
        call reflect_refract(dir, flatNormal, n1, n2, iseed, flag)

        ! intersect curved side and move to it
        flag = intersect(pos, dir, t, centre, R)
        if(.not. flag)error stop "Help"
        pos = pos + t * dir

        !refract
        curvedNormal = centre - pos
        curvedNormal = curvedNormal%magnitude()
        flag = .false.
        call reflect_refract(dir, curvedNormal, n2, n3, iseed, flag)
    end subroutine plano_convex


    subroutine achromatic_doublet(pos, dir, D1, D2, D3, R1, R2, R3, tc1, tc2, fb, lensRadius, u, iseed, skip)

        use vector_class

        implicit none

        real ,        intent(IN)    :: D1, D2, D3, R1, R2, R3, tc1, tc2, fb, lensRadius
        type(vector), intent(INOUT) :: pos, dir
        integer,      intent(INOUT) :: iseed, u
        logical,      intent(OUT)   :: skip

        type(vector) :: normal, centre
        logical      :: flag
        real         :: n1, n2, n3, n4, t, tc

        skip = .false.

        tc = tc1 + tc2
        n1 = 1.0 !air
        n2 = 1.64 !https://www.filmetrics.com/refractive-index-database/Schott+N-LAK22
        n3 = 1.79 !https://www.filmetrics.com/refractive-index-database/Schott+N-SF6
        n4 = n1

        !first sphere
        centre = vector(0., 0., D1 + D2 + D3 - tc - fb + R1)
        flag = intersect(pos, dir, t, centre, R1)
        if(.not. flag)then
            skip=.true.
            return
        end if
        pos = pos + t * dir
        ! make sure no rays get propagated that are outside lens radius
        ! this can double as an Iris
        if(sqrt(pos%x**2 + pos%y**2) > (lensRadius/1.))then
            skip=.true.
            return
        end if

        normal = pos - centre
        normal = normal%magnitude()

        flag = .false.
        call reflect_refract(dir, normal, n1, n2, iseed, flag)


        !second sphere
        centre = vector(0., 0., D1 + D2 + D3 - tc - fb + tc1 - R2)
        flag = intersect(pos, dir, t, centre, R2)
        if(.not. flag)then
            skip = .true.
            return
        end if
        pos = pos + t * dir
        normal = centre - pos
        normal = normal%magnitude()

        flag = .false.
        call reflect_refract(dir, normal, n2, n3, iseed, flag)


        !second sphere
        centre = vector(0., 0., D1 + D2 + D3 - fb - R3)
        flag = intersect(pos, dir, t, centre, R3)
        if(.not. flag)error stop "Help3"
        pos = pos + t * dir
        normal = centre - pos
        normal = normal%magnitude()

        flag = .false.
        call reflect_refract(dir, normal, n3, n4, iseed, flag)

    end subroutine achromatic_doublet


    logical function intersect(orig, dir, t, centre, radius)
    ! calculates where a line, with origin:orig and direction:dir hits a sphere, centre:centre and radius:radius
    ! returns true if intersection exists
    ! returns t, the paramertised parameter of the line equation
    ! adapted from scratchapixel
        
        use vector_class, only : vector

        implicit none

        type(vector), intent(IN)  :: dir, orig, centre
        real,         intent(OUT) :: t
        real,         intent(IN)  :: radius

        type(vector) :: L
        real         :: t0, t1, a, b, c, tmp

        intersect = .false.

        L = orig - centre
        a = dir .dot. dir
        b = 2.d0 * (dir .dot. L)
        c = (l .dot. l) - radius**2

        if(.not. solveQuadratic(a, b, c, t0, t1))return
        if(t0 > t1)then
            tmp = t1
            t1 = t0
            t0 = tmp
        end if
        if(t0 < 0.d0)then
            t0 = t1
            if(t0 < 0.)return
        end if

        t = t0
        intersect = .true.
        return

    end function intersect

    logical function solveQuadratic(a, b, c, x0, x1)
    ! solves quadratic equation given coeffs a, b, and c
    ! returns true if real soln
    ! returns x0 and x1
    ! adapted from scratchapixel

        implicit none

        real, intent(IN)  :: a, b, c
        real, intent(OUT) :: x0, x1

        real :: discrim, q

        solveQuadratic = .false.

        discrim = b**2 - 4.d0 * a * c
        if(discrim < 0.d0)then
            return
        elseif(discrim == 0.d0)then
            x0 = -0.5*b/a
            x1 = x0
        else
            if(b > 0.d0)then
                q = -0.5d0 * (b + sqrt(discrim))
            else
                q = -0.5d0 * (b - sqrt(discrim))
            end if
            x0 = q / a
            x1 = c / q
        end if
        solveQuadratic = .true.
        return

    end function solveQuadratic

    subroutine reflect_refract(I, N, n1, n2, iseed, rflag)

        use vector_class
        use random, only : ran2

        implicit none

        type(vector), intent(INOUT) :: I
        type(vector), intent(INOUT) :: N
        real,         intent(IN)    :: n1, n2
        integer,      intent(INOUT) :: iseed
        logical,      intent(OUT)   :: rflag

        rflag = .FALSE.
        ! print*,fresnel(I, N, n1, n2)
        ! if(ran2(iseed) <= fresnel(I, N, n1, n2))then
        !     ! call reflect(I, N)
        !     rflag = .true.
        ! else
            call refract(I, N, n1/n2)
        ! end if

    end subroutine reflect_refract


    subroutine reflect(I, N)
    !   get vector of reflected photon
    !
    !
        use vector_class

        implicit none

        type(vector), intent(INOUT) :: I
        type(vector), intent(IN)    :: N

        type(vector) :: R

        R = I - 2. * (N .dot. I) * N
        I = R

    end subroutine reflect


    subroutine refract(I, N, eta)
    !   get vector of refracted photon
    !
    !
        use vector_class

        implicit none

        type(vector), intent(INOUT) :: I
        type(vector), intent(IN)    :: N
        real,         intent(IN)    :: eta

        type(vector) :: T, Ntmp

        real :: c1, c2

        Ntmp = N

        c1 = (Ntmp .dot. I)
        if(c1 < 0.)then
            c1 = -c1
        else
            Ntmp = (-1.) * N
        end if

        c2 = sqrt(1.d0 - (eta)**2 * (1.d0 - c1**2))

        T = eta*I + (eta * c1 - c2) * Ntmp 

        I = T

    end subroutine refract


    function fresnel(I, N, n1, n2) result (tir)
    !   calculates the fresnel coefficents
    !
    !
        use vector_class
        use ieee_arithmetic, only : ieee_is_nan

        implicit none

        real, intent(IN)         :: n1, n2
        type(vector), intent(IN) :: I, N

        real             ::  costt, sintt, sint2, cost2, tir, f1, f2

        costt = abs(I .dot. N)

        sintt = sqrt(1. - costt * costt)
        sint2 = n1/n2 * sintt
        if(sint2 > 1.)then
            tir = 1.0
            return
        elseif(costt == 1.)then
            tir = 0.
            return
        else
            sint2 = (n1/n2)*sintt
            cost2 = sqrt(1. - sint2 * sint2)
            f1 = abs((n1*costt - n2*cost2) / (n1*costt + n2*cost2))**2
            f2 = abs((n1*cost2 - n2*costt) / (n1*cost2 + n2*costt))**2

            tir = 0.5 * (f1 + f2)
        if(ieee_is_nan(tir) .or. tir > 1. .or. tir < 0.)print*,'TIR: ', tir, f1, f2, costt,sintt,cost2,sint2
            return
        end if
    end function fresnel
end module lensMod