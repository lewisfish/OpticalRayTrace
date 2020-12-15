module source

    use random,    only : ranu, ran2
    use constants, only : Pi, twopi
    use vector_class

    implicit none
        
    contains
    
    subroutine point(pos, dir, cosThetaMax)
    ! emit isotropically from a point
    ! except bias towards lens
        use vector_class

        implicit none

        type(vector), intent(OUT)   :: pos, dir
        real,         intent(IN)    :: cosThetaMax

        real :: phi, cosp, cost, sinp, sint, nxp, nyp, nzp, ran

        phi = twopi * ran2()
        cosp = cos(phi)
        sinp = sin(phi)  
        
        ran = ran2()
        !sample cone taken from pbrt
        cost = (1.d0 - ran) + ran*cosThetaMax
        sint = sqrt(1.d0 - cost**2)

        nxp = sint * cosp
        nyp = sint * sinp
        nzp = cost

        dir = vector(nxp, nyp, nzp)
        pos = vector(0.d0, 0.d0, 0.d0)

    end subroutine point


    subroutine create_spot(pos, dir, cosThetaMax, nrays, n)

        use vector_class

        implicit none

        type(vector), intent(OUT)   :: pos, dir
        real,         intent(IN)    :: cosThetaMax
        integer,      intent(IN)    :: nrays, n

        real :: phi, cosp, cost, sinp, sint, nxp, nyp, nzp
        real :: nrays_sqrt, phimax, thetaMax, deltaPhi, deltaTheta, theta

        nrays_sqrt = sqrt(real(nrays))

        phimax = twopi
        thetaMax = acos(cosThetaMax)

        deltaPhi = phimax / nrays_sqrt
        deltaTheta = thetaMax / nrays_sqrt

        phi = deltaPhi * mod(n, 10)
        theta = deltaTheta * int(n / 10)

        sinp = sin(phi)
        cosp = cos(phi)

        cost = cos(theta)
        sint = sqrt(1. - cost*cost)

        nxp = sint * cosp
        nyp = sint * sinp
        nzp = cost

        dir = vector(nxp, nyp, nzp)
        pos = vector(0., 0, 0.)

    end subroutine create_spot

    subroutine ring(pos, dir, r1, r2, bottleRadius, bottleOffset)

        use vector_class

        implicit none

        type(vector), intent(OUT)   :: pos, dir
        real,         intent(IN)    :: r1, r2, bottleRadius, bottleOffset

        type(vector) :: lenspoint
        real         :: r, theta, posx, posy, posz, dist, nxp, nyp, nzp

        !sample on an annulus radii, r1, r2
        r = ranu(r1, r2)
        theta = ran2() * twopi
        posx = sqrt(r) * cos(theta)
        posy = sqrt(r) * sin(theta)

        ! sample on curved bottle r^2 - y
        ! b is the offset of the bottle centre on the z axis
        ! assuming a = 0
        ! (y-a)^2 + (z-b)^2 = r^2
        ! bottleOffset = -bottleRadius / 2.d0
        posz = bottleOffset + sqrt((bottleRadius)**2 - posy**2)
        pos = vector(posx, posy, posz)

        !create new z-axis
        r = ranu(0., 12.7d-3**2)!0.21787185!0.20427731
        theta = ran2() * twopi
        posx = sqrt(r) * cos(theta)
        posy = sqrt(r) * sin(theta)
        lenspoint = vector(posx, posy, 37.5d-3)

        dist = sqrt((lenspoint%x - pos%x)**2 + (lenspoint%y - pos%y)**2 + (lenspoint%z - pos%z)**2)

        nxp = (lenspoint%x - pos%x) / dist
        nyp = (lenspoint%y - pos%y) / dist
        nzp = (lenspoint%z - pos%z) / dist

        dir = vector(nxp, nyp, nzp)
        dir = dir%magnitude()

    end subroutine ring


    subroutine emit_image(img, pos, dir)

        implicit none

        type(vector), intent(OUT)   :: pos, dir
        integer,      intent(INOUT) :: img(:, :)

        integer :: i, j

        do i = 1, size(img, 2)
            do j = 1, size(img, 1)
                if(img(j, i) > 0)then
                    call emit(pos, dir, j, i)
!$OMP atomic
                    img(j , i) = img(j, i) - 1
                    return
                end if
            end do
        end do

    end subroutine emit_image

    subroutine emit(pos, dir, i, j)
        
        implicit none

        type(vector), intent(OUT) :: pos, dir 
        integer,      intent(IN)  :: i, j

        real :: x, y, z, dist, r, theta, posx, posy, nxp, nyp, nzp
        real :: dx, dy
        type(vector) :: lenspoint

        dx = 2d-2 / 101.
        dy = dx

        x = ranu((i-1.) * dx, i * dx) - 1d-2
        y = ranu((j-1.) * dx, j * dx) - 1d-2
        z = 0.d0
        pos = vector(x, y, z)


        r = ranu(0., 12.7d-3**2)!
        theta = ran2() * twopi
        posx = sqrt(r) * cos(theta)
        posy = sqrt(r) * sin(theta)
        lenspoint = vector(posx, posy, 37.5d-3)

        dist = sqrt((lenspoint%x - pos%x)**2 + (lenspoint%y - pos%y)**2 + (lenspoint%z - pos%z)**2)

        nxp = (lenspoint%x - pos%x) / dist
        nyp = (lenspoint%y - pos%y) / dist
        nzp = (lenspoint%z - pos%z) / dist

        dir = vector(nxp, nyp, nzp)
        dir = dir%magnitude()

    end subroutine emit
end module source