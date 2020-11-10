module source

    use random,    only : ranu, ran2
    use constants, only : Pi, twopi
    use vector_class

    implicit none
        
    contains
    
    subroutine point(pos, dir, cosThetaMax, iseed)
    ! emit isotropically from a point
    ! except bias towards lens
        use vector_class

        implicit none

        integer,      intent(INOUT) :: iseed
        type(vector), intent(OUT)   :: pos, dir
        real,         intent(IN)    :: cosThetaMax

        real :: phi, cosp, cost, sinp, sint, nxp, nyp, nzp, ran

        phi = twopi * ran2(iseed)
        cosp = cos(phi)
        sinp = sin(phi)  
        
        ran = ran2(iseed)
        !sample cone taken from pbrt
        cost = (1.d0 - ran) + ran*cosThetaMax
        sint = sqrt(1.d0 - cost**2)

        nxp = sint * cosp
        nyp = sint * sinp
        nzp = cost

        dir = vector(nxp, nyp, nzp)
        pos = vector(0.d0, 0.d0, 0.d0)

    end subroutine point


    subroutine ring(pos, dir, r1, r2, cosThetaMax, bottleRadius, bottleOffset, iseed)

        use vector_class

        implicit none

        type(vector), intent(OUT)   :: pos, dir
        integer,      intent(INOUT) :: iseed
        real,         intent(IN)    :: r1, r2, cosThetaMax, bottleRadius, bottleOffset

        type(vector) :: lenspoint, x, y, z
        real         :: r, theta, u(2), posx, posy, posz, dist, nxp, nyp, nzp

        !sample on an annulus radii, r1, r2
        r = ranu(r1, r2, iseed)
        theta = ran2(iseed) * twopi
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
        r = ranu(0., 12.7d-3**2, iseed)!0.21787185!0.20427731
        theta = ran2(iseed) * twopi
        posx = sqrt(r) * cos(theta)
        posy = sqrt(r) * sin(theta)
        lenspoint = vector(posx, posy, 37.5d-3)

        ! z = lenspoint - pos
        ! z = z%magnitude()

        dist = sqrt((lenspoint%x - pos%x)**2 + (lenspoint%y - pos%y)**2 + (lenspoint%z - pos%z)**2)
        ! print*,dist
        nxp = (lenspoint%x - pos%x) / dist
        nyp = (lenspoint%y - pos%y) / dist
        nzp = (lenspoint%z - pos%z) / dist
        dir = vector(nxp, nyp, nzp)
        dir = dir%magnitude()
        ! create new basis
        ! x = z .cross. vector(1.,0.,0.)
        ! x = x%magnitude()
        ! y = z .cross. x
        ! y = y%magnitude()

        ! ! sample onto lens
        ! u = [ran2(iseed), ran2(iseed)]
        ! dir = UniformSampleCone(u, 0.8708510322, y, x, z)
        ! dir = dir%magnitude()

    end subroutine ring

    type(vector) function UniformSampleCone(u, cosThetaMax, x, y, z)
    ! sample around a vector, e.g a in a cone
    ! takes three basis vectors for the coordinate system to be used where
    ! samples taken are with respect to the z axis of the given coordinate system.
    ! taken from PBRT
        use vector_class

        implicit none

        type(vector), intent(IN) :: x, y, z
        real,         intent(IN) :: u(2), cosThetaMax

        real :: sinTheta, cosTheta, phi

        cosTheta = lerp(u(1), cosThetaMax, 1.)
        sinTheta = sqrt(1. - cosTheta**2)
        phi = u(2) * twopi

        UniformSampleCone = cos(phi) * sinTheta * x + sin(phi) * sinTheta * y + cosTheta * z 

    end function UniformSampleCone


    real function Lerp(t, v1, v2)
    ! linear interpolate
    ! taken from PBRT

        implicit none

        real, intent(IN) :: t, v1, v2

        Lerp = (1. - t) * v1 + t * v2

    end function Lerp

end module source