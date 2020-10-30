module source

    use random, only : ranu, ran2
    use constants, only: Pi, twopi
    use vector_class

    implicit none
        
    contains
    
    subroutine point(pos, dir, cosThetaMax, iseed)
    ! emit isotropically from a point
    ! except bias towards lens
        use vector_class

        implicit none

        type(vector), intent(OUT) :: pos, dir
        integer, intent(INOUT) :: iseed
        real , intent(IN) :: cosThetaMax

        real :: phi, cosp, cost, sinp, sint, nxp, nyp, nzp, ran

        phi = twopi * ran2(iseed)
        cosp = cos(phi)
        sinp = sin(phi)  
        
        ran = ran2(iseed)
        !sample cone taken from pbrt
        cost = (1. - ran) + ran*cosThetaMax
        sint = sqrt(1. - cost**2)

        nxp = sint * cosp  
        nyp = sint * sinp
        nzp = cost

        dir = vector(nxp, nyp, nzp)
        pos = vector(0, 0., 0.)

    end subroutine point


    subroutine ring(pos, dir, r1, r2, cosThetaMax, iseed)

        use vector_class

        implicit none

        type(vector), intent(OUT)   :: pos, dir
        integer,      intent(INOUT) :: iseed
        real,         intent(IN)    :: r1, r2, cosThetaMax

        type(vector) :: lenspoint, x, y, z, angles, rot1, rot2, rot3
        real :: r, theta, nxp, nyp, nzp, u(2)
        real :: cosa, cosb, cosg, sina, sinb, sing

        !sample on an annulus radii, r1, r2
        r = ranu(r1, r2, iseed)
        theta = ran2(iseed) * twopi
        pos%x = sqrt(r) * cos(theta)
        pos%y = sqrt(r) * sin(theta)
        ! sample on curved bottle r^2 - y
        pos%z = sqrt((17.5d-3)**2 - pos%y**2)

        !create new z-axis
        lenspoint = vector(0., 0., 40d-3)! ring around point source
        z = lenspoint - pos
        z = z%magnitude()

        ! create new basis
        x = z .cross. vector(1.,0.,0.)
        x = x%magnitude()
        y = z .cross. x
        y = y%magnitude()

        ! sample onto lens
        u = [ran2(iseed), ran2(iseed)]
        dir = UniformSampleCone(u, 0.8708510322, x, y, z)

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