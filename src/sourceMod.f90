module source

    use constants, only : Pi, twopi
    use lensMod,   only : plano_convex
    use random,    only : ranu, ran2, rang
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


    subroutine point_on_bottle(pos, dir, cosThetaMax, bottle, sors_radius, spot_radius)
    ! emit isotropically from a point
    ! except bias towards lens
        use vector_class
        use surfaces, only : intersect_cylinder
        use lensMod,     only : glass_bottle

        implicit none

        type(vector),       intent(OUT) :: pos, dir
        real,               intent(IN)  :: cosThetaMax, sors_radius, spot_radius
        type(glass_bottle), intent(IN)  :: bottle

        real :: phi, cosp, cost, sinp, sint, nxp, nyp, nzp, ran, t, tmp1, tmp2
        logical :: flag

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

        tmp1 = ranu(sors_radius-spot_radius, sors_radius+spot_radius)
        tmp2 = ranu(sors_radius-spot_radius, sors_radius+spot_radius)
        call rang(tmp1, tmp2, 0., 100d-6)

        pos = vector(tmp1, tmp2, 1.d0)
        dir = vector(0., 0., -1.)

        flag = intersect_cylinder(pos, dir, t, bottle%centre, bottle%radiusa+bottle%thickness)
        pos = pos + dir*t

        dir = vector(nxp, nyp, nzp)

    end subroutine point_on_bottle


    subroutine cross(pos, dir)
    ! creates a cross for debugging purposes
    !
    !
        use vector_class

        implicit none
        
        type(vector), intent(OUT) :: pos, dir
        real :: x, y, VorH

        VorH = ran2()
        if(VorH > 0.5)then
            x = ranu(-.05d-2, .05d-2)
            y = ranu(-.25d-2, .25d-2)
        else
            y = ranu(-.05d-2, .05d-2)
            if(ran2() > 0.5)then
                x = ranu(-.25d-2, -.05d-2)
            else
                x = ranu(.05d-2, .25d-2)
            end if
        end if


        pos = vector(x, y, 0.)
        dir  = vector(0., 0., 1.)

    end subroutine cross

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


    subroutine iSORS(pos, dir, bottle, lens, seperation, beam_width, ring)

        use lensMod

        implicit none

        type(plano_convex), intent(IN)  :: lens
        type(vector),       intent(OUT) :: pos, dir
        type(glass_bottle), intent(IN)  :: bottle
        real,               intent(IN)  :: seperation, beam_width
        logical,            intent(IN)  :: ring

        logical :: flag, skip
        real :: alpha, axicon_n, radius, height, k, base_pos
        real :: posx, posy, posz, rad1, rad2, t, dist, nxp, nyp, nzp, r, theta
        type(vector) :: centre, normal, lenspoint

        axicon_n = 1.4
        radius = 12.7d-3
        height = 1.1d-3
        alpha = atan(height / radius)
        k = (radius / height)**2

        base_pos = (seperation + beam_width) / tan(alpha * (axicon_n -1.))
        centre = vector(0., 0., 0.)

        call rang(posx, posy, 0., beam_width)
        pos = centre + vector(posx, posy, 2*height)
        dir = vector(0., 0., -1.)

        flag = intersect_cone(pos, dir, t, centre, radius, height)
        if(flag)then
            pos = pos + t*dir

            normal = vector(2*(pos%x-centre%x) / k, 2*(pos%y-centre%y) / k, -2*(pos%z-centre%z)+2*height)
            normal = normal *(-1.)! upper cone so invert normals
            normal = normal%magnitude()
            call reflect_refract(dir, normal, axicon_n, 1., flag)
            ! move packet to required distance to get given seperation
            t = (base_pos) / dir%z
            pos = pos + t*dir
            ! move packet into proper frame of reference, ie just beside the bottle
            pos%z = bottle%radiusa + epsilon(1.)
            if(ring)then
                if(bottle%ellipse)then
                    !need to divide by 2 to get a,b for ellipse equation
                    rad1 = bottle%radiusa - bottle%thickness
                    rad2 = bottle%radiusb - bottle%thickness
                    flag = intersect_ellipse(pos, dir, t, bottle%centre, rad1, rad2)
                else
                    flag = intersect_cylinder(pos, dir, t, bottle%centre, bottle%radiusa - bottle%thickness)
                end if
                if(.not. flag)then
                    error stop "no intersection with bottle!"
                end if
                pos = pos + t * dir
            else
                call bottle%backward(pos, dir, skip)
                !move to middle of bottle
                t = (bottle%z - pos%z) / dir%z
                pos = pos + t*dir
            end if
        end if

        if(ring)then
            r = ranu(0., (lens%radius)**2)
        else
            r = ranu(0., (lens%radius+10d-3)**2)
        end if
        theta = ran2() * twopi
        posx = sqrt(r) * cos(theta)
        posy = sqrt(r) * sin(theta)
        lenspoint = vector(posx, posy, lens%fb)

        dist = sqrt((lenspoint%x - pos%x)**2 + (lenspoint%y - pos%y)**2 + (lenspoint%z - pos%z)**2)

        nxp = (lenspoint%x - pos%x) / dist
        nyp = (lenspoint%y - pos%y) / dist
        nzp = (lenspoint%z - pos%z) / dist

        dir = vector(nxp, nyp, nzp)
        dir = dir%magnitude()

    end subroutine iSORS


    subroutine ring(pos, dir, lens, r1, r2, bottleRadiusa, bottleRadiusb, ellipse, bottleOffset)

        use vector_class

        implicit none

        type(plano_convex), intent(IN)  :: lens
        type(vector),       intent(OUT) :: pos, dir
        real,               intent(IN)  :: r1, r2, bottleRadiusa, bottleOffset, bottleRadiusb
        logical,            intent(IN)  :: ellipse

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
        if(ellipse)then
            posz = bottleOffset + sqrt(bottleRadiusa**2 - (posy*bottleRadiusa/bottleRadiusb)**2)
        else
            posz = bottleOffset + sqrt((bottleRadiusa)**2 - posy**2)
        end if

        pos = vector(posx, posy, posz)

        !create new z-axis
        r = ranu(0., (lens%radius+10d-3)**2)
        theta = ran2() * twopi
        posx = sqrt(r) * cos(theta)
        posy = sqrt(r) * sin(theta)
        lenspoint = vector(posx, posy, lens%fb)

        dist = sqrt((lenspoint%x - pos%x)**2 + (lenspoint%y - pos%y)**2 + (lenspoint%z - pos%z)**2)

        nxp = (lenspoint%x - pos%x) / dist
        nyp = (lenspoint%y - pos%y) / dist
        nzp = (lenspoint%z - pos%z) / dist

        dir = vector(nxp, nyp, nzp)
        dir = dir%magnitude()

    end subroutine ring


    subroutine emit_image(img, pos, dir, lens)

        implicit none

        type(plano_convex), intent(IN)    :: lens
        type(vector),       intent(OUT)   :: pos, dir
        integer,            intent(INOUT) :: img(:, :)

        integer :: i, j

        do i = 1, size(img, 2)
            do j = 1, size(img, 1)
                if(img(j, i) > 0)then
                    call emit(pos, dir, lens, j, i)
                    img(j , i) = img(j, i) - 1
                    return
                end if
            end do
        end do

    end subroutine emit_image

    subroutine emit(pos, dir, lens, i, j)
        
        implicit none

        type(plano_convex), intent(IN)  :: lens
        type(vector),       intent(OUT) :: pos, dir 
        integer,            intent(IN)  :: i, j

        real :: x, y, z, dist, r, theta, posx, posy, nxp, nyp, nzp
        real :: dx, dy
        type(vector) :: lenspoint

        dx = 5000d-6 / 512.
        dy = dx

        x = ranu((i-1.) * dx, i * dx) - 2500d-6
        y = ranu((j-1.) * dx, j * dx) - 2500d-6
        z = 0.d0
        pos = vector(x, y, z)


        r = ranu(0., lens%radius**2)!
        theta = ran2() * twopi
        posx = sqrt(r) * cos(theta)
        posy = sqrt(r) * sin(theta)
        lenspoint = vector(posx, posy, lens%fb)

        dist = sqrt((lenspoint%x - pos%x)**2 + (lenspoint%y - pos%y)**2 + (lenspoint%z - pos%z)**2)

        nxp = (lenspoint%x - pos%x) / dist
        nyp = (lenspoint%y - pos%y) / dist
        nzp = (lenspoint%z - pos%z) / dist

        dir = vector(nxp, nyp, nzp)
        dir = dir%magnitude()

    end subroutine emit

    subroutine init_emit_image(filename, imgin, nphotons, nphotonsLocal)
        
        use omp_lib

        implicit none
        
        character(*), intent(IN)  :: filename
        integer,      intent(IN)  :: nphotons
        integer,      intent(OUT) :: nphotonsLocal
        integer,      intent(OUT) :: imgin(512, 512)

        real :: imgout(512, 512), tot, tmp, diff
        integer :: i, j, u
        
        imgout = 0

        open(newunit=u,file=filename, form="unformatted", access="stream", status="old")
        read(u)imgout
        close(u)

        ! array written out in wrong fashion
        imgout = transpose(imgout)

        !normalise
        tot = sum(imgout)
        imgin = 0
#ifdef _OPENMP
        nphotonsLocal = nphotons / omp_get_max_threads()
#else
        nphotonsLocal = nphotons
#endif
        do i = 1, 512
            do j = 1, 512
                tmp = (dble(nphotonsLocal) * imgout(i, j)) / dble(tot)

                diff = tmp - int(tmp)
                if(ran2() < diff .and. diff > 0)then
                    imgin(i,j) = imgin(i, j) + (int(tmp) + 1)
                else
                    imgin(i,j) = imgin(i, j) + int(tmp)
                end if
            end do
        end do
    end subroutine init_emit_image
end module source