module lensMod

    use vector_class

    implicit none

    type :: plano_convex

        real :: thickness, curve_radius, diameter, radius, fb, f
        real :: n1, n2
        type(vector) :: centre, flatNormal
        contains
        procedure :: forward  => plano_forward_sub !flat face first
    end type plano_convex

    interface plano_convex
        module procedure init_plano_convex
    end interface plano_convex

    type :: achromatic_doublet
        real :: thickness1, thickness2, thickness, R1, R2, R3
        real :: diameter, radius, fb, f
        real :: n1, n2, n3
        type(vector) :: centre1, centre2, centre3
        contains
        procedure :: forward => doublet_forward_sub!r1 face first
    end type achromatic_doublet

    interface achromatic_doublet
        module procedure init_achromatic_doublet
    end interface achromatic_doublet


    contains

    type(achromatic_doublet) function init_achromatic_doublet(file, wavelength, D1) result(this)

        implicit none

        character(*), intent(IN) :: file
        real,         intent(IN) :: D1, wavelength
        
        integer :: u
        real    :: b11, b21, b31, c11, c21, c31
        real    :: b12, b22, b32, c12, c22, c32

        open(newunit=u,file=trim(file),status='old')
            read(u,*) this%thickness1
            read(u,*) this%thickness2
            read(u,*) this%R1
            read(u,*) this%R2
            read(u,*) this%R3
            read(u,*) this%diameter
            read(u,*) this%f
            read(u,*) this%fb
            read(u,*) this%n1
            read(u,*) b11
            read(u,*) b21
            read(u,*) b31
            read(u,*) c11
            read(u,*) c21
            read(u,*) c31
            read(u,*) b12
            read(u,*) b22
            read(u,*) b32
            read(u,*) c12
            read(u,*) c22
            read(u,*) c32
        close(u)

        this%n2 = Sellmeier(wavelength, b11, b21, b31, c11, c21, c31)
        this%n3 = Sellmeier(wavelength, b12, b22, b32, c12, c22, c32)

        this%radius = this%diameter / 2.d0
        this%thickness = this%thickness1 + this%thickness2

        this%centre1 = vector(0., 0., D1 + (D1 + this%f) + this%f - this%thickness - this%fb + this%R1)
        this%centre2 = vector(0., 0., D1 + (D1 + this%f) + this%f - this%thickness - this%fb + this%thickness1 - this%R2)
        this%centre3 = vector(0., 0., D1 + (D1 + this%f) + this%f - this%fb - this%R3)

    end function init_achromatic_doublet


    type(plano_convex) function init_plano_convex(file, wavelength) result(this)

        implicit none

        character(*), intent(IN) :: file
        real,         intent(IN) :: wavelength

        integer :: u
        real    :: b1, b2, b3, c1, c2, c3

        open(newunit=u,file=trim(file),status='old')
            read(u,*) this%thickness
            read(u,*) this%curve_radius
            read(u,*) this%diameter
            read(u,*) this%f
            read(u,*) this%fb
            read(u,*) this%n1
            read(u,*) b1
            read(u,*) b2
            read(u,*) b3
            read(u,*) c1
            read(u,*) c2
            read(u,*) c3
        close(u)

        this%n2 = Sellmeier(wavelength, b1, b2, b3, c1, c2, c3)
        this%radius = this%diameter / 2.d0

        this%centre = vector(0., 0., (this%fb + this%thickness) - this%curve_radius)
        this%flatNormal = vector(0., 0., -1.)

    end function init_plano_convex

    real function Sellmeier(wave, b1, b2, b3, c1, c2, c3)
    ! Sellmeier equation
    ! wave is in units of nm

        implicit none

        real, intent(IN) :: wave, b1, b2, b3, c1, c2, c3
        real :: wave2, a, b ,c

        !convert to units of um
        wave2 = (wave*1d6)**2

        a = (b1 * wave2) / (wave2 - c1)
        b = (b2 * wave2) / (wave2 - c2)
        c = (b3 * wave2) / (wave2 - c3)

        Sellmeier = sqrt(1.d0 + (a + b + c))

    end function Sellmeier


    subroutine plano_forward_sub(this, pos, dir, u, iseed, skip)

        use stackMod, only : pointtype, stack

        implicit none

        class(plano_convex) :: this

        type(vector), intent(INOUT) :: pos, dir
        integer,      intent(INOUT) :: iseed
        type(stack),  intent(INOUT) :: u
        logical,      intent(OUT) :: skip

        type(vector) :: centre, flatNormal, curvedNormal
        real :: n1, n2, n3, d, t
        logical :: flag

        skip = .false.


        ! move to flat surface
        d = (this%fb - pos%z) / dir%z
        pos = pos + dir * d
        ! call u%push(pointtype(pos%x, pos%y, pos%z))

        flag = .false.
        call reflect_refract(dir, this%flatNormal, this%n1, this%n2, iseed, flag)

        ! intersect curved side and move to it
        flag = intersect_sphere(pos, dir, t, this%centre, this%curve_radius)
        !ray exits lens through top
        if(.not. flag)then
            skip = .true.
            return
        end if

        pos = pos + t * dir

        !refract
        curvedNormal = this%centre - pos
        curvedNormal = curvedNormal%magnitude()
        flag = .false.
        call reflect_refract(dir, curvedNormal, this%n2, this%n1, iseed, flag)
    end subroutine plano_forward_sub


    subroutine doublet_forward_sub(this, pos, dir, D1, D2, u, iseed, skip)

        use stackMod, only : pointtype, stack

        implicit none

        class(achromatic_doublet) :: this

        real,         intent(IN)    :: D1, D2
        type(vector), intent(INOUT) :: pos, dir
        integer,      intent(INOUT) :: iseed
        type(stack),  intent(INOUT) :: u
        logical,      intent(OUT)   :: skip

        type(vector) :: normal, origpos
        logical      :: flag
        real         :: t, d, r

        skip = .false.

        ! r = sqrt(pos%x**2 + pos%y**2)
        ! if(r < this%radius*(8./10.).or. r > this%radius*(10./10.))then
        !     skip=.true.
        !     return
        ! end if

        !first sphere
        flag = intersect_sphere(pos, dir, t, this%centre1, this%R1)
        if(.not. flag)then
            skip=.true.
            return
        end if
        pos = pos + t * dir
        ! make sure no rays get propagated that are outside lens radius
        ! this can double as an Iris
        r = sqrt(pos%x**2 + pos%y**2)
        if(r > (this%radius*(1.)))then
            skip=.true.
            return
        end if

        normal = pos - this%centre1
        normal = normal%magnitude()

        flag = .false.
        call reflect_refract(dir, normal, this%n1, this%n2, iseed, flag)


        !second sphere
        flag = intersect_sphere(pos, dir, t, this%centre2, this%R2)
        if(.not. flag)then
            skip = .true.
            return
        end if
        pos = pos + t * dir

        normal = this%centre2 - pos
        normal = normal%magnitude()

        flag = .false.
        call reflect_refract(dir, normal,this%n2, this%n3, iseed, flag)

        !third sphere
        flag = intersect_sphere(pos, dir, t, this%centre3, this%R3)
        if(.not. flag)error stop "Help3"
        pos = pos + t * dir

        normal = this%centre3 - pos
        normal = normal%magnitude()

        flag = .false.
        call reflect_refract(dir, normal, this%n3, this%n1, iseed, flag)

    end subroutine doublet_forward_sub


    logical function intersect_sphere(orig, dir, t, centre, radius)
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

        intersect_sphere = .false.

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
        intersect_sphere = .true.
        return

    end function intersect_sphere

    logical function intersect_cylinder(orig, dir, t, centre, radius)
    ! calculates where a line, with origin:orig and direction:dir hits a cylinder, centre:centre and radius:radius
    ! returns true if intersection exists
    ! returns t, the paramertised parameter of the line equation
    ! adapted from scratchapixel
    ! need to check z height after moving ray
    ! if not this is an infinte cylinder
        
        use vector_class, only : vector

        implicit none

        type(vector), intent(IN)  :: dir, orig, centre
        real,         intent(OUT) :: t
        real,         intent(IN)  :: radius

        type(vector) :: L
        real         :: t0, t1, a, b, c, tmp

        intersect_cylinder = .false.

        L = orig - centre
        a = dir%x**2 + dir%y**2
        b = 2 * (dir%x * L%x + dir%y * L%y)
        c = L%x**2 + L%y**2 - radius**2

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
        intersect_cylinder = .true.
        return

    end function intersect_cylinder

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