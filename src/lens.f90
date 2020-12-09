module lensMod

    use vector_class
    use surfaces

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


    subroutine glass_bottle(pos, dir, wave, bottleThickness, u, iseed, skip)

        use stackMod, only : stack

        implicit none

        type(vector), intent(INOUT) :: pos, dir
        type(stack),  intent(INOUT) :: u
        integer,      intent(INOUT) :: iseed
        logical,      intent(OUT)   :: skip
        real,         intent(IN)    :: wave, bottleThickness

        type(vector) :: centre, normal, orig
        real         :: t, radius, n, wave2
        logical      :: flag

        !centre of bottle
        centre = vector(0., 0., 0.)

        ! refractive index for clear container glass
        !https://refractiveindex.info/?shelf=glass&book=soda-lime&page=Rubin-clear
        ! wave2 = (wave*1d6 )**2
        n = 1.5163231388184517!1.5130 - 0.003169*wave2 + 0.003962/wave2
        !inner surface
        radius = 17.5d-3 - bottleThickness
        flag = intersect_cylinder(pos, dir, t, centre, radius)
        if(.not. flag)then
            skip = .true.
            return
        end if

        orig = pos
        pos = pos + t * dir
        ! call u%push(vector(pos%x, pos%y, pos%z))

        orig%z = pos%z
        normal = centre - orig
        normal = normal%magnitude()
        flag = .false.
        call reflect_refract(dir, normal, 1.d0, n, iseed, flag)
        if(flag)then
            skip = .true.
            return
        end if

        !outer surface
        radius = 17.5d-3
        flag = intersect_cylinder(pos, dir, t, centre, radius)
        if(.not. flag)then
            skip = .true.
            return
        end if

        orig = pos
        pos = pos + t * dir
        ! call u%push(vector(pos%x, pos%y, pos%z))

        orig%z = pos%z
        normal = centre - orig

        normal = normal%magnitude()

        flag = .false.
        call reflect_refract(dir, normal, n, 1.d0, iseed, flag)
        if(flag)then
            skip = .true.
            return
        end if

    end subroutine glass_bottle

    subroutine plano_forward_sub(this, pos, dir, u, iseed, skip)

        use stackMod, only : stack

        implicit none

        class(plano_convex) :: this

        type(vector), intent(INOUT) :: pos, dir
        integer,      intent(INOUT) :: iseed
        type(stack),  intent(INOUT) :: u
        logical,      intent(OUT) :: skip

        type(vector) :: curvedNormal
        real         :: d, t, r
        logical      :: flag

        skip = .false.


        ! move to flat surface
        d = (this%fb - pos%z) / dir%z
        pos = pos + dir * d
        r = sqrt(pos%x**2 + pos%y**2)
        if(r > this%radius)then
            skip=.true.
            return
        end if

        ! call u%push(vector(pos%x, pos%y, pos%z))

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
        if(flag)then
            skip = .true.
            return
        end if


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

        ! origpos = pos

        ! d = ((this%centre1%z - this%r1) - pos%z) / dir%z
        ! pos = pos + dir * d

        ! ! make sure no rays get propagated that are outside lens radius
        ! ! this can double as an Iris
        ! r = sqrt(pos%x**2 + pos%y**2)
        ! if(r > this%radius*(1./1.))then
        !     skip=.true.
        !     return
        ! end if
        ! pos = origpos

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
        if(flag)then
            skip = .true.
            return
        end if


        ! call u%push(vector(pos%x, pos%y, pos%z))


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
        if(flag)then
            skip = .true.
            return
        end if


        ! call u%push(vector(pos%x, pos%y, pos%z))

        !third sphere
        flag = intersect_sphere(pos, dir, t, this%centre3, this%R3)
        if(.not. flag)error stop "Help3"
        pos = pos + t * dir

        normal = this%centre3 - pos
        normal = normal%magnitude()

        flag = .false.
        call reflect_refract(dir, normal, this%n3, this%n1, iseed, flag)
        if(flag)then
            skip = .true.
            return
        end if


        ! call u%push(vector(pos%x, pos%y, pos%z))

        ! origpos = pos

        ! d = ((this%centre3%z + this%r3) - pos%z) / dir%z
        ! pos = pos + dir * d

        ! r = sqrt(pos%x**2 + pos%y**2)
        ! if(r > this%radius*(1./1.))then
        !     skip=.true.
        !     return
        ! end if
        ! pos = origpos



    end subroutine doublet_forward_sub
end module lensMod