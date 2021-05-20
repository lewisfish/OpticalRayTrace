module lensMod

    use vector_class
    use surfaces

    implicit none

    type :: lens
        real :: thickness, diameter, radius, fb, f, n1, n2
        contains
        procedure :: forward => fwdLens
    end type lens

    type, extends(lens) :: plano_convex
        real :: curve_radius
        type(vector) :: centre, flatNormal
        contains
        procedure :: forward  => plano_forward_sub !flat face first
        procedure :: backward  => plano_backward_sub
    end type plano_convex

    interface plano_convex
        module procedure init_plano_convex
    end interface plano_convex


    type, extends(lens) :: achromatic_doublet
        real :: thickness1, thickness2, R1, R2, R3
        real :: n3
        type(vector) :: centre1, centre2, centre3
        contains
        procedure :: forward => doublet_forward_sub!r1 face first
    end type achromatic_doublet

    interface achromatic_doublet
        module procedure init_achromatic_doublet
    end interface achromatic_doublet


    type :: glass_bottle
        real         :: nbottle, ncontents, thickness, radiusa, radiusb, mua, mus
        type(vector) :: centre
        logical      :: ellipse, beer_lambert
        contains
        procedure :: forward => bottle_forward_sub
        procedure :: backward => bottle_backward_sub
    end type glass_bottle

    interface glass_bottle
        module procedure init_bottle
    end interface glass_bottle


    contains

    subroutine fwdLens(this, pos, dir, u, skip, iris, iris_radius)

        use stackMod, only: stack

        implicit none
        
        class(lens) :: this

        logical, optional, intent(IN)    :: iris(2)
        real,    optional, intent(IN)    :: iris_radius
        logical,           intent(OUT)   :: skip
        type(stack),       intent(INOUT) :: u
        type(vector),      intent(INOUT) :: pos, dir

    end subroutine fwdLens

    type(achromatic_doublet) function init_achromatic_doublet(file, wavelength, offset_in) result(this)

        implicit none

        character(*),   intent(IN) :: file
        real,           intent(IN) :: wavelength
        real, optional, intent(IN) :: offset_in

        integer :: u
        real    :: b11, b21, b31, c11, c21, c31
        real    :: b12, b22, b32, c12, c22, c32, offset

        if(present(offset_in))then
            offset = offset_in
        else
            offset = 0.0d0
        end if


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

        this%centre1 = vector(0., 0., offset + this%fb + this%r1)
        this%centre2 = vector(0., 0., offset + this%fb + this%thickness1 - this%r2)
        this%centre3 = vector(0., 0., offset + this%fb + this%thickness - this%r3)

    end function init_achromatic_doublet


    type(plano_convex) function init_plano_convex(file, wavelength, offset_in) result(this)

        implicit none

        character(*),   intent(IN) :: file
        real,           intent(IN) :: wavelength
        real, optional, intent(IN) :: offset_in

        integer :: u
        real    :: b1, b2, b3, c1, c2, c3, offset

        if(present(offset_in))then
            offset = offset_in
        else
            offset = 0.0d0
        end if

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

        this%centre = vector(0., 0., offset + (this%fb + this%thickness) - this%curve_radius)
        this%flatNormal = vector(0., 0., -1.)

    end function init_plano_convex


    type(glass_bottle) function init_bottle(file, wavelength) result(this)

        use iso_fortran_env, Only : iostat_end

        implicit none
    
        real,         intent(IN) :: wavelength
        character(*), intent(IN) :: file

        integer :: u, io
        real    :: x, y, z, a1, b1, c1, a2, b2, c2, mua, mus

        open(newunit=u, file=file, status="old")
            read(u,*)this%thickness
            read(u,*)this%radiusa
            read(u,*)this%radiusb
            read(u,*)x
            read(u,*)y
            read(u,*)z
            read(u,*)a1
            read(u,*)b1
            read(u,*)c1
            read(u,*)a2
            read(u,*)b2
            read(u,*)c2
            read(u,*,iostat=io)mua
            if(io == iostat_end)then
                this%mua = 0.d0
            else
                this%mua = mua
            end if
            read(u,*,iostat=io)mus
            if(io == iostat_end)then
                this%mus = 0.d0
            else
                this%mus = mus
            end if
        close(u)

        this%centre = vector(x, y, z)
        this%nbottle = dispersion(wavelength, a1, b1, c1)
        this%ncontents = cauchy(wavelength, a2, b2, c2)

        this%beer_lambert = .false.
        if(this%mua /= 0.0)this%beer_lambert = .true.

        if(this%radiusa /= this%radiusb)then
            this%ellipse = .true.
        else
            this%ellipse = .false.
        end if

    end function init_bottle


    subroutine bottle_forward_sub(this, pos, dir, u, skip)

        use stackMod, only : stack
        use random,   only : ran2

        implicit none

        class(glass_bottle) :: this
        type(vector), intent(INOUT) :: pos, dir
        type(stack),  intent(INOUT) :: u
        logical,      intent(OUT)   :: skip

        type(vector) :: normal, orig, newpos
        real         :: t, rad1, rad2, dist, taudist, tau
        logical      :: flag

        !inner surface
        if(this%ellipse)then
            !need to divide by 2 to get a,b for ellipse equation
            rad1 = this%radiusa - this%thickness
            rad2 = this%radiusb - this%thickness
            flag = intersect_ellipse(pos, dir, t, this%centre, rad1, rad2)
        else
            flag = intersect_cylinder(pos, dir, t, this%centre, this%radiusa - this%thickness)
        end if
        if(.not. flag)then
            skip = .true.
            return
        end if
        pos = pos + t * dir
        orig = pos
        ! call u%push(pos)

        orig%x = this%centre%x
        normal = this%centre - orig
        normal = normal%magnitude()

        flag = .false.
        call reflect_refract(dir, normal, this%ncontents, this%nbottle, flag)
        if(flag)then
            skip = .true.
            return
        end if

        !outer surface
        if(this%ellipse)then
            flag = intersect_ellipse(pos, dir, t, this%centre, this%radiusa/2., this%radiusb/2.)
        else
            flag = intersect_cylinder(pos, dir, t, this%centre, this%radiusa)
        end if
        if(.not. flag)then
            skip = .true.
            return
        end if

        flag = .false.
        if(this%beer_lambert)then
            tau = -log(ran2())
            !calculate optical distance of medium
            newpos = pos + t*dir
            dist = e_dist(pos, newpos)
            taudist = dist * (this%mua + this%mus)
            if(tau < taudist)then!interact
                t = tau / (this%mua + this%mus)
                !packet absorbed
                flag = .true.
            end if
        end if

        pos = pos + t * dir
        if(flag)then
            skip = .true.
            return
        end if

        orig = pos
        ! call u%push(pos)
        orig%x = this%centre%x

        normal = this%centre - orig
        normal = normal%magnitude()

        flag = .false.
        call reflect_refract(dir, normal, this%nbottle, 1.d0, flag)
        if(flag)then
            skip = .true.
            return
        end if

    end subroutine bottle_forward_sub

    subroutine bottle_backward_sub(this, pos, dir, skip, u)

        use stackMod, only : stack
        use random,   only : ran2

        implicit none

        class(glass_bottle) :: this
        type(vector), intent(INOUT) :: pos, dir
        type(stack),  intent(INOUT), optional :: u
        logical,      intent(OUT)   :: skip

        type(vector) :: normal, orig, newpos
        real         :: t, rad1, rad2, dist, tau, taudist
        logical      :: flag

        !outer surface
        if(this%ellipse)then
            !need to divide by 2 to get a,b for ellipse equation
            rad1 = this%radiusa
            rad2 = this%radiusb
            flag = intersect_ellipse(pos, dir, t, this%centre, rad1, rad2)
        else
            flag = intersect_cylinder(pos, dir, t, this%centre, this%radiusa)
        end if
        if(.not. flag)then
            skip = .true.
            return
        end if
        pos = pos + t * dir
        orig = pos
        ! call u%push(pos)

        orig%x = this%centre%x
        normal = orig - this%centre
        normal = normal%magnitude()

        flag = .false.
        call reflect_refract(dir, normal, 1., this%nbottle, flag)
        if(flag)then
            skip = .true.
            return
        end if

        !outer surface
        if(this%ellipse)then
            flag = intersect_ellipse(pos, dir, t, this%centre, this%radiusa-this%thickness, this%radiusb-this%thickness)
        else
            flag = intersect_cylinder(pos, dir, t, this%centre, this%radiusa - this%thickness)
        end if

        if(.not. flag)then
            skip = .true.
            return
        end if

        if(this%beer_lambert)then
            tau = -log(ran2())
            !calculate optical distance of medium
            newpos = pos + t*dir
            dist = e_dist(pos, newpos)
            taudist = dist * (this%mua + this%mus)
            print*,tau,taudist
            if(tau < taudist)then!interact
                t = tau / (this%mua + this%mus)
            end if
        end if

        pos = pos + t * dir
        orig = pos
        ! call u%push(pos)
        orig%x = this%centre%x

        normal = orig - this%centre
        normal = normal%magnitude()

        flag = .false.
        call reflect_refract(dir, normal, this%nbottle, this%ncontents, flag)
        if(flag)then
            skip = .true.
            return
        end if

    end subroutine bottle_backward_sub

    subroutine plano_forward_sub(this, pos, dir, u, skip, iris, iris_radius)

        use stackMod, only : stack

        implicit none

        class(plano_convex) :: this

        type(vector), intent(INOUT) :: pos, dir
        type(stack),  intent(INOUT) :: u
        logical,      intent(OUT) :: skip
        logical, optional, intent(IN) :: iris(2)
        real, optional, intent(IN) :: iris_radius

        type(vector) :: curvedNormal
        real         :: d, t, r, a
        logical      :: flag

        skip = .false.


        ! move to flat surface
        a = this%centre%z + this%curve_radius - this%thickness
        d = (a - pos%z) / dir%z
        pos = pos + dir * d
        r = sqrt(pos%x**2 + pos%y**2)
        if(r > this%radius)then
            skip=.true.
            return
        end if

        ! call u%push(pos)

        flag = .false.
        call reflect_refract(dir, this%flatNormal, this%n1, this%n2, flag)

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
        call reflect_refract(dir, curvedNormal, this%n2, this%n1, flag)
        if(flag)then
            skip = .true.
            return
        end if

    end subroutine plano_forward_sub

    subroutine plano_backward_sub(this, pos, dir, u, skip)

        use stackMod, only : stack

        implicit none

        class(plano_convex) :: this

        type(vector), intent(INOUT) :: pos, dir
        type(stack),  intent(INOUT) :: u
        logical,      intent(OUT) :: skip

        type(vector) :: curvedNormal
        real         :: d, t, r
        logical      :: flag

        skip = .false.

        flag = intersect_sphere(pos, dir, t, this%centre, this%curve_radius)
        if(.not. flag)then
            skip = .true.
            return
        end if

        pos = pos + t * dir

        flag = .false.
        curvedNormal = pos - this%centre
        curvedNormal = curvedNormal%magnitude()
        call reflect_refract(dir, curvedNormal, this%n1, this%n2, flag)

        ! move to flat surface
        d = (this%thickness) / dir%z
        pos = pos + dir * d
        r = sqrt(pos%x**2 + pos%y**2)
        if(r > this%radius)then
            skip=.true.
            return
        end if

        ! call u%push(pos)

        flag = .false.
        call reflect_refract(dir, this%flatNormal, this%n1, this%n2, flag)

    end subroutine plano_backward_sub


    subroutine doublet_forward_sub(this, pos, dir, u, skip, iris, iris_radius)

        use stackMod, only : stack

        implicit none

        class(achromatic_doublet) :: this

        type(vector), intent(INOUT) :: pos, dir
        type(stack),  intent(INOUT) :: u
        logical,optional,      intent(IN)    :: iris(2)
        real,optional,         intent(IN)    :: iris_radius
        logical,      intent(OUT)   :: skip

        type(vector) :: normal, origpos
        logical      :: flag
        real         :: t, r

        skip = .false.

        if(iris(1))then
            origpos = pos

            t = ((this%centre1%z - this%r1) - pos%z) / dir%z
            pos = pos + dir * t

            ! make sure no rays get propagated that are outside lens radius
            ! this can double as an Iris
            r = sqrt(pos%x**2 + pos%y**2)
            if(r > this%radius*iris_radius)then
                skip=.true.
                return
            end if
            pos = origpos
        end if

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
        if(r > (this%radius*1.0))then
            skip=.true.
            return
        end if

        normal = pos - this%centre1
        normal = normal%magnitude()

        flag = .false.
        call reflect_refract(dir, normal, this%n1, this%n2, flag)
        if(flag)then
            skip = .true.
            return
        end if

        ! call u%push(pos)

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
        call reflect_refract(dir, normal,this%n2, this%n3, flag)
        if(flag)then
            skip = .true.
            return
        end if


        ! call u%push(pos)

        !third sphere
        flag = intersect_sphere(pos, dir, t, this%centre3, this%R3)
        if(.not. flag)error stop "Help3"
        pos = pos + t * dir

        normal = this%centre3 - pos
        normal = normal%magnitude()

        flag = .false.
        call reflect_refract(dir, normal, this%n3, this%n1, flag)
        if(flag)then
            skip = .true.
            return
        end if


        ! call u%push(pos)
        if(iris(2))then
            origpos = pos
    
            t = ((this%centre3%z + this%r3) - pos%z) / dir%z
            pos = pos + dir * t
    
            r = sqrt(pos%x**2 + pos%y**2)
            if(r > this%radius*iris_radius)then
                skip=.true.
                return
            end if
            pos = origpos
        end if
    end subroutine doublet_forward_sub

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

    real function cauchy(wave, a, b, c)
    !Cauchy equation for refractive index of alcohol
    !wave is in nm

        implicit none

        real, intent(IN) :: wave, a, b, c

        real :: wavetmp

        wavetmp = wave*1d6
        cauchy = a + b *wavetmp**(-2) + c*wavetmp**(-4)

    end function cauchy

    real function dispersion(wave, a, b, c)
    !dispersion equation for soda-lime glass
    !wave is in nm
        implicit none

        real, intent(IN) :: wave, a, b, c

        real :: wave2
        !convert to um and square
        wave2 = (wave*1d6)**2

        dispersion = a - b*wave2 + (c / wave2)

    end function dispersion
end module lensMod