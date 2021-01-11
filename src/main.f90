module imageMod

    implicit none

    contains

    subroutine makeImage(image, pos, diameter, layer)

        use vector_class

        implicit none

        integer,      intent(INOUT) :: image(-200:200,-200:200, 2)
        integer,      intent(IN)    :: layer
        real,         intent(IN)    :: diameter
        type(vector), intent(IN)    :: pos

        real :: binwid
        integer :: xp, yp

        binwid = diameter / 401.

        xp = floor(pos%x / binwid)
        yp = floor(pos%y / binwid)
        if(abs(xp) > 200 .or. abs(yp) > 200)then
            return
        end if
!$omp atomic
        image(xp,yp, layer) = image(xp,yp, layer) + 1

    end subroutine makeImage


    subroutine writeImage(image, name)

        implicit none
        
        integer,      intent(IN) :: image(-200:200, -200:200, 2)
        character(*), intent(IN) :: name

        integer :: u

        open(newunit=u,file=name//"ring.dat", access="stream", form="unformatted",status="replace")
        write(u)real(image(:,:,1))
        close(u)

        open(newunit=u,file=name//"point.dat", access="stream", form="unformatted",status="replace")
        write(u)real(image(:,:,2))
        close(u)

        open(newunit=u,file=name//"total.dat", access="stream", form="unformatted",status="replace")
        write(u)real(image(:,:,1)) + real(image(:,:,2))
        close(u)


    end subroutine writeImage

end module imageMod

program raytrace

    use constants,    only : pi
    use lensMod,      only : plano_convex, achromatic_doublet, glass_bottle
    use source,       only : ring, point, emit_image, create_spot
    use stackMod,     only : stack
    use utils,        only : str
    use vector_class, only : vector
    use imageMod
    use OMP_LIB
    use random, only : ran2
    use iso_fortran_env, only: int64

    implicit none
    

    type(vector) :: pos, dir
    type(stack)  :: tracker
    
    !lens and bottles
    type(plano_convex)       :: L2
    type(achromatic_doublet) :: L3
    type(glass_bottle)       :: bottle

    integer        :: image(-200:200, -200:200, 2), nphotons, i, uring, upoint, u
    integer(int64) :: rcount, pcount, nphotonsLocal
    real           :: d, angle, cosThetaMax, r1, r2, wavelength
    real           :: n, alpha, besselDiameter, distance, ringwidth
    logical        :: skip
    integer        :: imgin(101,101)
    character(len=256) :: filename

    image = 0
    rcount = 0_int64! counter for how many photons dont make from the ring source
    pcount = 0_int64! counter for how many photons dont make from the point source

    open(newunit=u, file="../res/settings.params", status="old")
        read(u,*) ringWidth
        read(u,*) wavelength
        read(u,*) nphotons
        read(u,*) alpha
        read(u,*) n
        !read in params files
        read(u,*) filename
        bottle = glass_bottle("../res/"//trim(filename), wavelength)
        read(u,*) filename
        L2 = plano_convex("../res/"//trim(filename), wavelength)
        read(u,*) filename
        L3 = achromatic_doublet("../res/"//trim(filename), wavelength, L2%fb, L2%thickness)
        read(u,*) filename
        call init_emit_image("../res/"//trim(filename), imgin, nphotons, nphotonsLocal)
    close(u)

    alpha = alpha * pi / 180.
    !max angle for point source to lens
    angle = atan(L2%radius / L2%fb)
    cosThetaMax = cos(angle)

    if(l2%fb <= bottle%radiusa + bottle%centre%z)error stop "Bottle offset too large!"
    ! distance to surface of bottle
    distance = L2%fb - (bottle%radius - bottle%centre%z)
    !https://en.wikipedia.org/wiki/Axicon#/media/File:Erzeugen_von_Besselstrahlen_durch_ein_Axicon.png
    besselDiameter = 2.* distance * tan(alpha*(n-1.))

    !annulus radii, squared as this is required for sampling
    r1 = besselDiameter/1.1d0 - ringWidth
    r2 = (besselDiameter / 2.d0)**2
    r1 = r1**2
    ! open(newunit=uring, file="test/ring-smallf-rays.dat", position="append")

!$OMP parallel default(none)&
!$OMP& shared(L2, L3, bottle, uring, upoint, nphotons, cosThetaMax, r1, r2, image, wavelength)&
!$omp& private(d, pos, dir, skip, tracker), firstprivate(imgin), reduction(+:rcount, pcount)
!$OMP do
    do i = 1, 1

        skip=.false.
        if(mod(i, 10000000) == 0)print*,i,"photons run from ring"
        
        call ring(pos, dir, r1, r2, bottle%Radiusa, bottle%radiusb, bottle%ellipse, bottle%centre%z)
        ! call tracker%push(pos)

        !propagate though lens 2
        call L2%forward(pos, dir, tracker, skip)
        ! call tracker%push(pos)
        
        if(skip)then
            rcount = rcount + 1_int64
            ! call tracker%zero()
            ! write(uring,*)" "
            ! write(uring,*)" "
            ! write(uring,*)" "
            cycle
        end if

        !propagate through lens 3
        call L3%forward(pos, dir, tracker, skip)
        ! call tracker%push(pos)

        if(skip)then
            rcount = rcount + 1_int64
            ! call tracker%zero()
            ! write(uring,*)" "
            ! write(uring,*)" "
            ! write(uring,*)" "
            cycle
        end if

        !move to image plane
        d = (L2%thickness+L2%fb +L2%fb+L3%fb+  L3%thickness+L3%fb- pos%z) / dir%z
        pos = pos + dir * d
        ! call tracker%push(pos)
        call makeImage(image, pos, 1d-2, 1)
        ! do while(.not. tracker%empty())
        !     write(uring,"(3(F10.7,1x))")tracker%pop()
        ! end do
        ! write(uring,*)" "
        ! write(uring,*)" "
        ! write(uring,*)" "
    end do
!$OMP end do
! close(uring)
! $OMP single
    ! wavelength = 843d-9
    ! L2 = plano_convex("../res/planoConvex.params", wavelength)
    ! L3 = achromatic_doublet("../res/achromaticDoublet.params", wavelength, L2%fb)
! $OMP end single
    !max angle for point source to lens
    ! angle = (13.*pi)/90.
    ! cosThetaMax = cos(angle)

    open(newunit=upoint, file="../data/spot-diags/point-large-bottle-thick-785nm.dat")
!$omp do
    do i=1, nphotons
        skip=.false.
        if(mod(i, 10000000) == 0)print*,i,"photons run from point"

        call create_spot(pos, dir, cosThetaMax, nphotons, i)
        ! call emit_image(imgin, pos, dir)
        ! call point(pos, dir, cosThetaMax)
        ! call tracker%push(pos)

        ! call bottle%forward(pos, dir, tracker, skip)
        ! call tracker%push(pos)

        call L2%forward(pos, dir, tracker, skip)
        ! call tracker%push(pos)
        
        if(skip)then
            pcount = pcount + 1_int64
            call tracker%zero()
            write(upoint,*)" "
            write(upoint,*)" "
            write(upoint,*)" "
            cycle
        end if

        call L3%forward(pos, dir, tracker, skip)
        ! call tracker%push(pos)

        if(skip)then
            pcount = pcount + 1_int64
            call tracker%zero()
            write(upoint,*)" "
            write(upoint,*)" "
            write(upoint,*)" "
            cycle
        end if

        !move to image plane
        d = (L2%thickness+L2%fb +L2%fb+L3%fb+  L3%thickness+L3%fb- pos%z) / dir%z
        pos = pos + dir * d
        ! call tracker%push(pos)
        call makeImage(image, pos, 1d-2, 2)
        do while(.not. tracker%empty())
            write(upoint,"(3(F10.7,1x))")tracker%pop()
        end do
        write(upoint,*)" "
        write(upoint,*)" "
        write(upoint,*)" "
    end do
!$OMP end do
!$omp end parallel

close(upoint)

print"(A,1X,f8.2,A)","Ring  transmitted: ",100.*(1.-(rcount/(real(nphotons)))),"%"
print"(A,1X,f8.2,A)","Point transmitted: ",100.*(1.-(pcount/(real(nphotons)))),"%"

! call writeImage(image, "../test/img-emit-")

end program raytrace