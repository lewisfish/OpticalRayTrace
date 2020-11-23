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
    use lensMod,      only : plano_convex, achromatic_doublet
    use source,       only : ring, point
    use stackMod,     only : pointtype, stack
    use utils,        only : str
    use vector_class, only : vector
    use imageMod
    use OMP_LIB

    implicit none
    
    type(vector) :: pos, dir
    type(stack)  :: tracker
    type(plano_convex) :: L2
    type(achromatic_doublet) :: L3
    integer :: image(-200:200, -200:200, 2)
    integer :: nphotons, i, iseed, u, count
    real    :: d, angle, cosThetaMax, r1, r2
    real    :: n, alpha, besselDiameter, distance, ringwidth, bottleRadius, bottleOffset
    logical :: skip

    iseed = 42
    nphotons = 2000000000
    image = 0
    count = 0

    bottleRadius = 17.5d-3

    L2 = plano_convex("planoConvex.params")
    L3 = achromatic_doublet("achromaticDoublet.params", L2%fb)

    !max angle for point source to lens
    angle = atan(L2%radius / L2%fb)
    cosThetaMax = cos(angle)

    !size of bessel annulus 
    ringWidth = 0.5d-3 ! half beam size incident on axicon -> guess is 1mm in diameter
    alpha = 5.d0 * pi / 180.
    n = 1.45d0
    bottleradius = 17.5d-3
    bottleOffset = 0.d0!17.4d-3!8.75d-3!17.4d-3!0.*3.5d-3!bottleradius / 2.d0
    ! distance to surface of bottle
    distance = L2%fb - (bottleradius - bottleOffset)
    !https://en.wikipedia.org/wiki/Axicon#/media/File:Erzeugen_von_Besselstrahlen_durch_ein_Axicon.png
    besselDiameter = 2.* distance * tan(alpha*(n-1.))
    !annulus radii, squared as this is required for sampling
    r1 = besselDiameter/1.1d0 - ringWidth
    r2 = (besselDiameter / 2.d0)**2
    r1 = r1**2
    ! open(newunit=u, file="tmp.dat", position="append")

!$OMP parallel default(none)&
!$OMP& shared(L2, L3, bottleradius, u, nphotons, cosThetaMax, bs1, bs2, imagel31)&
!$OMP& shared(r1, r2, image, l2image, L3image, bottleOffset),&
!$omp& private(iseed, d, pos, dir, skip, tracker), reduction(+:count)
!$OMP do
    do i = 1, nphotons
        skip=.false.
        if(mod(i, 1000000) == 0)print*,i,"photons run from ring"
        call ring(pos, dir, r1, r2, cosThetaMax, bottleRadius, bottleOffset, iseed)
        ! call tracker%push(pointtype(pos%x, pos%y, pos%z))

        !propagate though lens 2
        call L2%forward(pos, dir, tracker, iseed)
        ! call tracker%push(pointtype(pos%x, pos%y, pos%z))

        !propagate through lens 3
        call L3%forward(pos, dir, L2%fb, L2%fb+L3%f, tracker, iseed, skip)
        ! call tracker%push(pointtype(pos%x, pos%y, pos%z))

        if(skip)then
            count = count + 1
            cycle
        end if

        !move to image plane
        d = (L2%f+L2%fb+2.*L3%f - pos%z) / dir%z
        pos = pos + dir * d
        call makeImage(image, pos, 1d-2, 1)
    end do
!$OMP end do

    ! open(newunit=u, file="traj-annuli-4-point.dat", position="append")
!$omp do
    do i=1, nphotons
        skip=.false.
        if(mod(i, 1000000) == 0)print*,i,"photons run from point"

        call point(pos, dir, cosThetaMax, iseed)
        ! call tracker%push(pointtype(pos%x, pos%y, pos%z))

        call L2%forward(pos, dir, tracker, iseed)
        ! call tracker%push(pointtype(pos%x, pos%y, pos%z))

        call L3%forward(pos, dir, L2%fb, L2%fb+L3%f, tracker, iseed, skip)
        ! call tracker%push(pointtype(pos%x, pos%y, pos%z))

        if(skip)then
            cycle
        end if

        !move to image plane
        d = (L2%f+L2%fb+2.*L3%f - pos%z) / dir%z
        pos = pos + dir * d
        ! call tracker%push(pointtype(pos%x, pos%y, pos%z))
        call makeImage(image, pos, 1d-2, 2)
    end do
!$OMP end do
!$omp end parallel
! close(u)
print*,count/real(nphotons),count

call writeImage(image, "test-new-")

end program raytrace