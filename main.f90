module imageMod

    implicit none

    contains

    subroutine makeImage(image, pos, diameter, layer)

        use vector_class

        implicit none

        integer,      intent(INOUT) :: image(-100:100,-100:100, 2)
        integer,      intent(IN)    :: layer
        real,         intent(IN)    :: diameter
        type(vector), intent(IN)    :: pos

        real :: binwid
        integer :: xp, yp

        binwid = diameter / 201.

        xp = floor(pos%x / binwid)
        yp = floor(pos%y / binwid)
        if(abs(xp) > 100 .or. abs(yp) > 100)then
            return
        end if
        !$omp atomic
        image(xp,yp, layer) = image(xp,yp, layer) + 1

    end subroutine makeImage


    subroutine writeImage(image, name)

        implicit none
        
        integer,      intent(IN) :: image(-100:100, -100:100, 2)
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

    use lensMod,      only : plano_convex, achromatic_doublet
    use source,       only : ring, point
    use vector_class, only : vector
    use constants,    only : pi
    use imageMod
    use OMP_LIB

    implicit none
    

    type(vector) :: pos, dir
    integer :: image(-100:100, -100:100, 2), l2image(-100:100, -100:100, 2), l3image(-100:100, -100:100, 2)
    integer :: nphotons, i, iseed, u, count
    real    :: total_length, d, D1, D2, D3, angle, cosThetaMax, r1, r2
    real    :: L2Diameter, L2radius, L2Thickness, L2CurveRdaius, L3r1, L3r2, L3r3, L3tc1, L3tc2,L3radius
    real    :: n, alpha, besselDiameter, distance, ringwidth, bottleRadius, L3fb
    logical :: skip

    iseed = 42
    nphotons = 20000000
    image = 0
    l2image = 0
    l3image = 0
    count = 0
    D1 = 37.5d-3 ! distance to lens 2 from centre of bottle
    D2 = 90d-3 ! distance between lenses
    D3 = 50d-3 ! distance from fibre to lens 1
    total_length = D1 + D2 + D3
    bottleRadius = 17.5d-3

    !lens 1 properties - achromatic doublet
    L3r1 = 33.55d-3     !radius 1
    L3r2 = 27.05d-3     !radius 2
    L3r3 = 125.60d-3    !radius 3
    L3tc1 = 7.5d-3      !thickness 1
    L3tc2 = 1.8d-3      !thickness2
    L3radius = 12.7d-3  !lens aperture
    L3fb = 45.0d-3      !back focal length

    !lens 2 properties - plano-convex
    L2Diameter = 25.4d-3
    L2radius = L2Diameter / 2.
    L2Thickness = 6.4d-3
    L2CurveRdaius = 20.6d-3

    angle = atan(L2radius / D1)
    cosThetaMax = cos(angle)

    ringWidth = 0.5d-3 ! half beam size incident on axicon
    alpha = 5.d0 * pi / 180.
    n = 1.45d0
    distance = 40d-3 - 17.5d-3
    !https://en.wikipedia.org/wiki/Axicon#/media/File:Erzeugen_von_Besselstrahlen_durch_ein_Axicon.png
    besselDiameter = 2.* distance * tan(alpha*(n-1.))
    !annulus radii, squared as this is required for sampling
    r1 = besselDiameter/2.d0 - ringwidth
    r2 = (besselDiameter / 2.d0)**2
    r1 = r1**2

!$OMP parallel default(none)&
!$OMP& shared(L3r1, L3r2, L3r3, L3fb, bottleradius, u, L3radius, nphotons, cosThetaMax, d1, d2, d3)&
!$OMP& shared(L2Thickness, r1, r2, L2CurveRdaius, image, L3tc1, L3tc2, L2Diameter, l2image, L3image),&
!$omp& private(iseed, d, pos, dir, skip), reduction(+:count)
!$OMP do
    do i = 1, nphotons
        if(mod(i, 1000000) == 0)print*,i,"photons run from ring"
        call ring(pos, dir, r1, r2, cosThetaMax, bottleRadius, iseed)
        call plano_convex(pos, dir, D1, L2Thickness, L2CurveRdaius, u, iseed)
        call makeImage(l2image, pos, L2Diameter, 1) 
        call achromatic_doublet(pos, dir, D1, D2, D3, L3r1, L3r2, L3r3, L3tc1, L3tc2, L3fb, L3radius, u, iseed, skip)
        call makeImage(l3image, pos, L3radius*2, 1)

        if(skip)then
!$OMP atomic
            count = count + 1
            cycle
        end if

        d = (D1+D2+D3 - pos%z) / dir%z
        pos = pos + dir * d
        call makeImage(image, pos, 1d-2, 1)
    end do
!$OMP end do
!$omp do
    do i = 1, nphotons
        if(mod(i, 1000000) == 0)print*,i,"photons run from point"
        call point(pos, dir, cosThetaMax, iseed)
        call plano_convex(pos, dir, D1, L2Thickness, L2CurveRdaius, u, iseed)
        call makeImage(l2image, pos, L2Diameter, 2)
        call achromatic_doublet(pos, dir, D1, D2, D3, L3r1, L3r2, L3r3, L3tc1, L3tc2, L3fb, L3radius, u, iseed, skip)
        call makeImage(l3image, pos, L3radius*2, 2)
        !move to image plane
        d = (D1+D2+D3 - pos%z) / dir%z
        pos = pos + dir * d
        call makeImage(image, pos, 1d-2, 2)
    end do
!$omp end parallel

print*,count/real(nphotons),count


call writeImage(l2image, "L2-")
call writeImage(l3image, "L3-")
call writeImage(image, "Final-")


end program raytrace