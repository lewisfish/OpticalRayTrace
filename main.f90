module imageMod

    implicit none

    contains

    subroutine makeImage(image, pos, layer)

        use vector_class

        implicit none

        integer,      intent(INOUT) :: image(-100:100,-100:100, 2)
        integer,      intent(IN)    :: layer
        type(vector), intent(IN)    :: pos

        real :: binwid
        integer :: xp, yp

        binwid = 1.d-2 / 201.

        xp = floor(pos%x / binwid)
        yp = floor(pos%y / binwid)
        if(abs(xp) > 100 .or. abs(yp) > 100)then
            return
        end if
        !$omp atomic
        image(xp,yp, layer) = image(xp,yp, layer) + 1

    end subroutine makeImage

end module imageMod

program raytrace

    use lensMod, only : plano_convex, achromatic_doublet
    use source, only : ring, point
    use vector_class, only : vector
    use imageMod

    implicit none
    

    type(vector) :: pos, dir
    integer :: image(-100:100, -100:100, 2)
    integer :: nphotons, i, iseed, u, count
    real :: total_length, d, D1, D2, D3, angle, cosThetaMax, r1, r2
    real :: L2Diameter, L2radius, L2Thickness, L2CurveRdaius, L1r1, L1r2, L1r3, L1tc1, L1tc2,L1radius
    logical :: skip

    iseed = 42
    nphotons = 20000000
    image = 0
    count = 0
    D1 = 40d-3 ! distance to lens 2 from centre of bottle
    D2 = 90d-3 ! distance between lenses
    D3 = 50d-3 ! distance from fibre to lens 1
    total_length = D1 + D2 + D3

    !lens 1 properties
    L1r1 = 33.55d-3
    L1r2 = 27.05d-3
    L1r3 = 125.60d-3
    L1tc1 = 7.5d-3
    L1tc2 = 1.8d-3
    L1radius = 12.7d-3

    !lens 2 properties - plano-convex
    L2Diameter = 25.4d-3
    L2radius = L2Diameter / 2.
    L2Thickness = 6.4d-3
    L2CurveRdaius = 20.6d-3

    angle = atan(L2radius / D1)
    cosThetaMax = cos(angle)

    ! open(newunit=u,file="ring-rays.dat")
    !annulus radii, squared as this is required for sampling
    r1 = 1d-3
    r2 = (r1 + .5d-3)**2
    r1 = r1**2
!$OMP parallel shared(L1radius, nphotons, cosThetaMax, d1, d2, d3, L2Thickness, r1, r2, L2CurveRdaius,image),&
!$omp& private(iseed, d, pos, dir) 
!$OMP do
    do i = 1, nphotons
        if(mod(i, 100000) == 0)print*,i,"photons run from ring"
        call ring(pos, dir, r1, r2, cosThetaMax, iseed)
        ! write(u,*)pos
        call plano_convex(pos, dir, D1, L2Thickness, L2CurveRdaius, u, iseed)
        ! write(u,*)pos
        call achromatic_doublet(pos, dir, D1, D2, L1r1, L1r2, L1r3, L1tc1, L1tc2, L1radius, u, iseed, skip)
        if(skip)then
            count = count + 1
            ! write(u,*)
            ! write(u,*)
            ! write(u,*)
            cycle
        end if
        ! write(u,*)pos

        d = (D1+d2+D3 - pos%z) / dir%z
        pos = pos + dir * d
        ! write(u,*)pos
        ! write(u,*)
        ! write(u,*)
        ! write(u,*)
        call makeImage(image, pos, 1)
    end do
!$OMP end do
! close(u)
! open(newunit=u,file="point-rays.dat",status="unknown")
print*,count/real(nphotons),count
! stop
!$omp do
    do i = 1, nphotons
        if(mod(i, 100000) == 0)print*,i,"photons run from point"
        call point(pos, dir, cosThetaMax, iseed)
        ! write(u,*)pos
        call plano_convex(pos, dir, D1, L2Thickness, L2CurveRdaius, u, iseed)
        ! write(u,*)pos
        call achromatic_doublet(pos, dir, D1, D2, L1r1, L1r2, L1r3, L1tc1, L1tc2, L1radius, u, iseed, skip)
        ! write(u,*)pos
        d = (D1+D2+D3 - pos%z) / dir%z
        pos = pos + dir * d
        ! write(u,*)pos
        ! write(u,*)" "
        ! write(u,*)" "
        ! write(u,*)" "

        call makeImage(image, pos, 2)
    end do
!$omp end parallel
! print*,count/real(nphotons)
! close(u)


    open(newunit=u,file="image-ring.dat", access="stream", form="unformatted",status="replace")
    write(u)real(image(:,:,1))
    close(u)

    open(newunit=u,file="image-dot.dat", access="stream", form="unformatted",status="replace")
    write(u)real(image(:,:,2))
    close(u)

    open(newunit=u,file="image-total.dat", access="stream", form="unformatted",status="replace")
    write(u)real(image(:,:,1) + image(:,:,2))
    close(u)

end program raytrace