program raytrace

    use lensMod,      only : plano_convex, achromatic_doublet, glass_bottle
    use source,       only : ring, point, emit_image, create_spot, init_emit_image
    use stackMod,     only : stack
    use utils,        only : str, pbar
    use vector_class, only : vector
    use setup
    use imageMod
    use OMP_LIB
    use random, only : ran2
    use iso_fortran_env, only: int64

    implicit none
    
    type(vector) :: pos, dir
    type(stack)  :: tracker
    type(pbar)   :: bar

    !lens and bottles
    type(plano_convex)       :: L2
    type(achromatic_doublet) :: L3
    type(glass_bottle)       :: bottle

    integer        :: image(-200:200, -200:200, 2), i, uring, upoint
    integer(int64) :: rcount, pcount
    real           :: d, angle, cosThetaMax, r1, r2
    real           :: besselDiameter, distance
    logical        :: skip
    integer        :: imgin(512,512), nphotonsLocal
    character(len=:), allocatable :: filename

    image = 0
    rcount = 0_int64! counter for how many photons dont make from the ring source
    pcount = 0_int64! counter for how many photons dont make from the point source

    call setup_sim(L2, L3, bottle, imgin, nphotonsLocal)

    !max angle for point source to lens
    angle = atan(L2%radius / L2%fb)
    cosThetaMax = cos(angle)

    if(l2%fb <= bottle%radiusa + bottle%centre%z)error stop "Bottle offset too large!"
    ! distance to surface of bottle
    distance = L2%fb - (bottle%radiusa - bottle%centre%z)
    !https://en.wikipedia.org/wiki/Axicon#/media/File:Erzeugen_von_Besselstrahlen_durch_ein_Axicon.png
    besselDiameter = 2.* distance * tan(alpha*(n-1.))

    !annulus radii, squared as this is required for sampling
    r1 = besselDiameter/1.0d0 - ringWidth
    r2 = (besselDiameter / 2.d0)**2
    r1 = r1**2

    if(use_tracker)then
        open(newunit=uring, file=folder//"blah-delete.dat")
    end if
    
    bar = pbar(nphotons)

!$OMP parallel default(none)&
!$OMP& shared(L2, L3, bottle, uring, upoint, nphotons, cosThetaMax, r1, r2, image, wavelength)&
!$OMP& shared(use_tracker, use_bottle, image_source, point_source, spot_source, filename, folder)&
!$OMP& shared(source_type)&
!$omp& private(d, pos, dir, skip, tracker), firstprivate(imgin), reduction(+:rcount, pcount)
!$OMP do
    do i = 1, nphotons
        if(mod(i, 1000000) == 0)call bar%progress(i)
        skip=.false.
        ! if(mod(i, 10000000) == 0)print*,i,"photons run from ring"

        call ring(pos, dir, L2, r1, r2, bottle%Radiusa, bottle%radiusb, bottle%ellipse, bottle%centre%z)
        if(use_tracker)call tracker%push(pos)

        !propagate though lens 2
        call L2%forward(pos, dir, tracker, skip)
        if(use_tracker)call tracker%push(pos)
        
        if(skip)then
            rcount = rcount + 1_int64
            if(use_tracker)call tracker%write_empty(uring)
            cycle
        end if

        !propagate through lens 3
        call L3%forward(pos, dir, tracker, skip)
        if(use_tracker)call tracker%push(pos)

        if(skip)then
            rcount = rcount + 1_int64
            if(use_tracker)call tracker%write_empty(uring)
            cycle
        end if

        !move to image plane
        d = (L2%thickness+L2%fb +L2%fb+L3%fb+  L3%thickness+L3%fb- pos%z) / dir%z
        pos = pos + dir * d
        if(use_tracker)call tracker%push(pos)
        call makeImage(image, pos, 1d-2, 1)
        if(use_tracker)call tracker%write(uring)
    end do
!$OMP end do
if(use_tracker)close(uring)
! $OMP single
    ! wavelength = 843d-9
    ! L2 = plano_convex("../res/planoConvex.params", wavelength)
    ! L3 = achromatic_doublet("../res/achromaticDoublet.params", wavelength, L2%fb)
! $OMP end single
    !max angle for point source to lens
    ! angle = (13.*pi)/90.
    ! cosThetaMax = cos(angle)
    if(use_tracker)then
        filename = trim(adjustl(source_type))//"_track_bottle_"//str(use_bottle)//"_Ra_"// &
                   str(bottle%radiusa,7)//"_Rb_"//str(bottle%radiusb,7)//"_offset_"//&
                   str(bottle%centre%z,7)

        open(newunit=upoint, file=folder//filename//".dat")
        deallocate(filename)
    end if

!$omp do
    do i=1, nphotons
        skip=.false.

        ! if(mod(i, 10000000) == 0)print*,i,"photons run from point"

        if(image_source)then
            call emit_image(imgin, pos, dir, L2)
        elseif(point_source)then
            call point(pos, dir, cosThetaMax)
        elseif(spot_source)then
            call create_spot(pos, dir, cosThetaMax, nphotons, i)
        end if
        if(use_tracker)call tracker%push(pos)

        if(use_bottle)then
            call bottle%forward(pos, dir, tracker, skip)
            if(use_tracker)call tracker%push(pos)
        end if

        call L2%forward(pos, dir, tracker, skip)
        if(use_tracker)call tracker%push(pos)
        
        if(skip)then
            pcount = pcount + 1_int64
            if(use_tracker)call tracker%write_empty(upoint)
            cycle
        end if

        call L3%forward(pos, dir, tracker, skip)
        if(use_tracker)call tracker%push(pos)

        if(skip)then
            pcount = pcount + 1_int64
            if(use_tracker)call tracker%write_empty(upoint)
            cycle
        end if

        !move to image plane
        d = (L2%thickness+L2%fb +L2%fb+L3%fb+  L3%thickness+L3%fb- pos%z) / dir%z

        pos = pos + dir * d
        if(use_tracker)call tracker%push(pos)
        call makeImage(image, pos, 1d-2, 2)
        if(use_tracker)call tracker%write(upoint)
    end do
!$OMP end do
!$omp end parallel

if(use_tracker)close(upoint)

print"(A,1X,f8.2,A)","Ring  transmitted: ",100.*(1.-(rcount/(real(nphotons)))),"%"
print"(A,1X,f8.2,A)","Point transmitted: ",100.*(1.-(pcount/(real(nphotons)))),"%"

if(makeImages)then
    filename = trim(adjustl(source_type))//"_image_bottle_"//str(use_bottle)//"_Ra_"// &
                   str(bottle%radiusa,7)//"_Rb_"//str(bottle%radiusb,7)//"_offset_"//&
                   str(bottle%centre%z,7)
    call writeImage(image, filename)
end if

end program raytrace