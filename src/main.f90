program raytrace

    use lensMod,      only : plano_convex, achromatic_doublet, glass_bottle
    use source,       only : ring, point, emit_image, create_spot, init_emit_image, cross, point_on_bottle, iSORS
    use stackMod,     only : stack
    use utils,        only : str, pbar
    use vector_class, only : vector
    use constants,    only : pi
    use setup
    use imageMod
    use opticsystem, only : telescope
    use OMP_LIB
    use random, only : ran2, init_rng
    use iso_fortran_env, only: int64

    implicit none
    
    type(vector) :: pos, dir
    type(stack)  :: tracker
    type(pbar)   :: bar

    !lens and bottles
    type(plano_convex)       :: L2
    type(achromatic_doublet) :: L3
    type(glass_bottle)       :: bottle

    integer, allocatable :: image(:, :, :), imgin(:,:)
    integer(int64) :: rcount, pcount
    real           :: angle, cosThetaMax, r1, r2
    real           :: besselDiameter, distance, img_plane_1
    logical        :: skip, file_exists
    integer        :: nphotonsLocal, i, uring, upoint
    character(len=:), allocatable :: filename

    allocate(image(-200:200, -200:200, 2), imgin(512, 512))


    file_exists = .false.
    image = 0
    rcount = 0_int64! counter for how many photons dont make from the ring source
    pcount = 0_int64! counter for how many photons dont make from the point source

    call setup_sim(L2, L3, bottle, imgin, nphotonsLocal)

    filename = trim(adjustl(source_type))//"_bottle_"//str(use_bottle)//"_Ra_"// &
            str(bottle%radiusa,7)//"_Rb_"//str(bottle%radiusb,7)//"_offset_"//&
            str(bottle%centre%z,7)//"_"//str(iris)//"_"//str(iris_radius, 7)//"_L2f_"//str(L2%f,6)//"_L3f_"//str(L3%f,6)//&
            "_fo_"//str(fibre_offset, 7)//"_alp_"//str(alpha*180/pi,7)//"_bwidth_"//str(ringWidth,7)//"_sep_"//str(isors_offset,7)

    !max angle for point source to lens
    angle = atan(L2%radius / L2%fb)
    cosThetaMax = cos(angle)

    if(l2%fb <= bottle%radiusa + bottle%centre%z)then
        print*,"Bottle offset too large! Adjusting so that there is a minimum of 2mm offset from lens."
        bottle%centre%z = L2%fb - bottle%radiusa - 2d-3
        print*,"Now bottle set at z position:",bottle%centre%z
    end if
    ! distance to surface of bottle
    if(isors_source)then
        distance = bottle%radiusa + isors_offset
    else
        distance = (bottle%radiusa + bottle%centre%z)
    end if
    !https://en.wikipedia.org/wiki/Axicon#/media/File:Erzeugen_von_Besselstrahlen_durch_ein_Axicon.png
    besselDiameter = distance*97.3d-3*tan(alpha* (n - 1)) /(l2%fb) ! 97.3d-3 is the fb of L1
    !annulus radii, squared as this is required for sampling
    r1 = besselDiameter - ringWidth
    r2 = (besselDiameter / 2.d0)**2
    r1 = r1**2

    if(use_tracker)then
        open(newunit=uring, file=folder//filename//"-ringtrace.dat")
    end if
    
    bar = pbar(nphotons/ 1000000)

    !make repeatable
    call init_rng(123456789)

    img_plane_1 = 2.*(L2%fb + L3%fb) + L2%thickness + L3%thickness

!$OMP parallel default(none)&
!$OMP& shared(L2, L3, bottle, uring, upoint, nphotons, cosThetaMax, r1, r2, image, wavelength)&
!$OMP& shared(use_tracker, use_bottle, image_source, point_source, spot_source, filename, folder)&
!$OMP& shared(source_type, bar, iris, iris_radius, image_diameter, fibre_offset, L2file, L3file)&
!$OMP& shared(isors_source, isors_offset, spot_size, ringWidth, img_plane_1, crs_source)&
!$omp& private(pos, dir, skip, tracker), firstprivate(imgin), reduction(+:rcount, pcount)
!$OMP do
    do i = 1, nphotons
        skip=.false.

        if(mod(i, 1000000) == 0)call bar%progress()

        if(isors_source)then
            call iSORS(pos, dir, bottle, L2, isors_offset, ringWidth, .true.)
        elseif(crs_source)then
            call point_on_bottle(pos, dir, cosThetaMax, bottle, spot_size)
        else
            call ring(pos, dir, L2, r1, r2, bottle%Radiusa, bottle%radiusb, bottle%ellipse, bottle%centre%z)
        end if

        if(use_tracker)call tracker%push(pos)
        call telescope(pos, dir, L2, L3, img_plane_1, rcount, tracker, uring, skip)
        if(skip)cycle

        if(use_tracker)call tracker%write(uring)
        call makeImage(image, dir, pos, image_diameter, 1)
    end do
!$OMP end do

if(use_tracker)close(uring)
!$OMP single
    wavelength = 843d-9
    L2 = plano_convex("../res/"//trim(L2file), wavelength)
    L3 = achromatic_doublet("../res/"//trim(L3file), wavelength, 2.*L2%fb+ L2%thickness)
!$OMP end single

    bar = pbar(nphotons/ 1000000)

    if(use_tracker)then
        open(newunit=upoint, file=folder//filename//"-pointtrace.dat")
        deallocate(filename)
    end if

!$omp do
    do i=1, nphotons
        skip=.false.

        if(mod(i, 1000000) == 0)call bar%progress()

        if(image_source)then
            call emit_image(imgin, pos, dir, L2)
        elseif(point_source .or. crs_source)then
            ! call cross(pos, dir)
            call point(pos, dir, cosThetaMax)
        elseif(spot_source)then
            call create_spot(pos, dir, cosThetaMax, nphotons, i)
        elseif(isors_source)then
            call point(pos, dir, cosThetaMax, bottle%centre%z)
            ! call iSORS(pos, dir, bottle, L2, isors_offset, ringWidth, .false.)
        end if

        if(use_tracker)call tracker%push(pos)
        if(use_bottle)then
            call bottle%forward(pos, dir, tracker, skip)
            if(use_tracker)call tracker%push(pos)
        end if

        if(skip)then
            pcount = pcount + 1_int64
            if(use_tracker)call tracker%write(upoint)
            if(use_tracker)call tracker%write_empty(upoint)
            cycle
        end if
        
        call telescope(pos, dir, L2, L3, img_plane_1, pcount, tracker, upoint, skip)
        if(skip)cycle

        if(use_tracker)call tracker%write(upoint)
        call makeImage(image, dir, pos, image_diameter, 2)
    end do
!$OMP end do
!$omp end parallel

if(use_tracker)close(upoint)

inquire(file=folder//"trans-stats.dat", exist=file_exists)
if(.not. file_exists)then
    open(newunit=upoint, file=folder//"trans-stats.dat")
    write(upoint, *)"r/%, p/%, l2%f, l3%f, bottle?, radiusA, radiusB, iris_pos, iris_radius, offset, source_type, seperation"
else
    open(newunit=upoint, file=folder//"trans-stats.dat", position="append")
end if
write(upoint,*)100.*(1.-(rcount/(real(nphotons)))),",", 100.*(1.-(pcount/(real(nphotons)))),",", L2%f,",",&
          L3%f,",", use_bottle,",", bottle%radiusa,",", bottle%radiusb,",", iris,",", str(iris_radius,7), ",", bottle%centre%z, &
          ",", trim(adjustl(source_type))//",", isors_offset
close(upoint)

print"(A,1X,f8.2,A)","Ring  transmitted: ",100.*(1.-(rcount/(real(nphotons)))),"%"
print"(A,1X,f8.2,A)","Point transmitted: ",100.*(1.-(pcount/(real(nphotons)))),"%"

if(makeImages)then
    call writeImage(image, folder//filename//"_image")
end if

end program raytrace