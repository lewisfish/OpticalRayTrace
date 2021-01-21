module setup

    use lensMod, only : plano_convex, achromatic_doublet, glass_bottle

    implicit none
        
    real    :: alpha, ringwidth, wavelength, n
    integer :: nphotons
    logical :: use_tracker, use_bottle, point_source, spot_source, image_source
    character(len=256) :: source_type

    contains
    
    subroutine setup_sim(L2, L3, bottle, image, nphotonsLocal)
        
        implicit none

        type(plano_convex),       intent(OUT) :: L2
        type(achromatic_doublet), intent(OUT) :: L3
        type(glass_bottle),       intent(OUT) :: bottle
        integer,                  intent(OUT) :: nphotonsLocal, image(:, :)

        call read_settings(L2, L3, bottle, image, nphotonsLocal)

    end subroutine setup_sim

    subroutine read_settings(L2, L3, bottle, image, nphotonsLocal)
        
        use constants, only : pi
        use source,    only : init_emit_image

        implicit none

        type(plano_convex),       intent(OUT) :: L2
        type(achromatic_doublet), intent(OUT) :: L3
        type(glass_bottle),       intent(OUT) :: bottle
        integer,                  intent(OUT) :: nphotonsLocal, image(:, :)

        integer            :: i, u
        character(len=256) :: arg, filename

        point_source = .false.
        spot_source  = .false.
        image_source = .false.

        ! read in settings file from command line
        i = 1
        call get_command_argument(i, arg)
        print*,"Using "//trim(arg)//" settings."

        ! read in settingds from file
        open(newunit=u, file="../res/"//trim(arg), status="old")
            read(u,*) ringWidth
            read(u,*) wavelength
            read(u,*) nphotons
            read(u,*) alpha
            alpha = alpha * pi / 180.
            read(u,*) n
            read(u,*) use_bottle
            read(u,*) use_tracker
            if(nphotons > 10000 .and. use_tracker)error stop "Too many photons for tracker use!"
            read(u,*) source_type

            if(trim(source_type) == "image")then
                image_source = .true.
            elseif(trim(source_type) == "spot")then
                spot_source = .true.
            elseif(trim(source_type) == "point")then
                point_source = .true.        
            else
                error stop "No such source type!"
            end if

            !read in params files
            read(u,*) filename
            bottle = glass_bottle("../res/"//trim(filename), wavelength)
            read(u,*) filename
            L2 = plano_convex("../res/"//trim(filename), wavelength)
            read(u,*) filename
            L3 = achromatic_doublet("../res/"//trim(filename), wavelength, L2%fb, L2%thickness)
            read(u,*) filename
            call init_emit_image("../res/"//trim(filename), image, nphotons, nphotonsLocal)
        close(u)    

    end subroutine read_settings

end module setup