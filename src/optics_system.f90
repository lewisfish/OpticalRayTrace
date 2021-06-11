module opticsystem
    implicit none

    contains

        subroutine telescope(pos, dir, L1, L2, img_plane, count, tracker, u, skip)
            
            use lensMod, only : plano_convex, achromatic_doublet, lens
            use stackMod, only : stack
            use vector_class
            use setup, only : iris, iris_radius, use_tracker, fibre_offset

            use iso_fortran_env, only: int64

            implicit none

            class(lens),              intent(IN)    :: L1, L2 
            integer,                  intent(IN)    :: u
            real,                     intent(IN)    :: img_plane
            type(stack),              intent(INOUT) :: tracker
            type(vector),             intent(INOUT) :: pos, dir
            integer(int64),           intent(INOUT) :: count
            logical,                  intent(INOUT) :: skip
            
            real :: d

            !propagate though lens 1
            call L1%forward(pos, dir, tracker, skip, iris, iris_radius)
            if(use_tracker)call tracker%push(pos)

            if(skip)then
                count = count + 1_int64
                if(use_tracker)call tracker%write_empty(u)
                return
            end if

            !propagate through lens 2
            call L2%forward(pos, dir, tracker, skip, iris, iris_radius)
            if(use_tracker)call tracker%push(pos)

            if(skip)then
                count = count + 1_int64
                if(use_tracker)call tracker%write_empty(u)
                return
            end if

            !move to image plane
            d = ((img_plane + fibre_offset) - pos%z) / dir%z
            pos = pos + dir * d
            if(use_tracker)call tracker%push(pos)
            
        end subroutine telescope
    
end module opticsystem