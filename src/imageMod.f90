module imageMod

    use vector_class

    implicit none

    interface writeImage
        module procedure writeImage2D
        module procedure writeImage3D
    end interface

    interface makeImage
        module procedure makeImage2D
        module procedure makeImage3D
    end interface

    contains

    subroutine makeImage2D(image, dir, pos, diameter, layer)

        implicit none

        integer,      intent(INOUT) :: image(-200:200,-200:200, 2)
        integer,      intent(IN)    :: layer
        real,         intent(IN)    :: diameter
        type(vector), intent(IN)    :: pos, dir

        real :: binwid, angle, na, bottom, top
        type(vector) :: n, d
        integer :: xp, yp

        n = vector(0., 0., -1.)
        n = n%magnitude()
        d = dir%magnitude()
        d = (-1.)*d

        top = n .dot. d
        bottom = sqrt(d .dot. d) * sqrt(n .dot. n)
        angle = acos(top / bottom)
        na = asin(0.22)

        if(angle > na)then
            return
        end if
        binwid = diameter / 401.

        !if pos is large then return as wint hit detector
        if(pos%x > 1000 .or. pos%y > 1000)return

        xp = floor(pos%x / binwid)
        yp = floor(pos%y / binwid)
        if(abs(xp) > 200 .or. abs(yp) > 200)then
            return
        end if
!$omp atomic
        image(xp,yp, layer) = image(xp,yp, layer) + 1

    end subroutine makeImage2D


    subroutine makeImage3D(image, dir, pos, diameter, layer)
        
        implicit none
        
        integer,      intent(INOUT) :: image(-200:200, -200:200, 200, 2)
        integer,      intent(IN)    :: layer
        real,         intent(IN)    :: diameter
        type(vector), intent(IN)    :: pos, dir
        
        real :: binwid, dz
        integer :: xp, yp, zp, i
        type(vector) :: new_pos

        binwid = diameter / 401.
        dz = diameter / 200.

        do i = 0, 199

            new_pos = pos + (i*dz)*dir
            xp = floor(new_pos%x / binwid)
            yp = floor(new_pos%y / binwid)

            if(abs(xp) > 200 .or. abs(yp) > 200)then
                return
            end if
!$omp atomic
            image(xp, yp, i+1, layer) = image(xp, yp, i+1, layer) + 1
        end do

    end subroutine makeImage3D


    subroutine writeImage2D(image, name)

        implicit none
        
        integer,      intent(IN) :: image(-200:200, -200:200, 2)
        character(*), intent(IN) :: name

        integer :: u

        open(newunit=u,file=name//"-ring.dat", access="stream", form="unformatted",status="replace")
        write(u)real(image(:,:,1))
        close(u)

        open(newunit=u,file=name//"-point.dat", access="stream", form="unformatted",status="replace")
        write(u)real(image(:,:,2))
        close(u)

        open(newunit=u,file=name//"-total.dat", access="stream", form="unformatted",status="replace")
        write(u)real(image(:,:,1)) + real(image(:,:,2))
        close(u)

    end subroutine writeImage2D


    subroutine writeImage3D(image, name)

        implicit none
        
        integer,      intent(IN) :: image(-200:200, -200:200, 200, 2)
        character(*), intent(IN) :: name

        integer :: u

        open(newunit=u,file=name//"-vol-ring.dat", access="stream", form="unformatted",status="replace")
        write(u)real(image(:,:,:,1))
        close(u)

        open(newunit=u,file=name//"-vol-point.dat", access="stream", form="unformatted",status="replace")
        write(u)real(image(:,:,:,2))
        close(u)
    end subroutine writeImage3D

end module imageMod