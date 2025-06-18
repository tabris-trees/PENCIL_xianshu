!  Initial condition (density, magnetic field, velocity)
!  for a field from the PIC code(EPOCH).
!
!  04-sep-14/simon: coded
!  2025/05/06/xianshu: modified
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
module InitialCondition
!
    use Cparam
    use Cdata
    use General, only: keep_compiler_quiet
    use Messages
!
    implicit none
!
    include '../initial_condition.h'
!
    character(len=50) :: BxFile = 'bx2d_test.dat'
    character(len=50) :: ByFile = 'by2d_test.dat'
    character(len=50) :: BzFile = 'bz2d_test.dat'
    character(len=50) :: NFile = 'n2d_test.dat'
    character(len=50) :: UxFile = 'ux2d_test.dat'
    character(len=50) :: UyFile = 'uy2d_test.dat'
    character(len=50) :: UzFile = 'uz2d_test.dat'
    character(len=50) :: BxmFile = 'bx2d_m1_test.dat'
    character(len=50) :: BxpFile = 'bx2d_p1_test.dat'
!   real :: B_bkg = 0.0
!
    namelist /initial_condition_pars/ &
        BxFile, ByFile, BzFile, NFile
!
contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  07-oct-09/wlad: coded
!
!  Identify CVS/SVN version information.
!
        if (lroot) call svn_id( &
            "$Id: iucaa_logo.f90 19193 2012-06-30 12:55:46Z wdobler $")
!
    end subroutine register_initial_condition
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  07-may-09/wlad: coded
!
        real, dimension(mx, my, mz, mfarray), intent(inout) :: f

        ! real, dimension(mx, my, mz, mfarray), intent(inout) :: f
        integer :: l
        real, dimension(nxgrid, nygrid, 4) :: bb2d_pg  ! dat field
!
        call readdat(bb2d_pg)

        ! do n = 1, nz
        do n = 1, nzgrid
          ! do m = 1, ny
          do m = 1, nygrid
            ! do l = 1, nx
            do l = 1, nxgrid
              if ((l > nx + nx*ipx + nghost .or. m > ny + ny*ipy + nghost .or. n > nz + nz*ipz + nghost .or. &
                  l < 1 + nx*ipx - nghost .or. m < 1 + ny*ipy - nghost .or. n < 1 + nz*ipz - nghost) .eqv. .false.) then
                f(l + nghost - nx*ipx, m + nghost - ny*ipy, n + nghost - nz*ipz, iuu:iuu+2) = log(bb2d_pg(l, m, 7:9))
              end if
            end do
          end do
        end do
!
        ! call keep_compiler_quiet(f)
!
    end subroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density. init_lnrho
!  will take care of converting it to linear
!  density if you use ldensity_nolog
!
!  07-may-09/wlad: coded
!
        real, dimension(mx, my, mz, mfarray), intent(inout) :: f
        integer :: l
        real, dimension(nxgrid, nygrid, 4) :: bb2d_pg  ! dat field
!
        call readdat(bb2d_pg)

        ! do n = 1, nz
        do n = 1, nzgrid
          ! do m = 1, ny
          do m = 1, nygrid
            ! do l = 1, nx
            do l = 1, nxgrid
              if ((l > nx + nx*ipx + nghost .or. m > ny + ny*ipy + nghost .or. n > nz + nz*ipz + nghost .or. &
                  l < 1 + nx*ipx - nghost .or. m < 1 + ny*ipy - nghost .or. n < 1 + nz*ipz - nghost) .eqv. .false.) then
                f(l + nghost - nx*ipx, m + nghost - ny*ipy, n + nghost - nz*ipz, ilnrho) = log(bb2d_pg(l, m, 4))
              end if
            end do
          end do
        end do

        ! f(:,:,:,ilnrho) = log(f(:,:,:,irho))

!
    end subroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential.
!
!  04-sep-14/simon: coded
!  08-may-25/xianshu: modified
!
!  Load the magnetic field from three .dat file created by PIC code
!  (convert from .sdf file by MATLAB script)
!  The .dat file contains the magnetic field in the form of a 3D array
!  and compute the magnetic vector potential from it.
!
!  Created 2014-09-04 by Simon Candelaresi (Iomsn)
!
        use Poisson
        use Sub

        real, dimension(mx, my, mz, mfarray) :: f
        integer :: pos, l, j, ju
        ! The next 2 variables are used for the uncurling.
        real, dimension(nx, ny, nz, 3) :: jj, tmpJ  ! This is phi for poisson.f90
        ! real, dimension(nx, ny, 4) :: bb2d_pg  ! dat field
        real, dimension(nxgrid, nygrid, 4) :: bb2d_pg  ! dat field


! read the magnetic field
        call readdat(bb2d_pg)

        ! do n = 1, nz
        do n = 1, nzgrid
          ! do m = 1, ny
          do m = 1, nygrid
            ! do l = 1, nx
            do l = 1, nxgrid
              if ((l > nx + nx*ipx + nghost .or. m > ny + ny*ipy + nghost .or. n > nz + nz*ipz + nghost .or. &
                  l < 1 + nx*ipx - nghost .or. m < 1 + ny*ipy - nghost .or. n < 1 + nz*ipz - nghost) .eqv. .false.) then
                f(l + nghost - nx*ipx, m + nghost - ny*ipy, n + nghost - nz*ipz, iax:iaz) = bb2d_pg(l, m, 1:3)
              end if
            end do
          end do
        end do

        ! write (*, *) 'ibx = ', ibx
        ! write (*, *) 'test for f(1,1,1,ibx) = ', f(1+nghost, 1+nghost, 1+nghost, ibx)
        ! write (*, *) 'test for bb2d_pg(1,1,1) = ', bb2d_pg(1, 1, 1)
!
        ! deallocate (bb_pg)

!  Transform the magnetic field into a vector potential

!  Compute curl(B) = J for the Poisson solver
        do m = m1, m2
          do n = n1, n2
            call curl(f, iaa, jj(:, m - nghost, n - nghost, :))
          end do
        end do
        tmpJ = -jj
!  Use the Poisson solver to solve \nabla^2 A = -J for A
        do j = 1, 3
          call inverse_laplacian(tmpJ(:, :, :, j))
        end do

!  Overwrite the f-array with the correct vector potential A
        do j = 1, 3
          ju = iaa - 1 + j
          f(l1:l2, m1:m2, n1:n2, ju) = tmpJ(:, :, :, j)
        end do
!
!     Add a background field (no background field now)
        ! do l = 1, mx
        !   do m = 1, my
        !     f(l, m, :, iax) = f(l, m, :, iax) - y(m)*B_bkg/2.
        !     f(l, m, :, iay) = f(l, m, :, iay) + x(l)*B_bkg/2.
        !   end do
        ! end do
!
    end subroutine initial_condition_aa
!***********************************************************************
    subroutine read_initial_condition_pars(iostat)
!
        use File_io, only: parallel_unit
!
        integer, intent(out) :: iostat
!
        read (parallel_unit, NML=initial_condition_pars, IOSTAT=iostat)
!
    end subroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!
        integer, intent(in) :: unit
!
        write (unit, NML=initial_condition_pars)
!
    end subroutine write_initial_condition_pars
!***********************************************************************
    subroutine readdat(bb2d_pg)
! real, dimension(mx, my, mz, mfarray) :: f
        integer :: fileSize
        ! character(len=:), allocatable :: raw  ! contains the unformatted data
        ! real, allocatable :: bb(:, :, :, :)     ! magnetic field
        real, allocatable :: bb2d(:, :, :), bigmain(:, :)    ! temporary 2d magnetic field
        ! real*8, allocatable :: bb64(:,:,:,:)    ! for double precision
        ! real, dimension(nx, ny, 4) :: bb2d_pg
        real, dimension(nxgrid, nygrid, 4) :: bb2d_pg
        real :: valuemain, x1, y1, z1
        integer :: filenx, fileny, filenz
        integer :: unit, ios, ix, iy, iz, index, nx, ny, nz, index, idim
        ! integer :: nx, ny, nz
        ! integer :: pos, l, j, ju
        logical :: lexist
        ! The next 2 variables are used for the uncurling.
        ! real, dimension(nx, ny, nz, 3) :: jj, tmpJ  ! This is phi for poisson.f90
        character(len=100) :: outfile, infile
        integer :: unit_out

        !
        !  First one is the Bx field
        !
        ! check the length of the file and allocate array for the raw unformatted data
        ! inquire(file = BxFile, size = fileSize)
        ! nx = 512
        ! ny = 512
        ! nz = 512

        ! allocate (bb_pg(nx, ny, nz, 4)) ! has been allocated
        allocate (bigmain(filenx, fileny))

        do index = 1, 9

          select case (index)
          case (1)
            infile = ByFile
          case (2)
            infile = BzFile
          case (3)
            infile = BxFile
          case (4)
            infile = NFile
          case (5)
            infile = BxmFile
          case (6)
            infile = BxpFile
          case (7)
            infile = UxFile
          case (8)
            infile = UyFile
          case (9)
            infile = UzFile
          end select

          inquire (file=infile, exist=lexist)
          if (.not. lexist) then
            if (lroot) print *, 'ERROR: File ', trim(infile), ' not found.'
            ! call fatal_error('initial_condition_aa', 'Missing Bx2d.dat')
          end if

          unit = 10
          open (unit=unit, file=infile, status='old', action='read', iostat=ios)
          ! if (ios /= 0) call fatal_error('initial_condition_aa', 'Failed to open Bx2d.dat')

          read (unit, *) filenx, fileny
          filenz = filenx
          fileSize = filenx*fileny*filenz

          if (allocated(bb2d)) then
            continue
          else if (.not. allocated(bb2d)) then
            ! allocate (bb(filenx, fileny, filenz, 9))
            allocate (bb2d(filenx, fileny, 9))
          end if

          ! write (*, *) 'infile = ', infile, '  box = ', fileSize, '  ngrid', &
          !     '  filenx = ', filenx, '  fileny = ', fileny, '  filenz = ', filenz
          ! allocate(character(len = fileSize) :: raw)
          ! extract the parameters and data
          ! do iy = 1, fileny
          do ix = 1, filenx
            ! write (*, *) 'ix = ', ix
            read (unit, *) bb2d(ix, :, index)
          end do
          ! end do

          ! Map 2D Bx into 3D f-array (repeat in z)
          ! do iz = 1, filenx
          !   do iy = 1, fileny
          !     do ix = 1, filenz
          !       bb(ix, iy, iz, index) = bb2d(mod(ix - 1, filenx) + 1, mod(iy - 1, fileny) + 1, 1)
          !     end do
          !   end do
          ! end do

          ! write (*, *) 'bb', size(bb, dim=1), size(bb, dim=2), size(bb, dim=3), size(bb, dim=4)

          ! Interpolates to the grid specified in PENCIL
          bigmain = bb2d(:, :, index)

          ! do iz = 1, nz
          !   z1 = (real(iz) - 0.5)*real(filenz)/real(nz)
          ! do iy = 1, ny
          do iy = 1, nygrid
            ! y1 = (real(iy) - 0.5)*real(fileny)/real(ny)
            y1 = (real(iy) - 0.5)*real(fileny)/real(nygrid)
            ! do ix = 1, nx
            do ix = 1, nxgrid
              ! x1 = (real(ix) - 0.5)*real(filenx)/real(nx)
              x1 = (real(ix) - 0.5)*real(filenx)/real(nxgrid)
              call bilinear_interp(bigmain, real(x1), real(y1), valuemain)
              bb2d_pg(ix, iy, index) = valuemain
            end do
          end do
          ! end do
        end do

        ! write (*, *) 'bb2d_pg', size(bb2d_pg, dim=1), size(bb2d_pg, dim=2)

        close(unit)
        ! deallocate (bb)
        deallocate (bb2d)
        deallocate (bigmain)

    end subroutine readdat
!***********************************************************************
    subroutine trilinear_interp(big, xx, yy, zz, value)
        real, dimension(:, :, :) :: big
        real :: xx, yy, zz
        real :: value
        integer :: i, j, k
        real :: xd, yd, zd
        real :: c00, c01, c10, c11, c0, c1

        i = int(xx); xd = xx - i
        j = int(yy); yd = yy - j
        k = int(zz); zd = zz - k

        ! Border Protection
        if (i < 1) i = 1
        if (j < 1) j = 1
        if (k < 1) k = 1
        if (i >= size(big, 1)) i = size(big, 1) - 1
        if (j >= size(big, 2)) j = size(big, 2) - 1
        if (k >= size(big, 3)) k = size(big, 3) - 1

        ! Trilinear interpolation
        c00 = big(i, j, k)*(1 - xd) + big(i + 1, j, k)*xd
        c01 = big(i, j, k + 1)*(1 - xd) + big(i + 1, j, k + 1)*xd
        c10 = big(i, j + 1, k)*(1 - xd) + big(i + 1, j + 1, k)*xd
        c11 = big(i, j + 1, k + 1)*(1 - xd) + big(i + 1, j + 1, k + 1)*xd

        c0 = c00*(1 - yd) + c10*yd
        c1 = c01*(1 - yd) + c11*yd

        value = c0*(1 - zd) + c1*zd
    end subroutine trilinear_interp
!***********************************************************************
    subroutine bilinear_interp(big, xx, yy, value)
        real, dimension(:, :) :: big
        real :: xx, yy
        real :: value
        integer :: i, j
        real :: xd, yd
        real :: c0, c1

        ! Gets the integer index and decimal offset
        i = int(xx); xd = xx - i
        j = int(yy); yd = yy - j

        ! Border Protection
        if (i < 1) i = 1
        if (j < 1) j = 1
        if (i >= size(big, 1)) i = size(big, 1) - 1
        if (j >= size(big, 2)) j = size(big, 2) - 1

        ! Bilinear interpolation
        c0 = big(i, j)*(1 - xd) + big(i + 1, j)*xd
        c1 = big(i, j + 1)*(1 - xd) + big(i + 1, j + 1)*xd

        value = c0*(1 - yd) + c1*yd
    end subroutine bilinear_interp

!***********************************************************************
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noinitial_condition.f90 for any    **
!**  InitialCondition routines not implemented in this file        **
!**                                                                **
    include '../initial_condition_dummies.inc'
!********************************************************************
end module InitialCondition
