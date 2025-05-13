! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM integer, parameter :: nghost = 1
!
!***************************************************************
module Deriv
!
  use Messages
  use Cdata
  use General, only: keep_compiler_quiet
!
  implicit none
!
  private
!
  include 'deriv.h'
!
  contains
!
!***********************************************************************
    subroutine initialize_deriv()
!
!  Initialize stencil coefficients (dummy routine)
!
    endsubroutine initialize_deriv
!***********************************************************************
    subroutine calc_coeffs_1( grid, coeffs )
!
!  dummy
!
      !real, dimension(-2:3), intent(in ) :: grid
      !real, dimension(-3:3), intent(out) :: coeffs
      real, dimension(-0:1), intent(in ) :: grid
      real, dimension(-1:1), intent(out) :: coeffs
!
      if (lroot) print*,'calc_coeffs_1 is not evaluated'
!--   call fatal_error("calc_coeffs_1","not coded for deriv_2nd")
!
  endsubroutine calc_coeffs_1
!***********************************************************************
    subroutine der_main(f,k,df,j,ignoredx)
!
!  calculate derivative df_k/dx_j
!  accurate to 6th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!
!   1-oct-97/axel: coded
!  18-jul-98/axel: corrected mx -> my and mx -> mz in all y and z ders
!   1-apr-01/axel+wolf: pencil formulation
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!  21-feb-07/axel: added 1/r and 1/pomega factors for non-coord basis
!  25-aug-09/axel: adapted from deriv
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df,fac
      logical, intent(in), optional :: ignoredx
      integer :: j,k
!
      intent(in)  :: f,k,j
      intent(out) :: df
!
!debug      if (loptimise_ders) der_call_count(k,icount_der,j,1) = & !DERCOUNT
!debug                            der_call_count(k,icount_der,j,1)+1 !DERCOUNT
!
      if (present(ignoredx)) call fatal_error('der_main', 'optional argument ignoredx is not implemented. ')
!
      if (j==1) then
        if (nxgrid/=1) then
          fac=.5*dx_1(l1:l2)
          df=fac*(f(l1+1:l2+1,m,n,k)-f(l1-1:l2-1,m,n,k))
        else
          df=0.
          if (ip<=5) print*, 'der_main: Degenerate case in x-direction'
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          fac=.5*dy_1(m)
          df=fac*(f(l1:l2,m+1,n,k)-f(l1:l2,m-1,n,k))
          if (lspherical_coords)     df=df*r1_mn
          if (lcylindrical_coords)   df=df*rcyl_mn1
        else
          df=0.
          if (ip<=5) print*, 'der_main: Degenerate case in y-direction'
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          fac=.5*dz_1(n)
          df=fac*(f(l1:l2,m,n+1,k)-f(l1:l2,m,n-1,k))
          if (lspherical_coords) df=df*r1_mn*sin1th(m)
        else
          df=0.
          if (ip<=5) print*, 'der_main: Degenerate case in z-direction'
        endif
      endif
!
    endsubroutine der_main
!***********************************************************************
    subroutine der_other(f,df,j)
!
!  Along one pencil in NON f variable
!  calculate derivative of a scalar, get scalar
!  accurate to 6th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!
!  26-nov-02/tony: coded - duplicate der_main but without k subscript
!                          then overload the der interface.
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!  21-feb-07/axel: added 1/r and 1/pomega factors for non-coord basis
!  25-aug-09/axel: not yet adapted from deriv
!
      real, dimension (mx,my,mz) :: f
      real, dimension (nx) :: df,fac
      integer :: j
!
      intent(in)  :: f,j
      intent(out) :: df
!
!debug      if (loptimise_ders) der_call_count(1,icount_der_other,j,1) = &
!debug                          der_call_count(1,icount_der_other,j,1) + 1
!
      if (j==1) then
        if (nxgrid/=1) then
          fac=.5*dx_1(l1:l2)
          df=fac*(f(l1+1:l2+1,m,n)-f(l1-1:l2-1,m,n))
        else
          df=0.
          if (ip<=5) print*, 'der_other: Degenerate case in x-direction'
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          fac=.5*dy_1(m)
          df=fac*(f(l1:l2,m+1,n)-f(l1:l2,m-1,n))
          if (lspherical_coords)     df=df*r1_mn
          if (lcylindrical_coords)   df=df*rcyl_mn1
        else
          df=0.
          if (ip<=5) print*, 'der_other: Degenerate case in y-direction'
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          fac=.5*dz_1(n)
          df=fac*(f(l1:l2,m,n+1)-f(l1:l2,m,n-1))
          if (lspherical_coords) df=df*r1_mn*sin1th(m)
        else
          df=0.
          if (ip<=5) print*, 'der_other: Degenerate case in z-direction'
        endif
      endif
!
    endsubroutine der_other
!***********************************************************************
    subroutine der_pencil(j,pencil,df)
!
!  Calculate first derivative of any x, y or z pencil.
!
!  01-nov-07/anders: adapted from der
!  25-aug-09/axel: added fatal_error, because it is not adapted yet
!
      real, dimension (:) :: pencil, df
      integer :: j
!
      intent(in)  :: j, pencil
      intent(out) :: df
!
!  x-derivative
!
      if (j==1) then
        if (size(pencil)/=mx) then
          if (lroot) print*, 'der_pencil: pencil must be of size mx for x derivative'
          call fatal_error('der_pencil','')
        endif
        df(l1:l2)=0.5*dx_1(l1:l2)*(pencil(l1+1:l2+1)-pencil(l1-1:l2-1))
      else if (j==2) then
!
!  y-derivative
!
        if (size(pencil)/=my) then
          if (lroot) print*, 'der_pencil: pencil must be of size my for y derivative'
          call fatal_error('der_pencil','')
        endif
        df(m1:m2)=0.5*dy_1(m1:m2)*(pencil(m1+1:m2+1)-pencil(m1-1:m2-1))
      else if (j==3) then
!
!  z-derivative
!
        if (size(pencil)/=mz) then
          if (lroot) print*, 'der_pencil: pencil must be of size mz for z derivative'
          call fatal_error('der_pencil','')
        endif
        df(n1:n2)=0.5*dz_1(n1:n2)*(pencil(n1+1:n2+1)-pencil(n1-1:n2-1))
      else
        if (lroot) print*, 'der_pencil: no such direction j=', j
        call fatal_error('der_pencil','')
      endif
!
      if (lcylindrical_coords.or.lspherical_coords) &
           call fatal_error("der_pencil","Not implemented for non-cartesian")
!
    endsubroutine der_pencil
!***********************************************************************
    subroutine distr_der(arr,idir,der,order)
!
!  Dummy
!
    real, dimension(:,:), intent(in) :: arr
    integer, intent(in) :: idir
    real, dimension(:,:), intent(out) :: der
    integer, intent(in), optional :: order
!
    call not_implemented('distr_der','for 2nd order')
    call keep_compiler_quiet(arr)
    call keep_compiler_quiet(idir)
    call keep_compiler_quiet(der)
!
    endsubroutine distr_der
!***********************************************************************
    subroutine der2_main(f, k, df2, j, lwo_line_elem)
!
!  calculate 2nd derivative d^2f_k/dx_j^2
!  accurate to 6th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!
!   1-oct-97/axel: coded
!   1-apr-01/axel+wolf: pencil formulation
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      integer, intent(in) :: j, k
      real, dimension (nx), intent(out) :: df2
      logical, intent(in), optional :: lwo_line_elem
!
      real, parameter :: der2_coef0=-2., der2_coef1=1.
      real, dimension (nx) :: fac, df
!
!debug      if (loptimise_ders) der_call_count(k,icount_der2,j,1) = & !DERCOUNT
!debug                          der_call_count(k,icount_der2,j,1) + 1 !DERCOUNT
!
      if (present(lwo_line_elem)) then
        if (lwo_line_elem) call fatal_error("der2_main", "lwo_line_elem=T is not implemented")
      endif
!
      if (j==1) then
        if (nxgrid/=1) then
          fac=dx_1(l1:l2)**2
          df2=fac*(der2_coef0*f (l1  :l2  ,m,n,k) &
                  +der2_coef1*(f(l1+1:l2+1,m,n,k)+f(l1-1:l2-1,m,n,k)) )
          if (.not.lequidist(j)) then
            call der(f,k,df,j)
            df2=df2+dx_tilde(l1:l2)*df
          endif
        else
          df2=0.
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          fac=dy_1(m)**2
          df2=fac*(der2_coef0*f(l1:l2,m  ,n,k) &
                  +der2_coef1*(f(l1:l2,m+1,n,k)+f(l1:l2,m-1,n,k)) )
          if (lspherical_coords)     df2=df2*r2_mn
          if (lcylindrical_coords)   df2=df2*rcyl_mn2
          if (.not.lequidist(j)) then
            call der(f,k,df,j)
            df2=df2+dy_tilde(m)*df
          endif
        else
          df2=0.
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          fac=dz_1(n)**2
          df2=fac*(der2_coef0*f(l1:l2,m,n  ,k) &
                  +der2_coef1*(f(l1:l2,m,n+1,k)+f(l1:l2,m,n-1,k)) )
          if (lspherical_coords) df2=df2*r2_mn*sin2th(m)
          if (.not.lequidist(j)) then
            call der(f,k,df,j)
            df2=df2+dz_tilde(n)*df
          endif
        else
          df2=0.
        endif
      endif
!
!
    endsubroutine der2_main
!***********************************************************************
    subroutine der2_other(f,df2,j)
!
!  calculate 2nd derivative d^2f/dx_j^2 (of scalar f)
!  accurate to 6th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!
!   1-oct-97/axel: coded
!   1-apr-01/axel+wolf: pencil formulation
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!  25-aug-09/axel: added fatal_error, because it is not adapted yet
!
      real, dimension (mx,my,mz), intent(in) :: f
      real, dimension (nx), intent(out) :: df2
      integer, intent(in) :: j
!
      real, dimension (nx) :: fac, df
!
      if (j==1) then
        if (nxgrid/=1) then
          fac=dx_1(l1:l2)**2
          df2=fac*(-  2.0* f(l1:l2,m,n) &
                   +      (f(l1+1:l2+1,m,n)+f(l1-1:l2-1,m,n)) )
          if (.not.lequidist(j)) then
            call der(f,df,j)
            df2=df2+dx_tilde(l1:l2)*df
          endif
        else
          df2=0.
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          fac=dy_1(m)**2
          df2=fac*(-  2.0* f(l1:l2,m,n) &
                   +      (f(l1:l2,m+1,n)+f(l1:l2,m-1,n)) )
          if (lspherical_coords)     df2=df2*r2_mn
          if (lcylindrical_coords)   df2=df2*rcyl_mn2
          if (.not.lequidist(j)) then
            call der(f,df,j)
            df2=df2+dy_tilde(m)*df
          endif
        else
          df2=0.
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          fac=dz_1(n)**2
          df2=fac*(-  2.0* f(l1:l2,m,n) &
                   +      (f(l1:l2,m,n+1)+f(l1:l2,m,n-1)) )
          if (lspherical_coords) df2=df2*r2_mn*sin2th(m)
          if (.not.lequidist(j)) then
            call der(f,df,j)
            df2=df2+dz_tilde(n)*df
          endif
        else
          df2=0.
        endif
      endif
!
!
    endsubroutine der2_other
!***********************************************************************
    subroutine der2_pencil(j,pencil,df2)
!
!  Calculate 2nd derivative of any x, y or z pencil.
!
!  01-nov-07/anders: adapted from der2
!  25-aug-09/axel: added fatal_error, because it is not adapted yet
!
      real, dimension (:), intent(in) :: pencil
      real, dimension (:), intent(out) :: df2
      integer, intent(in) :: j
!
!  x-derivative
!
      if (j==1) then
        if (size(pencil)/=mx) then
          if (lroot) print*, 'der2_pencil: pencil must be of size mx for x derivative'
          call fatal_error('der2_pencil','')
        endif
        df2=dx_1(l1:l2)**2*(-  2.0*pencil(l1:l2) &
               +      (pencil(l1+1:l2+1)+pencil(l1-1:l2-1)) )
      else if (j==2) then
!
!  y-derivative
!
        if (size(pencil)/=my) then
          if (lroot) print*, 'der2_pencil: pencil must be of size my for y derivative'
          call fatal_error('der2_pencil','')
        endif
        df2=dy_1(m1:m2)**2*(-  2.0*pencil(m1:m2) &
               +      (pencil(m1+1:m2+1)+pencil(m1-1:m2-1)) )
      else if (j==3) then
!
!  z-derivative
!
        if (size(pencil)/=mz) then
          if (lroot) print*, 'der2_pencil: pencil must be of size mz for z derivative'
          call fatal_error('der2_pencil','')
        endif
        df2(n1:n2)=dz_1(n1:n2)**2*(-  2.0*pencil(n1:n2) &
               +      (pencil(n1+1:n2+1)+pencil(n1-1:n2-1)) )
      else
        if (lroot) print*, 'der2_pencil: no such direction j=', j
        call fatal_error('der2_pencil','')
      endif
!
    endsubroutine der2_pencil
!***********************************************************************
    subroutine der3(f,k,df,j,ignoredx)
!
!  Calculate 3rd derivative of a scalar, get scalar
!
!  10-feb-06/anders: adapted from der5
!  25-aug-09/axel: added fatal_error, because it is not adapted yet
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df
      integer :: j,k
      logical, optional :: ignoredx
!
      intent(in)  :: f,k,j,ignoredx
      intent(out) :: df
!
      call fatal_error('deriv_2nd','der3 not implemented yet')
      call keep_compiler_quiet(df)
!
    endsubroutine der3
!***********************************************************************
    subroutine der4(f,k,df,j,ignoredx,upwind)
!
!  Calculate 4th derivative of a scalar, get scalar
!    Used for hyperdiffusion that affects small wave numbers as little as
!  possible (useful for density).
!    The optional flag IGNOREDX is useful for numerical purposes, where
!  you want to affect the Nyquist scale in each direction, independent of
!  the ratios dx:dy:dz.
!
!   8-jul-02/wolf: coded
!   9-dec-03/nils: adapted from der6
!  10-feb-06/anders: corrected sign and factor
!  25-aug-09/axel: added fatal_error, because it is not adapted yet
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx), intent(out) :: df
      integer, intent(in) :: j, k
      logical, intent(in), optional :: ignoredx,upwind
!
      call fatal_error('deriv_2nd','der4 not implemented yet')
      call keep_compiler_quiet(df)
!
    endsubroutine der4
!***********************************************************************
    subroutine der5(f,k,df,j,ignoredx)
!
!  Calculate 5th derivative of a scalar, get scalar
!    Used for hyperdiffusion that affects small wave numbers as little as
!  possible (useful for density).
!    The optional flag IGNOREDX is useful for numerical purposes, where
!  you want to affect the Nyquist scale in each direction, independent of
!  the ratios dx:dy:dz.
!
!  29-oct-04/anders: adapted from der6
!  25-aug-09/axel: added fatal_error, because it is not adapted yet
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx), intent(out) :: df
      integer, intent(in) :: j, k
      logical, intent(in), optional :: ignoredx
!
      call fatal_error('deriv_2nd','der5 not implemented yet')
      call keep_compiler_quiet(df)
!
    endsubroutine der5
!***********************************************************************
    subroutine der6_main(f,k,df,j,ignoredx,upwind)
!
!  Calculate 6th derivative of a scalar, get scalar
!    Used for hyperdiffusion that affects small wave numbers as little as
!  possible (useful for density).
!    The optional flag IGNOREDX is useful for numerical purposes, where
!  you want to affect the Nyquist scale in each direction, independent of
!  the ratios dx:dy:dz.
!    The optional flag UPWIND is a variant thereof, which calculates
!  D^(6)*dx^5/60, which is the upwind correction of centered derivatives.
!
!   8-jul-02/wolf: coded
!  25-aug-09/axel: added fatal_error, because it is not adapted yet
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df,fac
      integer :: j,k
      logical, optional :: ignoredx,upwind
      logical :: igndx,upwnd
!
      intent(in)  :: f,k,j,ignoredx
      intent(out) :: df
!
      if (headtt) then
         call warning('deriv_2nd','der6 not implemented yet -- using 6th order')
      endif
!
!debug      if (loptimise_ders) der_call_count(k,icount_der6,j,1) = & !DERCOUNT
!debug                          der_call_count(k,icount_der6,j,1) + 1 !DERCOUNT
!
      if (present(ignoredx)) then
        igndx = ignoredx
      else
        igndx = .false.
      endif
      if (present(upwind)) then
        upwnd = upwind
      else
        upwnd = .false.
        if (.not. lequidist(j)) then
          call fatal_error('der6','NOT IMPLEMENTED for non-equidistant grid')
        endif
        if ((.not.lcartesian_coords).and.(.not.igndx)) then
          call fatal_error('der6','in non-cartesian coordinates '//&
               'only works if upwinding is used')
        endif
     endif
!
!
      if (j==1) then
        if (nxgrid/=1) then
          if (igndx) then
            fac=1.0
          else if (upwnd) then
            fac=(1.0/60)*dx_1(l1:l2)
          else
            fac=1/dx**6
          endif
          df=fac*(- 20.0* f(l1:l2,m,n,k) &
                  + 15.0*(f(l1+1:l2+1,m,n,k)+f(l1-1:l2-1,m,n,k)) &
                  -  6.0*(f(l1+2:l2+2,m,n,k)+f(l1-2:l2-2,m,n,k)) &
                  +      (f(l1+3:l2+3,m,n,k)+f(l1-3:l2-3,m,n,k)))
        else
          df=0.
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          if (igndx) then
            fac=1.0
          else if (upwnd) then
            fac=(1.0/60)*dy_1(m)
          else
            fac=1/dy**6
          endif
          df=fac*(- 20.0* f(l1:l2,m  ,n,k) &
                  + 15.0*(f(l1:l2,m+1,n,k)+f(l1:l2,m-1,n,k)) &
                  -  6.0*(f(l1:l2,m+2,n,k)+f(l1:l2,m-2,n,k)) &
                  +      (f(l1:l2,m+3,n,k)+f(l1:l2,m-3,n,k)))
         else
          df=0.
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          if (igndx) then
            fac=1.
          else if (upwnd) then
            fac=(1.0/60)*dz_1(n)
          else
            fac=1/dz**6
          endif
          df=fac*(- 20.0* f(l1:l2,m,n  ,k) &
                  + 15.0*(f(l1:l2,m,n+1,k)+f(l1:l2,m,n-1,k)) &
                  -  6.0*(f(l1:l2,m,n+2,k)+f(l1:l2,m,n-2,k)) &
                  +      (f(l1:l2,m,n+3,k)+f(l1:l2,m,n-3,k)))
         else
          df=0.
        endif
      endif
!
    endsubroutine der6_main
!***********************************************************************
    subroutine der6_other(f,df,j,ignoredx,upwind)
!
!  Calculate 6th derivative of a scalar, get scalar
!    Used for hyperdiffusion that affects small wave numbers as little as
!  possible (useful for density).
!    The optional flag IGNOREDX is useful for numerical purposes, where
!  you want to affect the Nyquist scale in each direction, independent of
!  the ratios dx:dy:dz.
!    The optional flag UPWIND is a variant thereof, which calculates
!  D^(6)*dx^5/60, which is the upwind correction of centered derivatives.
!
!   8-jul-02/wolf: coded
!  25-aug-09/axel: added fatal_error, because it is not adapted yet
!
      real, dimension (mx,my,mz), intent(in) :: f
      real, dimension (nx), intent(out) :: df
      integer, intent(in) :: j
      logical, intent(in), optional :: ignoredx, upwind
!
      call fatal_error('deriv_2nd','der6_other not implemented yet')
      call keep_compiler_quiet(df)
!
    endsubroutine der6_other
!***********************************************************************
    subroutine der6_pencil(j,pencil,df6,ignoredx,upwind)
!
!  Calculate 6th derivative of any x, y or z pencil.
!
!  20-jul-20/wlyra: adapted from der2_pencil
!
      integer :: j
      real, dimension (:) :: pencil
      real, dimension (:) :: df6
      logical, optional :: ignoredx, upwind
!
      real, dimension (nx) :: facx
      real, dimension (ny) :: facy
      real, dimension (nz) :: facz
      logical :: igndx,upwnd
!
      intent(in)  :: j, pencil,ignoredx,upwind
      intent(out) :: df6
!
      if (headtt) then
         call warning('deriv_2nd','der6_pencil not implemented yet -- using 6th order')
      endif
!
      if (present(ignoredx)) then
        igndx = ignoredx
      else
        igndx = .false.
      endif
!
      if (present(upwind)) then
        if (.not. lequidist(j).and..not.lignore_nonequi) then
          call fatal_error('der6','upwind cannot be used with '//&
              'non-equidistant grid.')
        endif
        upwnd = upwind
      else
        upwnd = .false.
      endif
!
!  x-derivative
!
      if (j==1) then
        if (size(pencil)/=mx) then
          if (lroot) print*, 'der6_pencil: pencil must be of size mx for x derivative'
          call fatal_error('der6_pencil','')
        endif
        if (igndx) then
          facx=1.
        elseif (upwnd) then
          facx=(1.0/60)*dx_1(l1:l2)
        else
          facx=dx_1(l1:l2)**6
        endif
        df6=facx*(- 20.0* pencil(l1:l2) &
                  + 15.0*(pencil(l1+1:l2+1)+pencil(l1-1:l2-1)) &
                  -  6.0*(pencil(l1+2:l2+2)+pencil(l1-2:l2-2)) &
                  +      (pencil(l1+3:l2+3)+pencil(l1-3:l2-3)))
      else if (j==2) then
!
!  y-derivative
!
        if (size(pencil)/=my) then
          if (lroot) &
              print*, 'der6_pencil: pencil must be of size my for y derivative'
          call fatal_error('der6_pencil','')
        endif
        if (igndx) then
          facy=1.
        else if (upwnd) then
          facy=(1.0/60)*dy_1(m1:m2)
        else
          facy=dy_1(m1:m2)**6
        endif
        df6=facy*(- 20.0* pencil(m1:m2) &
                  + 15.0*(pencil(m1+1:m2+1)+pencil(m1-1:m2-1)) &
                  -  6.0*(pencil(m1+2:m2+2)+pencil(m1-2:m2-2)) &
                  +      (pencil(m1+3:m2+3)+pencil(m1-3:m2-3)))
      else if (j==3) then
!
!  z-derivative
!
        if (size(pencil)/=mz) then
          if (lroot) &
              print*, 'der6_pencil: pencil must be of size mz for z derivative'
          call fatal_error('der6_pencil','')
        endif
        if (igndx) then
          facz=1.
        else if (upwnd) then
          facz=(1.0/60)*dz_1(n1:n2)
        else
          facz=dz_1(n1:n2)**6
        endif
        df6=facz*(- 20.0* pencil(n1:n2) &
                  + 15.0*(pencil(n1+1:n2+1)+pencil(n1-1:n2-1)) &
                  -  6.0*(pencil(n1+2:n2+2)+pencil(n1-2:n2-2)) &
                  +      (pencil(n1+3:n2+3)+pencil(n1-3:n2-3)))
      else
        if (lroot) print*, 'der6_pencil: no such direction j=', j
        call fatal_error('der6_pencil','')
      endif
!
    endsubroutine der6_pencil
!***********************************************************************
    real function der5_single(f,j,dc1)
!
!  computes 5th order derivative of function given by f at position j
!
!   3-oct-12/MR: coded
!
      real, dimension(:), intent(in) :: f, dc1
      integer, intent(in) :: j
!
      call fatal_error('deriv_2nd','der6_pencil not implemented yet')
      call keep_compiler_quiet(der5_single)
!
    endfunction der5_single
!***********************************************************************
    subroutine derij_main(f, k, df, i, j, lwo_line_elem)
!
!  calculate 2nd derivative with respect to two different directions
!  input: scalar, output: scalar
!  accurate to 6th order, explicit, periodic
!
!   8-sep-01/axel: coded
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!  14-nov-06/wolf: implemented bidiagonal scheme
!  25-aug-09/axel: added fatal_error, because it is not adapted yet
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx), intent(out) :: df
      integer, intent(in) :: i, j, k
      logical, intent(in), optional :: lwo_line_elem
!
      real, dimension (nx) :: fac
!
!debug      if (loptimise_ders) der_call_count(k,icount_derij,i,j) = & !DERCOUNT
!debug                          der_call_count(k,icount_derij,i,j) + 1 !DERCOUNT
!
      if (present(lwo_line_elem)) then
        if (lwo_line_elem) call fatal_error("derij_main", "lwo_line_elem=T is not implemented")
      endif
!
      if (lbidiagonal_derij) then
        !
        ! Use bidiagonal mixed-derivative operator, i.e.
        ! employ only the three neighbouring points on each of the four
        ! half-diagonals. This gives 6th-order mixed derivatives as the
        ! version below, but involves just 12 points instead of 36.
        !
        if ((i==1.and.j==2).or.(i==2.and.j==1)) then
          if (nxgrid/=1.and.nygrid/=1) then
            fac=.25*dx_1(l1:l2)*dy_1(m)
            df=fac*( &
                       ( f(l1+1:l2+1,m+1,n,k)-f(l1-1:l2-1,m+1,n,k)  &
                        +f(l1-1:l2-1,m-1,n,k)-f(l1+1:l2+1,m-1,n,k)) &
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij_main: Degenerate case in x- or y-direction'
          endif
        elseif ((i==2.and.j==3).or.(i==3.and.j==2)) then
          if (nygrid/=1.and.nzgrid/=1) then
            fac=.25*dy_1(m)*dz_1(n)
            df=fac*( &
                       ( f(l1:l2,m+1,n+1,k)-f(l1:l2,m+1,n-1,k)  &
                        +f(l1:l2,m-1,n-1,k)-f(l1:l2,m-1,n+1,k)) &
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij_main: Degenerate case in y- or z-direction'
          endif
        elseif ((i==3.and.j==1).or.(i==1.and.j==3)) then
          if (nzgrid/=1.and.nxgrid/=1) then
            fac=.25*dz_1(n)*dx_1(l1:l2)
            df=fac*( &
                       ( f(l1+1:l2+1,m,n+1,k)-f(l1-1:l2-1,m,n+1,k)  &
                        +f(l1-1:l2-1,m,n-1,k)-f(l1+1:l2+1,m,n-1,k)) &
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij_main: Degenerate case in x- or z-direction'
          endif
        endif
!
      else                      ! not using bidiagonal mixed derivatives
        !
        call fatal_error('derij_main','only implemented for lbidiagonal_derij=T')
        !
        ! This is the old, straight-forward scheme
        !
        if ((i==1.and.j==2).or.(i==2.and.j==1)) then
          if (nxgrid/=1.and.nygrid/=1) then
            fac=0.0625*dx_1(l1:l2)*dy_1(m)
            df=fac*( (f(l1+1:l2+1,m+1,n,k)-f(l1-1:l2-1,m+1,n,k)) &
                    -(f(l1+1:l2+1,m-1,n,k)-f(l1-1:l2-1,m-1,n,k)) )
          else
            df=0.
            if (ip<=5) print*, 'derij_main: Degenerate case in x- or y-direction'
          endif
        elseif ((i==2.and.j==3).or.(i==3.and.j==2)) then
          if (nygrid/=1.and.nzgrid/=1) then
            fac=0.0625*dy_1(m)*dz_1(n)
            df=fac*( (f(l1:l2,m+1,n+1,k)-f(l1:l2,m-1,n+1,k)) &
                    -(f(l1:l2,m+1,n-1,k)-f(l1:l2,m-1,n-1,k)) )
          else
            df=0.
            if (ip<=5) print*, 'derij_main: Degenerate case in y- or z-direction'
          endif
        elseif ((i==3.and.j==1).or.(i==1.and.j==3)) then
          if (nzgrid/=1.and.nxgrid/=1) then
            fac=0.0625*dz_1(n)*dx_1(l1:l2)
            df=fac*( (f(l1+1:l2+1,m,n+1,k)-f(l1-1:l2-1,m,n+1,k)) &
                    -(f(l1+1:l2+1,m,n-1,k)-f(l1-1:l2-1,m,n-1,k)) )
          else
            df=0.
            if (ip<=5) print*, 'derij_main: Degenerate case in x- or z-direction'
          endif
        endif
!
      endif
!
!  Spherical polars. The comments about "minus extra terms" refer to the
!  presence of extra terms that are being evaluated later in gij_etc.
!
      if (lspherical_coords) then
        if ((i==1.and.j==2)) df=df*r1_mn
        if ((i==2.and.j==1)) df=df*r1_mn !(minus extra terms)
        if ((i==1.and.j==3)) df=df*r1_mn*sin1th(m)
        if ((i==3.and.j==1)) df=df*r1_mn*sin1th(m) !(minus extra terms)
        if ((i==2.and.j==3)) df=df*r2_mn*sin1th(m)
        if ((i==3.and.j==2)) df=df*r2_mn*sin1th(m) !(minus extra terms)
      endif
!
      if (lcylindrical_coords) then
        if ((i==1.and.j==2)) df=df*rcyl_mn1
        if ((i==2.and.j==1)) df=df*rcyl_mn1
        if ((i==1.and.j==3)) df=df
        if ((i==3.and.j==1)) df=df
        if ((i==2.and.j==3)) df=df*rcyl_mn1
        if ((i==3.and.j==2)) df=df*rcyl_mn1
      endif
!
    endsubroutine derij_main
!***********************************************************************
    subroutine derij_other(f,df,i,j)
!
!  calculate 2nd derivative with respect to two different directions
!  input: scalar, output: scalar
!  accurate to 6th order, explicit, periodic
!
!   8-sep-01/axel: coded
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!  14-nov-06/wolf: implemented bidiagonal scheme
!  25-aug-09/axel: added fatal_error, because it is not adapted yet
!
      real, dimension (mx,my,mz), intent(in) :: f
      real, dimension (nx), intent(out) :: df
      integer, intent(in) :: i, j
!
      real, dimension (nx) :: fac
!
      call fatal_error('deriv_2nd','derij_other not implemented yet')
      call keep_compiler_quiet(df)
!
    endsubroutine derij_other
!***********************************************************************
    subroutine der5i1j(f,k,df,i,j)
!
!  Calculate 6th derivative with respect to two different directions.
!
!  05-dec-06/anders: adapted from derij
!  25-aug-09/axel: added fatal_error, because it is not adapted yet
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx), intent(out) :: df
      integer, intent(in) :: i, j, k
!
      call fatal_error('der5i1j','not implemented in deriv_2nd')
      call keep_compiler_quiet(df)
!
    endsubroutine der5i1j
!***********************************************************************
    subroutine der4i2j(f,k,df,i,j)
!
!  Calculate 6th derivative with respect to two different directions.
!
!  02-apr-17/wlyra: adapted from der5i1j
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx), intent(in) :: df
      integer, intent(in) :: i, j, k
!
      call fatal_error("der4i2j","not implemented in deriv_2nd")
      call keep_compiler_quiet(df)
!
    endsubroutine der4i2j
!***********************************************************************
    subroutine der2i2j2k(f,k,df)
!
!  Mixed 6th derivative of der2x(der2y(der2z(f))). Worked out symbolically
!  in python. Result as spit from the python routine.
!
!  02-apr-17/wlyra: coded
!
      real, dimension (mx,my,mz,mfarray),intent(in) :: f
      integer,intent(in) :: k
      real, dimension(nx), intent(out) :: df
!
      call fatal_error("der2i2j2k","not implemented in deriv_2nd")
      call keep_compiler_quiet(df)
!
    endsubroutine der2i2j2k
!***********************************************************************
    subroutine der3i3j(f,k,df,i,j)
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx), intent(out) :: df
      integer, intent(in) :: k,i,j
!
      call fatal_error("der3i3j","not implemented in deriv_2nd")
      call keep_compiler_quiet(df)
!
    endsubroutine der3i3j
!***********************************************************************
    subroutine der3i2j1k(f,ik,df,i,j,k)
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx), intent(out) :: df
      integer, intent(in) :: ik,i,j,k
!
      call fatal_error("der3i2j1k","not implemented in deriv_2nd")
      call keep_compiler_quiet(df)
!
    endsubroutine der3i2j1k
!***********************************************************************
    subroutine der4i1j1k(f,ik,df,i,j,k)
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx), intent(out) :: df
      integer, intent(in) :: ik,i,j,k
!
      call fatal_error("der4i1j1k","not implemented in deriv_10th")
      call keep_compiler_quiet(df)
!
    endsubroutine der4i1j1k
!***********************************************************************
    subroutine der_upwind1st(f,uu,k,df,j)
!
!  First order upwind derivative of variable
!  Useful for advecting non-logarithmic variables
!
!  25-aug-09/axel: added fatal_error, because it is not adapted yet
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx,3), intent(in) :: uu
      real, dimension (nx), intent(out) :: df
      integer, intent(in) :: j, k
!
      call fatal_error('deriv_2nd','der_upwind1st not implemented yet')
      call keep_compiler_quiet(df)
!
    endsubroutine der_upwind1st
!***********************************************************************
    subroutine der_onesided_4_slice_main(f,sgn,k,df,pos,j)
!
!   Calculate x/y/z-derivative on a yz/xz/xy-slice at gridpoint pos.
!   Uses a one-sided 4th order stencil.
!   sgn = +1 for forward difference, sgn = -1 for backwards difference.
!
!   Because of its original intended use in relation to solving
!   characteristic equations on boundaries (NSCBC), this sub should
!   return only PARTIAL derivatives, NOT COVARIANT. Applying the right
!   scaling factors and connection terms should instead be done when
!   solving the characteristic equations.
!
!   7-jul-08/arne: coded.
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (:,:), intent(out) :: df
      integer, intent(in) :: pos, k, sgn, j
!
      call fatal_error('deriv_2nd','der_onesided_4_slice_main not implemented yet')
      call keep_compiler_quiet(df)
!
    endsubroutine der_onesided_4_slice_main
!***********************************************************************
   subroutine der_onesided_4_slice_other(f,sgn,df,pos,j)
!
!   Calculate x/y/z-derivative on a yz/xz/xy-slice at gridpoint pos.
!   Uses a one-sided 4th order stencil.
!   sgn = +1 for forward difference, sgn = -1 for backwards difference.
!
!   Because of its original intended use in relation to solving
!   characteristic equations on boundaries (NSCBC), this sub should
!   return only PARTIAL derivatives, NOT COVARIANT. Applying the right
!   scaling factors and connection terms should instead be done when
!   solving the characteristic equations.
!
!   7-jul-08/arne: coded.
!
      real, dimension (mx,my,mz), intent(in) :: f
      real, dimension (:,:), intent(out) :: df
      integer, intent(in) :: pos, sgn, j
!
      call fatal_error('deriv_2nd','der_onesided_4_slice_other not implemented yet')
      call keep_compiler_quiet(df)
!
    endsubroutine der_onesided_4_slice_other
!***********************************************************************
    subroutine der_onesided_4_slice_main_pt(f,sgn,k,df,lll,mmm,nnn,j)
!
!  made using der_onesided_4_slice_main. One sided derivative is calculated
!  at one point (lll,mmm,nnn)
!
!  15-oct-09/Natalia: coded.
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, intent(out) :: df
      integer , intent(in):: lll, mmm, nnn, k, sgn, j
!
      call not_implemented('der_onesided_4_slice_main_pt','')
      call keep_compiler_quiet(df)
!
   endsubroutine der_onesided_4_slice_main_pt
!***********************************************************************
   subroutine der_onesided_4_slice_other_pt(f,sgn,df,lll,mmm,nnn,j)
!
!  One sided derivative is calculated
!  at one point (lll,mmm,nnn).
!
!  15-oct-09/Natalia: coded.
!  15-oct-09/axel: changed file name to shorter version
!
      real, dimension (mx,my,mz), intent(in) :: f
      real, intent(out) :: df
      integer, intent(in) :: lll, mmm, nnn, sgn, j
!
      call not_implemented('der_onesided_4_slice_other_pt','')
      call keep_compiler_quiet(df)
!
   endsubroutine der_onesided_4_slice_other_pt
!***********************************************************************
    subroutine der_z(f,df)
!
!  z derivative operating on a z-dependent 1-D array
!
!  19-may-11/bing: adapted from der_main; note that f is not the f array!
!
      real, dimension (mz), intent(in)  :: f
      real, dimension (nz), intent(out) :: df
!
      real, dimension (nz) :: fac
!
      if (nzgrid/=1) then
        fac=.5*dz_1(n1:n2)
        df=fac*(f(n1+1:n2+1)-f(n1-1:n2-1))
      else
        df=0.
        if (ip<=5) print*, 'der_z: Degenerate case in z-direction'
      endif
!
      if (lspherical_coords) &
          call fatal_error('der_z','not implemented for spherical coords')
!
    endsubroutine der_z
!***********************************************************************
    subroutine der2_z(f,df2)
!
!  z derivative operating on a z-dependent 1-D array
!
!  19-may-11/bing: adapted from der_z and der_main
!
      real, dimension (mz), intent(in)  :: f
      real, dimension (nz), intent(out) :: df2
!
      real, dimension (nz) :: fac,df
      real, parameter :: der2_coef0=-2., der2_coef1=1.
!
      if (nzgrid/=1) then
        fac=dz_1(n1:n2)**2
        df2=fac*(der2_coef0*f(n1:n2) &
            +der2_coef1*(f(n1+1:n2+1)+f(n1-1:n2-1)))
!
        if (.not.lequidist(3)) then
          call der_z(f,df)
          df2=df2+dz_tilde(n1:n2)*df
        endif
      else
        df2=0.
        if (ip<=5) print*, 'der2_z: Degenerate case in z-direction'
      endif
!
      if (lspherical_coords) &
          call fatal_error('der2_z','not implemented for spherical coords')
!
    endsubroutine der2_z
!***********************************************************************
    subroutine der_x(f,df)
!
! dummy routine
!
      use Cparam, only: mz, nz
      use Mpicomm, only: stop_it
!
      real, dimension (mx), intent(in)  :: f
      real, dimension (nx), intent(out) :: df
!
      call not_implemented("der_x","")
      call keep_compiler_quiet(df)
!
    endsubroutine der_x
!***********************************************************************
    subroutine der2_x(f,df2)
!
! dummy routine
!
      use Cparam, only: mz, nz
      use Mpicomm, only: stop_it
!
      real, dimension (mx), intent(in)  :: f
      real, dimension (nx), intent(out) :: df2
!
      call not_implemented("der2_x","")
      call keep_compiler_quiet(df2)
!
    endsubroutine der2_x
!***********************************************************************
    subroutine der2_minmod(f,j,delfk,delfkp1,delfkm1,k)
!
!  Dummy routine
!
!  09-Sep-2024/PABourdin: not yet implemented
!
      intent(in) :: f,k,j
      intent(out) :: delfk,delfkp1,delfkm1
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: delfk,delfkp1,delfkm1
      integer :: j,k
!
      call fatal_error('der2_minmod','Not implemented for deriv_2nd')
      call keep_compiler_quiet(delfk)
      call keep_compiler_quiet(delfkp1)
      call keep_compiler_quiet(delfkm1)
!
    endsubroutine der2_minmod
!***********************************************************************
    subroutine finalize_deriv()
!
!  Dummy
!
    endsubroutine finalize_deriv
!***********************************************************************
    subroutine deri_3d_inds(f,df,inds,j,lignored,lnometric)
!
!  dummy routine for compatibility
!
!  26-mar-12/MR: coded
!
      real, dimension (mx,my,mz)          :: f
      real, dimension (nx)                :: df
      integer                             :: j
      logical,                   optional :: lignored, lnometric
      integer, dimension(nx)              :: inds
!
      intent(in)  :: f,j,inds,lignored,lnometric
      intent(out) :: df
!
      call fatal_error('deri_3d_inds','Upwinding not implemented for nonuniform grids')
      call keep_compiler_quiet(df)
!
    endsubroutine deri_3d_inds
!************************************************************************
    logical function heatflux_deriv_x(f, inh, fac, topbot)
!
!   dummy routine
!
!  17-apr-12/MR: coded
!
      real, dimension(mx,my,mz,mfarray), intent(IN):: f
      real, dimension(my,mz)           , intent(IN):: inh
      real                             , intent(IN):: fac
      integer                          , intent(IN):: topbot
!
      heatflux_deriv_x = .false.
!
    endfunction heatflux_deriv_x
!***********************************************************************
    subroutine set_ghosts_for_onesided_ders(f,topbot,j,idir,l2nd_)
!
!  Dummy.
!
      real, dimension(mx,my,mz,*) :: f
      integer, intent(IN) :: topbot
      integer :: j,idir
      logical, optional :: l2nd_
!
      call fatal_error('set_ghosts_for_onesided_ders','Not implemented for 2nd order.')
!
    endsubroutine set_ghosts_for_onesided_ders
!***********************************************************************
    subroutine bval_from_neumann_scl(f,topbot,j,idir,val)
!
!  Dummy.
!
      real, dimension(mx,my,mz,*) :: f
      integer, intent(IN) :: topbot
      integer :: j,idir
      real :: val
!
      call fatal_error('bval_from_neumann_scl','Not implemented for 2nd order.')
!
    endsubroutine bval_from_neumann_scl
!***********************************************************************
    subroutine bval_from_neumann_arr(f,topbot,j,idir,val)
!
!  Dummy.
!
      real, dimension(mx,my,mz,*) :: f
      integer, intent(IN) :: topbot
      integer :: j,idir
      real, dimension(:,:) :: val
!
      call fatal_error('bval_from_neumann_arr','Not implemented for 2nd order.')
!
    endsubroutine bval_from_neumann_arr
!***********************************************************************
    subroutine bval_from_3rd_scl(f,topbot,j,idir,val)
!
!  Dummy.
!
      real, dimension(mx,my,mz,*) :: f
      integer, intent(IN) :: topbot
      integer :: j,idir
      real :: val
!
      call fatal_error('bval_from_3rd_scl','Not implemented for 2nd order.')
!
    endsubroutine bval_from_3rd_scl
!***********************************************************************
    subroutine bval_from_3rd_arr(f,topbot,j,idir,val,func)
!
!  Calculates the boundary value from the Neumann BC d f/d x_i = val employing
!  one-sided difference formulae. val depends on x,y.
!
!  30-sep-16/MR: coded
!
      real, dimension(mx,my,mz,*) :: f
      integer, intent(IN) :: topbot
      integer :: j,idir
      real, dimension(:,:) :: val
      external :: func
!
      call fatal_error('bval_from_3rd_arr','not implemented for 2nd order')
!
    endsubroutine bval_from_3rd_arr
!***********************************************************************
    subroutine bval_from_4th_scl(f,topbot,j,idir,val)
!
!  Dummy.
!
      real, dimension(mx,my,mz,*) :: f
      integer, intent(IN) :: topbot
      integer :: j,idir
      real :: val
!
      call fatal_error('bval_from_4th_scl','Not implemented for 2nd order.')
!
    endsubroutine bval_from_4th_scl
!***********************************************************************
    subroutine bval_from_4th_arr(f,topbot,j,idir,val)
!
!  Calculates the boundary value from the 4th kind BC d^2 f/d x_i^2 = val*f
!  employing one-sided difference formulae. val depends on x,y.
!
!  09-feb-17/Ivan: coded
!
      real, dimension(mx,my,mz,*) :: f
      integer, intent(IN) :: topbot
      integer :: j,idir
      real, dimension(:,:) :: val
!
      call not_implemented('bval_from_4th_arr','')
!
    endsubroutine bval_from_4th_arr
!***********************************************************************
endmodule Deriv
