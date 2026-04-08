!************************************************************************
!* This file is part of OpenMolcas.                                     *
!*                                                                      *
!* OpenMolcas is free software; you can redistribute it and/or modify   *
!* it under the terms of the GNU Lesser General Public License, v. 2.1. *
!* OpenMolcas is distributed in the hope that it will be useful, but it *
!* is provided "as is" and without any express or implied warranties.   *
!* For more details see the full text of the license in the file        *
!* LICENSE or in <http://www.gnu.org/licenses/>.                        *
!************************************************************************

!MODULE io_data
!
!PUBLIC :: read_magnetic_moment,       write_magnetic_moment,  &
!          read_electric_moment,       write_electric_moment,  &
!          read_spin_moment,           write_spin_moment,      &
!          read_angmom,                write_angmom,           &
!          read_edipmom,               write_edipmom,          &
!          read_amfi,                  write_amfi,             &
!          read_nss,                   write_nss,              &
!          read_nmult,                 write_nmult,            &
!          read_imult,                 write_imult,            &
!          read_format,                write_format,           &
!          read_nroot,                 write_nroot,            &
!          read_nstate,                write_nstate,           &
!          read_multiplicity,          write_multiplicity,     &
!          read_szproj,                write_szproj,           &
!          read_eso,                   write_eso,              &
!          read_esfs,                  write_esfs,             &
!          read_hso,                   write_hso,              &
!          read_eigen,                 write_eigen,            &
!          read_gtens,                 write_gtens,            &
!          read_stev_cfp,              write_stev_cfp,         &
!          read_susc,                  write_susc,             &
!          read_magn,                  write_magn,             &
!          read_complex_matrix,        write_complex_matrix,   &
!          open_datafile_read,         open_datafile_write,    &
!          close_datafile,                                     &
!          check_hermiticity_matrix,   check_commutation,      &
!          check_S_square
!
!PRIVATE :: file_advance_to_string,     inquire_key_presence,   &
!           read_INTEGER_scalar,        write_INTEGER_scalar,   &
!           read_1d_INTEGER_array,      write_1d_INTEGER_array, &
!           read_1d_real_array,         write_1d_real_array,    &
!           read_2d_real_array,         write_2d_real_array,    &
!           read_3d_real_array,         write_3d_real_array
!
!           write_string,                                       &
!           read_real_scalar,           write_real_scalar,      &
!           read_complex_scalar,        write_complex_scalar,   &
!           read_1d_size,                                       &
!           read_2d_size,                                       &
!           read_3d_size,                                       &
!           read_4d_size,                                       &
!           read_1d_complex_array,      write_1d_complex_array, &
!           read_2d_complex_array,      write_2d_complex_array, &
!           read_2d_INTEGER_array,      write_2d_INTEGER_array, &
!           read_3d_INTEGER_array,      write_3d_INTEGER_array, &
!           read_4d_INTEGER_array,      write_4d_INTEGER_array, &
!           read_4d_real_array,         write_4d_real_array,    &
!           read_3d_complex_array,      write_3d_complex_array, &
!           read_4d_complex_array,      write_4d_complex_array
!
! Keywords:
!    $dipm   --  electric dipole moment, SO basis  (6 components: _xr,_xi,  _yr,_yi,  _zr,_zi)
!    $magn   --  magnetic dipole moment, SO basis  (6 components: _xr,_xi,  _yr,_yi,  _zr,_zi)
!    $angmom --  orbital moment (L)    , SF basis  (3 compinents: _x, _y, _z)
!    $edip   --  electric dipole moment, SF basis  (3 components: _x, _y, _z)
!    $hso    --  spin-orbit hamiltonian, SO basis  (2 components: _r, _i)
!    $eigen  --  spin-orbit eigenstates, SO basis  (2 components: _r, _i)
!    $eso
!    $esfs
!    $multiplicity
!    $szvalue

!---- high level functions -------
!     >>> available <<<
! read_magnetic_moment
! read_electric_moment
! read_angmom
! read_edipmom
! read_eigenv
! read_hso
! read_nss
! read_nstate
! read_eso
! read_esfs
! read_szproj
! read_gtens
! xt,
! xt_field,
!    >>> verification functions <<<
! verify_commutation(L, S)
! verify_L2
! verify_S2
! verify_hermiticity
!---------------------------------
!   >>> yet_to_add <<<
! write_exch_ham
! write_exch_functions
! write_CFP
! write_magnetic_axes
! write_xyz
!
! x,
! x_field,
! x_minus_one,
! x_field,
! x_field_minus_one,
! magnetization,
!
! read_MOs
! read_CI_in_SD
! read_CI_in_CSF
! read_coeff_of_CSF
!---------------------------------

!------ low-level functions ----
! read_INTEGER_scalar
! read_real_scalar
! read_complex_scalar
! read_string
!
! read_string_length
!
! read_1d_size
! read_2d_size
! read_3d_size
! read_4d_size
!
! read_1d_INTEGER_array
! read_2d_INTEGER_array
! read_3d_INTEGER_array
! read_4d_INTEGER_array
!
! read_1d_real_array
! read_2d_real_array
! read_3d_real_array
! read_4d_real_array
!
! read_1d_complex_array
! read_2d_complex_array
! read_3d_complex_array
! read_4d_complex_array
!
!
! write_INTEGER_scalar
! write_real_scalar
! write_complex_scalar
! write_string
!
! write_1d_INTEGER_array
! write_2d_INTEGER_array
! write_3d_INTEGER_array
! write_4d_INTEGER_array
!
! write_1d_real_array
! write_2d_real_array
! write_3d_real_array
! write_4d_real_array
!
! write_1d_complex_array
! write_2d_complex_array
! write_3d_complex_array
! write_4d_complex_array
!
! open_datafile_write
! open_datafile_read
! close_datafile
! key_found
! file_advance_to_string
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
!           HIGH LEVEL VERIFICATION SUBROUTINES
!----------------------------------------------------------------------!

subroutine check_commutation(n,moment,dbg)
! valid for S, L, J in spin-orbit basis

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: cZero, cOne, Onei
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: N
complex(kind=wp), intent(in) :: moment(3,N,N)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, j, l
complex(kind=wp) :: tr
complex(kind=wp), allocatable :: X(:,:), Y(:,:), Z(:,:), XY(:,:), XZ(:,:), YX(:,:), YZ(:,:), ZX(:,:), ZY(:,:)

! verify the commutation relations for S
call mma_allocate(X,n,n,label='X')
call mma_allocate(Y,n,n,label='Y')
call mma_allocate(Z,n,n,label='Z')
call mma_allocate(XY,n,n,label='XY')
call mma_allocate(YX,n,n,label='YX')
call mma_allocate(YZ,n,n,label='YZ')
call mma_allocate(ZY,n,n,label='ZY')
call mma_allocate(ZX,n,n,label='ZX')
call mma_allocate(XZ,n,n,label='XZ')

X(:,:) = moment(1,:,:)
Y(:,:) = moment(2,:,:)
Z(:,:) = moment(3,:,:)

call zgemm_('n','n',n,n,n,cOne,X,n,Y,n,cZero,XY,n)
call zgemm_('n','n',n,n,n,cOne,Y,n,X,n,cZero,YX,n)

call zgemm_('n','n',n,n,n,cOne,Y,n,Z,n,cZero,YZ,n)
call zgemm_('n','n',n,n,n,cOne,Z,n,Y,n,cZero,ZY,n)

call zgemm_('n','n',n,n,n,cOne,Z,n,X,n,cZero,ZX,n)
call zgemm_('n','n',n,n,n,cOne,X,n,Z,n,cZero,XZ,n)

if (dbg) then
  do l=1,3
    write(u6,'(A,I2)') 'check_commutation:: moment, projection, L=',l
    do i=1,n
      write(u6,'(10(2F8.4,2x))') (moment(l,i,j),j=1,n)
    end do
  end do

  write(u6,'(A,I2)') 'check_commutation:: XY-YX'
  do i=1,n
    write(u6,'(10(2F8.4,2x))') (XY(i,j)-YX(i,j),j=1,n)
  end do

  write(u6,'(A,I2)') 'check_commutation:: i*Z'
  do i=1,n
    write(u6,'(10(2F8.4,2x))') (Onei*moment(3,i,j),j=1,n)
  end do

  write(u6,'(A,I2)') 'check_commutation:: (XY-YX) -i*Z'
  do i=1,n
    write(u6,'(10(2F8.4,2x))') (XY(i,j)-YX(i,j)-Onei*moment(3,i,j),j=1,n)
  end do
end if

tr = cZero
do j=1,n
  do i=1,n
    tr = tr+XY(i,j)-YX(i,j)-Onei*moment(3,i,j)+YZ(i,j)-ZY(i,j)-Onei*moment(1,i,j)+ZX(i,j)-XZ(i,j)-Onei*moment(2,i,j)
  end do
end do

if (dbg) write(u6,'(A,ES22.14)') 'check_commutation::  trace of [Sx,Sy]-iSz = ',abs(tr)
if (abs(tr) > 1.0e-6_wp) then
  call WarningMessage(1,'check_commutation:: trace of [Sx,Sy]-iSz  is larger than 1.0e-6. The input moment looks inacurate')
else
  write(u6,'(A,ES22.14)') 'check_commutation:  The input moment passes all three commutation tests.'
end if

call mma_deallocate(X)
call mma_deallocate(Y)
call mma_deallocate(Z)
call mma_deallocate(XY)
call mma_deallocate(YX)
call mma_deallocate(YZ)
call mma_deallocate(ZY)
call mma_deallocate(ZX)
call mma_deallocate(XZ)

return

end subroutine check_commutation

subroutine check_S_square(n,moment,dbg)
! valid also for J, L, S

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Quart, cZero, cOne
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n
complex(kind=wp), intent(in) :: moment(3,n,n)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i
real(kind=wp) :: S2_theoretic
complex(kind=wp) :: tr
complex(kind=wp), allocatable :: M(:,:), S2(:,:), X2(:,:), Y2(:,:), Z2(:,:)

call mma_allocate(X2,n,n,label='X2')
call mma_allocate(Y2,n,n,label='Y2')
call mma_allocate(Z2,n,n,label='Z2')
call mma_allocate(M,n,n,label='M')

M(:,:) = moment(1,:,:)
call zgemm_('c','n',n,n,n,cOne,M,n,M,n,cZero,X2,n)

M(:,:) = moment(2,:,:)
call zgemm_('c','n',n,n,n,cOne,M,n,M,n,cZero,Y2,n)

M(:,:) = moment(3,:,:)
call zgemm_('c','n',n,n,n,cOne,M,n,M,n,cZero,Z2,n)

call mma_deallocate(M)
call mma_allocate(S2,n,n,label='S2')

! matrix add:
! S2 = X2 + Y2 + Z2
S2(:,:) = X2(:,:)+Y2(:,:)+Z2(:,:)

tr = cZero
do i=1,n
  tr = tr+S2(i,i)
end do

S2_theoretic = real(n**3-n,kind=wp)*Quart

if (dbg) write(u6,'(A,ES22.14)') 'check_S_square::  trace of S2=(Sx).Sx+(Sy).Sy+(Sz).Sz = ',abs(tr)
if ((abs(tr)-S2_theoretic) > 1.0e-6_wp) then
  call WarningMessage(1,'check_S_square:: tr(S^2) - S2_theoretic is larger than 1.0e-6. The input moment looks inacurate')
else
  write(u6,'(A,ES22.14)') 'check_S_square:  The input moment passes the S^2 test.'
end if
call mma_deallocate(X2)
call mma_deallocate(Y2)
call mma_deallocate(Z2)
call mma_deallocate(S2)

return

end subroutine check_S_square

subroutine check_hermiticity_moment(n,moment,dbg)

use Constants, only: cZero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n
complex(kind=wp), intent(in) :: moment(3,n,n)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, j, l
complex(kind=wp) :: c

! build difference SUM ( M ATRIX(i,j) - CONJG(MATRIX(j,i)) )
c = cZero
do i=1,n
  do j=1,n
    if (i == j) cycle
    do l=1,3
      c = c+moment(l,i,j)-conjg(moment(l,j,i))
    end do
  end do
end do
if (dbg) write(u6,'(A,2ES22.14)') 'check_hermiticity_moment::  trace of A(i,j)-CONJG(A(j,i)) = ',c
if (abs(c) > 1.0e-6_wp) then
  call WarningMessage(1,'check_hermiticity_moment:: trace of M(:,i,j)-CONJG(A(:,j,i)) is larger than 1.0e-6. '// &
                      'The hermiticity of input moment is not quite fulfilled')
else
  write(u6,'(A,ES22.14)') 'check_hermiticity_moment:  The input moment passes the hermiticity test.'
end if

return

end subroutine check_hermiticity_moment

subroutine check_hermiticity_matrix(n,matrix,dbg)

use Constants, only: cZero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n
complex(kind=wp), intent(in) :: matrix(n,n)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, j
complex(kind=wp) :: c

! build difference ( MATRIX(i,j) - CONJG(MATRIX(j,i)) )
c = cZero
do i=1,n
  do j=i,n
    if (i == j) cycle
    c = c+(matrix(i,j)-conjg(matrix(j,i)))
  end do
end do
if (dbg) write(u6,'(A,2ES22.14)') 'check_hermiticity_matrix::  trace of A(i,j)-CONJG(A(j,i)) = ',c
if (abs(c) > 1.0e-6_wp) then
  call WarningMessage(1,'check_hermiticity_matrix:: trace of A(i,j)-CONJG(A(j,i)) is larger than 1.0e-6. '// &
                      'The hermiticity of input matrix is not quite fulfilled')
else
  write(u6,'(A,ES22.14)') 'check_hermiticity_matrix:  The input matrix passes the hermiticity test.'
end if

return

end subroutine check_hermiticity_matrix

!----------------------------------------------------------------------!
!           HIGH LEVEL READING SUBROUTINES
!----------------------------------------------------------------------!
! read_magnetic_moment
! read_electric_moment
! read_spin_moment
! read_amfi
! read_angmom
! read_edipmom
! read_eigen
! read_hso
! read_eso
! read_esfs
! read_szproj
! read_nss
! read_nstate
!----------------------------------------------------------------------!

subroutine read_magnetic_moment(DATA_FILE,n,moment,dbg)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Ten, cZero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: DATA_FILE, N
complex(kind=wp), intent(out) :: moment(3,N,N)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, j
real(kind=wp), allocatable :: ri(:,:), rr(:,:)
complex(kind=wp), allocatable :: M(:,:)
real(kind=wp), parameter :: MINIMAL_REAL = tiny(MINIMAL_REAL)*Ten
real(kind=wp), external :: dnrm2_, dznrm2_
logical(kind=iwp), external :: inquire_key_presence

moment(:,:,:) = cZero
call mma_allocate(rr,n,n,label='rr')
call mma_allocate(ri,n,n,label='ri')
! projection X
rr = Zero
ri = Zero
if (inquire_key_presence(DATA_FILE,'$magn_xr')) call read_2d_real_array(DATA_FILE,'$magn_xr',n,n,rr,dbg)
if (inquire_key_presence(DATA_FILE,'$magn_xi')) call read_2d_real_array(DATA_FILE,'$magn_xi',n,n,ri,dbg)
if (dbg) then
  write(u6,*) 'read_magnetic_moment::  norm of moment_xr=',dnrm2_(n*n,rr,1)
  write(u6,*) 'read_magnetic_moment::  norm of moment_xi=',dnrm2_(n*n,ri,1)
end if
moment(1,:,:) = cmplx(rr(:,:),ri(:,:),kind=wp)
if (dbg) then
  call mma_allocate(M,n,n,label='M')
  M(:,:) = moment(1,:,:)
  call check_hermiticity_matrix(n,M,dbg)
end if
! projection Y
rr = Zero
ri = Zero
if (inquire_key_presence(DATA_FILE,'$magn_yr')) call read_2d_real_array(DATA_FILE,'$magn_yr',n,n,rr,dbg)
if (inquire_key_presence(DATA_FILE,'$magn_yi')) call read_2d_real_array(DATA_FILE,'$magn_yi',n,n,ri,dbg)
if (dbg) then
  write(u6,*) 'read_magnetic_moment::  norm of moment_yr=',dnrm2_(n*n,rr,1)
  write(u6,*) 'read_magnetic_moment::  norm of moment_yi=',dnrm2_(n*n,ri,1)
end if
moment(2,:,:) = cmplx(rr(:,:),ri(:,:),kind=wp)
if (dbg) then
  M(:,:) = moment(2,:,:)
  call check_hermiticity_matrix(n,M,dbg)
end if
! projection Z
rr = Zero
ri = Zero
if (inquire_key_presence(DATA_FILE,'$magn_zr')) call read_2d_real_array(DATA_FILE,'$magn_zr',n,n,rr,dbg)
if (inquire_key_presence(DATA_FILE,'$magn_zi')) call read_2d_real_array(DATA_FILE,'$magn_zi',n,n,ri,dbg)
if (dbg) then
  write(u6,*) 'read_magnetic_moment::  norm of moment_zr=',dnrm2_(n*n,rr,1)
  write(u6,*) 'read_magnetic_moment::  norm of moment_zi=',dnrm2_(n*n,ri,1)
end if
moment(3,:,:) = cmplx(rr(:,:),ri(:,:),kind=wp)
if (dznrm2_(3*n*n,moment,1) <= MINIMAL_REAL) call WarningMessage(1,'read_magnetic_moment:: the norm of the read moment is zero!')
if (dbg) then
  M(:,:) = moment(3,:,:)
  call check_hermiticity_matrix(n,M,dbg)
  call mma_deallocate(M)
end if
call mma_deallocate(rr)
call mma_deallocate(ri)
if (dbg) then
  write(u6,*) 'read_magnetic_moment::  ELECTRIC MOMENT at the end of the suboutine'
  write(u6,*) 'projection X'
  do i=1,n
    write(u6,'(100(2F10.6,2x))') (moment(1,i,j),j=1,n)
  end do
  write(u6,*) 'projection Y'
  do i=1,n
    write(u6,'(100(2F10.6,2x))') (moment(2,i,j),j=1,n)
  end do
  write(u6,*) 'projection Z'
  do i=1,n
    write(u6,'(100(2F10.6,2x))') (moment(3,i,j),j=1,n)
  end do
end if

return

end subroutine read_magnetic_moment

subroutine read_electric_moment(DATA_FILE,n,moment,dbg)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Ten, cZero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: DATA_FILE, N
complex(kind=wp), intent(out) :: moment(3,N,N)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, j
real(kind=wp), allocatable :: ri(:,:), rr(:,:)
complex(kind=wp), allocatable :: M(:,:)
real(kind=wp), parameter :: MINIMAL_REAL = tiny(MINIMAL_REAL)*Ten
real(kind=wp), external :: dnrm2_, dznrm2_
logical(kind=iwp), external :: inquire_key_presence

moment(:,:,:) = cZero
call mma_allocate(rr,n,n,label='rr')
call mma_allocate(ri,n,n,label='ri')
! projection X
rr = Zero
ri = Zero
if (inquire_key_presence(DATA_FILE,'$edipm_xr')) call read_2d_real_array(DATA_FILE,'$edipm_xr',n,n,rr,dbg)
if (inquire_key_presence(DATA_FILE,'$edipm_xi')) call read_2d_real_array(DATA_FILE,'$edipm_xi',n,n,ri,dbg)
if (dbg) then
  write(u6,*) 'read_electric_moment::  norm of moment_xr=',dnrm2_(n*n,rr,1)
  write(u6,*) 'read_electric_moment::  norm of moment_xi=',dnrm2_(n*n,ri,1)
end if
moment(1,:,:) = cmplx(rr(:,:),ri(:,:),kind=wp)
if (dbg) then
  call mma_allocate(M,n,n,label='M')
  M(:,:) = moment(1,:,:)
  call check_hermiticity_matrix(n,M,dbg)
end if
! projection Y
rr = Zero
ri = Zero
if (inquire_key_presence(DATA_FILE,'$edipm_yr')) call read_2d_real_array(DATA_FILE,'$edipm_yr',n,n,rr,dbg)
if (inquire_key_presence(DATA_FILE,'$edipm_yi')) call read_2d_real_array(DATA_FILE,'$edipm_yi',n,n,ri,dbg)
if (dbg) then
  write(u6,*) 'read_electric_moment::  norm of moment_yr=',dnrm2_(n*n,rr,1)
  write(u6,*) 'read_electric_moment::  norm of moment_yi=',dnrm2_(n*n,ri,1)
end if
moment(2,:,:) = cmplx(rr(:,:),ri(:,:),kind=wp)
if (dbg) then
  M(:,:) = moment(2,:,:)
  call check_hermiticity_matrix(n,M,dbg)
end if
! projection Z
rr = Zero
ri = Zero
if (inquire_key_presence(DATA_FILE,'$edipm_zr')) call read_2d_real_array(DATA_FILE,'$edipm_zr',n,n,rr,dbg)
if (inquire_key_presence(DATA_FILE,'$edipm_zi')) call read_2d_real_array(DATA_FILE,'$edipm_zi',n,n,ri,dbg)
if (dbg) then
  write(u6,*) 'read_electric_moment::  norm of moment_zr=',dnrm2_(n*n,rr,1)
  write(u6,*) 'read_electric_moment::  norm of moment_zi=',dnrm2_(n*n,ri,1)
end if
moment(3,:,:) = cmplx(rr(:,:),ri(:,:),kind=wp)
if (dznrm2_(3*n*n,moment,1) <= MINIMAL_REAL) call WarningMessage(1,'read_electric_moment:: the norm of the read moment is zero!')
if (dbg) then
  M(:,:) = moment(3,:,:)
  call check_hermiticity_matrix(n,M,dbg)
  call mma_deallocate(M)
end if
call mma_deallocate(rr)
call mma_deallocate(ri)
if (dbg) then
  write(u6,*) 'read_electric_moment::  ELECTRIC MOMENT at the end of the suboutine'
  write(u6,*) 'projection X'
  do i=1,n
    write(u6,'(100(2F10.6,2x))') (moment(1,i,j),j=1,n)
  end do
  write(u6,*) 'projection Y'
  do i=1,n
    write(u6,'(100(2F10.6,2x))') (moment(2,i,j),j=1,n)
  end do
  write(u6,*) 'projection Z'
  do i=1,n
    write(u6,'(100(2F10.6,2x))') (moment(3,i,j),j=1,n)
  end do
end if

return

end subroutine read_electric_moment

subroutine read_spin_moment(DATA_FILE,n,moment,dbg)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Ten, cZero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: DATA_FILE, N
complex(kind=wp), intent(out) :: moment(3,N,N)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, j
real(kind=wp), allocatable :: ri(:,:), rr(:,:)
complex(kind=wp), allocatable :: M(:,:)
real(kind=wp), parameter :: MINIMAL_REAL = tiny(MINIMAL_REAL)*Ten
real(kind=wp), external :: dnrm2_, dznrm2_
logical(kind=iwp), external :: inquire_key_presence

if (dbg) write(u6,*) 'ENTER read_spin_moment'
moment(:,:,:) = cZero
call mma_allocate(rr,n,n,label='rr')
call mma_allocate(ri,n,n,label='ri')
! projection X
rr = Zero
ri = Zero
if (dbg) then
  write(u6,*) 'ENTER read_spin_moment p2  :spin_xr: ',inquire_key_presence(DATA_FILE,'$spin_xr')
  write(u6,*) 'ENTER read_spin_moment p2  :spin_xi: ',inquire_key_presence(DATA_FILE,'$spin_xi')
  write(u6,*) 'ENTER read_spin_moment p2  :spin_yr: ',inquire_key_presence(DATA_FILE,'$spin_yr')
  write(u6,*) 'ENTER read_spin_moment p2  :spin_yi: ',inquire_key_presence(DATA_FILE,'$spin_yi')
  write(u6,*) 'ENTER read_spin_moment p2  :spin_zr: ',inquire_key_presence(DATA_FILE,'$spin_zr')
  write(u6,*) 'ENTER read_spin_moment p2  :spin_zi: ',inquire_key_presence(DATA_FILE,'$spin_zi')
  call xFlush(u6)
end if
if (inquire_key_presence(DATA_FILE,'$spin_xr')) call read_2d_real_array(DATA_FILE,'$spin_xr',n,n,rr,dbg)
if (dbg) then
  write(u6,*) 'ENTER read_spin_moment p3'
  call xFlush(u6)
end if
if (inquire_key_presence(DATA_FILE,'$spin_xi')) call read_2d_real_array(DATA_FILE,'$spin_xi',n,n,ri,dbg)
if (dbg) then
  write(u6,*) 'ENTER read_spin_moment p4'
  write(u6,*) 'read_spin_moment::  norm of moment_zr=',dnrm2_(n*n,rr,1)
  write(u6,*) 'read_spin_moment::  norm of moment_zi=',dnrm2_(n*n,ri,1)
  write(u6,*) 'ENTER read_spin_moment p5'
  call xFlush(u6)
end if
moment(1,:,:) = cmplx(rr(:,:),ri(:,:),kind=wp)
if (dbg) then
  write(u6,*) 'ENTER read_spin_moment p6'
  call mma_allocate(M,n,n,label='M')
  M(:,:) = moment(1,:,:)
  call check_hermiticity_matrix(n,M,dbg)
  call xFlush(u6)
end if
! projection Y
rr = Zero
ri = Zero
if (inquire_key_presence(DATA_FILE,'$spin_yr')) call read_2d_real_array(DATA_FILE,'$spin_yr',n,n,rr,dbg)
if (inquire_key_presence(DATA_FILE,'$spin_yi')) call read_2d_real_array(DATA_FILE,'$spin_yi',n,n,ri,dbg)
if (dbg) then
  write(u6,*) 'read_spin_moment::  norm of moment_zr=',dnrm2_(n*n,rr,1)
  write(u6,*) 'read_spin_moment::  norm of moment_zi=',dnrm2_(n*n,ri,1)
  call xFlush(u6)
end if
moment(2,:,:) = cmplx(rr(:,:),ri(:,:),kind=wp)
if (dbg) then
  M(:,:) = moment(2,:,:)
  call check_hermiticity_matrix(n,M,dbg)
end if
! projection Z
rr = Zero
ri = Zero
if (inquire_key_presence(DATA_FILE,'$spin_zr')) call read_2d_real_array(DATA_FILE,'$spin_zr',n,n,rr,dbg)
if (inquire_key_presence(DATA_FILE,'$spin_zi')) call read_2d_real_array(DATA_FILE,'$spin_zi',n,n,ri,dbg)
if (dbg) then
  write(u6,*) 'read_spin_moment::  norm of moment_zr=',dnrm2_(n*n,rr,1)
  write(u6,*) 'read_spin_moment::  norm of moment_zi=',dnrm2_(n*n,ri,1)
  call xFlush(u6)
end if
moment(3,:,:) = cmplx(rr(:,:),ri(:,:),kind=wp)
if (dznrm2_(3*n*n,moment,1) <= MINIMAL_REAL) call WarningMessage(1,'read_spin:: the norm of the read moment is zero!')
if (dbg) then
  M(:,:) = moment(3,:,:)
  call check_hermiticity_matrix(n,M,dbg)
  call mma_deallocate(M)
end if
call mma_deallocate(rr)
call mma_deallocate(ri)
if (dbg) call check_commutation(n,moment,dbg)

if (dbg) then
  write(u6,*) 'read_spin_moment::  SPIN MOMENT at the end of the suboutine'
  write(u6,*) 'projection X'
  do i=1,n
    write(u6,'(100(2F10.6,2x))') (moment(1,i,j),j=1,n)
  end do
  write(u6,*) 'projection Y'
  do i=1,n
    write(u6,'(100(2F10.6,2x))') (moment(2,i,j),j=1,n)
  end do
  write(u6,*) 'projection Z'
  do i=1,n
    write(u6,'(100(2F10.6,2x))') (moment(3,i,j),j=1,n)
  end do
end if
if (dbg) write(u6,*) 'EXIT read_spin_moment'

return

end subroutine read_spin_moment

subroutine read_angmom(DATA_FILE,n,moment,dbg)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Ten
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: DATA_FILE, N
real(kind=wp), intent(out) :: moment(3,N,N)
logical(kind=iwp), intent(in) :: dbg
real(kind=wp), allocatable :: rr(:,:)
real(kind=wp), parameter :: MINIMAL_REAL = tiny(MINIMAL_REAL)*Ten
real(kind=wp), external :: dnrm2_
logical(kind=iwp), external :: inquire_key_presence

moment(:,:,:) = Zero
call mma_allocate(rr,n,n,label='rr')
! projection X
rr = Zero
if (inquire_key_presence(DATA_FILE,'$angmom_x')) call read_2d_real_array(DATA_FILE,'$angmom_x',n,n,rr,dbg)
if (dbg) write(u6,*) 'read_angmom::  norm of moment_x=',dnrm2_(n*n,rr,1)
moment(1,:,:) = rr(:,:)
! projection Y
rr = Zero
if (inquire_key_presence(DATA_FILE,'$angmom_y')) call read_2d_real_array(DATA_FILE,'$angmom_y',n,n,rr,dbg)
if (dbg) write(u6,*) 'read_angmom::  norm of moment_y=',dnrm2_(n*n,rr,1)
moment(2,:,:) = rr(:,:)
! projection Z
rr = Zero
if (inquire_key_presence(DATA_FILE,'$angmom_z')) call read_2d_real_array(DATA_FILE,'$angmom_z',n,n,rr,dbg)
if (dbg) write(u6,*) 'read_angmom::  norm of moment_z=',dnrm2_(n*n,rr,1)
moment(3,:,:) = rr(:,:)
call mma_deallocate(rr)
if (dnrm2_(3*n*n,moment,1) <= MINIMAL_REAL) call WarningMessage(1,'read_angmom:: the norm of the read moment is zero!')

return

end subroutine read_angmom

subroutine read_edipmom(DATA_FILE,n,moment,dbg)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Ten
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: DATA_FILE, N
real(kind=wp), intent(out) :: moment(3,N,N)
logical(kind=iwp), intent(in) :: dbg
real(kind=wp), allocatable :: rr(:,:)
real(kind=wp), parameter :: MINIMAL_REAL = tiny(MINIMAL_REAL)*Ten
real(kind=wp), external :: dnrm2_
logical(kind=iwp), external :: inquire_key_presence

moment(:,:,:) = Zero
call mma_allocate(rr,n,n,label='rr')
! projection X
rr = Zero
if (inquire_key_presence(DATA_FILE,'$edmom_x')) call read_2d_real_array(DATA_FILE,'$edmom_x',n,n,rr,dbg)
if (dbg) write(u6,*) 'read_edipmom::  norm of moment_x=',dnrm2_(n*n,rr,1)
moment(1,:,:) = rr(:,:)
! projection Y
rr = Zero
if (inquire_key_presence(DATA_FILE,'$edmom_y')) call read_2d_real_array(DATA_FILE,'$edmom_y',n,n,rr,dbg)
if (dbg) write(u6,*) 'read_edipmom::  norm of moment_y=',dnrm2_(n*n,rr,1)
moment(2,:,:) = rr(:,:)
! projection Z
rr = Zero
if (inquire_key_presence(DATA_FILE,'$edmom_z')) call read_2d_real_array(DATA_FILE,'$edmom_z',n,n,rr,dbg)
if (dbg) write(u6,*) 'read_edipmom::  norm of moment_z=',dnrm2_(n*n,rr,1)
moment(3,:,:) = rr(:,:)
call mma_deallocate(rr)
if (dnrm2_(3*n*n,moment,1) <= MINIMAL_REAL) call WarningMessage(1,'read_edipmom:: the norm of the read moment is zero!')

return

end subroutine read_edipmom

subroutine read_amfi(DATA_FILE,n,moment,dbg)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Ten
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: DATA_FILE, N
real(kind=wp), intent(out) :: moment(3,N,N)
logical(kind=iwp), intent(in) :: dbg
real(kind=wp), allocatable :: rr(:,:)
real(kind=wp), parameter :: MINIMAL_REAL = tiny(MINIMAL_REAL)*Ten
real(kind=wp), external :: dnrm2_
logical(kind=iwp), external :: inquire_key_presence

moment(:,:,:) = Zero
call mma_allocate(rr,n,n,label='rr')
! projection X
rr = Zero
if (inquire_key_presence(DATA_FILE,'$amfi_x')) call read_2d_real_array(DATA_FILE,'$amfi_x',n,n,rr,dbg)
if (dbg) write(u6,*) 'read_amfi::  norm of moment_x=',dnrm2_(n*n,rr,1)
moment(1,:,:) = rr(:,:)
! projection Y
rr = Zero
if (inquire_key_presence(DATA_FILE,'$amfi_y')) call read_2d_real_array(DATA_FILE,'$amfi_y',n,n,rr,dbg)
if (dbg) write(u6,*) 'read_amfi::  norm of moment_y=',dnrm2_(n*n,rr,1)
moment(2,:,:) = rr(:,:)
! projection Z
rr = Zero
if (inquire_key_presence(DATA_FILE,'$amfi_z')) call read_2d_real_array(DATA_FILE,'$amfi_z',n,n,rr,dbg)
if (dbg) write(u6,*) 'read_amfi::  norm of moment_z=',dnrm2_(n*n,rr,1)
moment(3,:,:) = rr(:,:)
call mma_deallocate(rr)
if (dnrm2_(3*n*n,moment,1) <= MINIMAL_REAL) call WarningMessage(1,'read_amfi:: the norm of the read moment is zero!')

return

end subroutine read_amfi

subroutine read_nss(DATA_FILE,n,dbg)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: DATA_FILE
integer(kind=iwp), intent(out) :: n
logical(kind=iwp), intent(in) :: dbg
logical(kind=iwp), external :: inquire_key_presence

n = 0
if (inquire_key_presence(DATA_FILE,'$nss')) call read_INTEGER_scalar(DATA_FILE,'$nss',n,dbg)
if (n <= 0) call WarningMessage(1,'read_nss:: nss value in DATA_FILE = 0. Is it really the case?')

return

end subroutine read_nss

subroutine read_nstate(DATA_FILE,n,dbg)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: DATA_FILE
integer(kind=iwp), intent(out) :: n
logical(kind=iwp), intent(in) :: dbg
logical(kind=iwp), external :: inquire_key_presence

n = 0
if (inquire_key_presence(DATA_FILE,'$nstate')) call read_INTEGER_scalar(DATA_FILE,'$nstate',n,dbg)
if (n <= 0) call WarningMessage(1,'read_nstate:: nstate value in DATA_FILE = 0. Is it really the case?')

return

end subroutine read_nstate

subroutine read_nmult(DATA_FILE,n,dbg)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: DATA_FILE
integer(kind=iwp), intent(out) :: n
logical(kind=iwp), intent(in) :: dbg
logical(kind=iwp), external :: inquire_key_presence

n = 0
if (inquire_key_presence(DATA_FILE,'$nmult')) call read_INTEGER_scalar(DATA_FILE,'$nmult',n,dbg)
if (n <= 0) call WarningMessage(1,'read_nmult:: nmult value in DATA_FILE = 0. Is it really the case?')

return

end subroutine read_nmult

subroutine read_multiplicity(DATA_FILE,n,array,dbg)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: DATA_FILE, n
integer(kind=iwp), intent(out) :: array(n)
logical(kind=iwp), intent(in) :: dbg
logical(kind=iwp), external :: inquire_key_presence

array(:) = 0
if (inquire_key_presence(DATA_FILE,'$multiplicity')) call read_1d_INTEGER_array(DATA_FILE,'$multiplicity',n,array,dbg)
if (sum(abs(array)) == 0) then
  call WarningMessage(1,'read_multiplicity:: it seems that all the multiplicities in DATA_FILE are 0. Is it really the case?')
  write(u6,*) 'read_multiplicity:: SUM(Sz) = ',sum(abs(array))
end if
if (sum(array) == 0) then
  call WarningMessage(1,'read_multiplicity:: it seems that all the multiplicities in DATA_FILE are 0. Is it really the case?')
  write(u6,*) 'read_szproj:: SUM(Sz) = ',sum(array)
end if

return

end subroutine read_multiplicity

subroutine read_imult(DATA_FILE,n,array,dbg)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: DATA_FILE, n
integer(kind=iwp), intent(out) :: array(n)
logical(kind=iwp), intent(in) :: dbg
logical(kind=iwp), external :: inquire_key_presence

array(:) = 0
if (inquire_key_presence(DATA_FILE,'$imult')) call read_1d_INTEGER_array(DATA_FILE,'$imult',n,array,dbg)
if (sum(array) == 0) then
  call WarningMessage(1,'read_imult:: it seems that all the multiplicities in DATA_FILE are 0. Is it really the case?')
  write(u6,*) 'read_imult:: SUM(mult()) = ',sum(array)
end if

return

end subroutine read_imult

subroutine read_format(DATA_FILE,n,dbg)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: DATA_FILE
integer(kind=iwp), intent(out) :: n
logical(kind=iwp), intent(in) :: dbg
logical(kind=iwp), external :: inquire_key_presence

n = 0
if (inquire_key_presence(DATA_FILE,'$format')) call read_INTEGER_scalar(DATA_FILE,'$format',n,dbg)
if (n <= 0) call WarningMessage(2,'read_format:: FORMAT value in DATA_FILE = 0. The FORMAT must be equal or larger than 2020. '// &
                                'Please check.')

return

end subroutine read_format

subroutine read_nroot(DATA_FILE,n,array,dbg)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: DATA_FILE, n
integer(kind=iwp), intent(out) :: array(n)
logical(kind=iwp), intent(in) :: dbg
logical(kind=iwp), external :: inquire_key_presence

array(:) = 0
if (inquire_key_presence(DATA_FILE,'$nroot')) call read_1d_INTEGER_array(DATA_FILE,'$nroot',n,array,dbg)
if (sum(array) == 0) then
  call WarningMessage(1,'read_nroot:: it seems that the number of roots included in spin-orbit interaction in DATA_FILE are 0. '// &
                      'Is it really the case?')
  write(u6,*) 'read_szproj:: SUM(array()) = ',sum(array)
end if

return

end subroutine read_nroot

subroutine read_szproj(DATA_FILE,n,array,dbg)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: DATA_FILE, n
integer(kind=iwp), intent(out) :: array(n)
logical(kind=iwp), intent(in) :: dbg
logical(kind=iwp), external :: inquire_key_presence

array(:) = 0
if (inquire_key_presence(DATA_FILE,'$szproj')) call read_1d_INTEGER_array(DATA_FILE,'$szproj',n,array,dbg)
if (sum(abs(array)) == 0) then
  call WarningMessage(1,'read_szproj:: it seems that SUM(ABS(Sz)) in DATA_FILE is 0. Is it really the case?')
  write(u6,*) 'read_szproj:: SUM(ABS(Sz)) = ',sum(abs(array))
end if
if (sum(array) /= 0) then
  call WarningMessage(1,'read_szproj:: it seems that SUM(Sz) in DATA_FILE is not 0. Is it really the case?')
  write(u6,*) 'read_szproj:: SUM(Sz) = ',sum(array)
end if

return

end subroutine read_szproj

subroutine read_eso(DATA_FILE,n,array,dbg)

use Constants, only: Zero, Ten
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: DATA_FILE, n
real(kind=wp), intent(out) :: array(n)
logical(kind=iwp), intent(in) :: dbg
real(kind=wp), parameter :: MINIMAL_REAL = tiny(MINIMAL_REAL)*Ten
real(kind=wp), external :: dnrm2_
logical(kind=iwp), external :: inquire_key_presence

array(:) = Zero
if (inquire_key_presence(DATA_FILE,'$eso')) call read_1d_real_array(DATA_FILE,'$eso',n,array,dbg)
if (dbg) write(u6,*) 'read_eso::  norm of eso=',dnrm2_(n,array,1)
if (dnrm2_(n,array,1) <= MINIMAL_REAL) then
  call WarningMessage(1,'read_eso:: it seems that the norm of ESO array in DATA_FILE is 0. Is it really the case?')
  write(u6,*) 'read_eso:: dnrm2_(eso) = ',dnrm2_(n,array,1)
end if

return

end subroutine read_eso

subroutine read_esfs(DATA_FILE,n,array,dbg)

use Constants, only: Zero, Ten
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: DATA_FILE, n
real(kind=wp), intent(out) :: array(n)
logical(kind=iwp), intent(in) :: dbg
real(kind=wp), parameter :: MINIMAL_REAL = tiny(MINIMAL_REAL)*Ten
real(kind=wp), external :: dnrm2_
logical(kind=iwp), external :: inquire_key_presence

array(:) = Zero
if (inquire_key_presence(DATA_FILE,'$esfs')) call read_1d_real_array(DATA_FILE,'$esfs',n,array,dbg)
if (dbg) write(u6,*) 'read_esfs::  norm of esfs=',dnrm2_(n,array,1)
if (dnrm2_(n,array,1) <= MINIMAL_REAL) then
  call WarningMessage(1,'read_esfs:: it seems that the norm of ESFS in DATA_FILE is 0. Is it really the case?')
  write(u6,*) 'read_esfs:: dnrm2_(esfs) = ',dnrm2_(n,array,1)
end if

return

end subroutine read_esfs

subroutine read_hso(DATA_FILE,n,array,dbg)

use Constants, only: Ten, cZero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: DATA_FILE, n
complex(kind=wp), intent(out) :: array(n)
logical(kind=iwp), intent(in) :: dbg
real(kind=wp), parameter :: MINIMAL_REAL = tiny(MINIMAL_REAL)*Ten
real(kind=wp), external :: dznrm2_
logical(kind=iwp), external :: inquire_key_presence

array(:) = cZero
if (inquire_key_presence(DATA_FILE,'$hso')) call read_complex_matrix(DATA_FILE,'$hso',n,array,dbg)
if (dbg) write(u6,*) 'read_hso::  norm of hso=',dznrm2_(n*n,array,1)
if (dznrm2_(n*n,array,1) <= MINIMAL_REAL) then
  call WarningMessage(1,'read_hso:: it seems that norm of HSO in DATA_FILE is 0. Is it really the case?')
  write(u6,*) 'read_hso:: dznrm2_(hso) = ',dznrm2_(n*n,array,1)
end if

return

end subroutine read_hso

subroutine read_eigen(DATA_FILE,n,array,dbg)

use Constants, only: Ten, cZero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: DATA_FILE, n
complex(kind=wp), intent(out) :: array(n)
logical(kind=iwp), intent(in) :: dbg
real(kind=wp), parameter :: MINIMAL_REAL = tiny(MINIMAL_REAL)*Ten
real(kind=wp), external :: dznrm2_
logical(kind=iwp), external :: inquire_key_presence

array(:) = cZero
if (inquire_key_presence(DATA_FILE,'$eigen')) call read_complex_matrix(DATA_FILE,'$eigen',n,array,dbg)
if (dbg) write(u6,*) 'read_eigen::  norm of eigenv=',dznrm2_(n*n,array,1)
if (dznrm2_(n*n,array,1) <= MINIMAL_REAL) then
  call WarningMessage(1,'read_eigen:: it seems that norm of EIGENV in DATA_FILE is 0. Is it really the case?')
  write(u6,*) 'read_eigen:: dznrm2_(array) = ',dznrm2_(n*n,array,1)
end if

return

end subroutine read_eigen

subroutine read_gtens(DATA_FILE,nmult,gtens,axes,dbg)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: DATA_FILE, nmult
real(kind=wp), intent(out) :: gtens(nmult,3), axes(nmult,3,3)
logical(kind=iwp), intent(in) :: dbg

call read_2d_real_array(DATA_FILE,'$gtens_main',nmult,3,gtens,dbg)
call read_3d_real_array(DATA_FILE,'$gtens_axes',nmult,3,3,axes,dbg)

return

end subroutine read_gtens

subroutine read_stev_cfp(DATA_FILE,s,n,cfp,dbg)

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: DATA_FILE, n   ! 2J+1, or 2L+1, i.e. the dimension of the J or L multiplet
character(len=*), intent(in) :: s
real(kind=wp), intent(out) :: cfp(n-1,-(n-1):n-1)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr, ik, iq, k, q
character(len=500) :: line

ierr = 0
if ((n <= 0)) then
  call WarningMessage(1,'read_stev_cfp_'//trim(s)//':: nothing to read. Array size = 0.')
  return
end if
CFP(:,:) = Zero

rewind(DATA_FILE)
call file_advance_to_string(DATA_FILE,'$stev_cfp_'//trim(s),line,ierr,dbg)

read(DATA_FILE,*,iostat=ierr) i
if (i /= n) call WarningMessage(2,'read_stev_cfp_'//trim(s)//':: size of the multiplet is not the same i/=n')

! if key is found, THEN read the data
if (ierr == 0) then
  do k=2,n-1,2
    do q=-k,k,2
      read(DATA_FILE,*,iostat=ierr) ik,iq,cfp(ik,iq)
      if (ierr /= 0) call WarningMessage(2,'read_stev_cfp_'//trim(s)//':: Something went wrong reading the array.')
      if (dbg) write(u6,*) 'read_stev_cfp_'//trim(s)//'::  k, q =',k,q
    end do
  end do
end if

return

end subroutine read_stev_cfp

subroutine read_susc(DATA_FILE,s,n,field,zj,t,x,x_tens,dbg)
! n     : number of temperature points
! zj    : intermolecular interaction
! field : applied field
! t     : temperature points
! x     : susceptibility X
! x_tens: susceptibility tensor, X_tens

use Constants, only: Zero, Ten
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: DATA_FILE
integer(kind=iwp), intent(inout) :: n
character(len=*), intent(in) :: s
real(kind=wp), intent(out) :: field, zj, t(n), x(n), x_tens(n,3,3)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr, j, k
character(len=500) :: line
real(kind=wp), parameter :: MINIMAL_REAL = tiny(MINIMAL_REAL)*Ten
real(kind=wp), external :: dnrm2_

t(:) = Zero
x(:) = Zero
x_tens(:,:,:) = Zero
! s takes the values:
! 'x'
! 'x_field'
! 'x_minus_one'
! 'x_field_minus_one'
! 'xt'
! 'xt_field'
! 'xt_minus_one'
! 'xt_field_minus_one'
! 'x_zj'
! 'x_field_zj'
! 'x_minus_one_zj'
! 'x_field_minus_one_zj'
! 'xt_zj'
! 'xt_field_zj'
! 'xt_minus_one_zj'
! 'xt_field_minus_one_zj'
!if ((trim(s) /= 'x')                    .or. &
!    (trim(s) /= 'x_field')              .or. &
!    (trim(s) /= 'x_minus_one')          .or. &
!    (trim(s) /= 'x_field_minus_one')    .or. &
!    (trim(s) /= 'xt')                   .or. &
!    (trim(s) /= 'xt_field')             .or. &
!    (trim(s) /= 'xt_minus_one')         .or. &
!    (trim(s) /= 'xt_field_minus_one')   .or. &
!    (trim(s) /= 'x_zj')                 .or. &
!    (trim(s) /= 'x_field_zj')           .or. &
!    (trim(s) /= 'x_minus_one_zj')       .or. &
!    (trim(s) /= 'x_field_minus_one_zj') .or. &
!    (trim(s) /= 'xt_zj')                .or. &
!    (trim(s) /= 'xt_field_zj')          .or. &
!    (trim(s) /= 'xt_minus_one_zj')      .or. &
!    (trim(s) /= 'xt_field_minus_one_zj')) then
!
!  call WarningMessage(1,'read_x '//trim(s)//' :: the parameter s='//trim(s)//'is not understood. RETURN without read of X data.')
!  return
!end if

ierr = 0
if ((n <= 0)) then
  call WarningMessage(1,'read_x '//trim(s)//' :: nothing to read. Array size = 0.')
  return
end if

rewind(DATA_FILE)
call file_advance_to_string(DATA_FILE,'$susceptibility_'//trim(s),line,ierr,dbg)

! if key is found, THEN read the data
if (ierr == 0) then
  read(DATA_FILE,*,iostat=ierr) n
  read(DATA_FILE,*,iostat=ierr) zj,field
  if (ierr /= 0) call WarningMessage(2,'read_x '//trim(s)//' :: Something went wrong reading the zJ and field values.')
  read(DATA_FILE,*,iostat=ierr) (T(i),i=1,n)
  if (dnrm2_(n,T,1) < MINIMAL_REAL) call WarningMessage(1,'read_x '//trim(s)//' :: all array T elements are zero = 0.')
  if (ierr /= 0) call WarningMessage(2,'read_x '//trim(s)//' :: Something went wrong reading the T array.')
  read(DATA_FILE,*,iostat=ierr) (X(i),i=1,n)
  if (dnrm2_(n,X,1) < MINIMAL_REAL) call WarningMessage(1,'read_x '//trim(s)//' :: all array X elements are zero = 0.')
  if (ierr /= 0) call WarningMessage(2,'read_x '//trim(s)//' :: Something went wrong reading the X array.')
  do j=1,3
    do k=1,3
      read(DATA_FILE,*,iostat=ierr) (X_tens(i,j,k),i=1,n)
      if (ierr /= 0) call WarningMessage(2,'read_x '//trim(s)//' :: Something went wrong reading the X_tens array.')
    end do
  end do
  if (dnrm2_(n*3*3,X_tens,1) < MINIMAL_REAL) call WarningMessage(1,'read_x '//trim(s)//' :: all array X_tens elements are zero.')
  if (dbg) call xFlush(u6)
else
  write(u6,*) 'keyword $susceptibility_'//trim(s)//' was not found in DATA_FILE'
end if

return

end subroutine read_susc

subroutine read_magn(DATA_FILE,nt,nh,nd,nss,zj,t,h,x,y,z,w,m,mav,energy,dbg)
! nt        :  number of temperature points
! nh        :  number of field points
! nd        :  number of directions of aplied field
! nss       :  number of spin-orbit states
! zj        :  inter-molecular parameter zJ
! t         :  temperature points
! h         :  field points
! x, y, z, w:  Lebedev grid directions and weight
! m         :  magnetisation vector
! mav       :  average magnetisation vector
! energy    :  Zeeman energy states

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: DATA_FILE, nt, nh, nd, nss
real(kind=wp), intent(out) :: zj, t(nt), h(nh), x(nd), y(nd), z(nd), w(nd), m(nd,3,nt,nh), mav(nt,nh), energy(nd,nh,nss)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, id, ierr, ih, iss, it, l
character(len=500) :: line

ierr = 0
if ((nt <= 0) .or. (nh <= 0) .or. (nd <= 0) .or. (nss <= 0)) then
  call WarningMessage(1,'read_magn :: nothing to read. Array size = 0.')
  return
end if

rewind(DATA_FILE)
call file_advance_to_string(DATA_FILE,'$magnetisation',line,ierr,dbg)

read(DATA_FILE,*,iostat=ierr) it,ih,id,iss
if (it /= nt) call WarningMessage(1,'read_magn :: nt read from DATA_FILE is not the same as the parameter used to CALL '// &
                                  'this function.')
if (ih /= nh) call WarningMessage(1,'read_magn :: nh read from DATA_FILE is not the same as the parameter used to CALL '// &
                                  'this function.')
if (id /= nd) call WarningMessage(1,'read_magn :: nd read from DATA_FILE is not the same as the parameter used to CALL '// &
                                  'this function.')
if (iss /= nss) call WarningMessage(1,'read_magn :: nss read from DATA_FILE is not the same as the parameter used to CALL '// &
                                    'this function.')
zJ = Zero
T(:) = Zero
H(:) = Zero
X(:) = Zero
Y(:) = Zero
Z(:) = Zero
W(:) = Zero
energy(:,:,:) = Zero
m(:,:,:,:) = Zero
mav(:,:) = Zero
read(DATA_FILE,*,iostat=ierr) zJ
read(DATA_FILE,*,iostat=ierr) (T(i),i=1,nt)
if (ierr /= 0) call WarningMessage(2,'read_magn :: Something went wrong reading the T array.')
! Magnetic field points
read(DATA_FILE,*,iostat=ierr) (H(i),i=1,nh)
if (ierr /= 0) call WarningMessage(2,'read_magn :: Something went wrong reading the H array.')
! Lebedev grid data (X, Y, Z, W)
read(DATA_FILE,*,iostat=ierr) (X(i),i=1,nd)
if (ierr /= 0) call WarningMessage(2,'read_magn :: Something went wrong reading the Lebedev grid X rray.')
read(DATA_FILE,*,iostat=ierr) (Y(i),i=1,nd)
if (ierr /= 0) call WarningMessage(2,'read_magn :: Something went wrong reading the Lebedev grid Y array.')
read(DATA_FILE,*,iostat=ierr) (Z(i),i=1,nd)
if (ierr /= 0) call WarningMessage(2,'read_magn :: Something went wrong reading the Lebedev grid Z array.')
read(DATA_FILE,*,iostat=ierr) (W(i),i=1,nd)
if (ierr /= 0) call WarningMessage(2,'read_magn :: Something went wrong reading the Lebedev grid W array.')
! Zeeman energy
do id=1,nd
  do ih=1,nh
    read(DATA_FILE,*,iostat=ierr) (energy(ih,id,iss),iss=1,nss)
    if (ierr /= 0) call WarningMessage(2,'read_magn :: Something went wrong reading the Zeeman energy data.')
  end do
end do
! magnetisation vector data
do id=1,nd
  do l=1,3
    do ih=1,nh
      read(DATA_FILE,*,iostat=ierr) (m(id,l,ih,it),it=1,nt)
      if (ierr /= 0) call WarningMessage(2,'read_magn :: Something went wrong reading the M data.')
    end do
  end do
end do
do ih=1,nh
  read(DATA_FILE,*,iostat=ierr) (mav(ih,it),it=1,nt)
  if (ierr /= 0) call WarningMessage(2,'read_magn :: Something went wrong reading the average M data.')
end do
if (dbg) call xFlush(u6)

return

end subroutine read_magn

subroutine read_complex_matrix(LU,key,n,matrix,dbg)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LU, N
character(len=*), intent(in) :: key
complex(kind=wp), intent(out) :: matrix(N,N)
logical(kind=iwp), intent(in) :: dbg
real(kind=wp), allocatable :: ri(:,:), rr(:,:)

call mma_allocate(rr,n,n,label='rr')
call mma_allocate(ri,n,n,label='ri')
call read_2d_real_array(LU,key//'r',n,n,rr,dbg)
call read_2d_real_array(LU,key//'i',n,n,ri,dbg)
matrix(:,:) = cmplx(rr(:,:),ri(:,:),kind=wp)
call mma_deallocate(rr)
call mma_deallocate(ri)

return

end subroutine read_complex_matrix

subroutine write_complex_matrix(LU,key,n,matrix,dbg)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LU, N
character(len=*), intent(in) :: key
complex(kind=wp), intent(in) :: matrix(N,N)
logical(kind=iwp), intent(in) :: dbg
real(kind=wp), allocatable :: ri(:,:), rr(:,:)

call mma_allocate(rr,n,n,label='rr')
call mma_allocate(ri,n,n,label='ri')
rr(:,:) = real(matrix(:,:))
ri(:,:) = aimag(matrix(:,:))
call write_2d_real_array(LU,key//'r',n,n,rr,dbg)
call write_2d_real_array(LU,key//'i',n,n,ri,dbg)
call mma_deallocate(rr)
call mma_deallocate(ri)

return

end subroutine write_complex_matrix

!----------------------------------------------------------------------!
!           HIGH LEVEL WRITING SUBROUTINES
!----------------------------------------------------------------------!
!
! write_magnetic_moment
! write_electric_moment
! write_amfi
! write_angmom
! write_edipmom
! write_eigenv
! write_hso
! write_eso
! write_esfs
! write_nss
! write_nmult
! write_nroot
! write_iroot
! write_nstate
! write_szproj
! write_multiplicty
! write_gtens
! write_stev_cfp
! write_susc
! write_magn
!
!----------------------------------------------------------------------!

subroutine write_magnetic_moment(ANISO_FILE,n,moment,dbg)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ANISO_FILE
integer(kind=iwp), intent(in) :: N
complex(kind=wp), intent(in) :: moment(3,N,N)
logical(kind=iwp), intent(in) :: dbg
real(kind=wp), allocatable :: ri(:,:), rr(:,:)

call mma_allocate(rr,n,n,label='rr')
call mma_allocate(ri,n,n,label='ri')
! projection X
rr(:,:) = real(moment(1,:,:))
ri(:,:) = aimag(moment(1,:,:))
call write_2d_real_array(ANISO_FILE,'$magn_xr',n,n,rr,dbg)
call write_2d_real_array(ANISO_FILE,'$magn_xi',n,n,ri,dbg)
! projection Y
rr(:,:) = real(moment(2,:,:))
ri(:,:) = aimag(moment(2,:,:))
call write_2d_real_array(ANISO_FILE,'$magn_yr',n,n,rr,dbg)
call write_2d_real_array(ANISO_FILE,'$magn_yi',n,n,ri,dbg)
! projection Z
rr(:,:) = real(moment(3,:,:))
ri(:,:) = aimag(moment(3,:,:))
call write_2d_real_array(ANISO_FILE,'$magn_zr',n,n,rr,dbg)
call write_2d_real_array(ANISO_FILE,'$magn_zi',n,n,ri,dbg)
call mma_deallocate(rr)
call mma_deallocate(ri)

return

end subroutine write_magnetic_moment

subroutine write_electric_moment(ANISO_FILE,n,moment,dbg)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ANISO_FILE
integer(kind=iwp), intent(in) :: N
complex(kind=wp), intent(in) :: moment(3,N,N)
logical(kind=iwp), intent(in) :: dbg
real(kind=wp), allocatable :: ri(:,:), rr(:,:)

call mma_allocate(rr,n,n,label='rr')
call mma_allocate(ri,n,n,label='ri')
! projection X
rr(:,:) = real(moment(1,:,:))
ri(:,:) = aimag(moment(1,:,:))
call write_2d_real_array(ANISO_FILE,'$edipm_xr',n,n,rr,dbg)
call write_2d_real_array(ANISO_FILE,'$edipm_xi',n,n,ri,dbg)
! projection Y
rr(:,:) = real(moment(2,:,:))
ri(:,:) = aimag(moment(2,:,:))
call write_2d_real_array(ANISO_FILE,'$edipm_yr',n,n,rr,dbg)
call write_2d_real_array(ANISO_FILE,'$edipm_yi',n,n,ri,dbg)
! projection Z
rr(:,:) = real(moment(3,:,:))
ri(:,:) = aimag(moment(3,:,:))
call write_2d_real_array(ANISO_FILE,'$edipm_zr',n,n,rr,dbg)
call write_2d_real_array(ANISO_FILE,'$edipm_zi',n,n,ri,dbg)
call mma_deallocate(rr)
call mma_deallocate(ri)

return

end subroutine write_electric_moment

subroutine write_spin_moment(ANISO_FILE,n,moment,dbg)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ANISO_FILE
integer(kind=iwp), intent(in) :: N
complex(kind=wp), intent(in) :: moment(3,N,N)
logical(kind=iwp), intent(in) :: dbg
real(kind=wp), allocatable :: ri(:,:), rr(:,:)

call mma_allocate(rr,n,n,label='rr')
call mma_allocate(ri,n,n,label='ri')
! projection X
rr(:,:) = real(moment(1,:,:))
ri(:,:) = aimag(moment(1,:,:))
call write_2d_real_array(ANISO_FILE,'$spin_xr',n,n,rr,dbg)
call write_2d_real_array(ANISO_FILE,'$spin_xi',n,n,ri,dbg)
! projection Y
rr(:,:) = real(moment(2,:,:))
ri(:,:) = aimag(moment(2,:,:))
call write_2d_real_array(ANISO_FILE,'$spin_yr',n,n,rr,dbg)
call write_2d_real_array(ANISO_FILE,'$spin_yi',n,n,ri,dbg)
! projection Z
rr(:,:) = real(moment(3,:,:))
ri(:,:) = aimag(moment(3,:,:))
call write_2d_real_array(ANISO_FILE,'$spin_zr',n,n,rr,dbg)
call write_2d_real_array(ANISO_FILE,'$spin_zi',n,n,ri,dbg)
call mma_deallocate(rr)
call mma_deallocate(ri)

return

end subroutine write_spin_moment

subroutine write_angmom(ANISO_FILE,n,moment,dbg)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ANISO_FILE
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(in) :: moment(3,N,N)
logical(kind=iwp), intent(in) :: dbg
real(kind=wp), allocatable :: tmp(:,:)

call mma_allocate(tmp,N,N,label='tmp')
tmp(:,:) = moment(1,:,:)
call write_2d_real_array(ANISO_FILE,'$angmom_x',n,n,tmp,dbg)
tmp(:,:) = moment(2,:,:)
call write_2d_real_array(ANISO_FILE,'$angmom_y',n,n,tmp,dbg)
tmp(:,:) = moment(3,:,:)
call write_2d_real_array(ANISO_FILE,'$angmom_z',n,n,tmp,dbg)
call mma_deallocate(tmp)

return

end subroutine write_angmom

subroutine write_edipmom(ANISO_FILE,n,moment,dbg)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ANISO_FILE
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(in) :: moment(3,N,N)
logical(kind=iwp), intent(in) :: dbg
real(kind=wp), allocatable :: tmp(:,:)

call mma_allocate(tmp,N,N,label='tmp')
tmp(:,:) = moment(1,:,:)
call write_2d_real_array(ANISO_FILE,'$edmom_x',n,n,tmp,dbg)
tmp(:,:) = moment(2,:,:)
call write_2d_real_array(ANISO_FILE,'$edmom_y',n,n,tmp,dbg)
tmp(:,:) = moment(3,:,:)
call write_2d_real_array(ANISO_FILE,'$edmom_z',n,n,tmp,dbg)
call mma_deallocate(tmp)

return

end subroutine write_edipmom

subroutine write_amfi(ANISO_FILE,n,moment,dbg)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ANISO_FILE
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(in) :: moment(3,N,N)
logical(kind=iwp), intent(in) :: dbg
real(kind=wp), allocatable :: tmp(:,:)

call mma_allocate(tmp,N,N,label='tmp')
tmp(:,:) = moment(1,:,:)
call write_2d_real_array(ANISO_FILE,'$amfi_x',n,n,tmp,dbg)
tmp(:,:) = moment(2,:,:)
call write_2d_real_array(ANISO_FILE,'$amfi_y',n,n,tmp,dbg)
tmp(:,:) = moment(3,:,:)
call write_2d_real_array(ANISO_FILE,'$amfi_z',n,n,tmp,dbg)
call mma_deallocate(tmp)

return

end subroutine write_amfi

subroutine write_format(ANISO_FILE,n,dbg)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ANISO_FILE
integer(kind=iwp), intent(in) :: n
logical(kind=iwp), intent(in) :: dbg

call write_INTEGER_scalar(ANISO_FILE,'$format',n,dbg)

return

end subroutine write_format

subroutine write_nss(ANISO_FILE,n,dbg)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ANISO_FILE, n
logical(kind=iwp), intent(in) :: dbg

call write_INTEGER_scalar(ANISO_FILE,'$nss',n,dbg)

return

end subroutine write_nss

subroutine write_nstate(ANISO_FILE,n,dbg)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ANISO_FILE, n
logical(kind=iwp), intent(in) :: dbg

call write_INTEGER_scalar(ANISO_FILE,'$nstate',n,dbg)

return

end subroutine write_nstate

subroutine write_nmult(ANISO_FILE,n,dbg)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ANISO_FILE, n
logical(kind=iwp), intent(in) :: dbg

call write_INTEGER_scalar(ANISO_FILE,'$nmult',n,dbg)

return

end subroutine write_nmult

subroutine write_multiplicity(ANISO_FILE,n,array,dbg)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ANISO_FILE, n, array(n)
logical(kind=iwp), intent(in) :: dbg

call write_1d_INTEGER_array(ANISO_FILE,'$multiplicity',n,array,dbg)

return

end subroutine write_multiplicity

subroutine write_imult(ANISO_FILE,n,array,dbg)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ANISO_FILE, n, array(n)
logical(kind=iwp), intent(in) :: dbg

call write_1d_INTEGER_array(ANISO_FILE,'$imult',n,array,dbg)

return

end subroutine write_imult

subroutine write_nroot(ANISO_FILE,n,array,dbg)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ANISO_FILE, n, array(n)
logical(kind=iwp), intent(in) :: dbg

call write_1d_INTEGER_array(ANISO_FILE,'$nroot',n,array,dbg)

return

end subroutine write_nroot

subroutine write_szproj(ANISO_FILE,n,array,dbg)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ANISO_FILE, n, array(n)
logical(kind=iwp), intent(in) :: dbg

call write_1d_INTEGER_array(ANISO_FILE,'$szproj',n,array,dbg)

return

end subroutine write_szproj

subroutine write_eso(ANISO_FILE,n,array,dbg)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ANISO_FILE, n
real(kind=wp), intent(in) :: array(n)
logical(kind=iwp), intent(in) :: dbg

if (dbg) write(u6,*) 'write_eso: '
call write_1d_real_array(ANISO_FILE,'$eso',n,array,dbg)

return

end subroutine write_eso

subroutine write_esfs(ANISO_FILE,n,array,dbg)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ANISO_FILE, n
real(kind=wp), intent(in) :: array(n)
logical(kind=iwp), intent(in) :: dbg

call write_1d_real_array(ANISO_FILE,'$esfs',n,array,dbg)

return

end subroutine write_esfs

subroutine write_hso(ANISO_FILE,n,array,dbg)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ANISO_FILE, n
complex(kind=wp), intent(in) :: array(n)
logical(kind=iwp), intent(in) :: dbg

call write_complex_matrix(ANISO_FILE,'$hso',n,array,dbg)

return

end subroutine write_hso

subroutine write_eigen(ANISO_FILE,n,array,dbg)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ANISO_FILE, n
complex(kind=wp), intent(in) :: array(n)
logical(kind=iwp), intent(in) :: dbg

call write_complex_matrix(ANISO_FILE,'$eigen',n,array,dbg)

return

end subroutine write_eigen

subroutine write_gtens(ANISO_FILE,nmult,gtens,axes,dbg)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ANISO_FILE, nmult
real(kind=wp), intent(in) :: gtens(nmult,3), axes(nmult,3,3)
logical(kind=iwp), intent(in) :: dbg

call write_2d_real_array(ANISO_FILE,'$gtens_main',nmult,3,gtens,dbg)
call write_3d_real_array(ANISO_FILE,'$gtens_axes',nmult,3,3,axes,dbg)

return

end subroutine write_gtens

subroutine write_stev_cfp(ANISO_FILE,s,n,cfp,dbg)

use Constants, only: Ten
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ANISO_FILE, n   ! 2J+1, or 2L+1, i.e. the dimension of the J or L multiplet
character(len=*), intent(in) :: s
real(kind=wp), intent(in) :: cfp(n-1,-(n-1):n-1)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: ierr, k, q
character(len=500) :: line
character(len=30) :: FMTCFP = '(2(I0,1x),ES22.14)'
real(kind=wp), parameter :: MINIMAL_REAL = tiny(MINIMAL_REAL)*Ten

! s takes the value:
! 'L' == CFP for a term
! 'J' == CFP for a spin-orbit multiplet

ierr = 0
if ((n <= 0)) then
  call WarningMessage(1,'write_stev_cfp_'//trim(s)//' :: nothing to write. Array size = 0.')
  return
end if

if ((trim(s) /= 'l') .or. (trim(s) /= 'j')) then
  call WarningMessage(1,'write_stev_cfp_'//trim(s)//' :: the parameter s='//trim(s)// &
                      'is not understood. RETURN without writing cfp')
  return
end if

if (sum(abs(cfp)) <= MINIMAL_REAL) call WarningMessage(1,'write_stev_cfp_'//trim(s)//':: all array elements are zero = 0.')

rewind(ANISO_FILE)
call file_advance_to_string(ANISO_FILE,'$stev_cfp_'//trim(s),line,ierr,dbg)

if (ierr /= 0) then
  write(ANISO_FILE,'(A)',iostat=ierr)
  write(ANISO_FILE,'(A)',iostat=ierr) '$stev_cfp_'//trim(s)
end if

! if key is found, THEN rewrite the data
write(ANISO_FILE,*,iostat=ierr) n
do k=2,n-1,2
  do q=-k,k,2
    if (abs(cfp(k,q)) > MINIMAL_REAL) write(ANISO_FILE,FMTCFP,iostat=ierr) k,q,cfp(k,q)
    if (ierr /= 0) call WarningMessage(2,'write_stev_cfp_'//trim(s)//':: Something went wrong writing the array.')
    if (dbg) write(u6,*) 'write_stev_cfp_'//trim(s)//'::  k, q =',k,q
  end do
end do
write(ANISO_FILE,'(A)',iostat=ierr)
call xFlush(ANISO_FILE)

return

end subroutine write_stev_cfp

subroutine write_susc(ANISO_FILE,s,n,field,zj,t,z,x,x_tens,dbg)
! n     :  number of temperature points
! field :  applied field
! zJ    :  intermolecular interaction
! t     :  temperature points
! z     :  partition function
! x     :  susceptibility X
! x_tens:  susceptibility tensor X_tens

use Constants, only: Ten
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ANISO_FILE, n
character(len=*), intent(in) :: s
real(kind=wp), intent(in) :: field, zJ, t(n), z(n), x(n), x_tens(n,3,3)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr, j, k
character(len=500) :: line
character(len=*), parameter :: FMTI = '(20(I0,1x))', FMTR = '(5ES22.14)'
real(kind=wp), parameter :: MINIMAL_REAL = tiny(MINIMAL_REAL)*Ten
real(kind=wp), external :: dnrm2_

! s takes the values:
! 'x'
! 'x_field'
! 'x_minus_one'
! 'x_field_minus_one'
! 'xt'
! 'xt_field'
! 'xt_minus_one'
! 'xt_field_minus_one'
! 'x_zj'
! 'x_field_zj'
! 'x_minus_one_zj'
! 'x_field_minus_one_zj'
! 'xt_zj'
! 'xt_field_zj'
! 'xt_minus_one_zj'
! 'xt_field_minus_one_zj'
!if ((trim(s) /= 'x')                    .or. &
!    (trim(s) /= 'x_field')              .or. &
!    (trim(s) /= 'x_minus_one')          .or. &
!    (trim(s) /= 'x_field_minus_one')    .or. &
!    (trim(s) /= 'xt')                   .or. &
!    (trim(s) /= 'xt_field')             .or. &
!    (trim(s) /= 'xt_minus_one')         .or. &
!    (trim(s) /= 'xt_field_minus_one')   .or. &
!    (trim(s) /= 'x_zj')                 .or. &
!    (trim(s) /= 'x_field_zj')           .or. &
!    (trim(s) /= 'x_minus_one_zj')       .or. &
!    (trim(s) /= 'x_field_minus_one_zj') .or. &
!    (trim(s) /= 'xt_zj')                .or. &
!    (trim(s) /= 'xt_field_zj')          .or. &
!    (trim(s) /= 'xt_minus_one_zj')      .or. &
!    (trim(s) /= 'xt_field_minus_one_zj')) then
!
!  call WarningMessage(1,'write_susc '//trim(s)//' :: the parameter s='//trim(s)// &
!                      ' is not understood. RETURN without write of X data.')
!  return
!end if

ierr = 0
if ((n <= 0)) then
  call WarningMessage(1,'write_susc '//trim(s)//' :: nothing to write. Array size = 0.')
  return
end if
if (dnrm2_(n,T,1) <= MINIMAL_REAL) call WarningMessage(1,'write_susc '//trim(s)//' :: all array T elements are zero = 0.')
if (dnrm2_(n,X,1) <= MINIMAL_REAL) call WarningMessage(1,'write_susc '//trim(s)//' :: all array X elements are zero = 0.')
if (dnrm2_(n,Z,1) <= MINIMAL_REAL) call WarningMessage(1,'write_susc '//trim(s)//' :: all array Z elements are zero = 0.')
if (dnrm2_(3*3*n,X_tens,1) <= MINIMAL_REAL) call WarningMessage(1,'write_susc '//trim(s)// &
                                                                ' :: all array X_tens elements are zero = 0.')

rewind(ANISO_FILE)
call file_advance_to_string(ANISO_FILE,'$susceptibility_'//trim(s),line,ierr,dbg)

if (ierr /= 0) then
  write(ANISO_FILE,'(A)',iostat=ierr)
  write(ANISO_FILE,'(A)',iostat=ierr) '$susceptibility_'//trim(s)
end if

write(ANISO_FILE,FMTI,iostat=ierr) n
write(ANISO_FILE,FMTR,iostat=ierr) zj,field
if (ierr /= 0) call WarningMessage(2,'write_susc '//trim(s)//' :: Something went wrong writing the zJ and field values.')
write(ANISO_FILE,FMTR,iostat=ierr) (T(i),i=1,n)
if (ierr /= 0) call WarningMessage(2,'write_susc '//trim(s)//' :: Something went wrong writing the T array.')
write(ANISO_FILE,FMTR,iostat=ierr) (Z(i),i=1,n)
if (ierr /= 0) call WarningMessage(2,'write_susc '//trim(s)//' :: Something went wrong writing the Z array.')
write(ANISO_FILE,FMTR,iostat=ierr) (X(i),i=1,n)
if (ierr /= 0) call WarningMessage(2,'write_susc '//trim(s)//' :: Something went wrong writing the X array.')
call xFlush(ANISO_FILE)
do j=1,3
  do k=1,3
    write(ANISO_FILE,FMTR,iostat=ierr) (X_tens(i,j,k),i=1,n)
    if (ierr /= 0) call WarningMessage(2,'write_susc '//trim(s)//' :: Something went wrong writing the X_tens array.')
  end do
end do
write(ANISO_FILE,*,iostat=ierr)
call xFlush(ANISO_FILE)
if (dbg) call xFlush(u6)

return

end subroutine write_susc

subroutine write_magn(ANISO_FILE,nt,nh,nd,nss,zj,t,h,x,y,z,w,m,mav,energy,dbg)
! nt        :  number of temperature points
! nh        :  number of field points
! nd        :  number of directions of aplied field
! nss       :  number of spin-orbit states included in the Zeeman interaction
! zj        :  inter-molecular parameter zJ
! t         :  temperature points
! h         :  field points
! x, y, z, w:  Lebedev grid directions and weight
! m         :  magnetisation vector
! mav       :  average magnetisation vector
! energy    :  Zeeman energy states

use Constants, only: Ten
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ANISO_FILE, nt, nh, nd, nss
real(kind=wp), intent(in) :: zj, t(nt), h(nh), x(nd), y(nd), z(nd), w(nd), m(nd,3,nh,nt), mav(nh,nt), energy(nd,nh,nss)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, id, ierr, ih, iss, it, l
character(len=500) :: line
real(kind=wp), parameter :: MINIMAL_REAL = tiny(MINIMAL_REAL)*Ten
character(len=*), parameter :: FMTI = '(20(I0,1x))', FMTR = '(5ES22.14)'
real(kind=wp), external :: dnrm2_

ierr = 0
if ((nt <= 0) .or. (nh <= 0) .or. (nd <= 0)) then
  call WarningMessage(1,'write_magn :: nothing to write. Array size = 0.')
  return
end if
if (dnrm2_(nt,t,1) <= MINIMAL_REAL) call WarningMessage(1,'write_magn :: all array T elements are zero = 0.')
if (dnrm2_(nh,h,1) <= MINIMAL_REAL) call WarningMessage(1,'write_magn :: all array H elements are zero = 0.')
if (dnrm2_(nd*3*nt*nh,m,1) <= MINIMAL_REAL) call WarningMessage(1,'write_magn :: all array M elements are zero = 0.')
if (dnrm2_(nt*nh,mav,1) <= MINIMAL_REAL) call WarningMessage(1,'write_magn :: all array MAV elements are zero = 0.')
if (dnrm2_(nd*nh*nss,energy,1) <= MINIMAL_REAL) call WarningMessage(1,'write_magn :: all array energy elements are zero = 0.')
if ((dnrm2_(nd,x,1) <= MINIMAL_REAL) .or. (dnrm2_(nd,y,1) <= MINIMAL_REAL) .or. &
    (dnrm2_(nd,z,1) <= MINIMAL_REAL) .or. (dnrm2_(nd,w,1) <= MINIMAL_REAL)) &
  call WarningMessage(1,'write_magn :: all array of X,Y,X,W elements are zero = 0.')

rewind(ANISO_FILE)
call file_advance_to_string(ANISO_FILE,'$magnetisation',line,ierr,dbg)

! if key is found, THEN rewrite the data
if (ierr /= 0) then
  write(ANISO_FILE,'(A)',iostat=ierr)
  write(ANISO_FILE,'(A)',iostat=ierr) '$magnetisation'
end if

write(ANISO_FILE,FMTI,iostat=ierr) nt,nh,nd,nss
write(ANISO_FILE,FMTR,iostat=ierr) zJ
write(ANISO_FILE,FMTR,iostat=ierr) (T(i),i=1,nt)
if (ierr /= 0) call WarningMessage(2,'write_magn :: Something went wrong writing the T array.')
! Magnetic field points
write(ANISO_FILE,FMTR,iostat=ierr) (H(i),i=1,nh)
if (ierr /= 0) call WarningMessage(2,'write_magn :: Something went wrong writing the H array.')
! Lebedev grid data (X, Y, Z, W)
write(ANISO_FILE,FMTR,iostat=ierr) (X(i),i=1,nd)
if (ierr /= 0) call WarningMessage(2,'write_magn :: Something went wrong writing the Lebedev grid X rray.')
write(ANISO_FILE,FMTR,iostat=ierr) (Y(i),i=1,nd)
if (ierr /= 0) call WarningMessage(2,'write_magn :: Something went wrong writing the Lebedev grid Y array.')
write(ANISO_FILE,FMTR,iostat=ierr) (Z(i),i=1,nd)
if (ierr /= 0) call WarningMessage(2,'write_magn :: Something went wrong writing the Lebedev grid Z array.')
write(ANISO_FILE,FMTR,iostat=ierr) (W(i),i=1,nd)
if (ierr /= 0) call WarningMessage(2,'write_magn :: Something went wrong writing the Lebedev grid W array.')
call xFlush(ANISO_FILE)
! Zeeman energy
do id=1,nd
  do ih=1,nh
    write(ANISO_FILE,FMTR,iostat=ierr) (energy(id,ih,iss),iss=1,nss)
    if (ierr /= 0) call WarningMessage(2,'write_magn :: Something went wrong writing the Zeeman energy data.')
  end do
end do
call xFlush(ANISO_FILE)
! magnetisation vector data
do id=1,nd
  do l=1,3
    do ih=1,nh
      write(ANISO_FILE,FMTR,iostat=ierr) (m(id,l,ih,it),it=1,nt)
      if (ierr /= 0) call WarningMessage(2,'write_magn :: Something went wrong writing the M data.')
    end do
  end do
end do
call xFlush(ANISO_FILE)
do ih=1,nh
  write(ANISO_FILE,FMTR,iostat=ierr) (mav(ih,it),it=1,nt)
  if (ierr /= 0) call WarningMessage(2,'write_magn :: Something went wrong writing the average M data.')
end do
write(ANISO_FILE,*,iostat=ierr)
call xFlush(ANISO_FILE)
if (dbg) call xFlush(u6)

return

end subroutine write_magn

!--------------------------------------------------------------------------------------------------!
!  BASIC  OPEN FILE / CLOSE FILE / SEARCH KEY INSIDE FILE etc. functions
!
!   open_datafile_write
!   open_datafile_read
!   close_datafile
!   key_found
!   file_advance_to_string
!
!--------------------------------------------------------------------------------------------------!

subroutine open_datafile_write(DATA_FILE,DATA_FILE_NAME)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: DATA_FILE
character(len=180), intent(in) :: DATA_FILE_NAME
!#ifdef _ALONE_
!integer :: ierr
!#endif

!#ifdef _ALONE_
!ierr = 0
!open(unit=DATA_FILE,file=DATA_FILE_NAME,form='formatted',status='old',access='sequential',action='write',iostat=ierr )
!if (ierr /= 0) call WarningMessage(2,'open_datafile_write:: Something went wrong opening DATA_FILE')
!#else
call molcas_open(DATA_FILE,DATA_FILE_NAME)
!#endif

return

end subroutine open_datafile_write

subroutine open_aniso_file(ANISO_FILE,ANISO_FILE_NAME)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ANISO_FILE
character(len=180), intent(in) :: ANISO_FILE_NAME
!#ifdef _ALONE_
!integer :: ierr
!#endif

!#ifdef _ALONE_
!ierr = 0
!open(unit=ANISO_FILE,form='formatted',status='old',access='sequential',action='write',iostat=ierr )
!if (ierr /= 0) call WarningMessage(2,'open_datafile_write:: Something went wrong opening ANISO_FILE')
!#else
call molcas_open(ANISO_FILE,ANISO_FILE_NAME)
!#endif

return

end subroutine open_aniso_file

subroutine open_datafile_read(DATA_FILE,DATA_FILE_NAME)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: DATA_FILE
character(len=180), intent(in) :: DATA_FILE_NAME
!#ifdef _ALONE_
!integer :: ierr
!#endif

!#ifdef _ALONE_
!ierr = 0
!open(unit=DATA_FILE,form='FORMATTED',status='old',access='sequential',action='read',iostat=ierr )
!if (ierr /= 0) call WarningMessage(2,'open_datafile_read:: Something went wrong opening DATA_FILE')
!#else
call molcas_open(DATA_FILE,DATA_FILE_NAME)
!#endif

return

end subroutine open_datafile_read

subroutine close_datafile(DATA_FILE)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: DATA_FILE
integer(kind=iwp) :: ierr

ierr = 0
close(unit=DATA_FILE,iostat=ierr)
if (ierr /= 0) call WarningMessage(2,'close_datafile:: Something went wrong closing DATA_FILE')

return

end subroutine close_datafile

subroutine close_anisofile(ANISO_FILE)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ANISO_FILE
integer(kind=iwp) :: ierr

ierr = 0
close(unit=ANISO_FILE,iostat=ierr)
if (ierr /= 0) call WarningMessage(2,'close_datafile:: Something went wrong closing ANISO_FILE')

return

end subroutine close_anisofile

function key_found(DATA_FILE,key,dbg)

use Definitions, only: iwp

implicit none
logical(kind=iwp) :: key_found
integer(kind=iwp), intent(in) :: DATA_FILE
character(len=*), intent(in) :: key
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: ierr
character(len=500) :: line

ierr = 0
key_found = .false.
rewind(DATA_FILE)
call file_advance_to_string(DATA_FILE,key,line,ierr,dbg)
if (index(line,trim(key)) /= 0) key_found = .true.

return

end function key_found

subroutine file_advance_to_string(LU,key,line,ierr,dbg)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU
character(len=*), intent(in) :: key
character(len=*), intent(out) :: line
integer(kind=iwp), intent(out) :: ierr
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: ios, num_read

ierr = 0
num_read = 0

rewind(LU)
do
  read(LU,'(a)',iostat=ios) line

  if (ios /= 0) then
    !write(u6,'(A,i0,2A)') 'file_advance_to_string:: error reading line: LU=',LU,' key=',key
    exit
  end if

  num_read = num_read+1

  if (index(line,trim(key)) /= 0) return
end do

line = ' '
ierr = 1

if (dbg) then
  write(u6,'(a)') ' '
  write(u6,'(a)') 'FILE_ADVANCE_TO_STRING - Warning!'
  write(u6,'(a)') '  Did not find the key:'
  write(u6,'(a)') '    '//trim(key)
  write(u6,'(a,i6)') '  Number of lines read was ',num_read
end if

return

end subroutine file_advance_to_string

function inquire_key_presence(LU,key)

use Definitions, only: iwp

implicit none
logical(kind=iwp) :: inquire_key_presence
integer(kind=iwp), intent(in) :: LU
character(len=*), intent(in) :: key
integer(kind=iwp) :: ios
character(len=500) :: line

inquire_key_presence = .false.
rewind(LU)
do
  read(LU,'(a)',iostat=ios) line
  if (ios /= 0) call WarningMessage(1,'inquire_key_presence:: error reading line')
  if (index(line,trim(key)) /= 0) then
    inquire_key_presence = .true.
    return
  end if
end do
line = ' '

return

end function inquire_key_presence

!----------------------------------------------------------------------!
!   READING  SUBROUTINES
!----------------------------------------------------------------------!

subroutine read_INTEGER_scalar(LU,key,i,dbg)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU
character(len=*), intent(in) :: key
integer(kind=iwp), intent(out) :: i
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: ierr
character(len=500) :: line

ierr = 0
i = 0
rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)
read(LU,*,iostat=ierr) i
if (ierr /= 0) call WarningMessage(2,'read_INTEGER_scalar:: Something went wrong reading key'//trim(key))
if (dbg) then
  write(u6,*) 'read_INTEGER_scalar:: key =',trim(key)
  write(u6,*) 'read_INTEGER_scalar::   i =',i
end if

return

end subroutine read_INTEGER_scalar

subroutine read_real_scalar(LU,key,r,dbg)

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU
character(len=*), intent(in) :: key
real(kind=wp), intent(out) :: r
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: ierr
character(len=500) :: line

ierr = 0
r = Zero
rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)
read(LU,*,iostat=ierr) r
if (ierr /= 0) call WarningMessage(2,'read_real_scalar:: Something went wrong reading key'//trim(key))
if (dbg) then
  write(u6,*) 'read_real_scalar:: key =',trim(key)
  write(u6,*) 'read_real_scalar::   r =',r
end if

return

end subroutine read_real_scalar

subroutine read_complex_scalar(LU,key,c,dbg)

use Constants, only: Zero, cZero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU
character(len=*), intent(in) :: key
complex(kind=wp), intent(out) :: c
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: ierr
real(kind=wp) :: ri, rr
character(len=500) :: line

ierr = 0
rr = Zero
ri = Zero
c = cZero
rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)
read(LU,*,iostat=ierr) rr,ri
if (ierr /= 0) call WarningMessage(2,'read_complex_scalar:: Something went wrong reading key'//trim(key))
if (dbg) then
  write(u6,*) 'read_complex_scalar::   key =',trim(key)
  write(u6,*) 'read_complex_scalar:: (r,i) =',rr,ri
  write(u6,*) 'read_complex_scalar::     c =',c
end if
c = cmplx(rr,ri,kind=wp)

return

end subroutine read_complex_scalar

subroutine read_string(LU,key,length,s,dbg)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU, length
character(len=*), intent(in) :: key
character(len=length), intent(out) :: s
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr
character(len=500) :: c, f, line

if (dbg) then
  write(u6,*) 'read_string::    key =',trim(key)
  write(u6,*) 'read_string:: length =',length
end if
write(f,'(A,i0,2A)') '"(A',length,')"'
write(u6,'(2A)') 'format =',trim(f)
call xFlush(u6)
rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)
read(LU,f,iostat=ierr) c
write(u6,'(2A)') 'c =',trim(c)
write(s,'(A)') trim(c)

do i=1,len(trim(LINE))
  read(LINE,'(A500)') c
  write(s,'(A)') trim(c)
  if (dbg) write(u6,*) 'read_string::   c =',trim(c)
  call xFlush(u6)
end do

return

end subroutine read_string

subroutine read_1d_size(LU,key,n,dbg)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU
character(len=*), intent(in) :: key
integer(kind=iwp), intent(out) :: n
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: ierr
character(len=500) :: line

ierr = 0
n = 0
rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)
read(LU,*,iostat=ierr) n
if (ierr /= 0) call WarningMessage(2,'read_1d_size:: Something went wrong reading key'//trim(key))
if (dbg) then
  write(u6,*) 'read_1d_size:: key =',trim(key)
  write(u6,*) 'read_1d_size::   n =',n
end if

return

end subroutine read_1d_size

subroutine read_2d_size(LU,key,n1,n2,dbg)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU
character(len=*), intent(in) :: key
integer(kind=iwp), intent(out) :: n1, n2
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: ierr
character(len=500) :: line

ierr = 0
n1 = 0
n2 = 0
rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)
read(LU,*,iostat=ierr) n1,n2
if (ierr /= 0) call WarningMessage(2,'read_2d_size:: Something went wrong reading key'//trim(key))
if (dbg) then
  write(u6,*) 'read_2d_size:: key =',trim(key)
  write(u6,*) 'read_2d_size::  n1 =',n1
  write(u6,*) 'read_2d_size::  n2 =',n2
end if

return

end subroutine read_2d_size

subroutine read_3d_size(LU,key,n1,n2,n3,dbg)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU
character(len=*), intent(in) :: key
integer(kind=iwp), intent(out) :: n1, n2, n3
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: ierr
character(len=500) :: line

ierr = 0
n1 = 0
n2 = 0
n3 = 0
rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)
read(LU,*,iostat=ierr) n1,n2,n3
if (ierr /= 0) call WarningMessage(2,'read_3d_size:: Something went wrong reading key'//trim(key))
if (dbg) then
  write(u6,*) 'read_3d_size:: key =',trim(key)
  write(u6,*) 'read_3d_size::  n1 =',n1
  write(u6,*) 'read_3d_size::  n2 =',n2
  write(u6,*) 'read_3d_size::  n3 =',n3
end if

return

end subroutine read_3d_size

subroutine read_4d_size(LU,key,n1,n2,n3,n4,dbg)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU
character(len=*), intent(in) :: key
integer(kind=iwp), intent(out) :: n1, n2, n3, n4
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: ierr
character(len=500) :: line

ierr = 0
n1 = 0
n2 = 0
n3 = 0
n4 = 0
rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)
read(LU,*,iostat=ierr) n1,n2,n3,n4
if (ierr /= 0) call WarningMessage(2,'read_4d_size:: Something went wrong reading key'//trim(key))
if (dbg) then
  write(u6,*) 'read_4d_size:: key =',trim(key)
  write(u6,*) 'read_4d_size::  n1 =',n1
  write(u6,*) 'read_4d_size::  n2 =',n2
  write(u6,*) 'read_4d_size::  n3 =',n3
  write(u6,*) 'read_4d_size::  n4 =',n4
end if

return

end subroutine read_4d_size

subroutine read_1d_INTEGER_array(LU,key,n,array,dbg)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU, n
character(len=*), intent(in) :: key
integer(kind=iwp), intent(out) :: array(n)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr
character(len=500) :: line

ierr = 0
array = 0
if (n <= 0) then
  call WarningMessage(1,'read_1d_INTEGER_array:: nothing to read. Array size = 0.')
  return
end if
rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)

read(LU,*,iostat=ierr) i
if (ierr /= 0) call WarningMessage(2,'read_1d_INTEGER_array:: Something went wrong reading key'//trim(key))
if (dbg) then
  write(u6,*) 'read_1d_INTEGER_array:: key =',trim(key)
  write(u6,*) 'read_1d_INTEGER_array::   n =',i
end if
if (i /= n) call WarningMessage(2,'read_1d_INTEGER_array:: sizes of the array are different from the ones used to CALL this '// &
                                'SUBROUTINE')

read(LU,*,iostat=ierr) (array(i),i=1,n)
if (ierr /= 0) call WarningMessage(2,'read_1d_INTEGER_array:: Something went wrong reading the array.')
if (dbg) write(u6,*) 'read_1d_INTEGER_array:: array =',(array(i),i=1,n)

return

end subroutine read_1d_INTEGER_array

subroutine read_2d_INTEGER_array(LU,key,n1,n2,array,dbg)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU, n1, n2
character(len=*), intent(in) :: key
integer(kind=iwp), intent(out) :: array(n1,n2)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr, j
character(len=500) :: line

ierr = 0
array = 0
if ((n1 <= 0) .or. (n2 <= 0)) then
  call WarningMessage(1,'read_2d_INTEGER_array:: nothing to read. Array size = 0.')
  if (dbg) then
    write(u6,*) 'read_2d_INTEGER_array::   n1 =',n1
    write(u6,*) 'read_2d_INTEGER_array::   n2 =',n2
  end if
  return
end if
rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)

read(LU,*,iostat=ierr) i,j
if (ierr /= 0) call WarningMessage(2,'read_2d_INTEGER_array:: Something went wrong reading key'//trim(key))
if (dbg) then
  write(u6,*) 'read_2d_INTEGER_array:: key =',trim(key)
  write(u6,*) 'read_2d_INTEGER_array::  n1 =',i
  write(u6,*) 'read_2d_INTEGER_array::  n2 =',j
end if
if ((i /= n1) .or. (j /= n2)) call WarningMessage(2,'read_2d_INTEGER_array:: sizes of the array are different from the ones '// &
                                                  'used to CALL this SUBROUTINE')
do i=1,n1
  read(LU,*,iostat=ierr) (array(i,j),j=1,n2)
  if (ierr /= 0) call WarningMessage(2,'read_2d_INTEGER_array:: Something went wrong reading the array.')
  if (dbg) write(u6,*) 'read_2d_INTEGER_array::  i =',i
end do

return

end subroutine read_2d_INTEGER_array

subroutine read_3d_INTEGER_array(LU,key,n1,n2,n3,array,dbg)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU, n1, n2, n3
character(len=*), intent(in) :: key
integer(kind=iwp), intent(out) :: array(n1,n2,n3)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr, j, k
character(len=500) :: line

ierr = 0
array = 0
if ((n1 <= 0) .or. (n2 <= 0) .or. (n3 <= 0)) then
  call WarningMessage(1,'read_3d_INTEGER_array:: nothing to read. Array size = 0.')
  if (dbg) then
    write(u6,*) 'read_3d_INTEGER_array::   n1 =',n1
    write(u6,*) 'read_3d_INTEGER_array::   n2 =',n2
    write(u6,*) 'read_3d_INTEGER_array::   n3 =',n3
  end if
  return
end if
rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)

read(LU,*,iostat=ierr) i,j,k
if (ierr /= 0) call WarningMessage(2,'read_3d_INTEGER_array:: Something went wrong reading key'//trim(key))
if (dbg) then
  write(u6,*) 'read_3d_INTEGER_array:: key =',trim(key)
  write(u6,*) 'read_3d_INTEGER_array::  n1 =',i
  write(u6,*) 'read_3d_INTEGER_array::  n2 =',j
  write(u6,*) 'read_3d_INTEGER_array::  n3 =',k
end if
if ((i /= n1) .or. (j /= n2) .or. (k /= n3)) call WarningMessage(2,'read_3d_INTEGER_array:: sizes of the array are different '// &
                                                                 'from the ones used to CALL this SUBROUTINE')
do i=1,n1
  do j=1,n2
    read(LU,*,iostat=ierr) (array(i,j,k),k=1,n3)
    if (ierr /= 0) call WarningMessage(2,'read_3d_INTEGER_array:: Something went wrong reading the array.')
    if (dbg) write(u6,*) 'read_3d_INTEGER_array::  i,j =',i,j
  end do
end do

return

end subroutine read_3d_INTEGER_array

subroutine read_4d_INTEGER_array(LU,key,n1,n2,n3,n4,array,dbg)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU, n1, n2, n3, n4
character(len=*), intent(in) :: key
integer(kind=iwp), intent(out) :: array(n1,n2,n3,n4)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr, j, k, l
character(len=500) :: line

ierr = 0
array = 0
if ((n1 <= 0) .or. (n2 <= 0) .or. (n3 <= 0) .or. (n4 <= 0)) then
  call WarningMessage(1,'read_4d_INTEGER_array:: nothing to read. Array size = 0.')
  if (dbg) then
    write(u6,*) 'read_4d_INTEGER_array::   n1 =',n1
    write(u6,*) 'read_4d_INTEGER_array::   n2 =',n2
    write(u6,*) 'read_4d_INTEGER_array::   n3 =',n3
    write(u6,*) 'read_4d_INTEGER_array::   n4 =',n4
  end if
  return
end if
rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)

read(LU,*,iostat=ierr) i,j,k,l
if (ierr /= 0) call WarningMessage(2,'read_4d_INTEGER_array:: Something went wrong reading key'//trim(key))
if (dbg) then
  write(u6,*) 'read_4d_INTEGER_array:: key =',trim(key)
  write(u6,*) 'read_4d_INTEGER_array::  n1 =',i
  write(u6,*) 'read_4d_INTEGER_array::  n2 =',j
  write(u6,*) 'read_4d_INTEGER_array::  n3 =',k
  write(u6,*) 'read_4d_INTEGER_array::  n4 =',l
end if
if ((i /= n1) .or. (j /= n2) .or. (k /= n3) .or. (l /= n4)) &
  call WarningMessage(2,'read_4d_INTEGER_array:: sizes of the array are different from the ones used to CALL this SUBROUTINE')
do i=1,n1
  do j=1,n2
    do k=1,n3
      read(LU,*,iostat=ierr) (array(i,j,k,l),l=1,n4)
      if (ierr /= 0) call WarningMessage(2,'read_4d_INTEGER_array:: Something went wrong reading the array.')
      if (dbg) write(u6,*) 'read_4d_INTEGER_array::  i,j,k =',i,j,k
    end do
  end do
end do

return

end subroutine read_4d_INTEGER_array

subroutine read_1d_real_array(LU,key,n,array,dbg)

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU, n
character(len=*), intent(in) :: key
real(kind=wp), intent(out) :: array(n)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr
character(len=500) :: line

ierr = 0
array(:) = Zero
if (n <= 0) then
  call WarningMessage(1,'read_1d_real_array:: nothing to read. Array size = 0.')
  return
end if
rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)

read(LU,*,iostat=ierr) i
if (ierr /= 0) call WarningMessage(2,'read_1d_real_array:: Something went wrong reading key'//trim(key))
if (dbg) then
  write(u6,*) 'read_1d_real_array:: key =',trim(key)
  write(u6,*) 'read_1d_real_array::   n =',i
end if
if ((i /= n)) call WarningMessage(2,'read_1d_real_array:: sizes of the array are different from the ones used to CALL this '// &
                                  'SUBROUTINE')

read(LU,*,iostat=ierr) (array(i),i=1,n)
if (ierr /= 0) call WarningMessage(2,'read_1d_real_array:: Something went wrong reading the array.')
if (dbg) write(u6,*) 'read_1d_real_array:: array =',(array(i),i=1,n)

return

end subroutine read_1d_real_array

subroutine read_2d_real_array(LU,key,n1,n2,array,dbg)

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU, n1, n2
character(len=*), intent(in) :: key
real(kind=wp), intent(out) :: array(n1,n2)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr, j
character(len=500) :: line

ierr = 0
array(:,:) = Zero
if ((n1 <= 0) .or. (n2 <= 0)) then
  call WarningMessage(1,'read_2d_real_array:: nothing to read. Array size = 0.')
  if (dbg) then
    write(u6,*) 'read_2d_real_array::   n1 =',n1
    write(u6,*) 'read_2d_real_array::   n2 =',n2
  end if
  return
end if

rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)

if (ierr /= 0) call WarningMessage(2,'read_2d_real_array:: Something went wrong reading key'//trim(key))
if (dbg) write(u6,*) 'read_2d_real_array:: key =',trim(key)

read(LU,*,iostat=ierr) i,j

if (dbg) then
  write(u6,*) 'read_2d_real_array::  n1 =',i
  write(u6,*) 'read_2d_real_array::  n2 =',j
end if
if ((i /= n1) .or. (j /= n2)) call WarningMessage(2,'read_2d_real_array:: sizes of the array are different from the ones used '// &
                                                  'to CALL this SUBROUTINE')

do i=1,n1
  read(LU,*,iostat=ierr) (array(i,j),j=1,n2)
  if (dbg) write(u6,*) (array(i,j),j=1,n2)
  if (ierr /= 0) call WarningMessage(2,'read_2d_real_array:: Something went wrong reading the array.')
  if (dbg) write(u6,*) 'read_2d_real_array::  i =',i
end do

return

end subroutine read_2d_real_array

subroutine read_3d_real_array(LU,key,n1,n2,n3,array,dbg)

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU, n1, n2, n3
character(len=*), intent(in) :: key
real(kind=wp), intent(out) :: array(n1,n2,n3)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr, j, k
character(len=500) :: line

ierr = 0
array(:,:,:) = Zero
if ((n1 <= 0) .or. (n2 <= 0) .or. (n3 <= 0)) then
  call WarningMessage(1,'read_3d_real_array:: nothing to read. Array size = 0.')
  if (dbg) then
    write(u6,*) 'read_3d_real_array::   n1 =',n1
    write(u6,*) 'read_3d_real_array::   n2 =',n2
    write(u6,*) 'read_3d_real_array::   n3 =',n3
  end if
  return
end if
rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)

read(LU,*,iostat=ierr) i,j,k
if (ierr /= 0) call WarningMessage(2,'read_3d_real_array:: Something went wrong reading key'//trim(key))
if (dbg) then
  write(u6,*) 'read_3d_real_array:: key =',trim(key)
  write(u6,*) 'read_3d_real_array::  n1 =',i
  write(u6,*) 'read_3d_real_array::  n2 =',j
  write(u6,*) 'read_3d_real_array::  n3 =',k
end if
if ((i /= n1) .or. (j /= n2) .or. (k /= n3)) call WarningMessage(2,'read_3d_real_array:: sizes of the array are different from '// &
                                                                 'the ones used to CALL this SUBROUTINE')
do i=1,n1
  do j=1,n2
    read(LU,*,iostat=ierr) (array(i,j,k),k=1,n3)
    if (ierr /= 0) call WarningMessage(2,'read_3d_real_array:: Something went wrong reading the array.')
    if (dbg) write(u6,*) 'read_3d_real_array::  i,j =',i,j
  end do
end do

return

end subroutine read_3d_real_array

subroutine read_4d_real_array(LU,key,n1,n2,n3,n4,array,dbg)

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU, n1, n2, n3, n4
character(len=*), intent(in) :: key
real(kind=wp), intent(out) :: array(n1,n2,n3,n4)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr, j, k, l
character(len=500) :: line

ierr = 0
array(:,:,:,:) = Zero
if ((n1 <= 0) .or. (n2 <= 0) .or. (n3 <= 0) .or. (n4 <= 0)) then
  call WarningMessage(1,'read_4d_real_array:: nothing to read. Array size = 0.')
  if (dbg) then
    write(u6,*) 'read_4d_real_array::   n1 =',n1
    write(u6,*) 'read_4d_real_array::   n2 =',n2
    write(u6,*) 'read_4d_real_array::   n3 =',n3
    write(u6,*) 'read_4d_real_array::   n4 =',n4
  end if
  return
end if
rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)

read(LU,*,iostat=ierr) i,j,k,l
if (ierr /= 0) call WarningMessage(2,'read_4d_real_array:: Something went wrong reading key'//trim(key))
if (dbg) then
  write(u6,*) 'read_4d_real_array:: key =',trim(key)
  write(u6,*) 'read_4d_real_array::  n1 =',i
  write(u6,*) 'read_4d_real_array::  n2 =',j
  write(u6,*) 'read_4d_real_array::  n3 =',k
  write(u6,*) 'read_4d_real_array::  n4 =',l
end if
if ((i /= n1) .or. (j /= n2) .or. (k /= n3) .or. (l /= n4)) &
  call WarningMessage(2,'read_4d_real_array:: sizes of the array are different from the ones used to CALL this SUBROUTINE')
do i=1,n1
  do j=1,n2
    do k=1,n3
      read(LU,*,iostat=ierr) (array(i,j,k,l),l=1,n4)
      if (ierr /= 0) call WarningMessage(2,'read_4d_real_array:: Something went wrong reading the array.')
      if (dbg) write(u6,*) 'read_4d_real_array::  i,j,k =',i,j,k
    end do
  end do
end do

return

end subroutine read_4d_real_array

subroutine read_1d_complex_array(LU,key,n,array,dbg)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU, n
character(len=*), intent(in) :: key
complex(kind=wp), intent(out) :: array(n)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr
character(len=500) :: line
real(kind=wp), allocatable :: ri(:), rr(:)

ierr = 0
if (n <= 0) then
  call WarningMessage(1,'read_1d_complex_array:: nothing to read. Array size = 0.')
  return
end if
rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)

read(LU,*,iostat=ierr) i
if (ierr /= 0) call WarningMessage(2,'read_1d_complex_array:: Something went wrong reading key'//trim(key))
if (dbg) then
  write(u6,*) 'read_1d_complex_array:: key =',trim(key)
  write(u6,*) 'read_1d_complex_array::   n =',i
end if
if ((i /= n)) call WarningMessage(2,'read_1d_complex_array:: sizes of the array are different from the ones used to CALL this '// &
                                  'SUBROUTINE')
call mma_allocate(rr,n,label='rr')
call mma_allocate(ri,n,label='ri')
rr(:) = Zero
ri(:) = Zero
read(LU,*,iostat=ierr) (rr(i),ri(i),i=1,n)
if (ierr /= 0) call WarningMessage(2,'read_1d_complex_array:: Something went wrong reading the array.')
if (dbg) write(u6,*) 'read_1d_complex_array:: array =',(rr(i),ri(i),i=1,n)
array(:) = cmplx(rr(:),ri(:),kind=wp)
call mma_deallocate(rr)
call mma_deallocate(ri)

return

end subroutine read_1d_complex_array

subroutine read_2d_complex_array(LU,key,n1,n2,array,dbg)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU, n1, n2
character(len=*), intent(in) :: key
complex(kind=wp), intent(out) :: array(n1,n2)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr, j
character(len=500) :: line
real(kind=wp), allocatable :: ri(:,:), rr(:,:)

ierr = 0
if ((n1 <= 0) .or. (n2 <= 0)) then
  call WarningMessage(1,'read_2d_complex_array:: nothing to read. Array size = 0.')
  if (dbg) then
    write(u6,*) 'read_2d_complex_array::   n1 =',n1
    write(u6,*) 'read_2d_complex_array::   n2 =',n2
  end if
  return
end if
rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)

read(LU,*,iostat=ierr) i,j
if (ierr /= 0) call WarningMessage(2,'read_2d_complex_array:: Something went wrong reading key'//trim(key))
if (dbg) then
  write(u6,*) 'read_2d_complex_array:: key =',trim(key)
  write(u6,*) 'read_2d_complex_array::  n1 =',i
  write(u6,*) 'read_2d_complex_array::  n2 =',j
end if
if ((i /= n1) .or. (j /= n2)) call WarningMessage(2,'read_2d_complex_array:: sizes of the array are different from the ones '// &
                                                  'used to CALL this SUBROUTINE')
call mma_allocate(rr,n1,n2,label='rr')
call mma_allocate(ri,n1,n2,label='ri')
rr(:,:) = Zero
ri(:,:) = Zero
do i=1,n1
  read(LU,*,iostat=ierr) (rr(i,j),ri(i,j),j=1,n2)
  if (ierr /= 0) call WarningMessage(2,'read_2d_complex_array:: Something went wrong reading the array.')
  if (dbg) write(u6,*) 'read_2d_complex_array::  i =',i
end do
array(:,:) = cmplx(rr(:,:),ri(:,:),kind=wp)
call mma_deallocate(rr)
call mma_deallocate(ri)

return

end subroutine read_2d_complex_array

subroutine read_3d_complex_array(LU,key,n1,n2,n3,array,dbg)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU, n1, n2, n3
character(len=*), intent(in) :: key
complex(kind=wp), intent(out) :: array(n1,n2,n3)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr, j, k
character(len=500) :: line
real(kind=wp), allocatable :: ri(:,:,:), rr(:,:,:)

ierr = 0
if ((n1 <= 0) .or. (n2 <= 0) .or. (n3 <= 0)) then
  call WarningMessage(1,'read_3d_complex_array:: nothing to read. Array size = 0.')
  if (dbg) then
    write(u6,*) 'read_3d_complex_array::   n1 =',n1
    write(u6,*) 'read_3d_complex_array::   n2 =',n2
    write(u6,*) 'read_3d_complex_array::   n3 =',n3
  end if
  return
end if
rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)

read(LU,*,iostat=ierr) i,j,k
if (ierr /= 0) call WarningMessage(2,'read_3d_complex_array:: Something went wrong reading key'//trim(key))
if (dbg) then
  write(u6,*) 'read_3d_complex_array:: key =',trim(key)
  write(u6,*) 'read_3d_complex_array::  n1 =',i
  write(u6,*) 'read_3d_complex_array::  n2 =',j
  write(u6,*) 'read_3d_complex_array::  n3 =',k
end if
if ((i /= n1) .or. (j /= n2) .or. (k /= n3)) call WarningMessage(2,'read_3d_complex_array:: sizes of the array are different '// &
                                                                 'from the ones used to CALL this SUBROUTINE')
call mma_allocate(rr,n1,n2,n3,label='rr')
call mma_allocate(ri,n1,n2,n3,label='ri')
rr(:,:,:) = Zero
ri(:,:,:) = Zero
do i=1,n1
  do j=1,n2
    read(LU,*,iostat=ierr) (rr(i,j,k),ri(i,j,k),k=1,n3)
    if (ierr /= 0) call WarningMessage(2,'read_3d_complex_array:: Something went wrong reading the array.')
    if (dbg) write(u6,*) 'read_3d_complex_array::  i,j =',i,j
  end do
end do
array(:,:,:) = cmplx(rr(:,:,:),ri(:,:,:),kind=wp)
call mma_deallocate(rr)
call mma_deallocate(ri)

return

end subroutine read_3d_complex_array

subroutine read_4d_complex_array(LU,key,n1,n2,n3,n4,array,dbg)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU, n1, n2, n3, n4
character(len=*), intent(in) :: key
complex(kind=wp), intent(out) :: array(n1,n2,n3,n4)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr, j, k, l
character(len=500) :: line
real(kind=wp), allocatable :: ri(:,:,:,:), rr(:,:,:,:)

ierr = 0
if ((n1 <= 0) .or. (n2 <= 0) .or. (n3 <= 0) .or. (n4 <= 0)) then
  call WarningMessage(1,'read_4d_complex_array:: nothing to read. Array size = 0.')
  if (dbg) then
    write(u6,*) 'read_4d_complex_array::   n1 =',n1
    write(u6,*) 'read_4d_complex_array::   n2 =',n2
    write(u6,*) 'read_4d_complex_array::   n3 =',n3
    write(u6,*) 'read_4d_complex_array::   n4 =',n4
  end if
  return
end if
rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)

read(LU,*,iostat=ierr) i,j,k,l
if (ierr /= 0) call WarningMessage(2,'read_4d_complex_array:: Something went wrong reading key'//trim(key))
if (dbg) then
  write(u6,*) 'read_4d_complex_array:: key =',trim(key)
  write(u6,*) 'read_4d_complex_array::  n1 =',i
  write(u6,*) 'read_4d_complex_array::  n2 =',j
  write(u6,*) 'read_4d_complex_array::  n3 =',k
  write(u6,*) 'read_4d_complex_array::  n4 =',l
end if
if ((i /= n1) .or. (j /= n2) .or. (k /= n3) .or. (l /= n4)) &
  call WarningMessage(2,'read_4d_complex_array:: sizes of the array are different from the ones used to CALL this SUBROUTINE')
call mma_allocate(rr,n1,n2,n3,n4,label='rr')
call mma_allocate(ri,n1,n2,n3,n4,label='ri')
rr(:,:,:,:) = Zero
ri(:,:,:,:) = Zero
do i=1,n1
  do j=1,n2
    do k=1,n3
      read(LU,*,iostat=ierr) (rr(i,j,k,l),ri(i,j,k,l),l=1,n4)
      if (ierr /= 0) call WarningMessage(2,'read_4d_real_array:: Something went wrong reading the array.')
      if (dbg) write(u6,*) 'read_4d_real_array::  i,j,k =',i,j,k
    end do
  end do
end do
array(:,:,:,:) = cmplx(rr(:,:,:,:),ri(:,:,:,:),kind=wp)
call mma_deallocate(rr)
call mma_deallocate(ri)

return

end subroutine read_4d_complex_array

!----------------------------------------------------------------------!
!  WRITING SUBROUTINES
!----------------------------------------------------------------------!

subroutine write_INTEGER_scalar(LU,key,i,dbg)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: LU, i
character(len=*), intent(in) :: key
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: ierr
character(len=500) :: line
character(len=*), parameter :: FMTI = '(20(I0,1x))'

ierr = 0
rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)
! if key is found, THEN rewrite the data
if (ierr == 0) then
  write(LU,FMTI,iostat=ierr) i
  ! if the keyword is not found, THEN append the new keyword and data
else if (ierr /= 0) then
  write(LU,'(A)',iostat=ierr)
  write(LU,'(A)',iostat=ierr) trim(key)
  if (ierr /= 0) call WarningMessage(1,'write_INTEGER_scalar:: Something went wrong writing key'//trim(key))
  write(LU,FMTI,iostat=ierr) i
  if (ierr /= 0) call WarningMessage(1,'write_INTEGER_scalar:: Something went wrong writing data')
end if
write(LU,*,iostat=ierr)
call xFlush(LU)

return

end subroutine write_INTEGER_scalar

subroutine write_real_scalar(LU,key,r,dbg)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LU
character(len=*), intent(in) :: key
real(kind=wp), intent(in) :: r
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: ierr
character(len=500) :: line
character(len=*), parameter :: FMTR = '(5ES22.14)'

ierr = 0
rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)
if (ierr == 0) then
  ! if key is found, THEN rewrite the data
  write(LU,FMTR,iostat=ierr) r
else if (ierr /= 0) then
  ! if the keyword is not found, THEN append the new keyword and data
  write(LU,'(A)',iostat=ierr)
  write(LU,'(A)',iostat=ierr) trim(key)
  if (ierr /= 0) call WarningMessage(1,'write_real_scalar:: Something went wrong writing key'//trim(key))
  write(LU,FMTR,iostat=ierr) r
  if (ierr /= 0) call WarningMessage(1,'write_real_scalar:: Something went wrong writing data')
end if
write(LU,*,iostat=ierr)
call xFlush(LU)

return

end subroutine write_real_scalar

subroutine write_complex_scalar(LU,key,c,dbg)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LU
character(len=*), intent(in) :: key
complex(kind=wp), intent(in) :: c
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: ierr
character(len=500) :: line
character(len=*), parameter :: FMTC = '(3(2ES22.14))'

ierr = 0
rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)
if (ierr == 0) then
  ! if key is found, THEN rewrite the data
  write(LU,FMTC,iostat=ierr) c
else if (ierr /= 0) then
  ! if the keyword is not found, THEN append the new keyword and data
  write(LU,'(A)',iostat=ierr)
  write(LU,'(A)',iostat=ierr) trim(key)
  if (ierr /= 0) call WarningMessage(1,'write_complex_scalar:: Something went wrong writing key'//trim(key))
  write(LU,FMTC,iostat=ierr) c
  if (ierr /= 0) call WarningMessage(1,'write_complex_scalar:: Something went wrong writing data')
end if
write(LU,*,iostat=ierr)
call xFlush(LU)

return

end subroutine write_complex_scalar

subroutine write_string(LU,key,s,dbg)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: LU
character(len=*), intent(in) :: key, s
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: ierr
character(len=500) :: line

ierr = 0
rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)
if (ierr == 0) then
  ! if key is found, THEN rewrite the data
  write(LU,'(A   )',iostat=ierr) trim(s)
else if (ierr /= 0) then
  ! if the keyword is not found, THEN append the new keyword and data
  write(LU,'(A)',iostat=ierr)
  write(LU,'(A   )',iostat=ierr) trim(key)
  if (ierr /= 0) call WarningMessage(1,'write_string:: Something went wrong writing key'//trim(key))
  write(LU,'(100A)',iostat=ierr) trim(s)
  if (ierr /= 0) call WarningMessage(1,'write_string:: Something went wrong writing data')
end if
write(LU,*,iostat=ierr)
call xFlush(LU)

return

end subroutine write_string

subroutine write_1d_INTEGER_array(LU,key,n,array,dbg)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: LU, n, array(n)
character(len=*), intent(in) :: key
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr
character(len=500) :: line
character(len=*), parameter :: FMTI = '(20(I0,1x))'

ierr = 0
if (n <= 0) then
  call WarningMessage(1,'write_1d_INTEGER_array:: nothing to write. Array size = 0.')
  return
end if
if (sum(abs(array)) == 0) call WarningMessage(1,'write_1d_INTEGER_array:: all array elements are zero = 0.')

rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)

if (ierr == 0) then
  ! if key is found, THEN rewrite the data
  write(LU,FMTI,iostat=ierr) n
  write(LU,FMTI,iostat=ierr) (array(i),i=1,n)
else if (ierr /= 0) then
  ! if the keyword is not found, THEN append the new keyword and data
  write(LU,'(A)',iostat=ierr)
  write(LU,'(A)',iostat=ierr) trim(key)
  if (ierr /= 0) call WarningMessage(1,'write_1d_INTEGER_array:: Something went wrong writing key'//trim(key))
  write(LU,FMTI,iostat=ierr) n
  write(LU,FMTI,iostat=ierr) (array(i),i=1,n)
  if (ierr /= 0) call WarningMessage(1,'write_1d_INTEGER_array:: Something went wrong writing data')
end if
write(LU,*,iostat=ierr)
call xFlush(LU)

return

end subroutine write_1d_INTEGER_array

subroutine write_2d_INTEGER_array(LU,key,n1,n2,array,dbg)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU, n1, n2, array(n1,n2)
character(len=*), intent(in) :: key
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr, j
character(len=500) :: line
character(len=*), parameter :: FMTI = '(20(I0,1x))'

ierr = 0
if ((n1 <= 0) .or. (n2 <= 0)) then
  call WarningMessage(1,'write_2d_INTEGER_array:: nothing to write. Array size = 0.')
  return
end if
if (sum(abs(array)) == 0) call WarningMessage(1,'write_2d_INTEGER_array:: all array elements are zero = 0.')

rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)

if (ierr == 0) then
  ! if key is found, THEN rewrite the data
  write(LU,FMTI,iostat=ierr) n1,n2
  do i=1,n1
    write(LU,FMTI,iostat=ierr) (array(i,j),j=1,n2)
    if (ierr /= 0) call WarningMessage(2,'write_2d_INTEGER_array:: Something went wrong writing the array.')
    if (dbg) write(u6,*) 'write_2d_INTEGER_array::  i =',i
  end do
else if (ierr /= 0) then
  ! if the keyword is not found, THEN append the new keyword and data
  write(LU,'(A)',iostat=ierr)
  write(LU,'(A)',iostat=ierr) trim(key)
  if (ierr /= 0) call WarningMessage(1,'write_2d_INTEGER_array:: Something went wrong writing key'//trim(key))
  write(LU,FMTI,iostat=ierr) n1,n2
  do i=1,n1
    write(LU,FMTI,iostat=ierr) (array(i,j),j=1,n2)
    if (ierr /= 0) call WarningMessage(2,'write_2d_INTEGER_array:: Something went wrong writing data.')
    if (dbg) write(u6,*) 'write_2d_INTEGER_array::  i =',i
  end do
end if
write(LU,*,iostat=ierr)
call xFlush(LU)

return

end subroutine write_2d_INTEGER_array

subroutine write_3d_INTEGER_array(LU,key,n1,n2,n3,array,dbg)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU, n1, n2, n3, array(n1,n2,n3)
character(len=*), intent(in) :: key
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr, j, k
character(len=500) :: line
character(len=*), parameter :: FMTI = '(20(I0,1x))'

ierr = 0
if ((n1 <= 0) .or. (n2 <= 0) .or. (n3 <= 0)) then
  call WarningMessage(1,'write_3d_INTEGER_array:: nothing to write. Array size = 0.')
  return
end if
if (sum(abs(array)) == 0) call WarningMessage(1,'write_3d_INTEGER_array:: all array elements are zero = 0.')

rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)

if (ierr == 0) then
  ! if key is found, THEN rewrite the data
  write(LU,FMTI,iostat=ierr) n1,n2,n3
  do i=1,n1
    do j=1,n2
      write(LU,FMTI,iostat=ierr) (array(i,j,k),k=1,n3)
      if (ierr /= 0) call WarningMessage(2,'write_3d_INTEGER_array:: Something went wrong writing the array.')
      if (dbg) write(u6,*) 'write_3d_INTEGER_array::  i,j =',i,j
    end do
  end do
else if (ierr /= 0) then
  ! if the keyword is not found, THEN append the new keyword and data
  write(LU,'(A)',iostat=ierr)
  write(LU,'(A)',iostat=ierr) trim(key)
  if (ierr /= 0) call WarningMessage(1,'write_3d_INTEGER_array:: Something went wrong writing key'//trim(key))
  write(LU,FMTI,iostat=ierr) n1,n2,n3
  do i=1,n1
    do j=1,n2
      write(LU,FMTI,iostat=ierr) (array(i,j,k),k=1,n3)
      if (ierr /= 0) call WarningMessage(2,'write_3d_INTEGER_array:: Something went wrong writing the array.')
      if (dbg) write(u6,*) 'write_3d_INTEGER_array::  i,j =',i,j
    end do
  end do
end if
write(LU,*,iostat=ierr)
call xFlush(LU)

return

end subroutine write_3d_INTEGER_array

subroutine write_4d_INTEGER_array(LU,key,n1,n2,n3,n4,array,dbg)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU, n1, n2, n3, n4, array(n1,n2,n3,n4)
character(len=*), intent(in) :: key
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr, j, k, l
character(len=500) :: line
character(len=*), parameter :: FMTI = '(20(I0,1x))'

ierr = 0
if ((n1 <= 0) .or. (n2 <= 0) .or. (n3 <= 0) .or. (n4 <= 0)) then
  call WarningMessage(1,'write_4d_INTEGER_array:: nothing to write. Array size = 0.')
  return
end if
if (sum(abs(array)) == 0) call WarningMessage(1,'write_4d_INTEGER_array:: all array elements are zero = 0.')

rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)

if (ierr == 0) then
  ! if key is found, THEN rewrite the data
  write(LU,FMTI,iostat=ierr) n1,n2,n3,n4
  do i=1,n1
    do j=1,n2
      do k=1,n3
        write(LU,FMTI,iostat=ierr) (array(i,j,k,l),l=1,n4)
        if (ierr /= 0) call WarningMessage(2,'write_4d_INTEGER_array:: Something went wrong reading the array.')
        if (dbg) write(u6,*) 'write_4d_INTEGER_array::  i,j,k =',i,j,k
      end do
    end do
  end do
else if (ierr /= 0) then
  ! if the keyword is not found, THEN append the new keyword and data
  write(LU,'(A)',iostat=ierr)
  write(LU,'(A)',iostat=ierr) trim(key)
  if (ierr /= 0) call WarningMessage(1,'write_4d_INTEGER_array:: Something went wrong writing key'//trim(key))
  write(LU,FMTI,iostat=ierr) n1,n2,n3,n4
  do i=1,n1
    do j=1,n2
      do k=1,n3
        write(LU,FMTI,iostat=ierr) (array(i,j,k,l),l=1,n4)
        if (ierr /= 0) call WarningMessage(2,'write_4d_INTEGER_array:: Something went wrong writting the array.')
        if (dbg) write(u6,*) 'write_4d_INTEGER_array::  i,j,k =',i,j,k
      end do
    end do
  end do
end if
write(LU,*,iostat=ierr)
call xFlush(LU)

return

end subroutine write_4d_INTEGER_array

subroutine write_1d_real_array(LU,key,n,array,dbg)

use Constants, only: Ten
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LU, n
character(len=*), intent(in) :: key
real(kind=wp), intent(in) :: array(n)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr
character(len=500) :: line
real(kind=wp), parameter :: MINIMAL_REAL = tiny(MINIMAL_REAL)*Ten
character(len=*), parameter :: FMTI = '(20(I0,1x))', FMTR = '(5ES22.14)'
real(kind=wp), external :: dnrm2_

ierr = 0
if (n <= 0) then
  call WarningMessage(1,'write_1d_real_array:: nothing to write. Array size = 0.')
  return
end if
if (dnrm2_(n,array,1) <= MINIMAL_REAL) call WarningMessage(1,'write_1d_real_array:: all array elements are zero = 0.0')

rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)

if (ierr == 0) then
  ! if key is found, THEN rewrite the data
  write(LU,FMTI,iostat=ierr) n
  write(LU,FMTR,iostat=ierr) (array(i),i=1,n)
else if (ierr /= 0) then
  ! if the keyword is not found, THEN append the new keyword and data
  write(LU,'(A)',iostat=ierr)
  write(LU,'(A)',iostat=ierr) trim(key)
  if (ierr /= 0) call WarningMessage(1,'write_1d_real_array:: Something went wrong writing key'//trim(key))
  write(LU,FMTI,iostat=ierr) n
  write(LU,FMTR,iostat=ierr) (array(i),i=1,n)
  if (ierr /= 0) call WarningMessage(1,'write_1d_real_array:: Something went wrong writing data')
end if
write(LU,*,iostat=ierr)
call xFlush(LU)

return

end subroutine write_1d_real_array

subroutine write_2d_real_array(LU,key,n1,n2,array,dbg)

use Constants, only: Ten
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU, n1, n2
character(len=*), intent(in) :: key
real(kind=wp), intent(in) :: array(n1,n2)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr, j
character(len=500) :: line
real(kind=wp), parameter :: MINIMAL_REAL = tiny(MINIMAL_REAL)*Ten
character(len=*), parameter :: FMTI = '(20(I0,1x))', FMTR = '(5ES22.14)'
real(kind=wp), external :: dnrm2_

ierr = 0
if ((n1 <= 0) .or. (n2 <= 0)) then
  call WarningMessage(1,'write_2d_real_array:: nothing to write. Array size = 0.')
  return
end if
if (dnrm2_(n1*n2,array,1) <= MINIMAL_REAL) call WarningMessage(1,'write_2d_real_array:: all array elements are zero = 0.')

rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)

if (ierr == 0) then
  ! if key is found, THEN rewrite the data
  write(LU,FMTI,iostat=ierr) n1,n2
  do i=1,n1
    write(LU,FMTR,iostat=ierr) (array(i,j),j=1,n2)
    if (ierr /= 0) call WarningMessage(2,'write_2d_real_array:: Something went wrong writing the array.')
    if (dbg) write(u6,*) 'write_2d_real_array::  i =',i
  end do
else if (ierr /= 0) then
  ! if the keyword is not found, THEN append the new keyword and data
  write(LU,'(A)',iostat=ierr)
  write(LU,'(A)',iostat=ierr) trim(key)
  if (ierr /= 0) call WarningMessage(1,'write_2d_real_array:: Something went wrong writing key'//trim(key))
  write(LU,FMTI,iostat=ierr) n1,n2
  do i=1,n1
    write(LU,FMTR,iostat=ierr) (array(i,j),j=1,n2)
    if (ierr /= 0) call WarningMessage(2,'write_2d_real_array:: Something went wrong writing data.')
    if (dbg) write(u6,*) 'write_2d_real_array::  i =',i
  end do
end if
write(LU,*,iostat=ierr)
call xFlush(LU)

return

end subroutine write_2d_real_array

subroutine write_3d_real_array(LU,key,n1,n2,n3,array,dbg)

use Constants, only: Ten
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU, n1, n2, n3
character(len=*), intent(in) :: key
real(kind=wp), intent(in) :: array(n1,n2,n3)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr, j, k
character(len=500) :: line
real(kind=wp), parameter :: MINIMAL_REAL = tiny(MINIMAL_REAL)*Ten
character(len=*), parameter :: FMTI = '(20(I0,1x))', FMTR = '(5ES22.14)'
real(kind=wp), external :: dnrm2_

ierr = 0
if ((n1 <= 0) .or. (n2 <= 0) .or. (n3 <= 0)) then
  call WarningMessage(1,'write_3d_real_array:: nothing to write. Array size = 0.')
  return
end if
if (dnrm2_(n1*n2*n3,array,1) <= MINIMAL_REAL) call WarningMessage(1,'write_3d_real_array:: all array elements are zero = 0.')

rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)

if (ierr == 0) then
  ! if key is found, THEN rewrite the data
  write(LU,FMTI,iostat=ierr) n1,n2,n3
  do i=1,n1
    do j=1,n2
      write(LU,FMTR,iostat=ierr) (array(i,j,k),k=1,n3)
      if (ierr /= 0) call WarningMessage(2,'write_3d_real_array:: Something went wrong writing the array.')
      if (dbg) write(u6,*) 'write_3d_real_array::  i,j =',i,j
    end do
  end do
else if (ierr /= 0) then
  ! if the keyword is not found, THEN append the new keyword and data
  write(LU,'(A)',iostat=ierr)
  write(LU,'(A)',iostat=ierr) trim(key)
  if (ierr /= 0) call WarningMessage(1,'write_3d_real_array:: Something went wrong writing key'//trim(key))
  write(LU,FMTI,iostat=ierr) n1,n2,n3
  do i=1,n1
    do j=1,n2
      write(LU,FMTR,iostat=ierr) (array(i,j,k),k=1,n3)
      if (ierr /= 0) call WarningMessage(2,'write_3d_real_array:: Something went wrong writing the array.')
      if (dbg) write(u6,*) 'write_3d_real_array::  i,j =',i,j
    end do
  end do
end if
write(LU,*,iostat=ierr)
call xFlush(LU)

return

end subroutine write_3d_real_array

subroutine write_4d_real_array(LU,key,n1,n2,n3,n4,array,dbg)

use Constants, only: Ten
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU, n1, n2, n3, n4
character(len=*), intent(in) :: key
real(kind=wp), intent(in) :: array(n1,n2,n3,n4)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr, j, k, l
character(len=500) :: line
real(kind=wp), parameter :: MINIMAL_REAL = tiny(MINIMAL_REAL)*Ten
character(len=*), parameter :: FMTI = '(20(I0,1x))', FMTR = '(5ES22.14)'
real(kind=wp), external :: dnrm2_

ierr = 0
if ((n1 <= 0) .or. (n2 <= 0) .or. (n3 <= 0) .or. (n4 <= 0)) then
  call WarningMessage(1,'write_4d_real_array:: nothing to write. Array size = 0.')
  return
end if
if (dnrm2_(n1*n2*n3*n4,array,1) <= MINIMAL_REAL) call WarningMessage(1,'write_4d_real_array:: all array elements are zero = 0.')

rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)

if (ierr == 0) then
  ! if key is found, THEN rewrite the data
  write(LU,FMTI,iostat=ierr) n1,n2,n3,n4
  do i=1,n1
    do j=1,n2
      do k=1,n3
        write(LU,FMTR,iostat=ierr) (array(i,j,k,l),l=1,n4)
        if (ierr /= 0) call WarningMessage(2,'write_4d_real_array:: Something went wrong reading the array.')
        if (dbg) write(u6,*) 'write_4d_real_array::  i,j,k =',i,j,k
      end do
    end do
  end do
else if (ierr /= 0) then
  ! if the keyword is not found, THEN append the new keyword and data
  write(LU,'(A)',iostat=ierr)
  write(LU,'(A)',iostat=ierr) trim(key)
  if (ierr /= 0) call WarningMessage(1,'write_4d_real_array:: Something went wrong writing key'//trim(key))
  write(LU,FMTI,iostat=ierr) n1,n2,n3,n4
  do i=1,n1
    do j=1,n2
      do k=1,n3
        write(LU,FMTR,iostat=ierr) (array(i,j,k,l),l=1,n4)
        if (ierr /= 0) call WarningMessage(2,'write_4d_real_array:: Something went wrong writting the array.')
        if (dbg) write(u6,*) 'write_4d_real_array::  i,j,k =',i,j,k
      end do
    end do
  end do
end if
write(LU,*,iostat=ierr)
call xFlush(LU)

return

end subroutine write_4d_real_array

subroutine write_1d_complex_array(LU,key,n,array,dbg)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LU, n
character(len=*), intent(in) :: key
complex(kind=wp), intent(in) :: array(n)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr
character(len=500) :: line
character(len=*), parameter :: FMTC = '(3(2ES22.14))', FMTI = '(20(I0,1x))'

ierr = 0
if (n <= 0) then
  call WarningMessage(1,'write_1d_complex_array:: nothing to write. Array size = 0.')
  return
end if

rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)

if (ierr == 0) then
  ! if key is found, THEN rewrite the data
  write(LU,FMTI,iostat=ierr) n
  write(LU,FMTC,iostat=ierr) (array(i),i=1,n)
else if (ierr /= 0) then
  ! if the keyword is not found, THEN append the new keyword and data
  write(LU,'(A)',iostat=ierr)
  write(LU,'(A)',iostat=ierr) trim(key)
  if (ierr /= 0) call WarningMessage(1,'write_1d_complex_array:: Something went wrong writing key'//trim(key))
  write(LU,FMTI,iostat=ierr) n
  write(LU,FMTC,iostat=ierr) (array(i),i=1,n)
  if (ierr /= 0) call WarningMessage(1,'write_1d_complex_array:: Something went wrong writing data')
end if
write(LU,*,iostat=ierr)
call xFlush(LU)

return

end subroutine write_1d_complex_array

subroutine write_2d_complex_array(LU,key,n1,n2,array,dbg)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU, n1, n2
character(len=*), intent(in) :: key
complex(kind=wp), intent(in) :: array(n1,n2)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr, j
character(len=500) :: line
character(len=*), parameter :: FMTC = '(3(2ES22.14))', FMTI = '(20(I0,1x))'

ierr = 0
if ((n1 <= 0) .or. (n2 <= 0)) then
  call WarningMessage(1,'write_2d_complex_array:: nothing to write. Array size = 0.')
  return
end if

rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)

if (ierr == 0) then
  ! if key is found, THEN rewrite the data
  write(LU,FMTI,iostat=ierr) n1,n2
  do i=1,n1
    write(LU,FMTC,iostat=ierr) (array(i,j),j=1,n2)
    if (ierr /= 0) call WarningMessage(2,'write_2d_complex_array:: Something went wrong writing the array.')
    if (dbg) write(u6,*) 'write_2d_complex_array::  i =',i
  end do
else if (ierr /= 0) then
  ! if the keyword is not found, THEN append the new keyword and data
  write(LU,'(A)',iostat=ierr)
  write(LU,'(A)',iostat=ierr) trim(key)
  if (ierr /= 0) call WarningMessage(1,'write_2d_complex_array:: Something went wrong writing key'//trim(key))
  write(LU,FMTI,iostat=ierr) n1,n2
  do i=1,n1
    write(LU,FMTC,iostat=ierr) (array(i,j),j=1,n2)
    if (ierr /= 0) call WarningMessage(2,'write_2d_complex_array:: Something went wrong writing data.')
    if (dbg) write(u6,*) 'write_2d_complex_array::  i =',i
  end do
end if
write(LU,*,iostat=ierr)
call xFlush(LU)

return

end subroutine write_2d_complex_array

subroutine write_3d_complex_array(LU,key,n1,n2,n3,array,dbg)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU, n1, n2, n3
character(len=*), intent(in) :: key
complex(kind=wp), intent(in) :: array(n1,n2,n3)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr, j, k
character(len=500) :: line
character(len=*), parameter :: FMTC = '(3(2ES22.14))', FMTI = '(20(I0,1x))'

ierr = 0
if ((n1 <= 0) .or. (n2 <= 0) .or. (n3 <= 0)) then
  call WarningMessage(1,'write_3d_complex_array:: nothing to write. Array size = 0.')
  return
end if

rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)

if (ierr == 0) then
  ! if key is found, THEN rewrite the data
  write(LU,FMTI,iostat=ierr) n1,n2,n3
  do i=1,n1
    do j=1,n2
      write(LU,FMTC,iostat=ierr) (array(i,j,k),k=1,n3)
      if (ierr /= 0) call WarningMessage(2,'write_3d_complex_array:: Something went wrong writing the array.')
      if (dbg) write(u6,*) 'write_3d_complex_array::  i,j =',i,j
    end do
  end do
else if (ierr /= 0) then
  ! if the keyword is not found, THEN append the new keyword and data
  write(LU,'(A)',iostat=ierr)
  write(LU,'(A)',iostat=ierr) trim(key)
  if (ierr /= 0) call WarningMessage(1,'write_3d_complex_array:: Something went wrong writing key'//trim(key))
  write(LU,FMTI,iostat=ierr) n1,n2,n3
  do i=1,n1
    do j=1,n2
      write(LU,FMTC,iostat=ierr) (array(i,j,k),k=1,n3)
      if (ierr /= 0) call WarningMessage(2,'write_3d_complex_array:: Something went wrong writing the array.')
      if (dbg) write(u6,*) 'write_3d_complex_array::  i,j =',i,j
    end do
  end do
end if
write(LU,*,iostat=ierr)
call xFlush(LU)

return

end subroutine write_3d_complex_array

subroutine write_4d_complex_array(LU,key,n1,n2,n3,n4,array,dbg)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LU, n1, n2, n3, n4
character(len=*), intent(in) :: key
complex(kind=wp), intent(in) :: array(n1,n2,n3,n4)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: i, ierr, j, k, l
character(len=500) :: line
character(len=*), parameter :: FMTC = '(3(2ES22.14))', FMTI = '(20(I0,1x))'

ierr = 0
if ((n1 <= 0) .or. (n2 <= 0) .or. (n3 <= 0) .or. (n4 <= 0)) then
  call WarningMessage(1,'write_4d_complex_array:: nothing to write. Array size = 0.')
  return
end if

rewind(LU)
call file_advance_to_string(LU,key,line,ierr,dbg)

if (ierr == 0) then
  ! if key is found, THEN rewrite the data
  write(LU,FMTI,iostat=ierr) n1,n2,n3,n4
  do i=1,n1
    do j=1,n2
      do k=1,n3
        write(LU,FMTC,iostat=ierr) (array(i,j,k,l),l=1,n4)
        if (ierr /= 0) call WarningMessage(2,'write_4d_complex_array:: Something went wrong reading the array.')
        if (dbg) write(u6,*) 'write_4d_complex_array::  i,j,k =',i,j,k
      end do
    end do
  end do
else if (ierr /= 0) then
  ! if the keyword is not found, THEN append the new keyword and data
  write(LU,'(A)',iostat=ierr)
  write(LU,'(A)',iostat=ierr) trim(key)
  if (ierr /= 0) call WarningMessage(1,'write_4d_complex_array:: Something went wrong writing key'//trim(key))
  write(LU,FMTI,iostat=ierr) n1,n2,n3,n4
  do i=1,n1
    do j=1,n2
      do k=1,n3
        write(LU,FMTC,iostat=ierr) (array(i,j,k,l),l=1,n4)
        if (ierr /= 0) call WarningMessage(2,'write_4d_complex_array:: Something went wrong writting the array.')
        if (dbg) write(u6,*) 'write_4d_complex_array::  i,j,k =',i,j,k
      end do
    end do
  end do
end if
write(LU,*,iostat=ierr)
call xFlush(LU)

return

end subroutine write_4d_complex_array

!END MODULE io_data
