MODULE io_data

IMPLICIT NONE

PUBLIC ::      read_magnetic_moment,       write_magnetic_moment,  &
               read_electric_moment,       write_electric_moment,  &
               read_spin_moment,           write_spin_moment,      &
               read_angmom,                write_angmom,           &
               read_edipmom,               write_edipmom,          &
               read_amfi,                  write_amfi,             &
               read_nss,                   write_nss,              &
               read_nmult,                 write_nmult,            &
               read_imult,                 write_imult,            &
               read_format,                write_format,           &
               read_nroot,                 write_nroot,            &
               read_nstate,                write_nstate,           &
               read_multiplicity,          write_multiplicity,     &
               read_szproj,                write_szproj,           &
               read_eso,                   write_eso,              &
               read_esfs,                  write_esfs,             &
               read_hso,                   write_hso,              &
               read_eigen,                 write_eigen,            &
               read_gtens,                 write_gtens,            &
               read_stev_cfp,              write_stev_cfp,         &
               read_susc,                  write_susc,             &
               read_magn,                  write_magn,             &
               read_complex_matrix,        write_complex_matrix,   &
               open_datafile_read,         open_datafile_write,    &
               close_datafile,                                     &
               check_hermiticity_matrix,   check_commutation,      &
               check_S_square


INTEGER, PRIVATE                       :: DATA_FILE
INTEGER, PRIVATE                       :: ANISO_FILE
INTEGER, PARAMETER, PRIVATE            :: FILE_SUS=91
INTEGER, PARAMETER, PRIVATE            :: FILE_MAG=92
INTEGER, PARAMETER, PRIVATE            :: FILE_CFP=93
INTEGER, PARAMETER, PRIVATE            :: FILE_GTE=94

CHARACTER (LEN=20), PARAMETER, PRIVATE :: FILE_NAME_SUS='FILE_SUS'
CHARACTER (LEN=20), PARAMETER, PRIVATE :: FILE_NAME_MAG='FILE_MAG'
CHARACTER (LEN=20), PARAMETER, PRIVATE :: FILE_NAME_CFP='FILE_CFP'
CHARACTER (LEN=20), PARAMETER, PRIVATE :: FILE_NAME_GTE='FILE_GTE'

INTEGER, PARAMETER, PRIVATE            :: wp = kind(0.d0)
REAL(wp), PARAMETER, PRIVATE           :: One = 1.0_wp
REAL(wp), PARAMETER, PRIVATE           :: Zero = 0.0_wp
REAL(wp), PARAMETER, PRIVATE           :: MINIMAL_REAL=TINY(0.0_wp)*10.0_wp
COMPLEX(wp), PARAMETER, PRIVATE        :: ZeroC = (0.0_wp,0.0_wp)
COMPLEX(wp), PARAMETER, PRIVATE        :: OneC = (1.0_wp,0.0_wp)
COMPLEX(wp), PARAMETER, PRIVATE        :: cI = (0.0_wp,1.0_wp)
LOGICAL, PRIVATE :: DBG=.false.
INTEGER, PRIVATE :: StdIn=5
INTEGER, PRIVATE :: StdOut=6

PRIVATE ::     file_advance_to_string,     inquire_key_presence,   &
               write_string,                                       &
               read_INTEGER_scalar,        write_INTEGER_scalar,   &
               read_real_scalar,           write_real_scalar,      &
               read_complex_scalar,        write_complex_scalar,   &
               read_1d_size,                                       &
               read_2d_size,                                       &
               read_3d_size,                                       &
               read_4d_size,                                       &
               read_1d_INTEGER_array,      write_1d_INTEGER_array, &
               read_2d_INTEGER_array,      write_2d_INTEGER_array, &
               read_3d_INTEGER_array,      write_3d_INTEGER_array, &
               read_4d_INTEGER_array,      write_4d_INTEGER_array, &
               read_1d_real_array,         write_1d_real_array,    &
               read_2d_real_array,         write_2d_real_array,    &
               read_3d_real_array,         write_3d_real_array,    &
               read_4d_real_array,         write_4d_real_array,    &
               read_1d_complex_array,      write_1d_complex_array, &
               read_2d_complex_array,      write_2d_complex_array, &
               read_3d_complex_array,      write_3d_complex_array, &
               read_4d_complex_array,      write_4d_complex_array

INTEGER, PRIVATE            :: ierr
CHARACTER (LEN=500), PRIVATE :: LINE
CHARACTER (LEN=20), PRIVATE  :: FMTR='(5ES22.14)', FMTI='(20(I0,1x))', FMTC='(3(2ES22.14))'
CHARACTER (LEN=30), PRIVATE  :: FMTCFP='(2(I0,1x),ES22.14)'


CONTAINS


! Keywords:
!    $dipm   --  electric dipole moment, SO basis  (6 components: _xr,_xi,  _yr,_yi,  _zr,_zi)
!    $magn   --  magnetic dipole moment, SO basis  (6 components: _xr,_xi,  _yr,_yi,  _zr,_zi)
!    $angmom --  orbital moment (L)    , SF basis  (3 compinents: _x, _y, _z)
!    $edip   --  electric dipole moment, SF basis  (3 components: _x, _y, _z)
!    $hso    --  spin-orbit hamiltonian, SO basis  (2 components: _r, _i)
!    $eigen  --  spin-orbit eigenstates, SO basis  (2 components: _r, _i)
!
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
! xt,               &
! xt_field,         &
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
! x,                &
! x_field,          &
! x_minus_one,      &
! x_field,          &
! x_field_minus_one,&
! magnetization,    &
!
!
!
! read_MOs
! read_CI_in_SD
! read_CI_in_CSF
! read_coeff_of_CSF
!--------------------------





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
!--------------------------------------------------------------------------------------------------!


!--------------------------------------------------------------------------------------------------!
!           HIGH LEVEL VERIFICATION SUBROUTINES
!--------------------------------------------------------------------------------------------------!
SUBROUTINE check_commutation(n,moment)
   ! valid for S, L, J in spin-orbit basis
   IMPLICIT NONE
   INTEGER, INTENT (IN)      :: N
   COMPLEX (wp), INTENT (IN) :: moment(3,N,N)
   COMPLEX (wp), ALLOCATABLE :: XY(:,:), YX(:,:), YZ(:,:), ZY(:,:), ZX(:,:), XZ(:,:)
   COMPLEX (wp) :: tr
   INTEGER     :: i, j

   ! verify the commutation relations for S
   ALLOCATE (XY(n,n))
   ALLOCATE (YX(n,n))
   ALLOCATE (YZ(n,n))
   ALLOCATE (ZY(n,n))
   ALLOCATE (ZX(n,n))
   ALLOCATE (XZ(n,n))
   XY=ZeroC
   YX=ZeroC
   CALL zgemm('c', 'n', n, n, n, OneC,      &
              moment(1,1:n,1:n), n,         &
              moment(2,1:n,1:n), n, ZeroC,  &
               XY, n )

   CALL zgemm('c', 'n', n, n, n, OneC,      &
              moment(2,1:n,1:n), n,         &
              moment(1,1:n,1:n), n, ZeroC,  &
               YX, n )

   YZ=ZeroC
   ZY=ZeroC
   CALL zgemm('c', 'n', n, n, n, OneC,      &
              moment(2,1:n,1:n), n,         &
              moment(3,1:n,1:n), n, ZeroC,  &
               YZ, n )

   CALL zgemm('c', 'n', n, n, n, OneC,      &
              moment(3,1:n,1:n), n,         &
              moment(2,1:n,1:n), n, ZeroC,  &
               ZY, n )

   ZX=ZeroC
   XZ=ZeroC
   CALL zgemm('c', 'n', n, n, n, OneC,      &
              moment(3,1:n,1:n), n,         &
              moment(1,1:n,1:n), n, ZeroC,  &
               ZX, n )

   CALL zgemm('c', 'n', n, n, n, OneC,      &
              moment(1,1:n,1:n), n,         &
              moment(3,1:n,1:n), n, ZeroC,  &
               XZ, n )

   tr=ZeroC
   DO i=1,n
     DO j=1,n
       tr = tr  + XY(i,j) - YX(i,j) - cI*moment(3,i,j) &
                + YZ(i,j) - ZY(i,j) - cI*moment(1,i,j) &
                + ZX(i,j) - XZ(i,j) - cI*moment(2,i,j)
     END DO
   END DO
   IF (DBG) WRITE (StdOut,'(A,ES22.14)') 'check_commutation::  trace of [Sx,Sy]-iSz = ', ABS(tr)
   IF (ABS(tr)>1.0e-6_wp) THEN
      CALL WarningMessage(1,'check_commutation:: trace of [Sx,Sy]-iSz  is larger than 1.0e-6. '//&
                            'The input moment looks inacurate')
   ELSE
      WRITE (StdOut,'(A,ES22.14)') 'check_commutation:  The input moment passes all three commutation tests.'
   END IF

   DEALLOCATE (XY)
   DEALLOCATE (YX)
   DEALLOCATE (YZ)
   DEALLOCATE (ZY)
   DEALLOCATE (ZX)
   DEALLOCATE (XZ)

   RETURN
END SUBROUTINE check_commutation


SUBROUTINE check_S_square (n, moment)
!  valid also for J, L, S
   IMPLICIT NONE
   INTEGER, INTENT (IN)      :: n
   COMPLEX (wp), INTENT (IN) :: moment(3,n,n)
   COMPLEX (wp), ALLOCATABLE :: X2(:,:), Y2(:,:), Z2(:,:), S2(:,:)
   COMPLEX (wp) :: tr
   REAL (wp)    :: S2_theoretic
   INTEGER     :: i

   ALLOCATE (X2(n,n))
   ALLOCATE (Y2(n,n))
   ALLOCATE (Z2(n,n))
   ALLOCATE (S2(n,n))
   X2=ZeroC
   Y2=ZeroC
   Z2=ZeroC
   S2=ZeroC

   CALL zgemm('c', 'n', n, n, n, OneC,      &
              moment(1,1:n,1:n), n,         &
              moment(1,1:n,1:n), n, ZeroC,  &
               X2, n )

   CALL zgemm('c', 'n', n, n, n, OneC,      &
              moment(2,1:n,1:n), n,         &
              moment(2,1:n,1:n), n, ZeroC,  &
               Y2, n )

   CALL zgemm('c', 'n', n, n, n, OneC,      &
              moment(3,1:n,1:n), n,         &
              moment(3,1:n,1:n), n, ZeroC,  &
               Z2, n )

   ! matrix add:
   S2 = X2 + Y2 + Z2

   tr=ZeroC
   DO i=1,n
     tr = tr + S2(i,i)
   END DO

   S2_theoretic=Zero
   S2_theoretic=DBLE(n**3-n)/4.0_wp

   IF (DBG) WRITE (StdOut,'(A,ES22.14)') 'check_S_square::  trace of S2=(Sx).Sx+(Sy).Sy+(Sz).Sz = ', ABS(tr)
   IF ( (ABS(tr) - S2_theoretic) > 1.0e-6_wp) THEN
      CALL WarningMessage(1,'check_S_square:: tr(S^2) - S2_theoretic is larger than 1.0e-6. '//&
                            'The input moment looks inacurate')
   ELSE
      WRITE (StdOut,'(A,ES22.14)') 'check_S_square:  The input moment passes the S^2 test.'
   END IF
   DEALLOCATE (X2)
   DEALLOCATE (Y2)
   DEALLOCATE (Z2)
   DEALLOCATE (S2)

   RETURN
END SUBROUTINE check_S_square

SUBROUTINE check_hermiticity_moment (n, moment)
   IMPLICIT NONE
   INTEGER, INTENT (IN)      :: n
   COMPLEX (wp), INTENT (IN)  :: moment(3,n,n)
   COMPLEX (wp) :: c
   INTEGER     :: i, j, l
   ! build difference SUM ( MATRIX(i,j) - DCONJG(MATRIX(j,i)) )
   c=ZeroC
   DO i=1,n
      DO j=1,n
         IF ( i==j ) CYCLE
         DO l=1,3
            c = c + moment(l,i,j) - DCONJG( moment(l,j,i) )
         END DO
      END DO
   END DO
   IF (DBG) WRITE (StdOut,'(A,2ES22.14)') 'check_hermiticity_moment::  trace of A(i,j)-DCONJG(A(j,i)) = ', c
   IF (ABS(c) > 1.0e-6_wp) THEN
      CALL WarningMessage(1,'check_hermiticity_moment:: trace of M(1:3,i,j)-DCONJG(A(1:3,j,i)) is larger than 1.0e-6. '//&
                            'The hermiticity of input moment is not quite fulfilled')
   ELSE
      WRITE (StdOut,'(A,ES22.14)') 'check_hermiticity_moment:  The input moment passes the hermiticity test.'
   END IF
   RETURN
END SUBROUTINE check_hermiticity_moment

SUBROUTINE check_hermiticity_matrix (n, matrix)
   IMPLICIT NONE
   INTEGER, INTENT (IN)      :: n
   COMPLEX (wp), INTENT (IN) :: matrix(n,n)
   COMPLEX (wp) :: c
   INTEGER     :: i, j
   ! build difference ( MATRIX(i,j) - DCONJG(MATRIX(j,i)) )
   c=ZeroC
   DO i=1,n
      DO j=i,n
         IF ( i==j ) CYCLE
         c = c + ( matrix(i,j) - DCONJG(matrix(j,i)) )
      END DO
   END DO
   IF (DBG) WRITE (StdOut,'(A,2ES22.14)') 'check_hermiticity_matrix::  trace of A(i,j)-DCONJG(A(j,i)) = ', c
   IF (ABS(c) > 1.0e-6_wp) THEN
      CALL WarningMessage(1,'check_hermiticity_matrix:: trace of A(i,j)-DCONJG(A(j,i)) is larger than 1.0e-6. '//&
                            'The hermiticity of input matrix is not quite fulfilled')
   ELSE
      WRITE (StdOut,'(A,ES22.14)') 'check_hermiticity_matrix:  The input matrix passes the hermiticity test.'
   END IF
   RETURN
END SUBROUTINE check_hermiticity_matrix





!--------------------------------------------------------------------------------------------------!
!           HIGH LEVEL READING SUBROUTINES
!--------------------------------------------------------------------------------------------------!
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
!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_magnetic_moment (DATA_FILE, n, moment)
   IMPLICIT NONE
   INTEGER, INTENT (IN)       :: DATA_FILE
   INTEGER, INTENT (IN)       :: N
   COMPLEX (wp), INTENT (OUT) :: moment(3,N,N)
   REAL (wp), ALLOCATABLE     :: rr(:,:), ri(:,:)
   INTEGER                    :: i, j
   REAL (wp), EXTERNAL        :: dnrm2_, dznrm2_

   moment=ZeroC
   ALLOCATE (rr(n,n))
   ALLOCATE (ri(n,n))
! projection X
   rr=Zero; ri=Zero;
   IF (inquire_key_presence( DATA_FILE, '$magn_xr')) CALL read_2d_real_array( DATA_FILE, '$magn_xr',n,n,rr)
   IF (inquire_key_presence( DATA_FILE, '$magn_xi')) CALL read_2d_real_array( DATA_FILE, '$magn_xi',n,n,ri)
   IF (dbg) WRITE (StdOut,*) 'read_magnetic_moment::  norm of moment_xr=',dnrm2_(n*n,rr,1)
   IF (dbg) WRITE (StdOut,*) 'read_magnetic_moment::  norm of moment_xi=',dnrm2_(n*n,ri,1)
   DO i=1,n
      DO j=1,n
         moment(1,i,j) = CMPLX(rr(i,j),ri(i,j),wp)
      END DO
   END DO
   IF (dbg) CALL check_hermiticity_matrix(n,moment(1,1:n,1:n))
! projection Y
   rr=Zero; ri=Zero;
   IF (inquire_key_presence( DATA_FILE, '$magn_yr')) CALL read_2d_real_array( DATA_FILE, '$magn_yr',n,n,rr)
   IF (inquire_key_presence( DATA_FILE, '$magn_yi')) CALL read_2d_real_array( DATA_FILE, '$magn_yi',n,n,ri)
   IF (dbg) WRITE (StdOut,*) 'read_magnetic_moment::  norm of moment_yr=',dnrm2_(n*n,rr,1)
   IF (dbg) WRITE (StdOut,*) 'read_magnetic_moment::  norm of moment_yi=',dnrm2_(n*n,ri,1)
   DO i=1,n
      DO j=1,n
         moment(2,i,j) = CMPLX(rr(i,j),ri(i,j),wp)
      END DO
   END DO
   IF (dbg) CALL check_hermiticity_matrix(n,moment(2,1:n,1:n))
! projection Z
   rr=Zero; ri=Zero;
   IF (inquire_key_presence( DATA_FILE, '$magn_zr')) CALL read_2d_real_array( DATA_FILE, '$magn_zr',n,n,rr)
   IF (inquire_key_presence( DATA_FILE, '$magn_zi')) CALL read_2d_real_array( DATA_FILE, '$magn_zi',n,n,ri)
   IF (dbg) WRITE (StdOut,*) 'read_magnetic_moment::  norm of moment_zr=',dnrm2_(n*n,rr,1)
   IF (dbg) WRITE (StdOut,*) 'read_magnetic_moment::  norm of moment_zi=',dnrm2_(n*n,ri,1)
   DO i=1,n
      DO j=1,n
         moment(3,i,j) = CMPLX(rr(i,j),ri(i,j),wp)
      END DO
   END DO
   IF ( dznrm2_(3*n*n,moment,1) <=TINY(Zero)) CALL WarningMessage(1,'read_magnetic_moment:: the norm of the read moment is zero!')
   IF (dbg) CALL check_hermiticity_matrix(n,moment(3,1:n,1:n))
   DEALLOCATE (rr)
   DEALLOCATE (ri)
   RETURN
END SUBROUTINE read_magnetic_moment

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_electric_moment (DATA_FILE, n, moment)
   IMPLICIT NONE
   INTEGER, INTENT (IN)       :: DATA_FILE
   INTEGER, INTENT (IN)       :: N
   COMPLEX (wp), INTENT (OUT) :: moment(3,N,N)
   REAL (wp), ALLOCATABLE     :: rr(:,:), ri(:,:)
   INTEGER                    :: i, j
   REAL (wp), EXTERNAL        :: dnrm2_, dznrm2_

   moment=ZeroC
   ALLOCATE (rr(n,n))
   ALLOCATE (ri(n,n))
! projection X
   rr=Zero; ri=Zero;
   IF (inquire_key_presence( DATA_FILE, '$edipm_xr')) CALL read_2d_real_array( DATA_FILE, '$edipm_xr',n,n,rr)
   IF (inquire_key_presence( DATA_FILE, '$edipm_xi')) CALL read_2d_real_array( DATA_FILE, '$edipm_xi',n,n,ri)
   IF (dbg) WRITE (StdOut,*) 'read_electric_moment::  norm of moment_xr=',dnrm2_(n*n,rr,1)
   IF (dbg) WRITE (StdOut,*) 'read_electric_moment::  norm of moment_xi=',dnrm2_(n*n,ri,1)
   DO i=1,n
      DO j=1,n
         moment(1,i,j) = CMPLX(rr(i,j),ri(i,j),wp)
      END DO
   END DO
   IF (dbg) CALL check_hermiticity_matrix(n,moment(1,:,:))
! projection Y
   rr=Zero; ri=Zero;
   IF (inquire_key_presence( DATA_FILE, '$edipm_yr')) CALL read_2d_real_array( DATA_FILE, '$edipm_yr',n,n,rr)
   IF (inquire_key_presence( DATA_FILE, '$edipm_yi')) CALL read_2d_real_array( DATA_FILE, '$edipm_yi',n,n,ri)
   IF (dbg) WRITE (StdOut,*) 'read_electric_moment::  norm of moment_yr=',dnrm2_(n*n,rr,1)
   IF (dbg) WRITE (StdOut,*) 'read_electric_moment::  norm of moment_yi=',dnrm2_(n*n,ri,1)
   DO i=1,n
      DO j=1,n
         moment(2,i,j) = CMPLX(rr(i,j),ri(i,j),wp)
      END DO
   END DO
   IF (dbg) CALL check_hermiticity_matrix(n,moment(2,:,:))
! projection Z
   rr=Zero; ri=Zero;
   IF (inquire_key_presence( DATA_FILE, '$edipm_zr')) CALL read_2d_real_array( DATA_FILE, '$edipm_zr',n,n,rr)
   IF (inquire_key_presence( DATA_FILE, '$edipm_zi')) CALL read_2d_real_array( DATA_FILE, '$edipm_zi',n,n,ri)
   IF (dbg) WRITE (StdOut,*) 'read_electric_moment::  norm of moment_zr=',dnrm2_(n*n,rr,1)
   IF (dbg) WRITE (StdOut,*) 'read_electric_moment::  norm of moment_zi=',dnrm2_(n*n,ri,1)
   DO i=1,n
      DO j=1,n
         moment(3,i,j) = CMPLX(rr(i,j),ri(i,j),wp)
      END DO
   END DO
   IF ( dznrm2_(3*n*n,moment,1) <=TINY(Zero)) CALL WarningMessage(1,'read_electric_moment:: the norm of the read moment is zero!')
   IF (dbg) CALL check_hermiticity_matrix(n,moment(3,:,:))
   DEALLOCATE (rr)
   DEALLOCATE (ri)
   RETURN
END SUBROUTINE read_electric_moment

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_spin_moment (DATA_FILE, n, moment)
   IMPLICIT NONE
   INTEGER, INTENT (IN)       :: DATA_FILE
   INTEGER, INTENT (IN)       :: N
   COMPLEX (wp), INTENT (OUT) :: moment(3,N,N)
   REAL (wp), ALLOCATABLE     :: rr(:,:), ri(:,:)
   INTEGER                    :: i, j
   REAL (wp), EXTERNAL        :: dnrm2_, dznrm2_

   moment=ZeroC
   ALLOCATE (rr(n,n))
   ALLOCATE (ri(n,n))
! projection X
   rr=Zero; ri=Zero;
   IF (inquire_key_presence( DATA_FILE, '$spin_xr')) CALL read_2d_real_array( DATA_FILE, '$spin_xr',n,n,rr)
   IF (inquire_key_presence( DATA_FILE, '$spin_xi')) CALL read_2d_real_array( DATA_FILE, '$spin_xi',n,n,ri)
   IF (dbg) WRITE (StdOut,*) 'read_spin_moment::  norm of moment_zr=',dnrm2_(n*n,rr,1)
   IF (dbg) WRITE (StdOut,*) 'read_spin_moment::  norm of moment_zi=',dnrm2_(n*n,ri,1)
   DO i=1,n
      DO j=1,n
         moment(1,i,j) = CMPLX(rr(i,j),ri(i,j),wp)
      END DO
   END DO
   IF (dbg) CALL check_hermiticity_matrix(n,moment(1,:,:))
! projection Y
   rr=Zero; ri=Zero;
   IF (inquire_key_presence( DATA_FILE, '$spin_yr')) CALL read_2d_real_array( DATA_FILE, '$spin_yr',n,n,rr)
   IF (inquire_key_presence( DATA_FILE, '$spin_yi')) CALL read_2d_real_array( DATA_FILE, '$spin_yi',n,n,ri)
   IF (dbg) WRITE (StdOut,*) 'read_spin_moment::  norm of moment_zr=',dnrm2_(n*n,rr,1)
   IF (dbg) WRITE (StdOut,*) 'read_spin_moment::  norm of moment_zi=',dnrm2_(n*n,ri,1)
   DO i=1,n
      DO j=1,n
         moment(2,i,j) = CMPLX(rr(i,j),ri(i,j),wp)
      END DO
   END DO
   IF (dbg) CALL check_hermiticity_matrix(n,moment(2,:,:))
! projection Z
   rr=Zero; ri=Zero;
   IF (inquire_key_presence( DATA_FILE, '$spin_zr')) CALL read_2d_real_array( DATA_FILE, '$spin_zr',n,n,rr)
   IF (inquire_key_presence( DATA_FILE, '$spin_zi')) CALL read_2d_real_array( DATA_FILE, '$spin_zi',n,n,ri)
   IF (dbg) WRITE (StdOut,*) 'read_spin_moment::  norm of moment_zr=',dnrm2_(n*n,rr,1)
   IF (dbg) WRITE (StdOut,*) 'read_spin_moment::  norm of moment_zi=',dnrm2_(n*n,ri,1)
   DO i=1,n
      DO j=1,n
         moment(3,i,j) = CMPLX(rr(i,j),ri(i,j),wp)
      END DO
   END DO
   IF ( dznrm2_(3*n*n,moment,1) <=TINY(Zero)) CALL WarningMessage(1,'read_spin:: the norm of the read moment is zero!')
   IF (dbg) CALL check_hermiticity_matrix(n,moment(3,:,:))
   DEALLOCATE (rr)
   DEALLOCATE (ri)
   IF (dbg) CALL check_commutation(n,moment)
   RETURN
END SUBROUTINE read_spin_moment

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_angmom (DATA_FILE, n, moment)
   IMPLICIT NONE
   INTEGER, INTENT (IN)    :: DATA_FILE
   INTEGER, INTENT (IN)    :: N
   REAL (wp), INTENT (OUT) :: moment(3,N,N)
   REAL (wp), ALLOCATABLE  :: rr(:,:)
   INTEGER                 :: i, j
   REAL (wp), EXTERNAL     :: dnrm2_

   moment=Zero
   ALLOCATE (rr(n,n))
! projection X
   rr=Zero
   IF (inquire_key_presence( DATA_FILE, '$angmom_x')) CALL read_2d_real_array( DATA_FILE, '$angmom_x',n,n,rr)
   IF (dbg) WRITE (StdOut,*) 'read_angmom::  norm of moment_x=',dnrm2_(n*n,rr,1)
   DO i=1,n
      DO j=1,n
         moment(1,i,j) = rr(i,j)
      END DO
   END DO
! projection Y
   rr=Zero
   IF (inquire_key_presence( DATA_FILE, '$angmom_y')) CALL read_2d_real_array( DATA_FILE, '$angmom_y',n,n,rr)
   IF (dbg) WRITE (StdOut,*) 'read_angmom::  norm of moment_y=',dnrm2_(n*n,rr,1)
   DO i=1,n
      DO j=1,n
         moment(2,i,j) = rr(i,j)
      END DO
   END DO
! projection Z
   rr=Zero
   IF (inquire_key_presence( DATA_FILE, '$angmom_z')) CALL read_2d_real_array( DATA_FILE, '$angmom_z',n,n,rr)
   IF (dbg) WRITE (StdOut,*) 'read_angmom::  norm of moment_z=',dnrm2_(n*n,rr,1)
   DO i=1,n
      DO j=1,n
         moment(3,i,j) = rr(i,j)
      END DO
   END DO
   DEALLOCATE (rr)
   IF (dnrm2_(3*n*n,moment,1) <=TINY(Zero)) CALL WarningMessage(1,'read_angmom:: the norm of the read moment is zero!')
   RETURN
END SUBROUTINE read_angmom

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_edipmom (DATA_FILE, n, moment)
   IMPLICIT NONE
   INTEGER, INTENT (IN)    :: DATA_FILE
   INTEGER, INTENT (IN)    :: N
   REAL (wp), INTENT (OUT) :: moment(3,N,N)
   REAL (wp), ALLOCATABLE  :: rr(:,:)
   INTEGER                 :: i, j
   REAL (wp), EXTERNAL     :: dnrm2_

   moment=Zero
   ALLOCATE (rr(n,n))
! projection X
   rr=Zero
   IF (inquire_key_presence( DATA_FILE, '$edmom_x')) CALL read_2d_real_array( DATA_FILE, '$edmom_x',n,n,rr)
   IF (dbg) WRITE (StdOut,*) 'read_edipmom::  norm of moment_x=',dnrm2_(n*n,rr,1)
   DO i=1,n
      DO j=1,n
         moment(1,i,j) = rr(i,j)
      END DO
   END DO
! projection Y
   rr=Zero
   IF (inquire_key_presence( DATA_FILE, '$edmom_y')) CALL read_2d_real_array( DATA_FILE, '$edmom_y',n,n,rr)
   IF (dbg) WRITE (StdOut,*) 'read_edipmom::  norm of moment_y=',dnrm2_(n*n,rr,1)
   DO i=1,n
      DO j=1,n
         moment(2,i,j) = rr(i,j)
      END DO
   END DO
! projection Z
   rr=Zero
   IF (inquire_key_presence( DATA_FILE, '$edmom_z')) CALL read_2d_real_array( DATA_FILE, '$edmom_z',n,n,rr)
   IF (dbg) WRITE (StdOut,*) 'read_edipmom::  norm of moment_z=',dnrm2_(n*n,rr,1)
   DO i=1,n
      DO j=1,n
         moment(3,i,j) = rr(i,j)
      END DO
   END DO
   DEALLOCATE (rr)
   IF (dnrm2_(3*n*n,moment,1) <=TINY(Zero)) CALL WarningMessage(1,'read_edipmom:: the norm of the read moment is zero!')
   RETURN
END SUBROUTINE read_edipmom

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_amfi (DATA_FILE, n, moment)
   IMPLICIT NONE
   INTEGER, INTENT (IN)    :: DATA_FILE
   INTEGER, INTENT (IN)    :: N
   REAL (wp), INTENT (OUT) :: moment(3,N,N)
   REAL (wp), ALLOCATABLE  :: rr(:,:)
   INTEGER                 :: i, j
   REAL (wp), EXTERNAL     :: dnrm2_

   moment=Zero
   ALLOCATE (rr(n,n))
! projection X
   rr=Zero
   IF (inquire_key_presence( DATA_FILE, '$amfi_x')) CALL read_2d_real_array( DATA_FILE, '$amfi_x',n,n,rr)
   IF (dbg) WRITE (StdOut,*) 'read_amfi::  norm of moment_x=',dnrm2_(n*n,rr,1)
   DO i=1,n
      DO j=1,n
         moment(1,i,j) = rr(i,j)
      END DO
   END DO
! projection Y
   rr=Zero
   IF (inquire_key_presence( DATA_FILE, '$amfi_y')) CALL read_2d_real_array( DATA_FILE, '$amfi_y',n,n,rr)
   IF (dbg) WRITE (StdOut,*) 'read_amfi::  norm of moment_y=',dnrm2_(n*n,rr,1)
   DO i=1,n
      DO j=1,n
         moment(2,i,j) = rr(i,j)
      END DO
   END DO
! projection Z
   rr=Zero
   IF (inquire_key_presence( DATA_FILE, '$amfi_z')) CALL read_2d_real_array( DATA_FILE, '$amfi_z',n,n,rr)
   IF (dbg) WRITE (StdOut,*) 'read_amfi::  norm of moment_z=',dnrm2_(n*n,rr,1)
   DO i=1,n
      DO j=1,n
         moment(3,i,j) = rr(i,j)
      END DO
   END DO
   DEALLOCATE (rr)
   IF (dnrm2_(3*n*n,moment,1) <=TINY(Zero)) CALL WarningMessage(1,'read_amfi:: the norm of the read moment is zero!')
   RETURN
END SUBROUTINE read_amfi

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_nss (DATA_FILE, n)
   IMPLICIT NONE
   INTEGER, INTENT (IN)   :: DATA_FILE
   INTEGER, INTENT (OUT)  :: n
   n=0
   IF (inquire_key_presence( DATA_FILE, '$nss')) CALL read_INTEGER_scalar ( DATA_FILE, '$nss', n )
   IF ( n<=0 ) THEN
      CALL WarningMessage(1,'read_nss:: nss value in DATA_FILE = 0. Is it really the case?')
   ENDIF
   RETURN
END SUBROUTINE read_nss

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_nstate (DATA_FILE, n)
   IMPLICIT NONE
   INTEGER, INTENT (IN)  :: DATA_FILE
   INTEGER, INTENT (OUT) :: n
   n=0
   IF (inquire_key_presence( DATA_FILE, '$nstate')) CALL read_INTEGER_scalar ( DATA_FILE, '$nstate', n )
   IF ( n<=0 ) THEN
      CALL WarningMessage(1,'read_nstate:: nstate value in DATA_FILE = 0. Is it really the case?')
   ENDIF
   RETURN
END SUBROUTINE read_nstate

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_nmult (DATA_FILE, n)
   IMPLICIT NONE
   INTEGER, INTENT (IN)   :: DATA_FILE
   INTEGER, INTENT (OUT)  :: n
   n=0
   IF (inquire_key_presence( DATA_FILE, '$nmult')) CALL read_INTEGER_scalar ( DATA_FILE, '$nmult', n )
   IF ( n<=0 ) THEN
      CALL WarningMessage(1,'read_nmult:: nmult value in DATA_FILE = 0. Is it really the case?')
   ENDIF
   RETURN
END SUBROUTINE read_nmult

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_multiplicity (DATA_FILE, n, array)
   IMPLICIT NONE
   INTEGER, INTENT (IN)  :: DATA_FILE
   INTEGER, INTENT (IN)  :: n
   INTEGER, INTENT (OUT) :: array(n)
   array=0
   IF (inquire_key_presence( DATA_FILE, '$multiplicity')) CALL read_1d_INTEGER_array( DATA_FILE, '$multiplicity', n, array )
   IF ( SUM(ABS(array(1:n)))==0 ) THEN
      CALL WarningMessage(1,'read_multiplicity:: it seems that all the multiplicities in DATA_FILE are 0. '// &
                            'Is it really the case?')
      WRITE (StdOut, *) 'read_multiplicity:: SUM(Sz) = ',SUM(ABS(array(1:n)))
   ENDIF
   IF ( SUM(array(1:n))==0 ) THEN
      CALL WarningMessage(1,'read_multiplicity:: it seems that all the multiplicities in DATA_FILE are 0. '// &
                            'Is it really the case?')
      WRITE (StdOut, *) 'read_szproj:: SUM(Sz) = ',SUM(array(1:n))
   ENDIF
   RETURN
END SUBROUTINE read_multiplicity

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_imult (DATA_FILE, n, array)
   IMPLICIT NONE
   INTEGER, INTENT (IN)  :: DATA_FILE
   INTEGER, INTENT (IN)  :: n
   INTEGER, INTENT (OUT) :: array(n)
   array=0
   IF (inquire_key_presence( DATA_FILE, '$imult')) CALL read_1d_INTEGER_array( DATA_FILE, '$imult', n, array )
   IF ( SUM(array(1:n))==0 ) THEN
      CALL WarningMessage(1,'read_imult:: it seems that all the multiplicities in DATA_FILE are 0. '// &
                            'Is it really the case?')
      WRITE (StdOut, *) 'read_imult:: SUM(mult()) = ',SUM(array(1:n))
   ENDIF
   RETURN
END SUBROUTINE read_imult

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_format (DATA_FILE, n)
   IMPLICIT NONE
   INTEGER, INTENT (IN)   :: DATA_FILE
   INTEGER, INTENT (OUT)  :: n
   n=0
   IF (inquire_key_presence( DATA_FILE, '$format')) CALL read_INTEGER_scalar ( DATA_FILE, '$format', n )
   IF ( n<=0 ) THEN
      CALL WarningMessage(2,'read_format:: FORMAT value in DATA_FILE = 0.'// &
                            ' The FORMAT must be equal or larger than 2020. Please check.' )
   ENDIF
   RETURN
END SUBROUTINE read_format

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_nroot (DATA_FILE, n, array)
   IMPLICIT NONE
   INTEGER, INTENT (IN)  :: DATA_FILE
   INTEGER, INTENT (IN)  :: n
   INTEGER, INTENT (OUT) :: array(n)
   array=0
   IF (inquire_key_presence( DATA_FILE, '$nroot')) CALL read_1d_INTEGER_array( DATA_FILE, '$nroot', n, array )
   IF ( SUM(array(1:n))==0 ) THEN
      CALL WarningMessage(1,'read_nroot:: it seems that the number of roots included in spin-orbit '// &
                            'interaction in DATA_FILE are 0.  Is it really the case?')
      WRITE (StdOut, *) 'read_szproj:: SUM(array()) = ',SUM(array(1:n))
   ENDIF
   RETURN
END SUBROUTINE read_nroot

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_szproj (DATA_FILE, n, array)
   IMPLICIT NONE
   INTEGER, INTENT (IN)  :: DATA_FILE
   INTEGER, INTENT (IN)  :: n
   INTEGER, INTENT (OUT) :: array(n)
   array=0
   IF (inquire_key_presence( DATA_FILE, '$szproj')) CALL read_1d_INTEGER_array( DATA_FILE, '$szproj', n, array )
   IF ( SUM(ABS(array(1:n))) == 0 ) THEN
      CALL WarningMessage(1,'read_szproj:: it seems that SUM(ABS(Sz)) in DATA_FILE is 0. '// &
                            'Is it really the case?')
      WRITE (StdOut, *) 'read_szproj:: SUM(ABS(Sz)) = ',SUM(ABS(array(1:n)))
   ENDIF
   IF ( SUM(array(1:n)) /= 0 ) THEN
      CALL WarningMessage(1,'read_szproj:: it seems that SUM(Sz) in DATA_FILE is not 0. '// &
                            'Is it really the case?')
      WRITE (StdOut, *) 'read_szproj:: SUM(Sz) = ',SUM(array(1:n))
   ENDIF
   RETURN
END SUBROUTINE read_szproj

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_eso (DATA_FILE, n, array)
   IMPLICIT NONE
   INTEGER, INTENT (IN)    :: DATA_FILE
   INTEGER, INTENT (IN)    :: n
   REAL (wp), INTENT (OUT) :: array(n)
   REAL (wp), EXTERNAL     :: dnrm2_
   array=Zero
   IF (inquire_key_presence( DATA_FILE, '$eso')) CALL read_1d_real_array( DATA_FILE, '$eso', n, array )
   IF (DBG) WRITE (StdOut,*) 'read_eso::  norm of eso=',dnrm2_(n,array,1)
   IF ( dnrm2_(n,array,1) <= MINIMAL_REAL ) THEN
      CALL WarningMessage(1,'read_eso:: it seems that the norm of ESO array in DATA_FILE is 0. '// &
                            'Is it really the case?')
      WRITE (StdOut, *) 'read_eso:: dnrm2_(eso) = ',dnrm2_(n,array,1)
   ENDIF
   RETURN
END SUBROUTINE read_eso

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_esfs (DATA_FILE, n, array)
   IMPLICIT NONE
   INTEGER, INTENT (IN)    :: DATA_FILE
   INTEGER, INTENT (IN)    :: n
   REAL (wp), INTENT (OUT) :: array(n)
   REAL (wp), EXTERNAL     :: dnrm2_
   array=Zero
   IF (inquire_key_presence( DATA_FILE, '$esfs')) CALL read_1d_real_array( DATA_FILE, '$esfs', n, array )
   IF (DBG) WRITE (StdOut,*) 'read_esfs::  norm of esfs=',dnrm2_(n,array,1)
   IF ( dnrm2_(n,array,1) <= MINIMAL_REAL ) THEN
      CALL WarningMessage(1,'read_esfs:: it seems that the norm of ESFS in DATA_FILE is 0. '// &
                            'Is it really the case?')
      WRITE (StdOut, *) 'read_esfs:: dnrm2_(esfs) = ',dnrm2_(n,array,1)
   ENDIF
   RETURN
END SUBROUTINE read_esfs

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_hso (DATA_FILE, n, array)
   IMPLICIT NONE
   INTEGER, INTENT (IN)       :: DATA_FILE
   INTEGER, INTENT (IN)       :: n
   COMPLEX (wp), INTENT (OUT) :: array(n)
   REAL (wp), EXTERNAL        :: dznrm2_
   array=ZeroC
   IF (inquire_key_presence( DATA_FILE, '$hso')) CALL read_complex_matrix( DATA_FILE, '$hso', n, array )
   IF (DBG) WRITE (StdOut,*) 'read_hso::  norm of hso=',dznrm2_(n*n,array,1)
   IF ( dznrm2_(n*n,array,1) <= MINIMAL_REAL ) THEN
      CALL WarningMessage(1,'read_hso:: it seems that norm of HSO in DATA_FILE is 0. '// &
                            'Is it really the case?')
      WRITE (StdOut, *) 'read_hso:: dznrm2_(hso) = ',dznrm2_(n*n,array,1)
   END IF
   RETURN
END SUBROUTINE read_hso

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_eigen (DATA_FILE, n, array)
   IMPLICIT NONE
   INTEGER, INTENT (IN)       :: DATA_FILE
   INTEGER, INTENT (IN)       :: n
   COMPLEX (wp), INTENT (OUT) :: array(n)
   REAL (wp), EXTERNAL        :: dznrm2_
   array=ZeroC
   IF (inquire_key_presence( DATA_FILE, '$eigen')) CALL read_complex_matrix( DATA_FILE, '$eigen', n, array )
   IF (DBG) WRITE (StdOut,*) 'read_eigen::  norm of eigenv=',dznrm2_(n*n,array,1)
   IF ( dznrm2_(n*n,array,1) <= MINIMAL_REAL ) THEN
      CALL WarningMessage(1,'read_eigen:: it seems that norm of EIGENV in DATA_FILE is 0. '// &
                            'Is it really the case?')
      WRITE (StdOut, *) 'read_eigen:: dznrm2_(array) = ',dznrm2_(n*n,array,1)
   END IF
   RETURN
END SUBROUTINE read_eigen

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_gtens (DATA_FILE, nmult, gtens, axes)
   IMPLICIT NONE
   INTEGER, INTENT (IN)    :: DATA_FILE
   INTEGER, INTENT (IN)    :: nmult
   REAL (wp), INTENT (OUT) :: gtens(nmult,3)
   REAL (wp), INTENT (OUT) ::  axes(nmult,3,3)
   gtens=Zero ; axes=Zero
   CALL read_2d_real_array( DATA_FILE, '$gtens_main', nmult, 3, gtens )
   CALL read_3d_real_array( DATA_FILE, '$gtens_axes', nmult, 3, 3, axes )
   RETURN
END SUBROUTINE read_gtens

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_stev_cfp (DATA_FILE, s, n, cfp)
   IMPLICIT NONE
   INTEGER, INTENT (IN)    :: DATA_FILE
   INTEGER, INTENT (IN)    :: n   ! 2J+1, or 2L+1, i.e. the dimension of the J or L multiplet
   REAL (wp), INTENT (OUT) :: cfp(n-1, -(n-1):(n-1) )
   CHARACTER (LEN=1)       :: s
   INTEGER                 :: k, q, i, ik,iq

   ierr=0
   IF ( (n<=0) ) THEN
      CALL WarningMessage(1,'read_stev_cfp_'//trim(s)//':: nothing to read. Array size = 0.')
      RETURN
   END IF
   CFP=Zero

   REWIND ( DATA_FILE )
   CALL file_advance_to_string( DATA_FILE, '$stev_cfp_'//trim(s), line )

   READ ( DATA_FILE, FMT=*, IOSTAT=ierr ) i
   IF ( i /= n ) CALL WarningMessage(2,'read_stev_cfp_'//trim(s)//':: size of the multiplet is not the same i/=n')

   ! if key is found, THEN read the data
   IF ( ierr == 0 ) THEN
      DO k=2,n-1,2
         DO q=-k,k,2
            ik=99999;
            iq=9999999;
            READ ( DATA_FILE, FMT=*, IOSTAT=ierr ) ik,iq, cfp(ik,iq)
            IF (ierr /= 0) THEN
               CALL WarningMessage(2,'read_stev_cfp_'//trim(s)//':: Something went wrong reading the array.')
            END IF
            IF (DBG) WRITE (StdOut,*) 'read_stev_cfp_'//trim(s)//'::  k, q =',k, q
            IF (DBG) FLUSH (StdOut)
         END DO
      END DO
   END IF
   RETURN
END SUBROUTINE read_stev_cfp

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_susc (DATA_FILE, s, n, field, zj, t, x, x_tens)
   IMPLICIT NONE
   INTEGER, INTENT (IN)           :: DATA_FILE
   INTEGER, INTENT (INOUT)        :: n                  ! number of temperature points
   REAL (wp), INTENT (OUT)        :: zj                 ! intermolecular interaction
   REAL (wp), INTENT (OUT)        :: field              ! applied field
   REAL (wp), INTENT (OUT)        :: t(n)               ! temperature points
   REAL (wp), INTENT (OUT)        :: x(n)               ! susceptibility X
   REAL (wp), INTENT (OUT)        :: x_tens(n,3,3)      ! susceptibility tensor, X_tens
   CHARACTER (LEN=*), INTENT (IN) :: s
   INTEGER                        :: i, j, k
   REAL (wp), EXTERNAL            :: dnrm2_

   t=Zero ; x=Zero ; x_tens=Zero
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
!   IF ( (trim(s).ne.'x')                      .OR.   &
!       (trim(s).ne.'x_field')                .OR.   &
!       (trim(s).ne.'x_minus_one')            .OR.   &
!       (trim(s).ne.'x_field_minus_one')      .OR.   &
!       (trim(s).ne.'xt')                     .OR.   &
!       (trim(s).ne.'xt_field')               .OR.   &
!       (trim(s).ne.'xt_minus_one')           .OR.   &
!       (trim(s).ne.'xt_field_minus_one')     .OR.   &
!       (trim(s).ne.'x_zj')                   .OR.   &
!       (trim(s).ne.'x_field_zj')             .OR.   &
!       (trim(s).ne.'x_minus_one_zj')         .OR.   &
!       (trim(s).ne.'x_field_minus_one_zj')   .OR.   &
!       (trim(s).ne.'xt_zj')                  .OR.   &
!       (trim(s).ne.'xt_field_zj')            .OR.   &
!       (trim(s).ne.'xt_minus_one_zj')        .OR.   &
!       (trim(s).ne.'xt_field_minus_one_zj') ) THEN
!
!       CALL WarningMessage(1,'read_x '//trim(s)//' :: the parameter s='//trim(s)//'is not understood. RETURN without'// &
!                             ' read of X data.')
!       RETURN
!   END IF

   ierr=0
   IF ( (n<=0) ) THEN
      CALL WarningMessage(1,'read_x '//trim(s)//' :: nothing to read. Array size = 0.')
      RETURN
   END IF

   REWIND ( DATA_FILE )
   CALL file_advance_to_string( DATA_FILE, '$susceptibility_'//trim(s), line )

   ! if key is found, THEN read the data
   IF ( ierr == 0 ) THEN
      READ ( DATA_FILE, FMT=*, IOSTAT=ierr ) n
      READ ( DATA_FILE, FMT=*, IOSTAT=ierr ) zj, field
      IF (ierr /= 0) THEN
         CALL WarningMessage(2,'read_x '//trim(s)//' :: Something went wrong reading the zJ and field values.')
      END IF
      READ ( DATA_FILE, FMT=*, IOSTAT=ierr ) (T(i),i=1,n)
      IF (dnrm2_(n,T,1)<TINY(Zero)) CALL WarningMessage(1,'read_x '//trim(s)//' :: all array T elements are zero = 0.')
      IF (ierr /= 0) THEN
         CALL WarningMessage(2,'read_x '//trim(s)//' :: Something went wrong reading the T array.')
      END IF
      READ ( DATA_FILE, FMT=*, IOSTAT=ierr ) (X(i),i=1,n)
      IF (dnrm2_(n,X,1)<TINY(Zero)) CALL WarningMessage(1,'read_x '//trim(s)//' :: all array X elements are zero = 0.')
      IF (ierr /= 0) THEN
         CALL WarningMessage(2,'read_x '//trim(s)//' :: Something went wrong reading the X array.')
      END IF
      DO j=1,3
         DO k=1,3
            READ ( DATA_FILE, FMT=*, IOSTAT=ierr ) ( X_tens(i,j,k), i=1,n )
            IF (ierr /= 0) THEN
               CALL WarningMessage(2,'read_x '//trim(s)//' :: Something went wrong reading the X_tens array.')
            END IF
         END DO
      END DO
      IF (dnrm2_(n*3*3,X_tens,1)<TINY(Zero)) CALL WarningMessage(1,'read_x '//trim(s)//' :: all array X_tens elements are zero.')
      IF (DBG) FLUSH (StdOut)
   ELSE
      WRITE (StdOut,*) 'keyword $susceptibility_'//trim(s)//' was not found in DATA_FILE'
   END IF
   RETURN
END SUBROUTINE read_susc

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_magn (DATA_FILE, nt, nh, nd, nss, zj, t, h, x, y, z, w, m, mav, energy)
   IMPLICIT NONE
   INTEGER, INTENT (IN)    :: DATA_FILE
   INTEGER, INTENT (IN)    :: nt                         ! number of temperature points
   INTEGER, INTENT (IN)    :: nh                         ! number of field points
   INTEGER, INTENT (IN)    :: nd                         ! number of directions of aplied field
   INTEGER, INTENT (IN)    :: nss                        ! number of spin-orbit states
   REAL (wp), INTENT (OUT) :: zj                         ! inter-molecular parameter zJ
   REAL (wp), INTENT (OUT) :: t(nt)                      ! temperature points
   REAL (wp), INTENT (OUT) :: h(nh)                      ! field points
   REAL (wp), INTENT (OUT) :: x(nd), y(nd), z(nd), w(nd) ! Lebedev grid directions and weight
   REAL (wp), INTENT (OUT) :: m(nd, 3, nt, nh)           ! magnetisation vector
   REAL (wp), INTENT (OUT) :: mav(nt, nh)                ! average magnetisation vector
   REAL (wp), INTENT (OUT) :: energy(nd, nh, nss)        ! Zeeman energy states
   INTEGER                 :: id, ih, it, i, iss, l
   REAL (wp), EXTERNAL     :: dnrm2_

   ierr=0
   IF ( (nt<=0) .OR. (nh<=0) .OR. (nd<=0) .OR. (nss<=0) ) THEN
      CALL WarningMessage(1,'read_magn :: nothing to read. Array size = 0.')
      RETURN
   END IF

   REWIND ( DATA_FILE )
   CALL file_advance_to_string( DATA_FILE, '$magnetisation', line )

   READ ( DATA_FILE, FMT=*, IOSTAT=ierr ) it, ih, id, iss
   IF (it/=nt) CALL WarningMessage(1,'read_magn :: nt read from DATA_FILE is not the same as the parameter used to CALL '//&
                                    'this function.')
   IF (ih/=nh) CALL WarningMessage(1,'read_magn :: nh read from DATA_FILE is not the same as the parameter used to CALL '//&
                                    'this function.')
   IF (id/=nd) CALL WarningMessage(1,'read_magn :: nd read from DATA_FILE is not the same as the parameter used to CALL '//&
                                    'this function.')
   IF (iss/=nss) CALL WarningMessage(1,'read_magn :: nss read from DATA_FILE is not the same as the parameter used to CALL '//&
                                    'this function.')
   zJ=Zero ; T=Zero ; H=Zero ; X=Zero ; Y=Zero ; Z=Zero ; W=Zero ;
   energy=Zero ; m=Zero ; mav=Zero ;

   READ ( DATA_FILE, FMT=*, IOSTAT=ierr ) zJ
   READ ( DATA_FILE, FMT=*, IOSTAT=ierr ) (T(i),i=1,nt)
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_magn :: Something went wrong reading the T array.')
   END IF
   ! Magnetic field points
   READ ( DATA_FILE, FMT=*, IOSTAT=ierr ) (H(i),i=1,nh)
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_magn :: Something went wrong reading the H array.')
   END IF
   ! Lebedev grid data (X, Y, Z, W)
   READ ( DATA_FILE, FMT=*, IOSTAT=ierr ) ( X(i), i=1,nd )
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_magn :: Something went wrong reading the Lebedev grid X rray.')
   END IF
   READ ( DATA_FILE, FMT=*, IOSTAT=ierr ) ( Y(i), i=1,nd )
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_magn :: Something went wrong reading the Lebedev grid Y array.')
   END IF
   READ ( DATA_FILE, FMT=*, IOSTAT=ierr ) ( Z(i), i=1,nd )
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_magn :: Something went wrong reading the Lebedev grid Z array.')
   END IF
   READ ( DATA_FILE, FMT=*, IOSTAT=ierr ) ( W(i), i=1,nd )
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_magn :: Something went wrong reading the Lebedev grid W array.')
   END IF
   ! Zeeman energy
   DO id=1,nd
      DO ih=1,nh
         READ ( DATA_FILE, FMT=*, IOSTAT=ierr )  (energy(ih,id,iss),iss=1,nss)
         IF (ierr /= 0) THEN
            CALL WarningMessage(2,'read_magn :: Something went wrong reading the Zeeman energy data.')
         END IF
      END DO
   END DO
   ! magnetisation vector data
   DO id=1,nd
      DO l=1,3
         DO ih=1,nh
            READ ( DATA_FILE, FMT=*, IOSTAT=ierr ) (m(id,l,ih,it),it=1,nt)
            IF (ierr /= 0) THEN
               CALL WarningMessage(2,'read_magn :: Something went wrong reading the M data.')
            END IF
         END DO
      END DO
   END DO
   DO ih=1,nh
      READ ( DATA_FILE, FMT=*, IOSTAT=ierr ) (mav(ih,it),it=1,nt)
      IF (ierr /= 0) THEN
         CALL WarningMessage(2,'read_magn :: Something went wrong reading the average M data.')
      END IF
   END DO
   IF (DBG) FLUSH (StdOut)
   RETURN
END SUBROUTINE read_magn











!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_complex_matrix ( LU,key, n, matrix)
   IMPLICIT NONE
   INTEGER, INTENT (IN)           :: LU
   INTEGER, INTENT (IN)           :: N
   CHARACTER (LEN=*), INTENT (IN) :: key
   COMPLEX (wp), INTENT (OUT)     :: matrix(N,N)
   REAL (wp), ALLOCATABLE        :: rr(:,:), ri(:,:)
   INTEGER                       :: i, j

   matrix=ZeroC
   ALLOCATE (rr(n,n))
   ALLOCATE (ri(n,n))
   rr=Zero ; ri=Zero ;
   CALL read_2d_real_array( LU,key//'r',n,n,rr)
   CALL read_2d_real_array( LU,key//'i',n,n,ri)
   DO i=1,n
      DO j=1,n
         matrix(i,j) = CMPLX( rr(i,j), ri(i,j), wp)
      END DO
   END DO
   DEALLOCATE (rr)
   DEALLOCATE (ri)
   RETURN
END SUBROUTINE read_complex_matrix

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_complex_matrix ( LU, key, n, matrix)
   IMPLICIT NONE
   INTEGER, INTENT (IN)           :: LU
   INTEGER, INTENT (IN)           :: N
   CHARACTER (LEN=*), INTENT (IN) :: key
   COMPLEX (wp), INTENT (IN)      :: matrix(N,N)
   REAL (wp), ALLOCATABLE        :: rr(:,:), ri(:,:)
   INTEGER                       :: i, j

   ALLOCATE (rr(n,n))
   ALLOCATE (ri(n,n))
   rr=Zero ; ri=Zero ;
   DO i=1,n
      DO j=1,n
         rr(i,j)= DBLE(matrix(i,j))
         ri(i,j)=DIMAG(matrix(i,j))
      END DO
   END DO
   CALL write_2d_real_array( LU,key//'r',n,n,rr)
   CALL write_2d_real_array( LU,key//'i',n,n,ri)
   DEALLOCATE (rr)
   DEALLOCATE (ri)
   RETURN
END SUBROUTINE write_complex_matrix

!--------------------------------------------------------------------------------------------------!















!--------------------------------------------------------------------------------------------------!
!           HIGH LEVEL WRITING SUBROUTINES
!--------------------------------------------------------------------------------------------------!
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
!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_magnetic_moment (ANISO_FILE, n, moment)
   IMPLICIT NONE
   INTEGER, INTENT (IN)      :: ANISO_FILE
   INTEGER, INTENT (IN)      :: N
   COMPLEX (wp), INTENT (IN) :: moment(3,N,N)
   REAL (wp), ALLOCATABLE    :: rr(:,:), ri(:,:)
   INTEGER                   :: i, j

   ALLOCATE (rr(n,n))
   ALLOCATE (ri(n,n))
! projection X
   rr=Zero; ri=Zero;
   DO i=1,n
      DO j=1,n
         rr(i,j)= DBLE(moment(1,i,j))
         ri(i,j)=DIMAG(moment(1,i,j))
      END DO
   END DO
   CALL write_2d_real_array( ANISO_FILE, '$magn_xr',n,n,rr)
   CALL write_2d_real_array( ANISO_FILE, '$magn_xi',n,n,ri)
! projection Y
   rr=Zero; ri=Zero;
   DO i=1,n
      DO j=1,n
         rr(i,j)= DBLE(moment(2,i,j))
         ri(i,j)=DIMAG(moment(2,i,j))
      END DO
   END DO
   CALL write_2d_real_array( ANISO_FILE, '$magn_yr',n,n,rr)
   CALL write_2d_real_array( ANISO_FILE, '$magn_yi',n,n,ri)
! projection Z
   rr=Zero; ri=Zero;
   DO i=1,n
      DO j=1,n
         rr(i,j)= DBLE(moment(3,i,j))
         ri(i,j)=DIMAG(moment(3,i,j))
      END DO
   END DO
   CALL write_2d_real_array( ANISO_FILE, '$magn_zr',n,n,rr)
   CALL write_2d_real_array( ANISO_FILE, '$magn_zi',n,n,ri)
   DEALLOCATE (rr)
   DEALLOCATE (ri)
   RETURN
END SUBROUTINE write_magnetic_moment

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_electric_moment (ANISO_FILE, n, moment)
   IMPLICIT NONE
   INTEGER, INTENT (IN)       :: ANISO_FILE
   INTEGER, INTENT (IN)       :: N
   COMPLEX (wp), INTENT (OUT) :: moment(3,N,N)
   REAL (wp), ALLOCATABLE     :: rr(:,:), ri(:,:)
   INTEGER                    :: i, j

   ALLOCATE (rr(n,n))
   ALLOCATE (ri(n,n))
! projection X
   rr=Zero; ri=Zero;
   DO i=1,n
      DO j=1,n
         rr(i,j)= DBLE(moment(1,i,j))
         ri(i,j)=DIMAG(moment(1,i,j))
      END DO
   END DO
   CALL write_2d_real_array( ANISO_FILE, '$edipm_xr',n,n,rr)
   CALL write_2d_real_array( ANISO_FILE, '$edipm_xi',n,n,ri)
! projection Y
   rr=Zero; ri=Zero;
   DO i=1,n
      DO j=1,n
         rr(i,j)= DBLE(moment(2,i,j))
         ri(i,j)=DIMAG(moment(2,i,j))
      END DO
   END DO
   CALL write_2d_real_array( ANISO_FILE, '$edipm_yr',n,n,rr)
   CALL write_2d_real_array( ANISO_FILE, '$edipm_yi',n,n,ri)
! projection Z
   rr=Zero; ri=Zero;
   DO i=1,n
      DO j=1,n
         rr(i,j)=DBLE(moment(3,i,j))
         ri(i,j)=DIMAG(moment(3,i,j))
      END DO
   END DO
   CALL write_2d_real_array( ANISO_FILE, '$edipm_zr',n,n,rr)
   CALL write_2d_real_array( ANISO_FILE, '$edipm_zi',n,n,ri)
   DEALLOCATE (rr)
   DEALLOCATE (ri)

   RETURN
END SUBROUTINE write_electric_moment

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_spin_moment (ANISO_FILE, n, moment)
   IMPLICIT NONE
   INTEGER, INTENT (IN)       :: ANISO_FILE
   INTEGER, INTENT (IN)       :: N
   COMPLEX (wp), INTENT (IN)  :: moment(3,N,N)
   REAL (wp), ALLOCATABLE     :: rr(:,:), ri(:,:)
   INTEGER                    :: i, j

   ALLOCATE (rr(n,n))
   ALLOCATE (ri(n,n))
! projection X
   rr=Zero; ri=Zero;
   DO i=1,n
      DO j=1,n
         rr(i,j)= DBLE(moment(1,i,j))
         ri(i,j)=DIMAG(moment(1,i,j))
      END DO
   END DO
   CALL write_2d_real_array( ANISO_FILE, '$spin_xr',n,n,rr)
   CALL write_2d_real_array( ANISO_FILE, '$spin_xi',n,n,ri)
! projection Y
   rr=Zero; ri=Zero;
   DO i=1,n
      DO j=1,n
         rr(i,j)= DBLE(moment(2,i,j))
         ri(i,j)=DIMAG(moment(2,i,j))
      END DO
   END DO
   CALL write_2d_real_array( ANISO_FILE, '$spin_yr',n,n,rr)
   CALL write_2d_real_array( ANISO_FILE, '$spin_yi',n,n,ri)
! projection Z
   rr=Zero; ri=Zero;
   DO i=1,n
      DO j=1,n
         rr(i,j)= DBLE(moment(3,i,j))
         ri(i,j)=DIMAG(moment(3,i,j))
      END DO
   END DO
   CALL write_2d_real_array( ANISO_FILE, '$spin_zr',n,n,rr)
   CALL write_2d_real_array( ANISO_FILE, '$spin_zi',n,n,ri)
   DEALLOCATE (rr)
   DEALLOCATE (ri)
   RETURN
END SUBROUTINE write_spin_moment

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_angmom (ANISO_FILE, n, moment)
   IMPLICIT NONE
   INTEGER, INTENT (IN)     :: ANISO_FILE
   INTEGER, INTENT (IN)     :: N
   REAL (wp), INTENT (IN)   :: moment(3,N,N)
   CALL write_2d_real_array( ANISO_FILE, '$angmom_x',n,n,moment(1,:,:) )
   CALL write_2d_real_array( ANISO_FILE, '$angmom_y',n,n,moment(2,:,:) )
   CALL write_2d_real_array( ANISO_FILE, '$angmom_z',n,n,moment(3,:,:) )
   RETURN
END SUBROUTINE write_angmom

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_edipmom (ANISO_FILE, n, moment)
   IMPLICIT NONE
   INTEGER, INTENT (IN)     :: ANISO_FILE
   INTEGER, INTENT (IN)     :: N
   REAL (wp), INTENT (IN)   :: moment(3,N,N)

   CALL write_2d_real_array( ANISO_FILE, '$edmom_x',n,n,moment(1,:,:) )
   CALL write_2d_real_array( ANISO_FILE, '$edmom_y',n,n,moment(2,:,:) )
   CALL write_2d_real_array( ANISO_FILE, '$edmom_z',n,n,moment(3,:,:) )
   RETURN
END SUBROUTINE write_edipmom

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_amfi (ANISO_FILE, n, moment)
   IMPLICIT NONE
   INTEGER, INTENT (IN)   :: ANISO_FILE
   INTEGER, INTENT (IN)   :: N
   REAL (wp), INTENT (IN) :: moment(3,N,N)

   CALL write_2d_real_array( ANISO_FILE, '$amfi_x',n,n,moment(1,:,:) )
   CALL write_2d_real_array( ANISO_FILE, '$amfi_y',n,n,moment(2,:,:) )
   CALL write_2d_real_array( ANISO_FILE, '$amfi_z',n,n,moment(3,:,:) )
   RETURN
END SUBROUTINE write_amfi

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_format(ANISO_FILE, n)
   IMPLICIT NONE
   INTEGER, INTENT (IN) :: ANISO_FILE
   INTEGER, INTENT (IN) :: n
   CALL write_INTEGER_scalar ( ANISO_FILE, '$format', n )
   RETURN
END SUBROUTINE write_format

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_nss (ANISO_FILE, n)
   IMPLICIT NONE
   INTEGER, INTENT (IN) :: ANISO_FILE
   INTEGER, INTENT (IN) :: n
   CALL write_INTEGER_scalar ( ANISO_FILE, '$nss', n )
   RETURN
END SUBROUTINE write_nss

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_nstate (ANISO_FILE, n)
   IMPLICIT NONE
   INTEGER, INTENT (IN) :: ANISO_FILE
   INTEGER, INTENT (IN) :: n
   CALL write_INTEGER_scalar ( ANISO_FILE, '$nstate', n )
   RETURN
END SUBROUTINE write_nstate

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_nmult (ANISO_FILE, n)
   IMPLICIT NONE
   INTEGER, INTENT (IN) :: ANISO_FILE
   INTEGER, INTENT (IN) :: n
   CALL write_INTEGER_scalar ( ANISO_FILE, '$nmult', n )
   RETURN
END SUBROUTINE write_nmult

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_multiplicity (ANISO_FILE, n, array)
   IMPLICIT NONE
   INTEGER, INTENT (IN) :: ANISO_FILE
   INTEGER, INTENT (IN) :: n
   INTEGER, INTENT (IN) :: array(n)
   CALL write_1d_INTEGER_array( ANISO_FILE, '$multiplicity', n, array )
   RETURN
END SUBROUTINE write_multiplicity

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_imult (ANISO_FILE, n, array)
   IMPLICIT NONE
   INTEGER, INTENT (IN) :: ANISO_FILE
   INTEGER, INTENT (IN) :: n
   INTEGER, INTENT (IN) :: array(n)
   CALL write_1d_INTEGER_array( ANISO_FILE, '$imult', n, array )
   RETURN
END SUBROUTINE write_imult

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_nroot (ANISO_FILE, n, array)
   IMPLICIT NONE
   INTEGER, INTENT (IN) :: ANISO_FILE
   INTEGER, INTENT (IN) :: n
   INTEGER, INTENT (IN) :: array(n)
   CALL write_1d_INTEGER_array( ANISO_FILE, '$nroot', n, array )
   RETURN
END SUBROUTINE write_nroot
!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_szproj (ANISO_FILE, n, array)
   IMPLICIT NONE
   INTEGER, INTENT (IN) :: ANISO_FILE
   INTEGER, INTENT (IN) :: n
   INTEGER, INTENT (IN) :: array(n)
   CALL write_1d_INTEGER_array( ANISO_FILE, '$szproj', n, array )
   RETURN
END SUBROUTINE write_szproj

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_eso (ANISO_FILE, n, array)
   IMPLICIT NONE
   INTEGER, INTENT (IN)   :: ANISO_FILE
   INTEGER, INTENT (IN)   :: n
   REAL (wp), INTENT (IN) :: array(n)
   INTEGER :: i
   IF (DBG) WRITE (StdOut,*) 'write_eso: '
   CALL write_1d_real_array( ANISO_FILE, '$eso', n, array )
   RETURN
END SUBROUTINE write_eso

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_esfs (ANISO_FILE, n, array)
   IMPLICIT NONE
   INTEGER, INTENT (IN)   :: ANISO_FILE
   INTEGER, INTENT (IN)   :: n
   REAL (wp), INTENT (IN) :: array(n)
   CALL write_1d_real_array( ANISO_FILE, '$esfs', n, array )
   RETURN
END SUBROUTINE write_esfs

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_hso (ANISO_FILE, n, array)
   IMPLICIT NONE
   INTEGER, INTENT (IN)      :: ANISO_FILE
   INTEGER, INTENT (IN)      :: n
   COMPLEX (wp), INTENT (IN) :: array(n)
   CALL write_complex_matrix( ANISO_FILE, '$hso', n, array )
   RETURN
END SUBROUTINE write_hso

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_eigen (ANISO_FILE, n, array)
   IMPLICIT NONE
   INTEGER, INTENT (IN)      :: ANISO_FILE
   INTEGER, INTENT (IN)      :: n
   COMPLEX (wp), INTENT (IN) :: array(n)
   CALL write_complex_matrix( ANISO_FILE, '$eigen', n, array )
   RETURN
END SUBROUTINE write_eigen

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_gtens (ANISO_FILE, nmult, gtens, axes)
   IMPLICIT NONE
   INTEGER, INTENT (IN)   :: ANISO_FILE
   INTEGER, INTENT (IN)   :: nmult
   REAL (wp), INTENT (IN) :: gtens(nmult,3)
   REAL (wp), INTENT (IN) ::  axes(nmult,3,3)
   CALL write_2d_real_array( ANISO_FILE, '$gtens_main', nmult, 3, gtens )
   CALL write_3d_real_array( ANISO_FILE, '$gtens_axes', nmult, 3, 3, axes )
   RETURN
END SUBROUTINE write_gtens

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_stev_cfp(ANISO_FILE, s, n, cfp)
   IMPLICIT NONE
   INTEGER, INTENT (IN)   :: ANISO_FILE
   INTEGER, INTENT (IN)   :: n   ! 2J+1, or 2L+1, i.e. the dimension of the J or L multiplet
   REAL (wp), INTENT (IN) :: cfp(n-1, -(n-1):(n-1) )
   CHARACTER (LEN=*)      :: s
   INTEGER                :: k, q

   ! s takes the value:
   ! 'L' == CFP for a term
   ! 'J' == CFP for a spin-orbit multiplet

   ierr=0
   IF ( (n<=0) ) THEN
      CALL WarningMessage(1,'write_stev_cfp_'//trim(s)//' :: nothing to write. Array size = 0.')
      RETURN
   END IF

   IF ( (trim(s).ne.'l')  .OR.  (trim(s).ne.'j') ) THEN
       CALL WarningMessage(1,'write_stev_cfp_'//trim(s)//' :: the parameter s='//trim(s)//'is not understood. '// &
                             'RETURN without writing cfp')
       RETURN
   END IF

   IF ( SUM(ABS( cfp(  1:(n-1), -(n-1):(n-1) ) )) <= TINY(Zero) ) THEN
      CALL WarningMessage(1,'write_stev_cfp_'//trim(s)//':: all array elements are zero = 0.')
   END IF

   REWIND ( ANISO_FILE )
   CALL file_advance_to_string( ANISO_FILE, '$stev_cfp_'//trim(s), line )

   IF (ierr /= 0 ) THEN
      WRITE ( ANISO_FILE, FMT='(A)', IOSTAT=ierr )
      WRITE ( ANISO_FILE, FMT='(A)', IOSTAT=ierr ) '$stev_cfp_'//trim(s)
   END IF

   ! if key is found, THEN rewrite the data
   WRITE ( ANISO_FILE, FMT=*, IOSTAT=ierr ) n
   DO k=2,n-1,2
      DO q=-k,k,2
         IF ( ABS(cfp(k,q)) > TINY(Zero) ) THEN
            WRITE ( ANISO_FILE, FMT=FMTCFP, IOSTAT=ierr ) k,q, cfp(k,q)
         END IF
         IF (ierr /= 0) THEN
            CALL WarningMessage(2,'write_stev_cfp_'//trim(s)//':: Something went wrong writing the array.')
         END IF
         IF (DBG) WRITE (StdOut,*) 'write_stev_cfp_'//trim(s)//'::  k, q =',k, q
         IF (DBG) FLUSH (StdOut)
      END DO
   END DO
   WRITE ( ANISO_FILE, FMT='(A)', IOSTAT=ierr )
   FLUSH ( ANISO_FILE )
RETURN
END SUBROUTINE write_stev_cfp

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_susc (ANISO_FILE, s, n, field, zj, t, z, x, x_tens)
   IMPLICIT NONE
   INTEGER, INTENT (IN)   :: ANISO_FILE
   INTEGER, INTENT (IN)   :: n                  ! number of temperature points
   REAL (wp), INTENT (IN) :: field              ! applied field
   REAL (wp), INTENT (IN) :: zJ                 ! intermolecular interaction
   REAL (wp), INTENT (IN) :: t(n)               ! temperature points
   REAL (wp), INTENT (IN) :: z(n)               ! partition function
   REAL (wp), INTENT (IN) :: x(n)               ! susceptibility X
   REAL (wp), INTENT (IN) :: x_tens(n,3,3)      ! susceptibility tensor  X_tens
   CHARACTER (LEN=*), INTENT (IN) :: s
   INTEGER :: i, j, k
   REAL (wp), EXTERNAL :: dnrm2_

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
!    IF ( (trim(s).ne.'x')                      .OR.   &
!         (trim(s).ne.'x_field')                .OR.   &
!         (trim(s).ne.'x_minus_one')            .OR.   &
!         (trim(s).ne.'x_field_minus_one')      .OR.   &
!         (trim(s).ne.'xt')                     .OR.   &
!         (trim(s).ne.'xt_field')               .OR.   &
!         (trim(s).ne.'xt_minus_one')           .OR.   &
!         (trim(s).ne.'xt_field_minus_one')     .OR.   &
!         (trim(s).ne.'x_zj')                   .OR.   &
!         (trim(s).ne.'x_field_zj')             .OR.   &
!         (trim(s).ne.'x_minus_one_zj')         .OR.   &
!         (trim(s).ne.'x_field_minus_one_zj')   .OR.   &
!         (trim(s).ne.'xt_zj')                  .OR.   &
!         (trim(s).ne.'xt_field_zj')            .OR.   &
!         (trim(s).ne.'xt_minus_one_zj')        .OR.   &
!         (trim(s).ne.'xt_field_minus_one_zj') ) THEN
!
!       CALL WarningMessage(1,'write_susc '//trim(s)//' :: the parameter s='//trim(s)//' is not understood. RETURN without'// &
!                             ' write of X data.')
!       RETURN
!   END IF

   ierr=0
   IF ( (n<=0) ) THEN
      CALL WarningMessage(1,'write_susc '//trim(s)//' :: nothing to write. Array size = 0.')
      RETURN
   END IF
   IF ( dnrm2_(n,T,1) <= TINY(Zero) ) THEN
      CALL WarningMessage(1,'write_susc '//trim(s)//' :: all array T elements are zero = 0.')
   END IF
   IF ( dnrm2_(n,X,1) <= TINY(Zero) ) THEN
      CALL WarningMessage(1,'write_susc '//trim(s)//' :: all array X elements are zero = 0.')
   END IF
   IF ( dnrm2_(n,Z,1) <= TINY(Zero) ) THEN
      CALL WarningMessage(1,'write_susc '//trim(s)//' :: all array Z elements are zero = 0.')
   END IF
   IF ( dnrm2_(3*3*n,X_tens,1) <= TINY(Zero) ) THEN
      CALL WarningMessage(1,'write_susc '//trim(s)//' :: all array X_tens elements are zero = 0.')
   END IF

   REWIND ( ANISO_FILE )
   CALL file_advance_to_string( ANISO_FILE, '$susceptibility_'//trim(s), line )

   IF (ierr /= 0 ) THEN
      WRITE ( ANISO_FILE, FMT='(A)', IOSTAT=ierr )
      WRITE ( ANISO_FILE, FMT='(A)', IOSTAT=ierr ) '$susceptibility_'//trim(s)
   END IF

   WRITE ( ANISO_FILE, FMT=FMTI, IOSTAT=ierr ) n
   WRITE ( ANISO_FILE, FMT=FMTR, IOSTAT=ierr ) zj, field
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'write_susc '//trim(s)//' :: Something went wrong writing the zJ and field values.')
   END IF
   WRITE ( ANISO_FILE, FMT=FMTR, IOSTAT=ierr ) (T(i),i=1,n)
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'write_susc '//trim(s)//' :: Something went wrong writing the T array.')
   END IF
   WRITE ( ANISO_FILE, FMT=FMTR, IOSTAT=ierr ) (Z(i),i=1,n)
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'write_susc '//trim(s)//' :: Something went wrong writing the Z array.')
   END IF
   WRITE ( ANISO_FILE, FMT=FMTR, IOSTAT=ierr ) (X(i),i=1,n)
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'write_susc '//trim(s)//' :: Something went wrong writing the X array.')
   END IF
   FLUSH ( ANISO_FILE )
   DO j=1,3
      DO k=1,3
         WRITE ( ANISO_FILE, FMT=FMTR, IOSTAT=ierr ) ( X_tens(i,j,k), i=1,n )
         IF (ierr /= 0) THEN
            CALL WarningMessage(2,'write_susc '//trim(s)//' :: Something went wrong writing the X_tens array.')
         END IF
      END DO
   END DO
   WRITE ( ANISO_FILE, FMT=*, IOSTAT=ierr )
   FLUSH ( ANISO_FILE )
   IF (DBG) FLUSH (StdOut)
   RETURN
END SUBROUTINE write_susc

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_magn (ANISO_FILE, nt, nh, nd, nss, zj, t, h, x, y, z, w, m, mav, energy)
   IMPLICIT NONE
   INTEGER, INTENT (IN)   :: ANISO_FILE
   INTEGER, INTENT (IN)   :: nt                         ! number of temperature points
   INTEGER, INTENT (IN)   :: nh                         ! number of field points
   INTEGER, INTENT (IN)   :: nd                         ! number of directions of aplied field
   INTEGER, INTENT (IN)   :: nss                        ! number of spin-orbit states included in the Zeeman interaction
   REAL (wp), INTENT (IN) :: zj                         ! inter-molecular parameter zJ
   REAL (wp), INTENT (IN) :: t(nt)                      ! temperature points
   REAL (wp), INTENT (IN) :: h(nh)                      ! field points
   REAL (wp), INTENT (IN) :: x(nd), y(nd), z(nd), w(nd) ! Lebedev grid directions and weight
   REAL (wp), INTENT (IN) :: m(nd, 3, nh, nt)           ! magnetisation vector
   REAL (wp), INTENT (IN) :: mav(nh, nt)                ! average magnetisation vector
   REAL (wp), INTENT (IN) :: energy(nd, nh, nss)        ! Zeeman energy states
   INTEGER              :: id, ih, it, i, iss, l
   REAL (wp), EXTERNAL   :: dnrm2_

   ierr=0
   IF ( (nt<=0) .OR. (nh<=0) .OR. (nd<=0) ) THEN
      CALL WarningMessage(1,'write_magn :: nothing to write. Array size = 0.')
      RETURN
   END IF
   IF ( dnrm2_(nt,t,1) <= TINY(Zero) ) THEN
      CALL WarningMessage(1,'write_magn :: all array T elements are zero = 0.')
   END IF
   IF ( dnrm2_(nh,h,1) <= TINY(Zero) ) THEN
      CALL WarningMessage(1,'write_magn :: all array H elements are zero = 0.')
   END IF
   IF ( dnrm2_(nd*3*nt*nh,m,1) <= TINY(Zero) ) THEN
      CALL WarningMessage(1,'write_magn :: all array M elements are zero = 0.')
   END IF
   IF ( dnrm2_(nt*nh,mav,1) <= TINY(Zero) ) THEN
      CALL WarningMessage(1,'write_magn :: all array MAV elements are zero = 0.')
   END IF
   IF ( dnrm2_(nd*nh*nss,energy,1) <= TINY(Zero) ) THEN
      CALL WarningMessage(1,'write_magn :: all array energy elements are zero = 0.')
   END IF
   IF ( (dnrm2_(nd,x,1) <= TINY(Zero)) .OR. (dnrm2_(nd,y,1) <= TINY(Zero)) .OR.  &
       (dnrm2_(nd,z,1) <= TINY(Zero)) .OR. (dnrm2_(nd,w,1) <= TINY(Zero)) )    THEN
      CALL WarningMessage(1,'write_magn :: all array of X,Y,X,W elements are zero = 0.')
   END IF

   REWIND ( ANISO_FILE )
   CALL file_advance_to_string( ANISO_FILE, '$magnetisation', line )

   ! if key is found, THEN rewrite the data
   IF ( ierr /= 0 ) THEN
      WRITE ( ANISO_FILE, FMT='(A)', IOSTAT=ierr )
      WRITE ( ANISO_FILE, FMT='(A)', IOSTAT=ierr ) '$magnetisation'
   END IF

   WRITE ( ANISO_FILE, FMT=FMTI, IOSTAT=ierr ) nt, nh, nd, nss
   WRITE ( ANISO_FILE, FMT=FMTR, IOSTAT=ierr ) zJ
   WRITE ( ANISO_FILE, FMT=FMTR, IOSTAT=ierr ) (T(i),i=1,nt)
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'write_magn :: Something went wrong writing the T array.')
   END IF
   ! Magnetic field points
   WRITE ( ANISO_FILE, FMT=FMTR, IOSTAT=ierr ) (H(i),i=1,nh)
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'write_magn :: Something went wrong writing the H array.')
   END IF
   ! Lebedev grid data (X, Y, Z, W)
   WRITE ( ANISO_FILE, FMT=FMTR, IOSTAT=ierr ) ( X(i), i=1,nd )
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'write_magn :: Something went wrong writing the Lebedev grid X rray.')
   END IF
   WRITE ( ANISO_FILE, FMT=FMTR, IOSTAT=ierr ) ( Y(i), i=1,nd )
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'write_magn :: Something went wrong writing the Lebedev grid Y array.')
   END IF
   WRITE ( ANISO_FILE, FMT=FMTR, IOSTAT=ierr ) ( Z(i), i=1,nd )
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'write_magn :: Something went wrong writing the Lebedev grid Z array.')
   END IF
   WRITE ( ANISO_FILE, FMT=FMTR, IOSTAT=ierr ) ( W(i), i=1,nd )
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'write_magn :: Something went wrong writing the Lebedev grid W array.')
   END IF
   FLUSH ( ANISO_FILE )
   ! Zeeman energy
   DO id=1,nd
      DO ih=1,nh
         WRITE ( ANISO_FILE, FMT=FMTR, IOSTAT=ierr )  (energy(id,ih,iss),iss=1,nss)
         IF (ierr /= 0) THEN
            CALL WarningMessage(2,'write_magn :: Something went wrong writing the Zeeman energy data.')
         END IF
      END DO
   END DO
   FLUSH ( ANISO_FILE )
   ! magnetisation vector data
   DO id=1,nd
      DO l=1,3
         DO ih=1,nh
            WRITE ( ANISO_FILE, FMT=FMTR, IOSTAT=ierr ) (m(id,l,ih,it),it=1,nt)
            IF (ierr /= 0) THEN
               CALL WarningMessage(2,'write_magn :: Something went wrong writing the M data.')
            END IF
         END DO
      END DO
   END DO
   FLUSH ( ANISO_FILE )
   DO ih=1,nh
      WRITE ( ANISO_FILE, FMT=FMTR, IOSTAT=ierr ) (mav(ih,it),it=1,nt)
      IF (ierr /= 0) THEN
         CALL WarningMessage(2,'write_magn :: Something went wrong writing the average M data.')
      END IF
   END DO
   WRITE ( ANISO_FILE, FMT=*, IOSTAT=ierr )
   FLUSH ( ANISO_FILE )
   IF (DBG) FLUSH (StdOut)
   RETURN
END SUBROUTINE write_magn














































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
SUBROUTINE open_datafile_write(DATA_FILE,DATA_FILE_NAME)
   IMPLICIT NONE
   INTEGER            :: DATA_FILE
   CHARACTER(LEN=180) :: DATA_FILE_NAME
!#ifdef _ALONE_
!   INTEGER :: ierr
!   ierr=0
!   OPEN (UNIT=DATA_FILE, FILE=DATA_FILE_NAME, FORM='formatted', STATUS='old', ACCESS='sequential', ACTION='write', IOSTAT=ierr )
!   IF (ierr /= 0) THEN
!      CALL WarningMessage(2,'open_datafile_write:: Something went wrong opening DATA_FILE')
!      STOP
!   END IF
!#else
   Call molcas_open(DATA_FILE,DATA_FILE_NAME)
!#endif

   RETURN
END SUBROUTINE open_datafile_write

!--------------------------------------------------------------------------------------------------!
SUBROUTINE open_aniso_file(ANISO_FILE,ANISO_FILE_NAME)
   IMPLICIT NONE
   INTEGER            :: ANISO_FILE
   CHARACTER(LEN=180) :: ANISO_FILE_NAME
!#ifdef _ALONE_
!   INTEGER :: ierr
!   ierr=0
!   OPEN (UNIT=ANISO_FILE, FORM='formatted', STATUS='old', ACCESS='sequential', ACTION='write', IOSTAT=ierr )
!   IF (ierr /= 0) THEN
!      CALL WarningMessage(2,'open_datafile_write:: Something went wrong opening ANISO_FILE')
!      STOP
!   END IF
!#else
   Call molcas_open(ANISO_FILE,ANISO_FILE_NAME)
!#endif

   RETURN
END SUBROUTINE open_aniso_file

!--------------------------------------------------------------------------------------------------!
SUBROUTINE open_datafile_read(DATA_FILE,DATA_FILE_NAME)
   IMPLICIT NONE
   INTEGER            :: DATA_FILE
   CHARACTER(LEN=180) :: DATA_FILE_NAME
!#ifdef _ALONE_
!   INTEGER :: ierr
!   ierr=0
!   OPEN (UNIT=DATA_FILE, FORM='FORMATTED', STATUS='old', ACCESS='sequential', ACTION='read', IOSTAT=ierr )
!   IF (ierr /= 0) THEN
!      CALL WarningMessage(2,'open_datafile_read:: Something went wrong opening DATA_FILE')
!      STOP
!   END IF
!#else
   Call molcas_open(DATA_FILE,DATA_FILE_NAME)
!#endif
   RETURN
END SUBROUTINE open_datafile_read

!--------------------------------------------------------------------------------------------------!
SUBROUTINE close_datafile(DATA_FILE)
   IMPLICIT NONE
   INTEGER            :: DATA_FILE
   INTEGER :: ierr
   ierr=0
   CLOSE (UNIT=DATA_FILE, IOSTAT=ierr )
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'close_datafile:: Something went wrong closing DATA_FILE')
      !STOP
   END IF
   RETURN
END SUBROUTINE close_datafile

!--------------------------------------------------------------------------------------------------!
SUBROUTINE close_anisofile(ANISO_FILE)
   IMPLICIT NONE
   INTEGER :: ANISO_FILE
   INTEGER :: ierr
   ierr=0
   CLOSE (UNIT=ANISO_FILE, IOSTAT=ierr )
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'close_datafile:: Something went wrong closing ANISO_FILE')
      !STOP
   END IF
   RETURN
END SUBROUTINE close_anisofile
!--------------------------------------------------------------------------------------------------!
LOGICAL FUNCTION key_found (key)
   IMPLICIT NONE
   CHARACTER (LEN=*), INTENT (IN) :: key

   key_found=.false.
   REWIND ( DATA_FILE )
   CALL file_advance_to_string( DATA_FILE, key, line )
   IF ( index (line, trim(key)) /= 0 ) THEN
      key_found=.true.
   END IF

   RETURN
END FUNCTION key_found

!--------------------------------------------------------------------------------------------------!
SUBROUTINE file_advance_to_string ( LU, key, line )
   IMPLICIT NONE
   INTEGER, INTENT (IN)           :: LU
   INTEGER                       :: ios
   INTEGER                       :: num_read
   CHARACTER (LEN=*), INTENT (IN) :: key
   CHARACTER (LEN=*)             :: line

   ierr = 0
   num_read = 0

   REWIND ( LU )
   DO
      READ ( LU, '(a)', IOSTAT = ios ) line

      IF ( ios /= 0 ) THEN
         !STOP 'file_advance_to_string:: error reading line'
         !print '(A,i0,2A)', 'file_advance_to_string:: error reading line: LU=',LU,' key=',key
         GO TO 1
      END IF

      num_read = num_read + 1

      IF ( index (line, trim(key)) /= 0 ) THEN
         RETURN
      END IF
   END DO

1  CONTINUE
   line = ' '
   ierr = 1

   IF ( DBG ) THEN
      WRITE (StdOut, '(a)' ) ' '
      WRITE (StdOut, '(a)' ) 'FILE_ADVANCE_TO_STRING - Warning!'
      WRITE (StdOut, '(a)' ) '  Did not find the key:'
      WRITE (StdOut, '(a)' ) '    ' // trim ( key )
      WRITE (StdOut, '(a,i6)' ) '  Number of lines read was ', num_read
   END IF

   RETURN
END SUBROUTINE file_advance_to_string
!--------------------------------------------------------------------------------------------------!


!--------------------------------------------------------------------------------------------------!
LOGICAL FUNCTION inquire_key_presence ( LU, key )
   IMPLICIT NONE
   INTEGER, INTENT (IN)   :: LU
   CHARACTER (LEN=*)     :: key
   INTEGER               :: ios

   inquire_key_presence=.false.
   REWIND ( LU )
   DO
      READ ( LU, '(a)', IOSTAT = ios ) line
      IF ( ios /= 0 ) THEN
         !STOP 'inquire_key_presence:: error reading line'
         CALL WarningMessage(1,'inquire_key_presence:: error reading line')
      END IF
      IF ( index (line, trim(key)) /= 0 ) THEN
         inquire_key_presence=.true.
         RETURN
      END IF
   END DO
   line = ' '
   RETURN
END FUNCTION inquire_key_presence
!--------------------------------------------------------------------------------------------------!




































!--------------------------------------------------------------------------------------------------!
!   READING  SUBROUTINES
!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_INTEGER_scalar( LU, key, i )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (OUT)         :: i
   ierr=0
   i=0
   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )
   READ ( LU, FMT=*, IOSTAT=ierr ) i
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_INTEGER_scalar:: Something went wrong reading key'//trim(key))
   END IF
   IF (DBG) WRITE (StdOut,*) 'read_INTEGER_scalar:: key =',trim(key)
   IF (DBG) WRITE (StdOut,*) 'read_INTEGER_scalar::   i =',i
   RETURN
END SUBROUTINE read_INTEGER_scalar

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_real_scalar( LU, key, r )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   REAL (wp), INTENT (OUT)        :: r
   ierr=0
   r=Zero
   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )
   READ ( LU, FMT=*, IOSTAT=ierr ) r
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_real_scalar:: Something went wrong reading key'//trim(key))
   END IF
   IF (DBG) WRITE (StdOut,*) 'read_real_scalar:: key =',trim(key)
   IF (DBG) WRITE (StdOut,*) 'read_real_scalar::   r =',r
   RETURN
END SUBROUTINE read_real_scalar

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_complex_scalar( LU, key, c )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   COMPLEX (wp), INTENT (OUT)     :: c
   REAL (wp)                     :: rr, ri

   ierr=0
   rr=Zero ; ri=Zero ; c=ZeroC
   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )
   READ ( LU, FMT=*, IOSTAT=ierr ) rr, ri
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_complex_scalar:: Something went wrong reading key'//trim(key))
   END IF
   IF (DBG) WRITE (StdOut,*) 'read_complex_scalar::   key =',trim(key)
   IF (DBG) WRITE (StdOut,*) 'read_complex_scalar:: (r,i) =', rr, ri
   IF (DBG) WRITE (StdOut,*) 'read_complex_scalar::     c =', c
   c=CMPLX(rr,ri,wp)
   RETURN
END SUBROUTINE read_complex_scalar

!--------------------------------------------------------------------------------------------------!
!SUBROUTINE read_string (key, length, s)
!result ( s )
!   IMPLICIT NONE
!   CHARACTER (LEN=*), INTENT (IN)      :: key
!   INTEGER, INTENT (IN)               :: length
!   CHARACTER (LEN=length),INTENT (OUT) :: s
!   CHARACTER (LEN=300)                :: c, f
!   INTEGER                           :: i
!
!   IF (DBG) WRITE (StdOut,*) 'read_string::    key =',trim(key)
!   IF (DBG) WRITE (StdOut,*) 'read_string:: length =',length
!   FLUSH (StdOut)
!   WRITE (f,'(A,i0,2A)') '"(A',length,')"'
!   WRITE (StdOut,'(2A)') 'format =', trim(f)
!   FLUSH (StdOut)
!   REWIND ( DATA_FILE )
!   CALL file_advance_to_string( key, line )
!   READ ( DATA_FILE, FMT=f, IOSTAT=ierr ) c
!   WRITE (StdOut,'(2A)') 'c =', trim(c)
!   WRITE (s,'(A)') trim(c)
!
!   DO i=1,LEN(TRIM(LINE))
!     READ (LINE,'(A300)') c
!     WRITE (s,'(A)') trim(c)
!     FLUSH (StdOut)
!     IF (DBG) WRITE (StdOut,*) 'read_string::   c =', trim(c)
!     FLUSH (StdOut)
!   END DO
!
!   RETURN
!END SUBROUTINE read_string

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_1d_size( LU, key, n )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (OUT)         :: n
   ierr=0
   n=0
   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )
   READ ( LU, FMT=*, IOSTAT=ierr ) n
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_1d_size:: Something went wrong reading key'//trim(key))
   END IF
   IF (DBG) WRITE (StdOut,*) 'read_1d_size:: key =',trim(key)
   IF (DBG) WRITE (StdOut,*) 'read_1d_size::   n =',n
   RETURN
END SUBROUTINE read_1d_size

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_2d_size( LU, key, n1, n2 )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (OUT)         :: n1, n2
   ierr=0
   n1=0 ; n2=0 ;
   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )
   READ ( LU, FMT=*, IOSTAT=ierr ) n1, n2
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_2d_size:: Something went wrong reading key'//trim(key))
   END IF
   IF (DBG) WRITE (StdOut,*) 'read_2d_size:: key =',trim(key)
   IF (DBG) WRITE (StdOut,*) 'read_2d_size::  n1 =',n1
   IF (DBG) WRITE (StdOut,*) 'read_2d_size::  n2 =',n2
   RETURN
END SUBROUTINE read_2d_size

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_3d_size( LU, key, n1, n2, n3 )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (OUT)         :: n1, n2, n3
   ierr=0
   n1=0 ; n2=0 ; n3=0 ;
   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )
   READ ( LU, FMT=*, IOSTAT=ierr ) n1, n2, n3
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_3d_size:: Something went wrong reading key'//trim(key))
   END IF
   IF (DBG) WRITE (StdOut,*) 'read_3d_size:: key =',trim(key)
   IF (DBG) WRITE (StdOut,*) 'read_3d_size::  n1 =',n1
   IF (DBG) WRITE (StdOut,*) 'read_3d_size::  n2 =',n2
   IF (DBG) WRITE (StdOut,*) 'read_3d_size::  n3 =',n3
   RETURN
END SUBROUTINE read_3d_size

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_4d_size( LU, key, n1, n2, n3, n4 )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (OUT)         :: n1, n2, n3, n4
   ierr=0
   n1=0 ; n2=0 ; n3=0 ; n4=0 ;
   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )
   READ ( LU, FMT=*, IOSTAT=ierr ) n1, n2, n3, n4
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_4d_size:: Something went wrong reading key'//trim(key))
   END IF
   IF (DBG) WRITE (StdOut,*) 'read_4d_size:: key =',trim(key)
   IF (DBG) WRITE (StdOut,*) 'read_4d_size::  n1 =',n1
   IF (DBG) WRITE (StdOut,*) 'read_4d_size::  n2 =',n2
   IF (DBG) WRITE (StdOut,*) 'read_4d_size::  n3 =',n3
   IF (DBG) WRITE (StdOut,*) 'read_4d_size::  n4 =',n4
   RETURN
END SUBROUTINE read_4d_size

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_1d_INTEGER_array( LU, key, n, array )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (IN)          :: n
   INTEGER, INTENT (OUT)         :: array(n)
   INTEGER                      :: i
   ierr=0
   array=0
   IF (n<=0) THEN
      CALL WarningMessage(1,'read_1d_INTEGER_array:: nothing to read. Array size = 0.')
      RETURN
   END IF
   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )

   READ ( LU, FMT=*, IOSTAT=ierr ) i
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_1d_INTEGER_array:: Something went wrong reading key'//trim(key))
   END IF
   IF (DBG) WRITE (StdOut,*) 'read_1d_INTEGER_array:: key =',trim(key)
   IF (DBG) WRITE (StdOut,*) 'read_1d_INTEGER_array::   n =',i
   IF ( i/=n ) THEN
      CALL WarningMessage(2,'read_1d_INTEGER_array:: sizes of the array are different from the ones used '// &
                            'to CALL this SUBROUTINE')
   END IF

   READ ( LU, FMT=*, IOSTAT=ierr ) (array(i),i=1,n)
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_1d_INTEGER_array:: Something went wrong reading the array.')
   END IF
   IF (DBG) WRITE (StdOut,*) 'read_1d_INTEGER_array:: array =',(array(i),i=1,n)

   RETURN
END SUBROUTINE read_1d_INTEGER_array

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_2d_INTEGER_array( LU, key, n1, n2, array )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (IN)          :: n1, n2
   INTEGER, INTENT (OUT)         :: array(n1,n2)
   INTEGER                      :: i,j
   ierr=0
   array=0
   IF ( (n1<=0).OR.(n2<=0) ) THEN
      CALL WarningMessage(1,'read_2d_INTEGER_array:: nothing to read. Array size = 0.')
      IF (DBG) WRITE (StdOut,*) 'read_2d_INTEGER_array::   n1 =',n1
      IF (DBG) WRITE (StdOut,*) 'read_2d_INTEGER_array::   n2 =',n2
      IF (DBG) FLUSH (StdOut)
      RETURN
   END IF
   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )

   READ ( LU, FMT=*, IOSTAT=ierr ) i,j
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_2d_INTEGER_array:: Something went wrong reading key'//trim(key))
   END IF
   IF (DBG) WRITE (StdOut,*) 'read_2d_INTEGER_array:: key =',trim(key)
   IF (DBG) WRITE (StdOut,*) 'read_2d_INTEGER_array::  n1 =',i
   IF (DBG) WRITE (StdOut,*) 'read_2d_INTEGER_array::  n2 =',j
   IF ( (i/=n1).OR.(j/=n2) ) THEN
      CALL WarningMessage(2,'read_2d_INTEGER_array:: sizes of the array are different from the ones used '// &
                            'to CALL this SUBROUTINE')
   END IF
   DO i=1, n1
      READ ( LU, FMT=*, IOSTAT=ierr ) (array(i,j),j=1,n2)
      IF (ierr /= 0) THEN
         CALL WarningMessage(2,'read_2d_INTEGER_array:: Something went wrong reading the array.')
      END IF
      IF (DBG) WRITE (StdOut,*) 'read_2d_INTEGER_array::  i =',i
      IF (DBG) FLUSH (StdOut)
   END DO

   RETURN
END SUBROUTINE read_2d_INTEGER_array

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_3d_INTEGER_array( LU, key, n1, n2, n3, array )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (IN)          :: n1, n2, n3
   INTEGER, INTENT (OUT)         :: array(n1,n2,n3)
   INTEGER                      :: i,j,k
   ierr=0
   array=0
   IF ( (n1<=0) .OR. (n2<=0) .OR. (n3<=0) ) THEN
      CALL WarningMessage(1,'read_3d_INTEGER_array:: nothing to read. Array size = 0.')
      IF (DBG) WRITE (StdOut,*) 'read_3d_INTEGER_array::   n1 =',n1
      IF (DBG) WRITE (StdOut,*) 'read_3d_INTEGER_array::   n2 =',n2
      IF (DBG) WRITE (StdOut,*) 'read_3d_INTEGER_array::   n3 =',n3
      IF (DBG) FLUSH (StdOut)
      RETURN
   END IF
   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )

   READ ( LU, FMT=*, IOSTAT=ierr ) i,j,k
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_3d_INTEGER_array:: Something went wrong reading key'//trim(key))
   END IF
   IF (DBG) WRITE (StdOut,*) 'read_3d_INTEGER_array:: key =',trim(key)
   IF (DBG) WRITE (StdOut,*) 'read_3d_INTEGER_array::  n1 =',i
   IF (DBG) WRITE (StdOut,*) 'read_3d_INTEGER_array::  n2 =',j
   IF (DBG) WRITE (StdOut,*) 'read_3d_INTEGER_array::  n3 =',k
   IF ( (i/=n1) .OR. (j/=n2) .OR. (k/=n3) ) THEN
      CALL WarningMessage(2,'read_3d_INTEGER_array:: sizes of the array are different from the ones used '// &
                            'to CALL this SUBROUTINE')
   END IF
   DO i=1, n1
      DO j=1, n2
         READ ( LU, FMT=*, IOSTAT=ierr ) (array(i,j,k),k=1,n3)
         IF (ierr /= 0) THEN
            CALL WarningMessage(2,'read_3d_INTEGER_array:: Something went wrong reading the array.')
         END IF
         IF (DBG) WRITE (StdOut,*) 'read_3d_INTEGER_array::  i,j =',i,j
         IF (DBG) FLUSH (StdOut)
      END DO
   END DO

   RETURN
END SUBROUTINE read_3d_INTEGER_array

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_4d_INTEGER_array( LU, key, n1, n2, n3, n4, array )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (IN)          :: n1, n2, n3, n4
   INTEGER, INTENT (OUT)         :: array(n1,n2,n3,n4)
   INTEGER                      :: i,j,k,l
   ierr=0
   array=0
   IF ( (n1<=0) .OR. (n2<=0) .OR. (n3<=0) .OR. (n4<=0) ) THEN
      CALL WarningMessage(1,'read_4d_INTEGER_array:: nothing to read. Array size = 0.')
      IF (DBG) WRITE (StdOut,*) 'read_4d_INTEGER_array::   n1 =',n1
      IF (DBG) WRITE (StdOut,*) 'read_4d_INTEGER_array::   n2 =',n2
      IF (DBG) WRITE (StdOut,*) 'read_4d_INTEGER_array::   n3 =',n3
      IF (DBG) WRITE (StdOut,*) 'read_4d_INTEGER_array::   n4 =',n4
      IF (DBG) FLUSH (StdOut)
      RETURN
   END IF
   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )

   READ ( LU, FMT=*, IOSTAT=ierr ) i,j,k,l
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_4d_INTEGER_array:: Something went wrong reading key'//trim(key))
   END IF
   IF (DBG) WRITE (StdOut,*) 'read_4d_INTEGER_array:: key =',trim(key)
   IF (DBG) WRITE (StdOut,*) 'read_4d_INTEGER_array::  n1 =',i
   IF (DBG) WRITE (StdOut,*) 'read_4d_INTEGER_array::  n2 =',j
   IF (DBG) WRITE (StdOut,*) 'read_4d_INTEGER_array::  n3 =',k
   IF (DBG) WRITE (StdOut,*) 'read_4d_INTEGER_array::  n4 =',l
   IF ( (i/=n1) .OR. (j/=n2) .OR. (k/=n3) .OR. (l/=n4) ) THEN
      CALL WarningMessage(2,'read_4d_INTEGER_array:: sizes of the array are different from the ones used '// &
                            'to CALL this SUBROUTINE')
   END IF
   DO i=1, n1
      DO j=1, n2
         DO k=1, n3
            READ ( LU, FMT=*, IOSTAT=ierr ) (array(i,j,k,l),l=1,n4)
            IF (ierr /= 0) THEN
               CALL WarningMessage(2,'read_4d_INTEGER_array:: Something went wrong reading the array.')
            END IF
            IF (DBG) WRITE (StdOut,*) 'read_4d_INTEGER_array::  i,j,k =',i,j,k
            IF (DBG) FLUSH (StdOut)
         END DO
      END DO
   END DO

   RETURN
END SUBROUTINE read_4d_INTEGER_array

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_1d_real_array( LU, key, n, array )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (IN)          :: n
   REAL (wp), INTENT (OUT)        :: array(n)
   INTEGER                      :: i
   ierr=0
   array=Zero
   IF (n<=0) THEN
      CALL WarningMessage(1,'read_1d_real_array:: nothing to read. Array size = 0.')
      RETURN
   END IF
   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )

   READ ( LU, FMT=*, IOSTAT=ierr ) i
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_1d_real_array:: Something went wrong reading key'//trim(key))
   END IF
   IF (DBG) WRITE (StdOut,*) 'read_1d_real_array:: key =',trim(key)
   IF (DBG) WRITE (StdOut,*) 'read_1d_real_array::   n =',i
   IF ( (i/=n) ) THEN
      CALL WarningMessage(2,'read_1d_real_array:: sizes of the array are different from the ones used '// &
                            'to CALL this SUBROUTINE')
   END IF

   READ ( LU, FMT=*, IOSTAT=ierr ) (array(i),i=1,n)
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_1d_real_array:: Something went wrong reading the array.')
   END IF
   IF (DBG) WRITE (StdOut,*) 'read_1d_real_array:: array =',(array(i),i=1,n)

   RETURN
END SUBROUTINE read_1d_real_array

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_2d_real_array( LU, key, n1, n2, array )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (IN)          :: n1, n2
   REAL (wp), INTENT (OUT)        :: array(n1,n2)
   INTEGER                      :: i,j
   ierr=0
   array=Zero
   IF ( (n1<=0) .OR. (n2<=0) ) THEN
      CALL WarningMessage(1,'read_2d_real_array:: nothing to read. Array size = 0.')
      IF (DBG) WRITE (StdOut,*) 'read_2d_real_array::   n1 =',n1
      IF (DBG) WRITE (StdOut,*) 'read_2d_real_array::   n2 =',n2
      IF (DBG) FLUSH (StdOut)
      RETURN
   END IF
   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )

   READ ( LU, FMT=*, IOSTAT=ierr ) i,j
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_2d_real_array:: Something went wrong reading key'//trim(key))
   END IF
   IF (DBG) WRITE (StdOut,*) 'read_2d_real_array:: key =',trim(key)
   IF (DBG) WRITE (StdOut,*) 'read_2d_real_array::  n1 =',i
   IF (DBG) WRITE (StdOut,*) 'read_2d_real_array::  n2 =',j
   IF ( (i/=n1) .OR. (j/=n2) ) THEN
      CALL WarningMessage(2,'read_2d_real_array:: sizes of the array are different from the ones used '// &
                            'to CALL this SUBROUTINE')
   END IF
   DO i=1, n1
      READ ( LU, FMT=*, IOSTAT=ierr ) (array(i,j),j=1,n2)
      IF (ierr /= 0) THEN
         CALL WarningMessage(2,'read_2d_real_array:: Something went wrong reading the array.')
      END IF
      IF (DBG) WRITE (StdOut,*) 'read_2d_real_array::  i =',i
      IF (DBG) FLUSH (StdOut)
   END DO

   RETURN
END SUBROUTINE read_2d_real_array

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_3d_real_array( LU, key, n1, n2, n3, array )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (IN)          :: n1, n2, n3
   REAL (wp), INTENT (OUT)        :: array(n1,n2,n3)
   INTEGER                      :: i,j,k
   ierr=0
   array=Zero
   IF ( (n1<=0) .OR. (n2<=0) .OR. (n3<=0) ) THEN
      CALL WarningMessage(1,'read_3d_real_array:: nothing to read. Array size = 0.')
      IF (DBG) WRITE (StdOut,*) 'read_3d_real_array::   n1 =',n1
      IF (DBG) WRITE (StdOut,*) 'read_3d_real_array::   n2 =',n2
      IF (DBG) WRITE (StdOut,*) 'read_3d_real_array::   n3 =',n3
      IF (DBG) FLUSH (StdOut)
      RETURN
   END IF
   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )

   READ ( LU, FMT=*, IOSTAT=ierr ) i,j,k
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_3d_real_array:: Something went wrong reading key'//trim(key))
   END IF
   IF (DBG) WRITE (StdOut,*) 'read_3d_real_array:: key =',trim(key)
   IF (DBG) WRITE (StdOut,*) 'read_3d_real_array::  n1 =',i
   IF (DBG) WRITE (StdOut,*) 'read_3d_real_array::  n2 =',j
   IF (DBG) WRITE (StdOut,*) 'read_3d_real_array::  n3 =',k
   IF ( (i/=n1) .OR. (j/=n2) .OR. (k/=n3) ) THEN
      CALL WarningMessage(2,'read_3d_real_array:: sizes of the array are different from the ones used '// &
                            'to CALL this SUBROUTINE')
   END IF
   DO i=1, n1
      DO j=1, n2
         READ ( LU, FMT=*, IOSTAT=ierr ) (array(i,j,k),k=1,n3)
         IF (ierr /= 0) THEN
            CALL WarningMessage(2,'read_3d_real_array:: Something went wrong reading the array.')
         END IF
         IF (DBG) WRITE (StdOut,*) 'read_3d_real_array::  i,j =',i,j
         IF (DBG) FLUSH (StdOut)
      END DO
   END DO

   RETURN
END SUBROUTINE read_3d_real_array

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_4d_real_array( LU, key, n1, n2, n3, n4, array )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (IN)          :: n1, n2, n3, n4
   REAL (wp), INTENT (OUT)        :: array(n1,n2,n3,n4)
   INTEGER                      :: i,j,k,l
   ierr=0
   array=Zero
   IF ( (n1<=0) .OR. (n2<=0) .OR. (n3<=0) .OR. (n4<=0) ) THEN
      CALL WarningMessage(1,'read_4d_real_array:: nothing to read. Array size = 0.')
      IF (DBG) WRITE (StdOut,*) 'read_4d_real_array::   n1 =',n1
      IF (DBG) WRITE (StdOut,*) 'read_4d_real_array::   n2 =',n2
      IF (DBG) WRITE (StdOut,*) 'read_4d_real_array::   n3 =',n3
      IF (DBG) WRITE (StdOut,*) 'read_4d_real_array::   n4 =',n4
      IF (DBG) FLUSH (StdOut)
      RETURN
   END IF
   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )

   READ ( LU, FMT=*, IOSTAT=ierr ) i,j,k,l
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_4d_real_array:: Something went wrong reading key'//trim(key))
   END IF
   IF (DBG) WRITE (StdOut,*) 'read_4d_real_array:: key =',trim(key)
   IF (DBG) WRITE (StdOut,*) 'read_4d_real_array::  n1 =',i
   IF (DBG) WRITE (StdOut,*) 'read_4d_real_array::  n2 =',j
   IF (DBG) WRITE (StdOut,*) 'read_4d_real_array::  n3 =',k
   IF (DBG) WRITE (StdOut,*) 'read_4d_real_array::  n4 =',l
   IF ( (i/=n1) .OR. (j/=n2) .OR. (k/=n3) .OR. (l/=n4) ) THEN
      CALL WarningMessage(2,'read_4d_real_array:: sizes of the array are different from the ones used '// &
                            'to CALL this SUBROUTINE')
   END IF
   DO i=1, n1
      DO j=1, n2
         DO k=1, n3
            READ ( LU, FMT=*, IOSTAT=ierr ) (array(i,j,k,l),l=1,n4)
            IF (ierr /= 0) THEN
               CALL WarningMessage(2,'read_4d_real_array:: Something went wrong reading the array.')
            END IF
            IF (DBG) WRITE (StdOut,*) 'read_4d_real_array::  i,j,k =',i,j,k
            IF (DBG) FLUSH (StdOut)
         END DO
      END DO
   END DO

   RETURN
END SUBROUTINE read_4d_real_array

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_1d_complex_array( LU, key, n, array )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (IN)          :: n
   COMPLEX (wp), INTENT (OUT)     :: array(n)
   REAL (wp), ALLOCATABLE        :: rr(:), ri(:)
   INTEGER                      :: i
   ierr=0
   array=ZeroC
   IF (n<=0) THEN
      CALL WarningMessage(1,'read_1d_complex_array:: nothing to read. Array size = 0.')
      RETURN
   END IF
   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )

   READ ( LU, FMT=*, IOSTAT=ierr ) i
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_1d_complex_array:: Something went wrong reading key'//trim(key))
   END IF
   IF (DBG) WRITE (StdOut,*) 'read_1d_complex_array:: key =',trim(key)
   IF (DBG) WRITE (StdOut,*) 'read_1d_complex_array::   n =',i
   IF ( (i/=n) ) THEN
      CALL WarningMessage(2,'read_1d_complex_array:: sizes of the array are different from the ones used '// &
                            'to CALL this SUBROUTINE')
   END IF
   ALLOCATE (rr(n))
   ALLOCATE (ri(n))
   rr=Zero ; ri=Zero
   READ ( LU, FMT=*, IOSTAT=ierr ) ( rr(i), ri(i), i=1,n )
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_1d_complex_array:: Something went wrong reading the array.')
   END IF
   IF (DBG) WRITE (StdOut,*) 'read_1d_complex_array:: array =',(rr(i),ri(i),i=1,n)
   DO i=1, n
     array(i)=CMPLX(rr(i), ri(i), wp)
   END DO
   DEALLOCATE (rr)
   DEALLOCATE (ri)
   RETURN
END SUBROUTINE read_1d_complex_array

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_2d_complex_array( LU, key, n1, n2, array )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (IN)          :: n1, n2
   COMPLEX (wp), INTENT (OUT)     :: array(n1,n2)
   REAL (wp), ALLOCATABLE        :: rr(:,:), ri(:,:)
   INTEGER                      :: i,j
   ierr=0
   array=ZeroC
   IF ( (n1<=0) .OR. (n2<=0) ) THEN
      CALL WarningMessage(1,'read_2d_complex_array:: nothing to read. Array size = 0.')
      IF (DBG) WRITE (StdOut,*) 'read_2d_complex_array::   n1 =',n1
      IF (DBG) WRITE (StdOut,*) 'read_2d_complex_array::   n2 =',n2
      IF (DBG) FLUSH (StdOut)
      RETURN
   END IF
   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )

   READ ( LU, FMT=*, IOSTAT=ierr ) i,j
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_2d_complex_array:: Something went wrong reading key'//trim(key))
   END IF
   IF (DBG) WRITE (StdOut,*) 'read_2d_complex_array:: key =',trim(key)
   IF (DBG) WRITE (StdOut,*) 'read_2d_complex_array::  n1 =',i
   IF (DBG) WRITE (StdOut,*) 'read_2d_complex_array::  n2 =',j
   IF ( (i/=n1) .OR. (j/=n2) ) THEN
      CALL WarningMessage(2,'read_2d_complex_array:: sizes of the array are different from the ones used '// &
                            'to CALL this SUBROUTINE')
   END IF
   ALLOCATE (rr(n1,n2))
   ALLOCATE (ri(n1,n2))
   rr=Zero ; ri=Zero
   DO i=1, n1
      READ ( LU, FMT=*, IOSTAT=ierr ) (rr(i,j), ri(i,j), j=1,n2)
      IF (ierr /= 0) THEN
         CALL WarningMessage(2,'read_2d_complex_array:: Something went wrong reading the array.')
      END IF
      IF (DBG) WRITE (StdOut,*) 'read_2d_complex_array::  i =',i
      IF (DBG) FLUSH (StdOut)
   END DO
   DO i=1, n1
      DO j=1, n2
         array(i,j)=CMPLX(rr(i,j), ri(i,j), wp)
      END DO
   END DO
   DEALLOCATE (rr)
   DEALLOCATE (ri)
   RETURN
END SUBROUTINE read_2d_complex_array

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_3d_complex_array( LU, key, n1, n2, n3, array )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (IN)          :: n1, n2, n3
   COMPLEX (wp), INTENT (OUT)     :: array(n1,n2,n3)
   REAL (wp), ALLOCATABLE        :: rr(:,:,:), ri(:,:,:)
   INTEGER                      :: i,j,k
   ierr=0
   array=ZeroC
   IF ( (n1<=0) .OR. (n2<=0) .OR. (n3<=0) ) THEN
      CALL WarningMessage(1,'read_3d_complex_array:: nothing to read. Array size = 0.')
      IF (DBG) WRITE (StdOut,*) 'read_3d_complex_array::   n1 =',n1
      IF (DBG) WRITE (StdOut,*) 'read_3d_complex_array::   n2 =',n2
      IF (DBG) WRITE (StdOut,*) 'read_3d_complex_array::   n3 =',n3
      IF (DBG) FLUSH (StdOut)
      RETURN
   END IF
   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )

   READ ( LU, FMT=*, IOSTAT=ierr ) i,j,k
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_3d_complex_array:: Something went wrong reading key'//trim(key))
   END IF
   IF (DBG) WRITE (StdOut,*) 'read_3d_complex_array:: key =',trim(key)
   IF (DBG) WRITE (StdOut,*) 'read_3d_complex_array::  n1 =',i
   IF (DBG) WRITE (StdOut,*) 'read_3d_complex_array::  n2 =',j
   IF (DBG) WRITE (StdOut,*) 'read_3d_complex_array::  n3 =',k
   IF ( (i/=n1) .OR. (j/=n2) .OR. (k/=n3) ) THEN
      CALL WarningMessage(2,'read_3d_complex_array:: sizes of the array are different from the ones used '// &
                            'to CALL this SUBROUTINE')
   END IF
   ALLOCATE (rr(n1,n2,n3))
   ALLOCATE (ri(n1,n2,n3))
   rr=Zero ; ri=Zero
   DO i=1, n1
      DO j=1, n2
         READ ( LU, FMT=*, IOSTAT=ierr ) ( rr(i,j,k), ri(i,j,k), k=1,n3 )
         IF (ierr /= 0) THEN
            CALL WarningMessage(2,'read_3d_complex_array:: Something went wrong reading the array.')
         END IF
         IF (DBG) WRITE (StdOut,*) 'read_3d_complex_array::  i,j =',i,j
         IF (DBG) FLUSH (StdOut)
      END DO
   END DO
   DO i=1, n1
      DO j=1, n2
         DO k=1, n3
            array(i,j,k)=CMPLX(rr(i,j,k), ri(i,j,k), wp)
         END DO
      END DO
   END DO
   DEALLOCATE (rr)
   DEALLOCATE (ri)
   RETURN
END SUBROUTINE read_3d_complex_array

!--------------------------------------------------------------------------------------------------!
SUBROUTINE read_4d_complex_array( LU, key, n1, n2, n3, n4, array )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (IN)          :: n1, n2, n3, n4
   COMPLEX (wp), INTENT (OUT)     :: array(n1,n2,n3,n4)
   REAL (wp), ALLOCATABLE        :: rr(:,:,:,:), ri(:,:,:,:)
   INTEGER                      :: i,j,k,l
   ierr=0
   array=ZeroC
   IF ( (n1<=0) .OR. (n2<=0) .OR. (n3<=0) .OR. (n4<=0) ) THEN
      CALL WarningMessage(1,'read_4d_complex_array:: nothing to read. Array size = 0.')
      IF (DBG) WRITE (StdOut,*) 'read_4d_complex_array::   n1 =',n1
      IF (DBG) WRITE (StdOut,*) 'read_4d_complex_array::   n2 =',n2
      IF (DBG) WRITE (StdOut,*) 'read_4d_complex_array::   n3 =',n3
      IF (DBG) WRITE (StdOut,*) 'read_4d_complex_array::   n4 =',n4
      IF (DBG) FLUSH (StdOut)
      RETURN
   END IF
   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )

   READ ( LU, FMT=*, IOSTAT=ierr ) i,j,k,l
   IF (ierr /= 0) THEN
      CALL WarningMessage(2,'read_4d_complex_array:: Something went wrong reading key'//trim(key))
   END IF
   IF (DBG) WRITE (StdOut,*) 'read_4d_complex_array:: key =',trim(key)
   IF (DBG) WRITE (StdOut,*) 'read_4d_complex_array::  n1 =',i
   IF (DBG) WRITE (StdOut,*) 'read_4d_complex_array::  n2 =',j
   IF (DBG) WRITE (StdOut,*) 'read_4d_complex_array::  n3 =',k
   IF (DBG) WRITE (StdOut,*) 'read_4d_complex_array::  n4 =',l
   IF ( (i/=n1) .OR. (j/=n2) .OR. (k/=n3) .OR. (l/=n4) ) THEN
      CALL WarningMessage(2,'read_4d_complex_array:: sizes of the array are different from the ones used '// &
                            'to CALL this SUBROUTINE')
   END IF
   ALLOCATE (rr(n1,n2,n3,n4))
   ALLOCATE (ri(n1,n2,n3,n4))
   rr=Zero; ri=Zero;
   DO i=1, n1
      DO j=1, n2
         DO k=1, n3
            READ ( LU, FMT=*, IOSTAT=ierr ) (rr(i,j,k,l), ri(i,j,k,l), l=1,n4)
            IF (ierr /= 0) THEN
               CALL WarningMessage(2,'read_4d_real_array:: Something went wrong reading the array.')
            END IF
            IF (DBG) WRITE (StdOut,*) 'read_4d_real_array::  i,j,k =',i,j,k
            IF (DBG) FLUSH (StdOut)
         END DO
      END DO
   END DO
   DO i=1, n1
      DO j=1, n2
         DO k=1, n3
            DO l=1, n4
               array(i,j,k,l)=CMPLX(rr(i,j,k,l), ri(i,j,k,l), wp)
            END DO
         END DO
      END DO
   END DO
   DEALLOCATE (rr)
   DEALLOCATE (ri)
   RETURN
END SUBROUTINE read_4d_complex_array
!--------------------------------------------------------------------------------------------------!


























!--------------------------------------------------------------------------------------------------!
!  WRITING SUBROUTINES
!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_INTEGER_scalar( LU, key, i )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (IN)          :: i
   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )
   ! if key is found, THEN rewrite the data
   IF ( ierr == 0 ) THEN
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) i
   ! if the keyword is not found, THEN append the new keyword and data
   ELSE IF (ierr /= 0 ) THEN
      WRITE ( LU, FMT='(A)', IOSTAT=ierr )
      WRITE ( LU, FMT='(A)', IOSTAT=ierr ) trim(key)
      IF (ierr /= 0) CALL WarningMessage(1,'write_INTEGER_scalar:: Something went wrong writing key'//trim(key))
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) i
      IF (ierr /= 0) CALL WarningMessage(1,'write_INTEGER_scalar:: Something went wrong writing data')
   END IF
   WRITE ( LU, FMT=*, IOSTAT=ierr )
   FLUSH ( LU )
   RETURN
END SUBROUTINE write_INTEGER_scalar

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_real_scalar( LU, key, r )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   REAL (wp), INTENT (IN)         :: r
   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )
   ! if key is found, THEN rewrite the data
   IF ( ierr == 0 ) THEN
      WRITE ( LU, FMT=FMTR, IOSTAT=ierr ) r
   ! if the keyword is not found, THEN append the new keyword and data
   ELSE IF (ierr /= 0 ) THEN
      WRITE ( LU, FMT='(A)', IOSTAT=ierr )
      WRITE ( LU, FMT='(A)', IOSTAT=ierr ) trim(key)
      IF (ierr /= 0) CALL WarningMessage(1,'write_real_scalar:: Something went wrong writing key'//trim(key))
      WRITE ( LU, FMT=FMTR, IOSTAT=ierr ) r
      IF (ierr /= 0) CALL WarningMessage(1,'write_real_scalar:: Something went wrong writing data')
   END IF
   WRITE ( LU, FMT=*, IOSTAT=ierr )
   FLUSH ( LU )
   RETURN
END SUBROUTINE write_real_scalar

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_complex_scalar( LU, key, c )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   COMPLEX (wp), INTENT (IN)      :: c
   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )
   ! if key is found, THEN rewrite the data
   IF ( ierr == 0 ) THEN
      WRITE ( LU, FMT=FMTC, IOSTAT=ierr ) c
   ! if the keyword is not found, THEN append the new keyword and data
   ELSE IF (ierr /= 0 ) THEN
      WRITE ( LU, FMT='(A)', IOSTAT=ierr )
      WRITE ( LU, FMT='(A)', IOSTAT=ierr ) trim(key)
      IF (ierr /= 0) CALL WarningMessage(1,'write_complex_scalar:: Something went wrong writing key'//trim(key))
      WRITE ( LU, FMT=FMTC, IOSTAT=ierr ) c
      IF (ierr /= 0) CALL WarningMessage(1,'write_complex_scalar:: Something went wrong writing data')
   END IF
   WRITE ( LU, FMT=*, IOSTAT=ierr )
   FLUSH ( LU )
   RETURN
END SUBROUTINE write_complex_scalar

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_string( LU, key, s )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   CHARACTER (LEN=*), INTENT (IN) :: s
   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )
   ! if key is found, THEN rewrite the data
   IF ( ierr == 0 ) THEN
      WRITE ( LU, FMT='(A   )', IOSTAT=ierr ) trim(s)
   ! if the keyword is not found, THEN append the new keyword and data
   ELSE IF (ierr /= 0 ) THEN
      WRITE ( LU, FMT='(A)', IOSTAT=ierr )
      WRITE ( LU, FMT='(A   )', IOSTAT=ierr ) trim(key)
      IF (ierr /= 0) CALL WarningMessage(1,'write_string:: Something went wrong writing key'//trim(key))
      WRITE ( LU, FMT='(100A)', IOSTAT=ierr ) trim(s)
      IF (ierr /= 0) CALL WarningMessage(1,'write_string:: Something went wrong writing data')
   END IF
   WRITE ( LU, FMT=*, IOSTAT=ierr )
   FLUSH ( LU )
   RETURN
END SUBROUTINE write_string

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_1d_INTEGER_array( LU, key, n, array )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (IN)          :: n
   INTEGER, INTENT (IN)          :: array(n)
   INTEGER                      :: i
   ierr=0
   IF (n<=0) THEN
      CALL WarningMessage(1,'write_1d_INTEGER_array:: nothing to write. Array size = 0.')
      RETURN
   END IF
   IF (SUM(ABS(array(1:n)))==0) THEN
      CALL WarningMessage(1,'write_1d_INTEGER_array:: all array elements are zero = 0.')
   END IF

   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )

   ! if key is found, THEN rewrite the data
   IF ( ierr == 0 ) THEN
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) n
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) (array(i),i=1,n)
   ! if the keyword is not found, THEN append the new keyword and data
   ELSE IF (ierr /= 0 ) THEN
      WRITE ( LU, FMT='(A)', IOSTAT=ierr )
      WRITE ( LU, FMT='(A)', IOSTAT=ierr ) trim(key)
      IF (ierr /= 0) CALL WarningMessage(1,'write_1d_INTEGER_array:: Something went wrong writing key'//trim(key))
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) n
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) (array(i),i=1,n)
      IF (ierr /= 0) CALL WarningMessage(1,'write_1d_INTEGER_array:: Something went wrong writing data')
   END IF
   WRITE ( LU, FMT=*, IOSTAT=ierr )
   FLUSH ( LU )
   RETURN
END SUBROUTINE write_1d_INTEGER_array

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_2d_INTEGER_array( LU, key, n1, n2, array )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (IN)          :: n1, n2
   INTEGER, INTENT (IN)          :: array(n1,n2)
   INTEGER                      :: i, j
   ierr=0
   IF ( (n1<=0) .OR. (n2<=0) ) THEN
      CALL WarningMessage(1,'write_2d_INTEGER_array:: nothing to write. Array size = 0.')
      RETURN
   END IF
   IF (SUM(ABS(array(1:n1,1:n2)))==0) THEN
      CALL WarningMessage(1,'write_2d_INTEGER_array:: all array elements are zero = 0.')
   END IF

   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )

   ! if key is found, THEN rewrite the data
   IF ( ierr == 0 ) THEN
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) n1, n2
      DO i=1,n1
         WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) (array(i,j),j=1,n2)
         IF (ierr /= 0) THEN
            CALL WarningMessage(2,'write_2d_INTEGER_array:: Something went wrong writing the array.')
         END IF
         IF (DBG) WRITE (StdOut,*) 'write_2d_INTEGER_array::  i =',i
         IF (DBG) FLUSH (StdOut)
      END DO
   ! if the keyword is not found, THEN append the new keyword and data
   ELSE IF (ierr /= 0 ) THEN
      WRITE ( LU, FMT='(A)', IOSTAT=ierr )
      WRITE ( LU, FMT='(A)', IOSTAT=ierr ) trim(key)
      IF (ierr /= 0) CALL WarningMessage(1,'write_2d_INTEGER_array:: Something went wrong writing key'//trim(key))
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) n1, n2
      DO i=1,n1
         WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) (array(i,j),j=1,n2)
         IF (ierr /= 0) THEN
            CALL WarningMessage(2,'write_2d_INTEGER_array:: Something went wrong writing data.')
         END IF
         IF (DBG) WRITE (StdOut,*) 'write_2d_INTEGER_array::  i =',i
         IF (DBG) FLUSH (StdOut)
      END DO
   END IF
   WRITE ( LU, FMT=*, IOSTAT=ierr )
   FLUSH ( LU )
   RETURN
END SUBROUTINE write_2d_INTEGER_array

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_3d_INTEGER_array( LU, key, n1, n2, n3, array )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (IN)          :: n1, n2, n3
   INTEGER, INTENT (IN)          :: array(n1,n2,n3)
   INTEGER                      :: i, j, k
   ierr=0
   IF ( (n1<=0) .OR. (n2<=0) .OR. (n3<=0) ) THEN
      CALL WarningMessage(1,'write_3d_INTEGER_array:: nothing to write. Array size = 0.')
      RETURN
   END IF
   IF (SUM(ABS(array(1:n1,1:n2,1:n3)))==0) THEN
      CALL WarningMessage(1,'write_3d_INTEGER_array:: all array elements are zero = 0.')
   END IF

   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )

   ! if key is found, THEN rewrite the data
   IF ( ierr == 0 ) THEN
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) n1, n2, n3
      DO i=1,n1
         DO j=1,n2
            WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) (array(i,j,k),k=1,n3)
            IF (ierr /= 0) THEN
               CALL WarningMessage(2,'write_3d_INTEGER_array:: Something went wrong writing the array.')
            END IF
            IF (DBG) WRITE (StdOut,*) 'write_3d_INTEGER_array::  i,j =',i,j
            IF (DBG) FLUSH (StdOut)
         END DO
      END DO
   ! if the keyword is not found, THEN append the new keyword and data
   ELSE IF (ierr /= 0 ) THEN
      WRITE ( LU, FMT='(A)', IOSTAT=ierr )
      WRITE ( LU, FMT='(A)', IOSTAT=ierr ) trim(key)
      IF (ierr /= 0) CALL WarningMessage(1,'write_3d_INTEGER_array:: Something went wrong writing key'//trim(key))
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) n1, n2, n3
      DO i=1,n1
         DO j=1,n2
            WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) (array(i,j,k),k=1,n3)
            IF (ierr /= 0) THEN
               CALL WarningMessage(2,'write_3d_INTEGER_array:: Something went wrong writing the array.')
            END IF
            IF (DBG) WRITE (StdOut,*) 'write_3d_INTEGER_array::  i,j =',i,j
            IF (DBG) FLUSH (StdOut)
         END DO
      END DO
   END IF
   WRITE ( LU, FMT=*, IOSTAT=ierr )
   FLUSH ( LU )
   RETURN
END SUBROUTINE write_3d_INTEGER_array

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_4d_INTEGER_array( LU, key, n1, n2, n3, n4, array )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (IN)          :: n1, n2, n3, n4
   INTEGER, INTENT (IN)          :: array(n1,n2,n3,n4)
   INTEGER                      :: i, j, k, l
   ierr=0
   IF ( (n1<=0) .OR. (n2<=0) .OR. (n3<=0) .OR. (n4<=0) ) THEN
      CALL WarningMessage(1,'write_4d_INTEGER_array:: nothing to write. Array size = 0.')
      RETURN
   END IF
   IF (SUM(ABS(array(1:n1,1:n2,1:n3,1:n4)))==0) THEN
      CALL WarningMessage(1,'write_4d_INTEGER_array:: all array elements are zero = 0.')
   END IF

   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )

   ! if key is found, THEN rewrite the data
   IF ( ierr == 0 ) THEN
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) n1, n2, n3, n4
      DO i=1,n1
         DO j=1,n2
            DO k=1,n3
               WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) (array(i,j,k,l),l=1,n4)
               IF (ierr /= 0) THEN
                  CALL WarningMessage(2,'write_4d_INTEGER_array:: Something went wrong reading the array.')
               END IF
               IF (DBG) WRITE (StdOut,*) 'write_4d_INTEGER_array::  i,j,k =',i,j,k
               IF (DBG) FLUSH (StdOut)
            END DO
         END DO
      END DO
   ! if the keyword is not found, THEN append the new keyword and data
   ELSE IF (ierr /= 0 ) THEN
      WRITE ( LU, FMT='(A)', IOSTAT=ierr )
      WRITE ( LU, FMT='(A)', IOSTAT=ierr ) trim(key)
      IF (ierr /= 0) CALL WarningMessage(1,'write_4d_INTEGER_array:: Something went wrong writing key'//trim(key))
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) n1, n2, n3, n4
      DO i=1,n1
         DO j=1,n2
            DO k=1,n3
               WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) (array(i,j,k,l),l=1,n4)
               IF (ierr /= 0) THEN
                  CALL WarningMessage(2,'write_4d_INTEGER_array:: Something went wrong writting the array.')
               END IF
               IF (DBG) WRITE (StdOut,*) 'write_4d_INTEGER_array::  i,j,k =',i,j,k
               IF (DBG) FLUSH (StdOut)
            END DO
         END DO
      END DO
   END IF
   WRITE ( LU, FMT=*, IOSTAT=ierr )
   FLUSH ( LU )
   RETURN
END SUBROUTINE write_4d_INTEGER_array

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_1d_real_array( LU, key, n, array )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (IN)          :: n
   REAL (wp), INTENT (IN)         :: array(n)
   INTEGER                      :: i
   REAL (wp), EXTERNAL           :: dnrm2_
   ierr=0
   IF (n<=0) THEN
      CALL WarningMessage(1,'write_1d_real_array:: nothing to write. Array size = 0.')
      RETURN
   END IF
   IF (dnrm2_(n,array,1)<=MINIMAL_REAL) THEN
      CALL WarningMessage(1,'write_1d_real_array:: all array elements are zero = 0.0')
   END IF

   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )

   ! if key is found, THEN rewrite the data
   IF ( ierr == 0 ) THEN
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) n
      WRITE ( LU, FMT=FMTR, IOSTAT=ierr ) (array(i),i=1,n)
   ! if the keyword is not found, THEN append the new keyword and data
   ELSE IF (ierr /= 0 ) THEN
      WRITE ( LU, FMT='(A)', IOSTAT=ierr )
      WRITE ( LU, FMT='(A)', IOSTAT=ierr ) trim(key)
      IF (ierr /= 0) CALL WarningMessage(1,'write_1d_real_array:: Something went wrong writing key'//trim(key))
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) n
      WRITE ( LU, FMT=FMTR, IOSTAT=ierr ) (array(i),i=1,n)
      IF (ierr /= 0) CALL WarningMessage(1,'write_1d_real_array:: Something went wrong writing data')
   END IF
   WRITE ( LU, FMT=*, IOSTAT=ierr )
   FLUSH ( LU )
   RETURN
END SUBROUTINE write_1d_real_array

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_2d_real_array( LU, key, n1, n2, array )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (IN)          :: n1, n2
   REAL (wp), INTENT (IN)         :: array(n1,n2)
   INTEGER                      :: i, j
   REAL (wp), EXTERNAL           :: dnrm2_
   ierr=0
   IF ( (n1<=0) .OR. (n2<=0) ) THEN
      CALL WarningMessage(1,'write_2d_real_array:: nothing to write. Array size = 0.')
      RETURN
   END IF
   IF (dnrm2_(n1*n2,array,1)<=MINIMAL_REAL) THEN
      CALL WarningMessage(1,'write_2d_real_array:: all array elements are zero = 0.')
   END IF

   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )

   ! if key is found, THEN rewrite the data
   IF ( ierr == 0 ) THEN
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) n1, n2
      DO i=1,n1
         WRITE ( LU, FMT=FMTR, IOSTAT=ierr ) (array(i,j),j=1,n2)
         IF (ierr /= 0) THEN
            CALL WarningMessage(2,'write_2d_real_array:: Something went wrong writing the array.')
         END IF
         IF (DBG) WRITE (StdOut,*) 'write_2d_real_array::  i =',i
         IF (DBG) FLUSH (StdOut)
      END DO
   ! if the keyword is not found, THEN append the new keyword and data
   ELSE IF (ierr /= 0 ) THEN
      WRITE ( LU, FMT='(A)', IOSTAT=ierr )
      WRITE ( LU, FMT='(A)', IOSTAT=ierr ) trim(key)
      IF (ierr /= 0) CALL WarningMessage(1,'write_2d_real_array:: Something went wrong writing key'//trim(key))
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) n1, n2
      DO i=1,n1
         WRITE ( LU, FMT=FMTR, IOSTAT=ierr ) (array(i,j),j=1,n2)
         IF (ierr /= 0) THEN
            CALL WarningMessage(2,'write_2d_real_array:: Something went wrong writing data.')
         END IF
         IF (DBG) WRITE (StdOut,*) 'write_2d_real_array::  i =',i
         IF (DBG) FLUSH (StdOut)
      END DO
   END IF
   WRITE ( LU, FMT=*, IOSTAT=ierr )
   FLUSH ( LU )
   RETURN
END SUBROUTINE write_2d_real_array

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_3d_real_array( LU, key, n1, n2, n3, array )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (IN)          :: n1, n2, n3
   REAL (wp), INTENT (IN)         :: array(n1,n2,n3)
   INTEGER                      :: i, j, k
   REAL (wp), EXTERNAL           :: dnrm2_
   ierr=0
   IF ( (n1<=0) .OR. (n2<=0) .OR. (n3<=0) ) THEN
      CALL WarningMessage(1,'write_3d_real_array:: nothing to write. Array size = 0.')
      RETURN
   END IF
   IF (dnrm2_(n1*n2*n3,array,1) <= MINIMAL_REAL) THEN
      CALL WarningMessage(1,'write_3d_real_array:: all array elements are zero = 0.')
   END IF

   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )

   ! if key is found, THEN rewrite the data
   IF ( ierr == 0 ) THEN
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) n1, n2, n3
      DO i=1,n1
         DO j=1,n2
            WRITE ( LU, FMT=FMTR, IOSTAT=ierr ) (array(i,j,k),k=1,n3)
            IF (ierr /= 0) THEN
               CALL WarningMessage(2,'write_3d_real_array:: Something went wrong writing the array.')
            END IF
            IF (DBG) WRITE (StdOut,*) 'write_3d_real_array::  i,j =',i,j
            IF (DBG) FLUSH (StdOut)
         END DO
      END DO
   ! if the keyword is not found, THEN append the new keyword and data
   ELSE IF (ierr /= 0 ) THEN
      WRITE ( LU, FMT='(A)', IOSTAT=ierr )
      WRITE ( LU, FMT='(A)', IOSTAT=ierr ) trim(key)
      IF (ierr /= 0) CALL WarningMessage(1,'write_3d_real_array:: Something went wrong writing key'//trim(key))
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) n1, n2, n3
      DO i=1,n1
         DO j=1,n2
            WRITE ( LU, FMT=FMTR, IOSTAT=ierr ) (array(i,j,k),k=1,n3)
            IF (ierr /= 0) THEN
               CALL WarningMessage(2,'write_3d_real_array:: Something went wrong writing the array.')
            END IF
            IF (DBG) WRITE (StdOut,*) 'write_3d_real_array::  i,j =',i,j
            IF (DBG) FLUSH (StdOut)
         END DO
      END DO
   END IF
   WRITE ( LU, FMT=*, IOSTAT=ierr )
   FLUSH ( LU )
   RETURN
END SUBROUTINE write_3d_real_array

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_4d_real_array( LU, key, n1, n2, n3, n4, array )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (IN)          :: n1, n2, n3, n4
   REAL (wp), INTENT (IN)         :: array(n1,n2,n3,n4)
   INTEGER                      :: i, j, k, l
   REAL (wp), EXTERNAL           :: dnrm2_
   ierr=0
   IF ( (n1<=0) .OR. (n2<=0) .OR. (n3<=0) .OR. (n4<=0) ) THEN
      CALL WarningMessage(1,'write_4d_real_array:: nothing to write. Array size = 0.')
      RETURN
   END IF
   IF (dnrm2_(n1*n2*n3*n4,array,1) <= MINIMAL_REAL) THEN
      CALL WarningMessage(1,'write_4d_real_array:: all array elements are zero = 0.')
   END IF

   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )

   ! if key is found, THEN rewrite the data
   IF ( ierr == 0 ) THEN
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) n1, n2, n3, n4
      DO i=1,n1
         DO j=1,n2
            DO k=1,n3
               WRITE ( LU, FMT=FMTR, IOSTAT=ierr ) (array(i,j,k,l),l=1,n4)
               IF (ierr /= 0) THEN
                  CALL WarningMessage(2,'write_4d_real_array:: Something went wrong reading the array.')
               END IF
               IF (DBG) WRITE (StdOut,*) 'write_4d_real_array::  i,j,k =',i,j,k
               IF (DBG) FLUSH (StdOut)
            END DO
         END DO
      END DO
   ! if the keyword is not found, THEN append the new keyword and data
   ELSE IF (ierr /= 0 ) THEN
      WRITE ( LU, FMT='(A)', IOSTAT=ierr )
      WRITE ( LU, FMT='(A)', IOSTAT=ierr ) trim(key)
      IF (ierr /= 0) CALL WarningMessage(1,'write_4d_real_array:: Something went wrong writing key'//trim(key))
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) n1, n2, n3, n4
      DO i=1,n1
         DO j=1,n2
            DO k=1,n3
               WRITE ( LU, FMT=FMTR, IOSTAT=ierr ) (array(i,j,k,l),l=1,n4)
               IF (ierr /= 0) THEN
                  CALL WarningMessage(2,'write_4d_real_array:: Something went wrong writting the array.')
               END IF
               IF (DBG) WRITE (StdOut,*) 'write_4d_real_array::  i,j,k =',i,j,k
               IF (DBG) FLUSH (StdOut)
            END DO
         END DO
      END DO
   END IF
   WRITE ( LU, FMT=*, IOSTAT=ierr )
   FLUSH ( LU )
   RETURN
END SUBROUTINE write_4d_real_array

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_1d_complex_array( LU, key, n, array )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (IN)          :: n
   COMPLEX (wp), INTENT (IN)      :: array(n)
   INTEGER                      :: i
   ierr=0
   IF (n<=0) THEN
      CALL WarningMessage(1,'write_1d_complex_array:: nothing to write. Array size = 0.')
      RETURN
   END IF

   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )

   ! if key is found, THEN rewrite the data
   IF ( ierr == 0 ) THEN
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) n
      WRITE ( LU, FMT=FMTC, IOSTAT=ierr ) (array(i),i=1,n)
   ! if the keyword is not found, THEN append the new keyword and data
   ELSE IF (ierr /= 0 ) THEN
      WRITE ( LU, FMT='(A)', IOSTAT=ierr )
      WRITE ( LU, FMT='(A)', IOSTAT=ierr ) trim(key)
      IF (ierr /= 0) CALL WarningMessage(1,'write_1d_complex_array:: Something went wrong writing key'//trim(key))
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) n
      WRITE ( LU, FMT=FMTC, IOSTAT=ierr ) (array(i),i=1,n)
      IF (ierr /= 0) CALL WarningMessage(1,'write_1d_complex_array:: Something went wrong writing data')
   END IF
   WRITE ( LU, FMT=*, IOSTAT=ierr )
   FLUSH ( LU )
   RETURN
END SUBROUTINE write_1d_complex_array

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_2d_complex_array( LU, key, n1, n2, array )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (IN)          :: n1, n2
   COMPLEX (wp), INTENT (IN)      :: array(n1,n2)
   INTEGER                      :: i, j
   ierr=0
   IF ( (n1<=0) .OR. (n2<=0) ) THEN
      CALL WarningMessage(1,'write_2d_complex_array:: nothing to write. Array size = 0.')
      RETURN
   END IF

   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )

   ! if key is found, THEN rewrite the data
   IF ( ierr == 0 ) THEN
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) n1, n2
      DO i=1,n1
         WRITE ( LU, FMT=FMTC, IOSTAT=ierr ) (array(i,j),j=1,n2)
         IF (ierr /= 0) THEN
            CALL WarningMessage(2,'write_2d_complex_array:: Something went wrong writing the array.')
         END IF
         IF (DBG) WRITE (StdOut,*) 'write_2d_complex_array::  i =',i
         IF (DBG) FLUSH (StdOut)
      END DO
   ! if the keyword is not found, THEN append the new keyword and data
   ELSE IF (ierr /= 0 ) THEN
      WRITE ( LU, FMT='(A)', IOSTAT=ierr )
      WRITE ( LU, FMT='(A)', IOSTAT=ierr ) trim(key)
      IF (ierr /= 0) CALL WarningMessage(1,'write_2d_complex_array:: Something went wrong writing key'//trim(key))
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) n1, n2
      DO i=1,n1
         WRITE ( LU, FMT=FMTC, IOSTAT=ierr ) (array(i,j),j=1,n2)
         IF (ierr /= 0) THEN
            CALL WarningMessage(2,'write_2d_complex_array:: Something went wrong writing data.')
         END IF
         IF (DBG) WRITE (StdOut,*) 'write_2d_complex_array::  i =',i
         IF (DBG) FLUSH (StdOut)
      END DO
   END IF
   WRITE ( LU, FMT=*, IOSTAT=ierr )
   FLUSH ( LU )
   RETURN
END SUBROUTINE write_2d_complex_array

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_3d_complex_array( LU, key, n1, n2, n3, array )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (IN)          :: n1, n2, n3
   COMPLEX (wp), INTENT (IN)      :: array(n1,n2,n3)
   INTEGER                      :: i, j, k
   ierr=0
   IF ( (n1<=0) .OR. (n2<=0) .OR. (n3<=0) ) THEN
      CALL WarningMessage(1,'write_3d_complex_array:: nothing to write. Array size = 0.')
      RETURN
   END IF

   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )

   ! if key is found, THEN rewrite the data
   IF ( ierr == 0 ) THEN
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) n1, n2, n3
      DO i=1,n1
         DO j=1,n2
            WRITE ( LU, FMT=FMTC, IOSTAT=ierr ) (array(i,j,k),k=1,n3)
            IF (ierr /= 0) THEN
               CALL WarningMessage(2,'write_3d_complex_array:: Something went wrong writing the array.')
            END IF
            IF (DBG) WRITE (StdOut,*) 'write_3d_complex_array::  i,j =',i,j
            IF (DBG) FLUSH (StdOut)
         END DO
      END DO
   ! if the keyword is not found, THEN append the new keyword and data
   ELSE IF (ierr /= 0 ) THEN
      WRITE ( LU, FMT='(A)', IOSTAT=ierr )
      WRITE ( LU, FMT='(A)', IOSTAT=ierr ) trim(key)
      IF (ierr /= 0) CALL WarningMessage(1,'write_3d_complex_array:: Something went wrong writing key'//trim(key))
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) n1, n2, n3
      DO i=1,n1
         DO j=1,n2
            WRITE ( LU, FMT=FMTC, IOSTAT=ierr ) (array(i,j,k),k=1,n3)
            IF (ierr /= 0) THEN
               CALL WarningMessage(2,'write_3d_complex_array:: Something went wrong writing the array.')
            END IF
            IF (DBG) WRITE (StdOut,*) 'write_3d_complex_array::  i,j =',i,j
            IF (DBG) FLUSH (StdOut)
         END DO
      END DO
   END IF
   WRITE ( LU, FMT=*, IOSTAT=ierr )
   FLUSH ( LU )
   RETURN
END SUBROUTINE write_3d_complex_array

!--------------------------------------------------------------------------------------------------!
SUBROUTINE write_4d_complex_array( LU, key, n1, n2, n3, n4, array )
   IMPLICIT NONE
   INTEGER, INTENT (IN)          :: LU
   CHARACTER (LEN=*), INTENT (IN) :: key
   INTEGER, INTENT (IN)          :: n1, n2, n3, n4
   COMPLEX (wp), INTENT (IN)      :: array(n1,n2,n3,n4)
   INTEGER                      :: i, j, k, l
   ierr=0
   IF ( (n1<=0) .OR. (n2<=0) .OR. (n3<=0) .OR. (n4<=0) ) THEN
      CALL WarningMessage(1,'write_4d_complex_array:: nothing to write. Array size = 0.')
      RETURN
   END IF

   REWIND ( LU )
   CALL file_advance_to_string( LU, key, line )

   ! if key is found, THEN rewrite the data
   IF ( ierr == 0 ) THEN
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) n1, n2, n3, n4
      DO i=1,n1
         DO j=1,n2
            DO k=1,n3
               WRITE ( LU, FMT=FMTC, IOSTAT=ierr ) (array(i,j,k,l),l=1,n4)
               IF (ierr /= 0) THEN
                  CALL WarningMessage(2,'write_4d_complex_array:: Something went wrong reading the array.')
               END IF
               IF (DBG) WRITE (StdOut,*) 'write_4d_complex_array::  i,j,k =',i,j,k
               IF (DBG) FLUSH (StdOut)
            END DO
         END DO
      END DO
   ! if the keyword is not found, THEN append the new keyword and data
   ELSE IF (ierr /= 0 ) THEN
      WRITE ( LU, FMT='(A)', IOSTAT=ierr )
      WRITE ( LU, FMT='(A)', IOSTAT=ierr ) trim(key)
      IF (ierr /= 0) CALL WarningMessage(1,'write_4d_complex_array:: Something went wrong writing key'//trim(key))
      WRITE ( LU, FMT=FMTI, IOSTAT=ierr ) n1, n2, n3, n4
      DO i=1,n1
         DO j=1,n2
            DO k=1,n3
               WRITE ( LU, FMT=FMTC, IOSTAT=ierr ) (array(i,j,k,l),l=1,n4)
               IF (ierr /= 0) THEN
                  CALL WarningMessage(2,'write_4d_complex_array:: Something went wrong writting the array.')
               END IF
               IF (DBG) WRITE (StdOut,*) 'write_4d_complex_array::  i,j,k =',i,j,k
               IF (DBG) FLUSH (StdOut)
            END DO
         END DO
      END DO
   END IF
   WRITE ( LU, FMT=*, IOSTAT=ierr )
   FLUSH ( LU )
   RETURN
END SUBROUTINE write_4d_complex_array

END MODULE io_data
