!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2003, Jeppe Olsen                                      *
!***********************************************************************

subroutine MATML7(C,A,B,NCROW,NCCOL,NAROW,NACOL,NBROW,NBCOL,FACTORC,FACTORAB,ITRNSP)
! MULTIPLY A AND B TO GIVE C
!
!     C =  FACTORC*C + FACTORAB* A * B       FOR ITRNSP = 0
!
!     C =  FACTORC*C + FACTORAB* A(T) * B    FOR ITRNSP = 1
!
!     C =  FACTORC*C + FACTORAB* A * B(T)    FOR ITRNSP = 2
!
!     C =  FACTORC*C + FACTORAB* A(T) * B(T) FOR ITRNSP = 3
!
! Warning ITRNSP = 3 should only be used for small matrices,
! as this path involves notunit strides in inner loops.
! As Lasse points out, it is better to calculate C(T) = BA
! and then transpose C
!
! JEPPE OLSEN,
!
! ITRNSP = 3 added, march 2003
!
! Notice : If the summation index has dimension zero nothing
!          is performed

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NCROW, NCCOL, NAROW, NACOL, NBROW, NBCOL, ITRNSP
real(kind=wp), intent(inout) :: C(NCROW,NCCOL)
real(kind=wp), intent(in) :: A(NAROW,NACOL), B(NBROW,NBCOL), FACTORC, FACTORAB
integer(kind=iwp) :: I, IZERO, J, K, LDA, LDB, LDC
real(kind=wp) :: AKI, B1J, BJ1, BJK, BKJ, T

if ((NAROW == 0) .or. (NACOL == 0) .or. (NBROW == 0) .or. (NBCOL == 0) .or. (NCROW == 0) .or. (NCCOL == 0)) then
  IZERO = 1
else
  IZERO = 0
end if

if ((IZERO == 1) .and. (NCROW*NCCOL /= 0)) then
  if (FACTORC /= Zero) then
    C(:,:) = FACTORC*C(:,:)
  else
    C(:,:) = Zero
  end if
end if

if (IZERO == 0) then
  ! DGEMM from CONVEX/ESSL  lib
  LDA = max(1,NAROW)
  LDB = max(1,NBROW)

  LDC = max(1,NCROW)
  if (ITRNSP == 0) then
    call DGEMM_('N','N',NAROW,NBCOL,NACOL,FACTORAB,A,LDA,B,LDB,FACTORC,C,LDC)
  else if (ITRNSP == 1) then
    call DGEMM_('T','N',NACOL,NBCOL,NAROW,FACTORAB,A,LDA,B,LDB,FACTORC,C,LDC)
  else if (ITRNSP == 2) then
    call DGEMM_('N','T',NAROW,NBROW,NACOL,FACTORAB,A,LDA,B,LDB,FACTORC,C,LDC)
  end if

else
  ! here if iZERO == 1
  ! Use Jeppe's version (it should be working)
  if (ITRNSP == 0) then
    ! ======
    ! C=A*B
    ! ======

    !C(:,:) = FACTORC*C(:,:)
    do J=1,NCCOL
      ! Initialize with FACTORC*C(I,J) + FACTORAB*A(I,1)*B(1,J)
      if (NBROW >= 1) then
        B1J = FACTORAB*B(1,J)
        C(:,J) = FACTORC*C(:,J)+B1J*A(1:NCROW,1)
      end if
      ! and the major part
      do K=2,NBROW
        BKJ = FACTORAB*B(K,J)
        C(:,J) = C(:,J)+BKJ*A(1:NCROW,K)
      end do
    end do

  end if
  if (ITRNSP == 1) then

    ! =========
    ! C=A(T)*B
    ! =========

    do J=1,NCCOL
      do I=1,NCROW
        T = sum(A(1:NBROW,I)*B(1:NBROW,J))
        C(I,J) = FACTORC*C(I,J)+FACTORAB*T
      end do
    end do
  end if

  if (ITRNSP == 2) then
    !===========
    ! C = A*B(T)
    !===========
    do J=1,NCCOL
      ! Initialization
      if (NBCOL >= 1) then
        BJ1 = FACTORAB*B(J,1)
        C(:,J) = FACTORC*C(:,J)+BJ1*A(1:NCROW,1)
      end if
      ! And the rest
      do K=2,NBCOL
        BJK = FACTORAB*B(J,K)
        C(:,J) = C(:,J)+BJK*A(1:NCROW,K)
      end do
    end do
  end if
  ! end of iZero ==1
end if

if (ITRNSP == 3) then
  !================
  ! C = A(T)*B(T)
  !================
  ! C(I,J) = FACTORC*C(I,J) + FACTORAB*sum(K) A(K,I)*B(J,K)
  C(:,:) = FACTORC*C(:,:)
  do I=1,NCROW
    do K=1,NAROW
      AKI = FACTORAB*A(K,I)
      C(I,1:NBROW) = C(I,1:NBROW)+AKI*B(:,K)
    end do
  end do
end if

!#ifdef _DEBUGPRINT_
!write(u6,*)
!write(u6,*) ' C MATRIX FROM MATML7'
!write(u6,*)
!call WRTMAT(C,NCROW,NCCOL,NCROW,NCCOL)
!#endif

end subroutine MATML7
