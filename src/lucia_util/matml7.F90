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
integer(kind=iwp) :: NCROW, NCCOL, NAROW, NACOL, NBROW, NBCOL, ITRNSP
real(kind=wp) :: C(NCROW,NCCOL), A(NAROW,NACOL), B(NBROW,NBCOL), FACTORC, FACTORAB
integer(kind=iwp) :: I, IZERO, J, K, LDA, LDB, LDC
real(kind=wp) :: AKI, B1J, BJ1, BJK, BKJ, T

if ((NAROW == 0) .or. (NACOL == 0) .or. (NBROW == 0) .or. (NBCOL == 0) .or. (NCROW == 0) .or. (NCCOL == 0)) then
  IZERO = 1
else
  IZERO = 0
end if

if ((IZERO == 1) .and. (NCROW*NCCOL /= 0)) then
  if (FACTORC /= Zero) then
    call SCALVE(C,FACTORC,NCROW*NCCOL)
  else
    call SETVEC(C,Zero,NCROW*NCCOL)
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
  ! Use Jeppes version ( it should be working )
  if (ITRNSP == 0) then
    ! ======
    ! C=A*B
    ! ======

    !call SCALVE(C,FACTORC,NCROW*NCCOL)
    do J=1,NCCOL
      ! Initialize with FACTORC*C(I,J) + FACTORAB*A(I,1)*B(1,J)
      if (NBROW >= 1) then
        B1J = FACTORAB*B(1,J)
        do I=1,NCROW
          C(I,J) = FACTORC*C(I,J)+B1J*A(I,1)
        end do
      end if
      ! and the major part
      do K=2,NBROW
        BKJ = FACTORAB*B(K,J)
        do I=1,NCROW
          C(I,J) = C(I,J)+BKJ*A(I,K)
        end do
      end do
    end do

  end if
  if (ITRNSP == 1) then

    ! =========
    ! C=A(T)*B
    ! =========

    do J=1,NCCOL
      do I=1,NCROW
        T = Zero
        do K=1,NBROW
          T = T+A(K,I)*B(K,J)
        end do
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
        do I=1,NCROW
          C(I,J) = FACTORC*C(I,J)+BJ1*A(I,1)
        end do
      end if
      ! And the rest
      do K=2,NBCOL
        BJK = FACTORAB*B(J,K)
        do I=1,NCROW
          C(I,J) = C(I,J)+BJK*A(I,K)
        end do
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
  call SCALVE(C,FACTORC,NCROW*NCCOL)
  do I=1,NCROW
    do K=1,NAROW
      AKI = FACTORAB*A(K,I)
      do J=1,NBROW
        C(I,J) = C(I,J)+AKI*B(J,K)
      end do
    end do
  end do
end if

!if (NTEST == 0) then
!  write(u6,*)
!  write(u6,*) ' C MATRIX FROM MATML7'
!  write(u6,*)
!  call WRTMAT(C,NCROW,NCCOL,NCROW,NCCOL)
!end if

end subroutine MATML7
