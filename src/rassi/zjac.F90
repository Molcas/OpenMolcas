!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine ZJAC(NDIMEN,ARRRE,ARRIM,LDV,VECRE,VECIM)

use Constants, only: Zero, One, Two, Four
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: NDIMEN, LDV
real(kind=wp) :: ARRRE(NDIMEN,NDIMEN), ARRIM(NDIMEN,NDIMEN), VECRE(LDV,*), VECIM(LDV,*)
integer(kind=iwp) :: I, IFERR, J, K, NR, NROT, NSWEEP
real(kind=wp) :: AAIJ, AIIJ, AIIK, AIJK, AIKI, AIKJ, ARII, ARIJ, ARIK, ARJJ, ARJK, ARKI, ARKJ, CS, DIFF, DUM, EIM, ERE, ERRIM, &
                 ERRRE, SBDMAX, SGN, SN, TN, VIKI, VIKJ, VNOLD, VNSUM, VRKI, VRKJ
integer(kind=iwp), parameter :: IFTEST = 0
real(kind=wp), parameter :: EPS = 1.0e-12_wp

! von Neumanns sum should be ever decreasing. Check this:
VNSUM = 1.0e99_wp
! Max sub-diagonal element:
NROT = 0
NSWEEP = 0
do
  NSWEEP = NSWEEP+1
  NR = 0
  SBDMAX = EPS
  VNOLD = VNSUM
  VNSUM = Zero
  do I=2,NDIMEN
    do J=1,I-1
      ARII = ARRRE(I,I)
      ARJJ = ARRRE(J,J)
      ARIJ = ARRRE(I,J)
      AIIJ = ARRIM(I,J)
      AAIJ = sqrt(ARIJ**2+AIIJ**2)
      VNSUM = VNSUM+AAIJ**2
      SBDMAX = max(SBDMAX,AAIJ)
      if (Two*AAIJ < SBDMAX) cycle
      NR = NR+1
      ERE = ARIJ/AAIJ
      EIM = AIIJ/AAIJ
      DIFF = ARII-ARJJ
      SGN = One
      if (DIFF < Zero) then
        DIFF = -DIFF
        SGN = -SGN
      end if
      DUM = DIFF+sqrt(DIFF**2+Four*AAIJ**2)
      TN = Two*SGN*AAIJ/DUM
      CS = One/sqrt(One+TN**2)
      SN = CS*TN
      ! TN,CS,SN=TAN,COS AND SIN OF ROTATION ANGLE.
      ! Rotate vectors:
      do K=1,LDV
        VRKJ = CS*VECRE(K,J)-SN*(ERE*VECRE(K,I)-EIM*VECIM(K,I))
        VIKJ = CS*VECIM(K,J)-SN*(EIM*VECRE(K,I)+ERE*VECIM(K,I))
        VRKI = SN*VECRE(K,J)+CS*(ERE*VECRE(K,I)-EIM*VECIM(K,I))
        VIKI = SN*VECIM(K,J)+CS*(EIM*VECRE(K,I)+ERE*VECIM(K,I))
        VECRE(K,J) = VRKJ
        VECIM(K,J) = VIKJ
        VECRE(K,I) = VRKI
        VECIM(K,I) = VIKI
      end do
      ! Rotate matrix columns:
      do K=1,NDIMEN
        ARKJ = CS*ARRRE(K,J)-SN*(ERE*ARRRE(K,I)-EIM*ARRIM(K,I))
        AIKJ = CS*ARRIM(K,J)-SN*(EIM*ARRRE(K,I)+ERE*ARRIM(K,I))
        ARKI = SN*ARRRE(K,J)+CS*(ERE*ARRRE(K,I)-EIM*ARRIM(K,I))
        AIKI = SN*ARRIM(K,J)+CS*(EIM*ARRRE(K,I)+ERE*ARRIM(K,I))
        ARRRE(K,J) = ARKJ
        ARRIM(K,J) = AIKJ
        ARRRE(K,I) = ARKI
        ARRIM(K,I) = AIKI
      end do
      ! Rotate matrix rows:
      do K=1,NDIMEN
        ARJK = CS*ARRRE(J,K)-SN*(ERE*ARRRE(I,K)+EIM*ARRIM(I,K))
        AIJK = CS*ARRIM(J,K)-SN*(-EIM*ARRRE(I,K)+ERE*ARRIM(I,K))
        ARIK = SN*ARRRE(J,K)+CS*(ERE*ARRRE(I,K)+EIM*ARRIM(I,K))
        AIIK = SN*ARRIM(J,K)+CS*(-EIM*ARRRE(I,K)+ERE*ARRIM(I,K))
        ARRRE(J,K) = ARJK
        ARRIM(J,K) = AIJK
        ARRRE(I,K) = ARIK
        ARRIM(I,K) = AIIK
      end do
    end do
  end do
  NROT = NROT+NR

  if (IFTEST > 0) then
    ! --- CHECK IF DIVERGING (This should never happen):
    if (VNSUM >= VNOLD) then
      call WarningMessage(2,'ZJAC got increasing von-Neumann-sum.')
      write(u6,*) ' Panic exit from Jacobi iteration loop.'
      call Error()
    end if
    ! --- CHECK IF IDLE LOOPS (This should never happen):
    if ((NR == 0) .and. (SBDMAX > EPS)) then
      call WarningMessage(2,'ZJAC detected infinite idle loops.')
      write(u6,*) ' Panic exit from Jacobi iteration loop.'
      call Error()
    end if
  end if

  ! PAM 2008 If unchecked, then idle loop just generates warning and return:
  IFERR = 0
  if (IFTEST == 0) then
    if (VNSUM >= VNOLD) then
      call WarningMessage(1,'ZJAC got increasing von-Neumann-sum.')
      IFERR = 1
    end if
    if ((NR == 0) .and. (SBDMAX > EPS)) then
      call WarningMessage(1,'ZJAC detected infinite idle loops.')
      IFERR = 1
    end if
    if (IFERR /= 0) then
      write(u6,*) ' Probably, the convergence criteria of the ZJAC'
      write(u6,*) ' code are too strict. You may report this as a bug'
      write(u6,*) ' but since you probably simply hit the limits to'
      write(u6,*) ' accuracy caused by round-off, the program will'
      write(u6,*) ' continue.'
      write(u6,*) ' At this point, the largest subdiagonal element'
      write(u6,*) ' is SBDMAX=',SBDMAX
    end if
  end if

  ! --- CHECK IF CONVERGED:
  if ((IFERR /= 0) .or. (SBDMAX <= EPS)) exit
end do

!> order the eigenvalues in increasing sequence
call zorder(ndimen,ldv,vecre,vecim,arrre,1)

contains

subroutine Error()

  ! Jump here on error.
  call WarningMessage(2,'ZJAC abend diagnostics:')
  write(u6,*) ' Nr of sweeps      NSWEEP=',NSWEEP
  write(u6,*) ' Nr of 2-rotations   NROT=',NROT
  write(u6,*) ' Rotations this sweep  NR=',NR
  write(u6,*) ' von Neumann sum    VNSUM=',VNSUM
  write(u6,*) ' Previous value     VNOLD=',VNOLD
  write(u6,*) ' Largest subdiag   SBDMAX=',SBDMAX
  call HRMCHK(NDIMEN,ARRRE,ARRIM,ERRRE,ERRIM)
  write(u6,*) ' Antisymm. part of ARRRE:',ERRRE
  write(u6,*) ' Symmetric part of ARRIM:',ERRIM
  call ABEND()

end subroutine Error

end subroutine ZJAC
