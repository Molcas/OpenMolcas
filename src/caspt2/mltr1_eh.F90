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

subroutine MLTR1_EH(IMLTOP,LST1,X,NAS1,NIS1,JXOFF,F,NFT,NFA,Y,NAS2,jYLo,jYHi)
! this routine is adapted to use chunks of a distributed array (lg_Y)
! for the H case.  The chunks span all rows (NAS2) and columns jYlo to
! jYHi. The array Y points to the beginning of a chunk.

use definitions, only: iwp, wp
use Sigma_data, only: NLST1, INCX3, VAL1

implicit none
integer(kind=iwp), intent(in) :: IMLTOP, NAS1, NIS1, JXOFF, NFT, NFA, NAS2, jYLo, jYHi
real(kind=wp), intent(inout) :: X(NAS1,NIS1), Y(NAS2,jYHi-jYLo+1)
real(kind=wp), intent(inout) :: F(NFT,NFA)
integer(kind=iwp), intent(in) :: LST1(4,NLST1)
integer(kind=iwp) NA, ILST, L1, L2, L3, L4, JX, I, J, NI
real(kind=wp) V, A
real(kind=wp), external :: DDot_

select case (IMLTOP)
  case (0)
    NA = INCX3/NAS1
    do ILST=1,NLST1
      L1 = LST1(1,ILST)
      L2 = LST1(2,ILST)
      L3 = LST1(3,ILST)
      L4 = LST1(4,ILST)
      V = VAL1(L4)
      JX = JXOFF+L1
      do I=1,NAS1
        A = V*F(I,L2)
        ! X(L1,p,q):= Add V*F(L2,p)*Y(L3,q), p=1..LEN1, q=1..LEN2
        call DAXPY_(jYHi-jYLo+1,A,Y(L3,1),NAS2,X(I,JX+NA*(jYLo-1)),INCX3)
      end do
    end do
  case (1)
    NA = INCX3/NAS1
    do ILST=1,NLST1
      L1 = LST1(1,ILST)
      L2 = LST1(2,ILST)
      L3 = LST1(3,ILST)
      L4 = LST1(4,ILST)
      V = VAL1(L4)
      JX = JXOFF+L1
      do J=jYLo,jYHi
        ! Y(L3,q):= Add V*F(L2,p)*X(L1,p,q), p=1..LEN1, q=1..LEN2
        Y(L3,J-jYLo+1) = Y(L3,J-jYLo+1)+V*DDOT_(NAS1,F(1,L2),1,X(1,JX+NA*(J-1)),1)
      end do
    end do
  case DEFAULT
    NI = INCX3/NAS1
    do ILST=1,NLST1
      L1 = LST1(1,ILST)
      L2 = LST1(2,ILST)
      L3 = LST1(3,ILST)
      L4 = LST1(4,ILST)
      V = VAL1(L4)
      ! F(L2,p) := Add V*X(L1,p,q)*Y(L3,q)
      JX = JXOFF+L1
      do J=jYLo,jYHi
        A = V*Y(L3,J-jYLo+1)
        call DAXPY_(NAS1,A,X(1,JX+(J-1)*NI),1,F(1,L2),1)
      end do
    end do
end select

!NFR1 = NFR1+2*NLST1*LEN1*LEN2

end subroutine MLTR1_EH
