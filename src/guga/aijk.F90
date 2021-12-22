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
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!***********************************************************************

subroutine AIJK(ITAI,L0,L1,L2,L3)

use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(_OUT_) :: ITAI(*)
integer(kind=iwp), intent(in) :: L0(*), L1(*), L2(*), L3(*)
#include "SysDef.fh"
#include "files_addr.fh"
#include "real_guga.fh"
#include "integ.fh"
integer(kind=iwp) :: I, II, IID, IND, IT1, IT2, ITT1, ITT2, ITURN, J, JJ, JJD, K, L, NI, NJ, NK

IOUT = 0
NMAT = 0
L = 0
do NI=1,LN
  do NJ=1,NI
    I = ICH(NI)
    J = ICH(NJ)
    if (I > J) GO TO 29
    I = ICH(NJ)
    J = ICH(NI)
29  do NK=1,LN
      IOUT = IOUT+1
      ICOP1(IOUT) = 0
      if (IOUT < NBUF) GO TO 460
      ICOP1(nCOP+1) = NBUF
      call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT = NMAT+NBUF
      IOUT = 0
460   K = ICH(NK)
      IOUT = IOUT+1
      !IND = I+2**10*J
      IND = ior(I,ishft(J,10))
      !ICOP1(IOUT) = IND+2**20*K
      ICOP1(IOUT) = ior(IND,ishft(K,20))
      if (IOUT < NBUF) GO TO 31
      ICOP1(nCOP+1) = NBUF
      call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT = NMAT+NBUF
      IOUT = 0
      ! DOUBLET-VALENCE INTERACTIONS
31    ITURN = 1
      ITT1 = 1
      ITT2 = 0
150   IT1 = ITT1*MXVERT
      IT2 = ITT2*MXVERT
      JJ = 0
      if (ITT1 /= 0) JJ = IRC(ITT1)
      JJD = 0
      if (ITT1 /= 0) JJD = JRC(ITT1)
      II = 0
      if (ITT2 /= 0) II = IRC(ITT2)
      IID = 0
      if (ITT2 /= 0) IID = JRC(ITT2)
      if (I /= K) GO TO 520
      if (J == K) GO TO 524
      call INT61(L,J,I,IT1,IT2,II,IID,JJ,JJD,ITURN,ITAI,L0,L1,L2,L3)
      call INT62(L,J,I,IT1,IT2,II,IID,JJ,JJD,ITURN,ITAI,L0,L1,L2,L3)
      GO TO 35
524   call INT9(L,K,I,IT1,IT2,II,IID,JJ,JJD,ITURN,ITAI,L0,L1,L2,L3)
      GO TO 35
520   if (J /= K) GO TO 535
      call INT4(L,K,I,IT1,IT2,II,IID,JJ,JJD,ITURN,ITAI,L0,L1,L2,L3)
      GO TO 35
535   if (I /= J) GO TO 546
      call INT8(L,K,I,IT1,IT2,II,IID,JJ,JJD,ITURN,ITAI,L0,L1,L2,L3)
      GO TO 35
546   if (K < I) GO TO 540
      call INT3(J,I,L,K,IT1,IT2,II,IID,JJ,JJD,ITURN,ITAI,L0,L1,L2,L3)
      GO TO 35
540   if (K < J) GO TO 545
      call INT2(L,K,J,I,IT1,IT2,II,IID,JJ,JJD,ITURN,ITAI,L0,L1,L2,L3)
      GO TO 35
545   call INT1(L,K,J,I,IT1,IT2,II,IID,JJ,JJD,ITURN,ITAI,L0,L1,L2,L3)
35    GO TO(101,102,103,104,105,30),ITURN
      ! TRIPLET-DOUBLET INTERACTIONS
101   ITURN = 2
      ITT1 = 2
      ITT2 = 1
      GO TO 150
      ! SINGLET-DOUBLET INTERACTIONS
102   ITURN = 3
      ITT1 = 3
      ITT2 = 1
      GO TO 150
      ! VALENCE-DOUBLET INTERACTIONS
103   ITURN = 4
      ITT1 = 0
      ITT2 = 1
      GO TO 150
      ! DOUBLET-TRIPLET INTERACTIONS
104   ITURN = 5
      ITT1 = 1
      ITT2 = 2
      GO TO 150
      ! DOUBLET-SINGLET INTERACTIONS
105   ITURN = 6
      ITT1 = 1
      ITT2 = 3
      GO TO 150
30    continue
    end do
  end do
end do
ICOP1(nCOP+1) = IOUT
call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
NMAT = NMAT+IOUT
ICOP1(nCOP+1) = -1
call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
write(IW,600) NMAT
600 format(/6X,'COEFFICIENTS FOR AIJK',I9)

return

end subroutine AIJK
