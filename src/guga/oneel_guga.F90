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

subroutine ONEEL_GUGA()

use Definitions, only: iwp

implicit none
#include "SysDef.fh"
#include "real_guga.fh"
#include "integ.fh"
#include "files_addr.fh"
integer(kind=iwp) :: I, ISTOP, IT1, IT2, ITT, ITYP, J, K, KJL, KJS, KM, NI, NK, NSI, NSK

IOUT = 0
NMAT = 0
ITYP = 0
do NK=1,LN
  do NI=1,NK
    K = ICH(NK)
    I = ICH(NI)
    if (K > I) GO TO 19
    K = ICH(NI)
    I = ICH(NK)
19  NSK = NSM(K)
    KJS = IJ(K+1)+1
    KJL = IJ(K)
    NSI = NSM(I)
    if (NSI /= NSK) GO TO 20
    IOUT = IOUT+1
    ICOP1(IOUT) = 0
    if (IOUT < NBUF) GO TO 460
    ICOP1(nCOP+1) = NBUF
    call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
    call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
    NMAT = NMAT+NBUF
    IOUT = 0
460 IOUT = IOUT+1
    !ICOP1(IOUT) = I+2**10*K
    ICOP1(IOUT) = ior(I,ishft(K,10))
    if (IOUT < NBUF) GO TO 21
    ICOP1(nCOP+1) = NBUF
    call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
    call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
    NMAT = NMAT+NBUF
    IOUT = 0
21  if (I == K) GO TO 20
    do ITT=1,ILIM
      IT1 = (ITT-1)*MXVERT
      IT2 = IT1
      do J=KJS,KJL
        IWAY(K) = 1
32      KM = K
        J2(KM+1) = J
        J1(KM+1) = J
        call LOOP1(KM,ISTOP,IT1,IT2)
        if (ISTOP == 1) GO TO 30
41      KM = KM-1
        IWAY(KM) = 1
        if (KM == I) GO TO 51
42      call LOOP5(KM,ISTOP,IT1,IT2)
        if (ISTOP == 0) GO TO 41
        KM = KM+1
        if (KM == K) GO TO 32
        GO TO 42
51      IWAY(I) = 1
52      KM = I
        call LOOP3(KM,ISTOP,IT1,IT2)
        if (ISTOP == 1) GO TO 53
        call COMP(I,J,ITYP,I,IT1,IT2)
        GO TO 52
53      KM = KM+1
        if (KM == K) GO TO 32
        GO TO 42
30      continue
      end do
    end do
20  continue
  end do
end do
ICOP1(nCOP+1) = IOUT
call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
NMAT = NMAT+IOUT
ICOP1(nCOP+1) = -1
call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
write(IW,100) NMAT
100 format(/6X,'COEFFICIENTS FOR IJ',I11)

return

end subroutine ONEEL_GUGA
