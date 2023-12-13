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
!               2021, Ignacio Fdez. Galvan                             *
!***********************************************************************
! 2021: Remove GOTOs

subroutine ONEEL_GUGA()

use guga_global, only: IADD10, ICH, IJ, ILIM, IOUT, IWAY, J1, J2, LN, Lu_10, MXVERT, NBUF, NMAT, NSM
use guga_util_global, only: COP, ICOP1, nCOP
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: I, ISTOP, IT1, IT2, ITT, ITYP, J, K, KJL, KJS, KM, NI, NK, NSI, NSK
logical(kind=iwp) :: first

IOUT = 0
NMAT = 0
ITYP = 0
do NK=1,LN
  do NI=1,NK
    K = ICH(NK)
    I = ICH(NI)
    if (K <= I) then
      K = ICH(NI)
      I = ICH(NK)
    end if
    NSK = NSM(K)
    KJS = IJ(K+1)+1
    KJL = IJ(K)
    NSI = NSM(I)
    if (NSI /= NSK) cycle
    IOUT = IOUT+1
    ICOP1(IOUT) = 0
    if (IOUT >= NBUF) then
      ICOP1(nCOP+1) = NBUF
      call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT = NMAT+NBUF
      IOUT = 0
    end if
    IOUT = IOUT+1
    !ICOP1(IOUT) = I+2**10*K
    ICOP1(IOUT) = ior(I,ishft(K,10))
    if (IOUT >= NBUF) then
      ICOP1(nCOP+1) = NBUF
      call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT = NMAT+NBUF
      IOUT = 0
    end if
    if (I == K) cycle
    do ITT=1,ILIM
      IT1 = (ITT-1)*MXVERT
      IT2 = IT1
      do J=KJS,KJL
        IWAY(K) = 1
        do
          KM = K
          J2(KM+1) = J
          J1(KM+1) = J
          call LOOP1(KM,ISTOP,IT1,IT2)
          if (ISTOP == 1) exit
          first = .true.
          do
            if (first) then
              KM = KM-1
              IWAY(KM) = 1
              first = .false.
            end if
            if (KM /= I) then
              call LOOP5(KM,ISTOP,IT1,IT2)
              if (ISTOP == 0) then
                first = .true.
                cycle
              end if
            else
              IWAY(I) = 1
              do
                KM = I
                call LOOP3(KM,ISTOP,IT1,IT2)
                if (ISTOP == 1) exit
                call COMP(I,J,ITYP,I,IT1,IT2)
              end do
            end if
            KM = KM+1
            if (KM == K) exit
          end do
        end do
      end do
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
write(u6,100) NMAT

return

100 format(/6X,'COEFFICIENTS FOR IJ',I11)

end subroutine ONEEL_GUGA
