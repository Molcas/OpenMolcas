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
!               1986, Margareta R. A. Blomberg                         *
!***********************************************************************

subroutine IJIJ_CPF(JSY,HDIAG,FJI)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"
dimension JSY(*), HDIAG(*), FJI(*)
dimension HCOUT(nCOP)
! Statement function
JSYM(L) = JSUNP_CPF(JSY,L)

ICOUP = 0 ! dummy initialize
IVL = 0 ! dummy initialize
NSS = 0 ! dummy initialize

IADD25 = IAD25S
IAD27 = 0
if (IREF0 > nCOP) then
  write(6,*) 'IJIJ_CPF Error: IREF0>nCOP (See code.)'
end if
call dDAFILE(Lu_27,2,HDIAG,IRC(1),IAD27)

IFS = 0
TERM = 0.0d0
IVSAVE = 0
ICOUPS = 0
IOUT = 0
ICHK = 0
IADD10 = IAD10(3)

300 continue
! Read a new COP buffer:
301 continue
call dDAFILE(Lu_CIGuga,2,COP,nCOP,IADD10)
call iDAFILE(Lu_CIGuga,2,iCOP1,nCOP+1,IADD10)
LENGTH = ICOP1(nCOP+1)
if (LENGTH == 0) GO TO 301
if (LENGTH < 0) GO TO 350
! A long loop over the COP buffer:
do II=1,LENGTH
  IND = ICOP1(II)
  if (ICHK == 0) then
    if (IND /= 0) GO TO 361
    ICHK = 1
    GO TO 360
  end if

  ! Here, if ICHK is 1.
  ICHK = 0
  INDI = IND
  !PAM97 ICOUP = iand(INDI,65535)
  !PAM97 IVL = iand(ishft(INDI,-16),255)
  !ICOUP = mod(INDI,65536)
  !IVL = mod(INDI/65536,256)
  ICOUP = ibits(INDI,0,16)
  IVL = ibits(INDI,16,8)
  ICHK = 0
  INS = 1
  if (IVSAVE == IV0) then
    INS = ICOUPS
    INB = ICOUPS
  end if

  if (INB > 0) then
    ! Transfer HDIAG via buffer HCOUT, write it to unit 25:
    do J=INS,INB
      IOUT = IOUT+1
      HCOUT(IOUT) = HDIAG(J)
      if (IOUT >= nCOP) then
        ! Write out the filled HCOUT buffer:
        if (IFS /= 1) then
          POTNUC = HCOUT(IREF0)
          IFS = 1
        end if
        do KK=1,nCOP
          HCOUT(KK) = HCOUT(KK)-POTNUC
        end do
        call dDAFILE(Lu_25,1,HCOUT,nCOP,IADD25)
        IOUT = 0
      end if
    end do
  end if

  if (IVL /= IV0) then
    JJ = IRC(IVL)+ICOUP
    NSS = MUL(JSYM(JJ),LSYM)
    if (IVL == 1) INB = NVIR(NSS)
    if (IVL > 1) INB = NNS(NSS)
    if (INB > 0) call dDAFILE(Lu_27,2,HDIAG,INB,IAD27)
  end if
  IVSAVE = IVL
  ICOUPS = ICOUP
  GO TO 360

361 continue
  ! Here, if ICHK == 0 and IND /= 0
  !PAM97 ITYP = iand(IND,1)
  !PAM97 IJJ = iand(ishft(IND,-1),2047)
  !ITYP = mod(IND,2)
  !IJJ = mod(IND/2,2048)
  ITYP = ibits(IND,0,1)
  IJJ = ibits(IND,1,11)
  if (ITYP == 0) TERM = COP(II)*FJI(IJJ)
  if (IVL /= IV0) GO TO 362

  ! IVL == IV0, Valence:
  INB = ICOUP
  HDIAG(INB) = HDIAG(INB)+TERM
  GO TO 360

362 if (IVL /= IV1) GO TO 363
  ! IVL == IV1, Singles:
  INB = 0
  NA1 = NSYS(NSS)+1
  NA2 = NSYS(NSS+1)
  if (NA2 < NA1) GO TO 360
  do NA=NA1,NA2
    INB = INB+1
    if (ITYP /= 0) then
      IIJ = IROW(LN+NA)+IJJ
      TERM = COP(II)*FJI(IIJ)
    end if
    HDIAG(INB) = HDIAG(INB)+TERM
  end do
  GO TO 360

363 INB = 0
  ! Doubles:
  do NA=1,NVIRT
    NSA = MUL(NSS,NSM(LN+NA))
    NB1 = NSYS(NSA)+1
    NB2 = NSYS(NSA+1)
    if (NB2 > NA) NB2 = NA
    if (NB2 >= NB1) then
      IIJ1 = IROW(LN+NA)+IJJ
      do NB=NB1,NB2
        INB = INB+1
        if (ITYP /= 0) then
          IIJ2 = IROW(LN+NB)+IJJ
          TERM = COP(II)*(FJI(IIJ1)+FJI(IIJ2))
        end if
        HDIAG(INB) = HDIAG(INB)+TERM
      end do
    end if
  end do

360 continue
end do
GO TO 300

! Transfer remaining HDIAG elements to 25 via buffer HCOUT:
350 if (INB == 0) GO TO 21

do J=1,INB
  IOUT = IOUT+1
  HCOUT(IOUT) = HDIAG(J)
  if (IOUT >= nCOP) then
    ! Write out the filled HCOUT buffer:
    if (IFS /= 1) then
      POTNUC = HCOUT(IREF0)
      IFS = 1
    end if
    do KK=1,nCOP
      HCOUT(KK) = HCOUT(KK)-POTNUC
    end do
    call dDAFILE(Lu_25,1,HCOUT,nCOP,IADD25)
    IOUT = 0
  end if
end do

21 continue
! One last write of the HCOUT buffer:
if (IFS /= 1) then
  POTNUC = HCOUT(IREF0)
  IFS = 1
end if
do KK=1,IOUT
  HCOUT(KK) = HCOUT(KK)-POTNUC
end do
call dDAFILE(Lu_25,1,HCOUT,nCOP,IADD25)
write(6,50) POTNUC
call XFLUSH(6)
50 format(/,6X,'REFERENCE ENERGY',F18.8)

return

end subroutine IJIJ_CPF
