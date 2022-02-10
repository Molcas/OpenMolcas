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

subroutine IJIJ(INTSYM,HDIAG,FC,FIJIJ)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
dimension INTSYM(*), HDIAG(*), FC(*), FIJIJ(*)
dimension HCOUT(nCOP)
!Statement function
JSYM(L) = JSUNP(INTSYM,L)

!------
! POW: Unnecessary but warning stopping initializations
inb = -1234567
!------
IADD25 = IAD25S
IAD27 = 0
IREF0 = 1
call dDAFILE(Lu_27,2,HDIAG,IRC(1),IAD27)

!write(6,*) ' Hdiag'
!write(6,*) ( Hdiag(i),i=1,IRC(1) )

IFS = 0
IVL = 0
IVSAVE = 0
ICOUPS = 0
ICOUP = 0
NSS = 1
IOUT = 0
ICHK = 0
IADD10 = IAD10(3)
TERM = 0.0d00
300 continue
! READ A COP BUFFER:
call dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
call iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
LENGTH = ICOP1(nCOP+1)
if (LENGTH == 0) GO TO 300
if (LENGTH < 0) GO TO 350
! LOOP OVER THE COP BUFFER:
do II=1,LENGTH
  IND = ICOP1(II)
  if (ICHK /= 0) GO TO 460
  if (IND /= 0) GO TO 361
  ICHK = 1
  GO TO 360
460 ICHK = 0
  INDI = IND
  !ICOUP = mod(INDI,2**16)
  !IVL = mod(INDI/2**16,2**8)
  ICOUP = ibits(INDI,0,16)
  IVL = ibits(INDI,16,8)
  ICHK = 0
  INS = 1
  if (IVSAVE == IVVER) then
    INS = ICOUPS
    INB = ICOUPS
  end if
  if (INB /= 0) then
    do J=INS,INB
      IOUT = IOUT+1
      HCOUT(IOUT) = HDIAG(J)
      if (IOUT < nCOP) GO TO 10
      if (IFS == 0) then
        POTNUC = HCOUT(IREF0)
        IFS = 1
      end if
      do KK=1,nCOP
        HCOUT(KK) = HCOUT(KK)-POTNUC
      end do
      call dDAFILE(Lu_25,1,HCOUT,nCOP,IADD25)
      IOUT = 0
10    continue
    end do
  end if
  if (IVL /= IVVER) then
    JJ = IRC(IVL)+ICOUP
    NSS = MUL(JSYM(JJ),LSYM)
    if (IVL == IDVER) then
      INB = NVIR(NSS)
    else
      INB = NVPAIR(NSS)
    end if
    if (INB > 0) call dDAFILE(Lu_27,2,HDIAG,INB,IAD27)
  end if
  IVSAVE = IVL
  ICOUPS = ICOUP
  GO TO 360
361 continue
  !ITYP = mod(IND,2)
  !IJJ = mod(IND/2,2**11)
  ITYP = ibits(IND,0,1)
  IJJ = ibits(IND,1,11)
  if (ITYP == 0) TERM = COP(II)*FIJIJ(IJJ)
  if (IVL == IVVER) then
    INB = ICOUP
    HDIAG(INB) = HDIAG(INB)+TERM
  else if (IVL == IDVER) then
    INB = 0
    NA1 = NVIRP(NSS)+1
    NA2 = NVIRP(NSS)+NVIR(NSS)
    if (NA2 < NA1) GO TO 360
    do NA=NA1,NA2
      INB = INB+1
      if (ITYP == 1) then
        IIJ = IROW(LN+NA)+IJJ
        TERM = COP(II)*FIJIJ(IIJ)
      end if
      HDIAG(INB) = HDIAG(INB)+TERM
    end do
  else
    INB = 0
    do NA=1,NVIRT
      NSA = MUL(NSS,NSM(LN+NA))
      NB1 = NVIRP(NSA)+1
      NB2 = NVIRP(NSA)+NVIR(NSA)
      if (NB2 > NA) NB2 = NA
      if (NB2 < NB1) GO TO 375
      IIJ1 = IROW(LN+NA)+IJJ
      do NB=NB1,NB2
        INB = INB+1
        if (ITYP == 1) then
          IIJ2 = IROW(LN+NB)+IJJ
          TERM = COP(II)*(FIJIJ(IIJ1)+FIJIJ(IIJ2))
        end if
        HDIAG(INB) = HDIAG(INB)+TERM
      end do
375   continue
    end do
  end if
360 continue
end do
GO TO 300
! EMPTY LAST BUFFER
350 continue
do J=1,INB
  IOUT = IOUT+1
  HCOUT(IOUT) = HDIAG(J)
  if (IOUT < nCOP) GO TO 20
  if (IFS == 0) then
    POTNUC = HCOUT(IREF0)
    IFS = 1
  end if
  do KK=1,nCOP
    HCOUT(KK) = HCOUT(KK)-POTNUC
  end do
  call dDAFILE(Lu_25,1,HCOUT,nCOP,IADD25)
  IOUT = 0
20 continue
end do
if (IFS == 0) then
  POTNUC = HCOUT(IREF0)
  IFS = 1
end if
do KK=1,IOUT
  HCOUT(KK) = HCOUT(KK)-POTNUC
end do
call dDAFILE(Lu_25,1,HCOUT,nCOP,IADD25)

return
! Avoid unused argument warnings
if (.false.) call Unused_real_array(FC)

end subroutine IJIJ
