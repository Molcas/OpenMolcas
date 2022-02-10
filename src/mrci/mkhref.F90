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

!pgi$g opt=1
subroutine MKHREF(HREF,FC,FIJKL,JREFX)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
dimension HREF(*), FC(*), FIJKL(*), JREFX(NCVAL)

NHREF = (NREF*(NREF+1))/2
call FZERO(HREF,NHREF)
ICHK = 0
IK = 0
FINI = 0.0d00
IADD25 = 0
call dDAFILE(Lu_25,2,FC,NBTRI,IADD25)
IADD10 = IAD10(8)
100 call dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
call iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
LEN = ICOP1(nCOP+1)
if (LEN == 0) GO TO 100
if (LEN < 0) goto 200
do IN=1,LEN
  IND = ICOP1(IN)
  if (ICHK /= 0) then
    ICHK = 0
    INDI = IND
    !NI = mod(INDI,2**10)
    !NK = mod(INDI/2**10,2**10)
    NI = ibits(INDI,0,10)
    NK = ibits(INDI,10,10)
    IK = IROW(NK)+NI
    FINI = FC(IK)
    GO TO 10
  end if
  if (IND == 0) then
    ICHK = 1
    GO TO 10
  end if
  !IVL = mod(IND,2**6)
  IVL = ibits(IND,0,6)
  if (IVL /= IVVER) GO TO 10
  !IC2 = mod(IND/2**6,2**13)
  IC2 = ibits(IND,6,13)
  NA = JREFX(IC2)
  if (NA == 0) GO TO 10
  !IC1 = mod(IND/2**19,2**13)
  IC1 = ibits(IND,19,13)
  NB = JREFX(IC1)
  if (NB == 0) GO TO 10
  if (NA < NB) then
    NAT = NA
    NA = NB
    NB = NAT
  end if
  IVEC = (NA*(NA-1))/2+NB
  HREF(IVEC) = HREF(IVEC)+COP(IN)*FINI
10 continue
end do
GO TO 100
200 continue
ICHK = 0
NIJ = IROW(LN+1)
NIJKL = NIJ*(NIJ+1)/2
call FZERO(FIJKL,NIJKL)
FINI = 0.0d00
IADR = LASTAD(1)
201 continue
call dDAFILE(Lu_70,2,VALSRT,NSRTMX,IADR)
call iDAFILE(Lu_70,2,INDSRT,NSRTMX+2,IADR)
LENGTH = INDSRT(NSRTMX+1)
IADR = INDSRT(NSRTMX+2)
!if (LENGTH > 0) call SCATTER(LENGTH,FIJKL,INDSRT,VALSRT)
do i=1,length
  FIJKL(INDSRT(i)) = VALSRT(i)
end do
if (IADR /= -1) GO TO 201
IADD10 = IAD10(5)
300 continue
call dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
call iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
LEN = ICOP1(nCOP+1)
if (LEN == 0) GO TO 300
if (LEN < 0) goto 400
do IN=1,LEN
  IND = ICOP1(IN)
  if (ICHK /= 0) then
    ICHK = 0
    INDI = IND
    !PAM96 IP = iand(INDI,255)
    !PAM96 JP = iand(ishft(INDI,-8),255)
    !PAM96 KP = iand(ishft(INDI,-16),255)
    !PAM96 LP = iand(ishft(INDI,-24),255)
    !IP = mod(INDI,2**8)
    !JP = mod(INDI/2**8,2**8)
    !KP = mod(INDI/2**16,2**8)
    !LP = mod(INDI/2**24,2**8)
    IP = ibits(INDI,0,8)
    JP = ibits(INDI,8,8)
    KP = ibits(INDI,16,8)
    LP = ibits(INDI,24,8)
    NIJ = IROW(IP)+JP
    NKL = IROW(KP)+LP
    IND = NIJ*(NIJ-1)/2+NKL
    FINI = FIJKL(IND)
    goto 310
  end if
  if (IND == 0) then
    ICHK = 1
    goto 310
  end if
  !IVL = mod(IND,2**6)
  IVL = ibits(IND,0,6)
  if (IVL /= 0) GO TO 310
  !IC2 = mod(IND/2**6,2**13)
  IC2 = ibits(IND,6,13)
  NA = JREFX(IC2)
  if (NA == 0) GO TO 310
  !IC1 = mod(IND/2**19,2**13)
  IC1 = ibits(IND,19,13)
  NB = JREFX(IC1)
  if (NB == 0) GO TO 310
  if (NA < NB) then
    NAT = NA
    NA = NB
    NB = NAT
  end if
  IVEC = (NA*(NA-1))/2+NB
  HREF(IVEC) = HREF(IVEC)+COP(IN)*FINI
310 continue
end do
GO TO 300
400 continue
IADD25 = IAD25S
IBUF = nCOP
do I=1,IRC(1)
  IBUF = IBUF+1
  if (IBUF > nCOP) then
    call dDAFILE(Lu_25,2,COP,nCOP,IADD25)
    IBUF = 1
  end if
  IR = JREFX(I)
  if (IR > 0) then
    IIR = (IR*(IR+1))/2
    HREF(IIR) = HREF(IIR)+COP(IBUF)+POTNUC
  end if
end do

return

end subroutine MKHREF
