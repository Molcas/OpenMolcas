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

subroutine AID(INTSYM,INDX,C,DMO,A,B,FK)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
dimension INTSYM(*), INDX(*), C(*), DMO(*), A(*), B(*), FK(*)
dimension IPOB(9)
!Statement function
JSYM(L) = JSUNP(INTSYM,L)

! SCRATCH AREAS: A(),B() AND FK().
call CSCALE(INDX,INTSYM,C,SQ2)
ICHK = 0
IJOLD = 0
NK = 0
NSK = 1
IADD10 = IAD10(9)
! READ A COP BUFFER
100 continue
call dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
call iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
LEN = ICOP1(nCOP+1)
if (LEN == 0) GO TO 100
if (LEN < 0) GO TO 200
! LOOP THROUGH THE COP BUFFER:
do II=1,LEN
  IND = ICOP1(II)
  if (ICHK /= 0) GO TO 460
  if (IND /= 0) GO TO 11
  ! IND=0 INDICATES END OF THIS BLOCK OF COUPLING COEFFS.
  ICHK = 1
  GO TO 10
460 continue
  ! ICHK=1 INDICATES BEGINNING OF A NEW BLOCK OF COUPLING COEFFS.
  ICHK = 0
  if (IJOLD /= 0) then
    ! PUT AWAY FK INTO DMO
    NA1 = NVIRP(NSK)+1
    NA2 = NVIRP(NSK)+NVIR(NSK)
    INK = 0
    if (NA2 < NA1) GO TO 10
    do NA=NA1,NA2
      INK = INK+1
      NAK = IROW(LN+NA)+NK
      DMO(NAK) = FK(INK)
    end do
  end if
  NK = IND
  IJOLD = NK
  NSK = NSM(NK)
  ! PICK OUT ELEMENTS FROM DMO AND PUT INTO FK:
  NA1 = NVIRP(NSK)+1
  NA2 = NVIRP(NSK)+NVIR(NSK)
  INK = 0
  if (NA2 < NA1) GO TO 10
  do NA=NA1,NA2
    INK = INK+1
    NAK = IROW(LN+NA)+NK
    FK(INK) = DMO(NAK)
  end do
  GO TO 10
11 if (INK == 0) GO TO 10
  !ITYP = mod(IND,2**6)
  !ICP2 = mod(IND/2**6,2**13)
  !ICP1 = mod(IND/2**19,2**13)
  ITYP = ibits(IND,0,6)
  ICP2 = ibits(IND,6,13)
  ICP1 = ibits(IND,19,13)
  if (ITYP > 1) GO TO 12
  INDA = ICP1
  INDB = IRC(1)+ICP2
  INNY = INDX(INDB)+1
  COPI = C(INDA)*COP(II)/ENP
  call DAXPY_(INK,COPI,C(INNY),1,FK,1)
  GO TO 10
12 if ((ITER == 1) .and. (IREST == 0)) GO TO 10
  INDA = IRC(1)+ICP1
  INDB = IRC(ITYP)+ICP2
  INMY = INDX(INDA)+1
  INNY = INDX(INDB)+1
  MYSYM = JSYM(INDA)
  NYSYM = MUL(MYSYM,NSK)
  MYL = MUL(MYSYM,LSYM)
  NYL = MUL(NYSYM,LSYM)
  IFT = 0
  if (ITYP == 2) IFT = 1
  call IPO(IPOB,NVIR,MUL,NSYM,NYL,IFT)
  NVM = NVIR(MYL)
  call FZERO(B,INK)
  COPI = COP(II)/ENP
  if (NYL /= 1) then
    if (NSK > MYL) then
      call FMMM(C(INMY),C(INNY+IPOB(NSK)),B,1,INK,NVM)
      call VSMA(B,1,COPI,FK,1,FK,1,INK)
    else
      call FMMM(C(INNY+IPOB(MYL)),C(INMY),B,INK,1,NVM)
      if (IFT == 1) COPI = -COPI
      call VSMA(B,1,COPI,FK,1,FK,1,INK)
    end if
  else
    if (IFT == 0) call SQUAR(C(INNY+IPOB(MYL)),A,NVM)
    if (IFT == 1) call SQUARN(C(INNY+IPOB(MYL)),A,NVM)
    call FMMM(C(INMY),A,B,1,INK,NVM)
    call VSMA(B,1,COPI,FK,1,FK,1,INK)
  end if
10 continue
end do
GO TO 100
200 continue
NA1 = NVIRP(NSK)+1
NA2 = NVIRP(NSK)+NVIR(NSK)
INK = 0
do NA=NA1,NA2
  INK = INK+1
  NAK = IROW(LN+NA)+NK
  DMO(NAK) = FK(INK)
end do
call CSCALE(INDX,INTSYM,C,SQ2INV)

return

end subroutine AID
