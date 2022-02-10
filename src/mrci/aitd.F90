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

subroutine AITD(INTSYM,INDX,C1,C2,TDMO,A,FAK,FKA)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
dimension INTSYM(*), INDX(*), C1(*), C2(*), TDMO(NBAST,NBAST), A(*), FAK(*), FKA(*)
dimension IPOB(9)
!Statement function
JSYM(L) = JSUNP(INTSYM,L)

! CALCULATE TRANSITION DENSITY ELEMENTS TDMO(K,A) AND TDMO(A,K),
! WHERE K IS INTERNAL, A IS EXTERNAL ORBITAL.
! SCRATCH AREAS ARE: A(), SIZE NEEDED IS NVMAX**2
!       AND FAK(), FKA(), SIZE NEEDED IS NVMAX
call CSCALE(INDX,INTSYM,C1,SQ2)
call CSCALE(INDX,INTSYM,C2,SQ2)
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
  ICHK = 1
  GO TO 10
460 ICHK = 0
  if (IJOLD /= 0) then
    ! PUT FAK,FKA BACK INTO TDMO.
    NA1 = NVIRP(NSK)+1
    NA2 = NVIRP(NSK)+NVIR(NSK)
    INK = 0
    if (NA2 < NA1) GO TO 10
    do NA=NA1,NA2
      INK = INK+1
      TDMO(LN+NA,NK) = FAK(INK)
      TDMO(NK,LN+NA) = FKA(INK)
    end do
  end if
  NK = IND
  IJOLD = NK
  NSK = NSM(NK)
  ! PUT TDMO ELEMENTS INTO ARRAYS FAK, FKA.
  NA1 = NVIRP(NSK)+1
  NA2 = NVIRP(NSK)+NVIR(NSK)
  INK = 0
  if (NA2 < NA1) GO TO 10
  do NA=NA1,NA2
    INK = INK+1
    FAK(INK) = TDMO(LN+NA,NK)
    FKA(INK) = TDMO(NK,LN+NA)
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
  COPI = C1(INDA)*COP(II)
  call DAXPY_(INK,COPI,C2(INNY),1,FAK,1)
  COPI = C2(INDA)*COP(II)
  call DAXPY_(INK,COPI,C1(INNY),1,FKA,1)
  GO TO 10
12 continue
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
  COPI = COP(II)
  if (NYL /= 1) then
    if (NSK > MYL) then
      !call DGEMTX (NVM,INK,COPI,C1(INNY+IPOB(NSK)),NVM,C2(INMY),1,FAK,1)
      call DGEMV_('T',NVM,INK,COPI,C1(INNY+IPOB(NSK)),NVM,C2(INMY),1,1.0d0,FAK,1)
      !call DGEMTX (NVM,INK,COPI,C2(INNY+IPOB(NSK)),NVM,C1(INMY),1,FKA,1)
      call DGEMV_('T',NVM,INK,COPI,C2(INNY+IPOB(NSK)),NVM,C1(INMY),1,1.0d0,FKA,1)
    else
      if (IFT == 1) COPI = -COPI
      !call DGEMX (INK,NVM,COPI,C1(INNY+IPOB(MYL)),INK,C2(INMY),1,FAK,1)
      call DGEMV_('N',INK,NVM,COPI,C1(INNY+IPOB(MYL)),INK,C2(INMY),1,1.0d0,FAK,1)
      !call DGEMX (INK,NVM,COPI,C2(INNY+IPOB(MYL)),INK,C1(INMY),1,FKA,1)
      call DGEMV_('N',INK,NVM,COPI,C2(INNY+IPOB(MYL)),INK,C1(INMY),1,1.0d0,FKA,1)
    end if
  else
    if (IFT == 0) call SQUAR(C1(INNY+IPOB(MYL)),A,NVM)
    if (IFT == 1) call SQUARN(C1(INNY+IPOB(MYL)),A,NVM)
    !call DGEMTX (NVM,INK,COPI,A,NVM,C2(INMY),1,FAK,1)
    call DGEMV_('T',NVM,INK,COPI,A,NVM,C2(INMY),1,1.0d0,FAK,1)
    if (IFT == 0) call SQUAR(C2(INNY+IPOB(MYL)),A,NVM)
    if (IFT == 1) call SQUARN(C2(INNY+IPOB(MYL)),A,NVM)
    !call DGEMTX (NVM,INK,COPI,A,NVM,C1(INMY),1,FKA,1)
    call DGEMV_('T',NVM,INK,COPI,A,NVM,C1(INMY),1,1.0d0,FKA,1)
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
  TDMO(LN+NA,NK) = FAK(INK)
  TDMO(NK,LN+NA) = FKA(INK)
end do
call CSCALE(INDX,INTSYM,C1,SQ2INV)
call CSCALE(INDX,INTSYM,C2,SQ2INV)

return

end subroutine AITD
