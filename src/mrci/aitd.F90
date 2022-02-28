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

use mrci_global, only: IRC, LN, LSYM, LUSYMB, NBAST, NSM, NSYM, NVIR, NVIRP, SQ2, SQ2INV
use guga_util_global, only: COP, IAD10, ICOP1, nCOP
use Symmetry_Info, only: Mul
use Constants, only: One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: INTSYM(*), INDX(*)
real(kind=wp), intent(inout) :: C1(*), C2(*), TDMO(NBAST,NBAST), FAK(*), FKA(*)
real(kind=wp), intent(_OUT_) :: A(*)
integer(kind=iwp) :: IADD10, ICHK, ICP1, ICP2, IFT, II, IJOLD, ILEN, IND, INDA, INDB, INK, INMY, INNY, IPOB(9), ITYP, MYL, MYSYM, &
                     NA, NA1, NA2, NK, NSK, NVM, NYL, NYSYM
real(kind=wp) :: COPI
integer(kind=iwp), external :: JSUNP

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
do
  call dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
  call iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
  ILEN = ICOP1(nCOP+1)
  if (ILEN < 0) exit
  ! LOOP THROUGH THE COP BUFFER:
  do II=1,ILEN
    IND = ICOP1(II)
    if (ICHK == 0) then
      if (IND == 0) then
        ICHK = 1
      else if (INK /= 0) then
        ITYP = ibits(IND,0,6)
        ICP2 = ibits(IND,6,13)
        ICP1 = ibits(IND,19,13)
        if (ITYP <= 1) then
          INDA = ICP1
          INDB = IRC(1)+ICP2
          INNY = INDX(INDB)+1
          COPI = C1(INDA)*COP(II)
          FAK(1:INK) = FAK(1:INK)+COPI*C2(INNY:INNY+INK-1)
          COPI = C2(INDA)*COP(II)
          FKA(1:INK) = FKA(1:INK)+COPI*C1(INNY:INNY+INK-1)
        else
          INDA = IRC(1)+ICP1
          INDB = IRC(ITYP)+ICP2
          INMY = INDX(INDA)+1
          INNY = INDX(INDB)+1
          MYSYM = JSUNP(INTSYM,INDA)
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
              call DGEMV_('T',NVM,INK,COPI,C1(INNY+IPOB(NSK)),NVM,C2(INMY),1,One,FAK,1)
              call DGEMV_('T',NVM,INK,COPI,C2(INNY+IPOB(NSK)),NVM,C1(INMY),1,One,FKA,1)
            else
              if (IFT == 1) COPI = -COPI
              call DGEMV_('N',INK,NVM,COPI,C1(INNY+IPOB(MYL)),INK,C2(INMY),1,One,FAK,1)
              call DGEMV_('N',INK,NVM,COPI,C2(INNY+IPOB(MYL)),INK,C1(INMY),1,One,FKA,1)
            end if
          else
            if (IFT == 0) call SQUAR(C1(INNY+IPOB(MYL)),A,NVM)
            if (IFT == 1) call SQUARN(C1(INNY+IPOB(MYL)),A,NVM)
            call DGEMV_('T',NVM,INK,COPI,A,NVM,C2(INMY),1,One,FAK,1)
            if (IFT == 0) call SQUAR(C2(INNY+IPOB(MYL)),A,NVM)
            if (IFT == 1) call SQUARN(C2(INNY+IPOB(MYL)),A,NVM)
            call DGEMV_('T',NVM,INK,COPI,A,NVM,C1(INMY),1,One,FKA,1)
          end if
        end if
      end if
    else
      ICHK = 0
      if (IJOLD /= 0) then
        ! PUT FAK,FKA BACK INTO TDMO.
        NA1 = NVIRP(NSK)+1
        NA2 = NVIRP(NSK)+NVIR(NSK)
        INK = 0
        if (NA2 < NA1) cycle
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
      do NA=NA1,NA2
        INK = INK+1
        FAK(INK) = TDMO(LN+NA,NK)
        FKA(INK) = TDMO(NK,LN+NA)
      end do
    end if
  end do
end do
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
