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

use mrci_global, only: ENP, IRC, IREST, IROW, ITER, LN, LSYM, LUSYMB, NSM, NSYM, NVIR, NVIRP, SQ2, SQ2INV
use guga_util_global, only: COP, IAD10, ICOP1, nCOP
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: INTSYM(*), INDX(*)
real(kind=wp), intent(inout) :: C(*), DMO(*), FK(*)
real(kind=wp), intent(_OUT_) :: A(*), B(*)
integer(kind=iwp) :: IADD10, ICHK, ICP1, ICP2, IFT, II, IJOLD, ILEN, IND, INDA, INDB, INK, INMY, INNY, IPOB(9), ITYP, MYL, MYSYM, &
                     NA, NA1, NA2, NAK, NK, NSK, NVM, NYL, NYSYM
real(kind=wp) :: COPI
integer(kind=iwp), external :: JSUNP

! SCRATCH AREAS: A(),B() AND FK().
call CSCALE(INDX,INTSYM,C,SQ2)
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
  if (ILEN == 0) cycle
  if (ILEN < 0) exit
  ! LOOP THROUGH THE COP BUFFER:
  do II=1,ILEN
    IND = ICOP1(II)
    if (ICHK == 0) then
      if (IND == 0) then
        ! IND=0 INDICATES END OF THIS BLOCK OF COUPLING COEFFS.
        ICHK = 1
      else
        if (INK == 0) cycle
        ITYP = ibits(IND,0,6)
        ICP2 = ibits(IND,6,13)
        ICP1 = ibits(IND,19,13)
        if (ITYP > 1) then
          if ((ITER == 1) .and. (IREST == 0)) cycle
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
          B(1:INK) = Zero
          COPI = COP(II)/ENP
          if (NYL /= 1) then
            if (NSK > MYL) then
              call FMMM(C(INMY),C(INNY+IPOB(NSK)),B,1,INK,NVM)
            else
              call FMMM(C(INNY+IPOB(MYL)),C(INMY),B,INK,1,NVM)
              if (IFT == 1) COPI = -COPI
            end if
          else
            if (IFT == 0) call SQUAR(C(INNY+IPOB(MYL)),A,NVM)
            if (IFT == 1) call SQUARN(C(INNY+IPOB(MYL)),A,NVM)
            call FMMM(C(INMY),A,B,1,INK,NVM)
          end if
          FK(1:INK) = FK(1:INK)+COPI*B(1:INK)
        else
          INDA = ICP1
          INDB = IRC(1)+ICP2
          INNY = INDX(INDB)+1
          COPI = C(INDA)*COP(II)/ENP
          FK(1:INK) = FK(1:INK)+COPI*C(INNY:INNY+INK-1)
        end if
      end if
    else
      ! ICHK=1 INDICATES BEGINNING OF A NEW BLOCK OF COUPLING COEFFS.
      ICHK = 0
      if (IJOLD /= 0) then
        ! PUT AWAY FK INTO DMO
        NA1 = NVIRP(NSK)+1
        NA2 = NVIRP(NSK)+NVIR(NSK)
        INK = 0
        if (NA2 < NA1) cycle
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
      if (NA2 < NA1) cycle
      do NA=NA1,NA2
        INK = INK+1
        NAK = IROW(LN+NA)+NK
        FK(INK) = DMO(NAK)
      end do
    end if
  end do
end do
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
