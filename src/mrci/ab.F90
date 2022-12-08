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

subroutine AB(INTSYM,INDX,C,S,FC,A,B,FK)

use mrci_global, only: IFIRST, IRC, IROW, LN, LSYM, NSYM, NVIR, NVIRP, SQ2, SQ2INV
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: INTSYM(*), INDX(*)
real(kind=wp), intent(inout) :: C(*), S(*)
real(kind=wp), intent(_OUT_) :: FC(*), A(*), B(*), FK(*)
integer(kind=iwp) :: IAB, IASYM, ICSYM, IFT, INDA, INMY, IPOA(9), IPOF(9), ITAIL, J, MYL, MYSYM, NA, NA1, NA2, NAA, NAB, NAC, NB, &
                     NCLIM
integer(kind=iwp), external :: JSUNP

call CSCALE(INDX,INTSYM,C,SQ2)
call CSCALE(INDX,INTSYM,S,SQ2INV)
NCLIM = 4
if (IFIRST /= 0) NCLIM = 2
! MOVE FOCK MATRIX TO FK IN SYMMETRY BLOCKS
call IPO(IPOF,NVIR,MUL,NSYM,1,-1)
do IASYM=1,NSYM
  IAB = IPOF(IASYM)
  NA1 = NVIRP(IASYM)+1
  NA2 = NVIRP(IASYM)+NVIR(IASYM)
  do NA=NA1,NA2
    do NB=NA1,NA2
      IAB = IAB+1
      NAB = IROW(LN+NA)+LN+NB
      if (NB > NA) NAB = IROW(LN+NB)+LN+NA
      FK(IAB) = FC(NAB)
      if (NA == NB) FK(IAB) = Zero
    end do
  end do
end do
ITAIL = IRC(NCLIM)
do INDA=1,ITAIL
  if (INDA <= IRC(1)) cycle
  MYSYM = JSUNP(INTSYM,INDA)
  MYL = MUL(MYSYM,LSYM)
  INMY = INDX(INDA)+1
  if (INDA <= IRC(2)) then
    ! DOUBLET-DOUBLET INTERACTIONS
    if (NVIR(MYL) /= 0) then
      A(1:NVIR(MYL)) = Zero
      call FMMM(FK(IPOF(MYL)+1),C(INMY),A,NVIR(MYL),1,NVIR(MYL))
      S(INMY:INMY+NVIR(MYL)-1) = S(INMY:INMY+NVIR(MYL)-1)+A(1:NVIR(MYL))
    end if
  else
    ! TRIPLET-TRIPLET AND SINGLET-SINGLET INTERACTIONS
    IFT = 1
    if (INDA > IRC(3)) IFT = 0
    call IPO(IPOA,NVIR,MUL,NSYM,MYL,IFT)
    do IASYM=1,NSYM
      IAB = IPOF(IASYM+1)-IPOF(IASYM)
      if (IAB == 0) cycle
      ICSYM = MUL(MYL,IASYM)
      if (NVIR(ICSYM) == 0) cycle
      if (MYL == 1) then
        if (IFT == 0) call SQUAR(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
        if (IFT == 1) call SQUARM(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
        NAA = NVIR(IASYM)*NVIR(IASYM)
        B(1:NAA) = Zero
        call FMMM(FK(IPOF(IASYM)+1),A,B,NVIR(IASYM),NVIR(IASYM),NVIR(IASYM))
        A(1:NAA) = B(1:NAA)
        if (IFT /= 1) then
          call SIADD(A,S(INMY+IPOA(IASYM)),NVIR(IASYM))
        else
          call TRADD(A,S(INMY+IPOA(IASYM)),NVIR(IASYM))
        end if
        A(1:NAA) = Zero
      else
        NAC = NVIR(IASYM)*NVIR(ICSYM)
        A(1:NAC) = Zero
        if (IASYM <= ICSYM) then
          call FMMM(FK(IPOF(IASYM)+1),C(INMY+IPOA(ICSYM)),A,NVIR(IASYM),NVIR(ICSYM),NVIR(IASYM))
          J = INMY+IPOA(ICSYM)
          S(J:J+NAC-1) = S(J:J+NAC-1)+A(1:NAC)
        else
          call FMMM(C(INMY+IPOA(IASYM)),FK(IPOF(IASYM)+1),A,NVIR(ICSYM),NVIR(IASYM),NVIR(IASYM))
          J = INMY+IPOA(IASYM)
          S(J:J+NAC-1) = S(J:J+NAC-1)+A(1:NAC)
        end if
      end if
    end do
  end if
end do
call CSCALE(INDX,INTSYM,C,SQ2INV)
call CSCALE(INDX,INTSYM,S,SQ2)

return

end subroutine AB
