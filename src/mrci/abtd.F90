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

subroutine ABTD(ICSPCK,INTSYM,INDX,C1,C2,TDMO,A1,A2,F)

use mrci_global, only: IFIRST, IRC, LN, LSYM, NBAST, NSYM, NVIR, NVIRP, SQ2, SQ2INV
use Symmetry_Info, only: Mul
use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ICSPCK(*), INTSYM(*), INDX(*)
real(kind=wp), intent(inout) :: C1(*), C2(*), TDMO(NBAST,NBAST)
real(kind=wp), intent(_OUT_) :: A1(*), A2(*), F(*)
integer(kind=iwp) :: I, IA, IAB, IASYM, IC, ICSYM, IFT, II1, INDA, INMY, INN, IOC(55), IPF, IPOA(9), IPOF(9), ITAIL, LNA, LNC, &
                     MYL, MYSYM, NA, NA1, NA2, NAC, NB, NCLIM, NVIRA, NVIRC
real(kind=wp) :: TERM, TSUM
integer(kind=iwp), external :: ICUNP, JSUNP
real(kind=wp), external :: DDOT_

! CALCULATE A) TRANSITION DENSITY ELEMENTS OF TYPE TDMO(A,B)
!           B) DIAGONAL ELEMENTS TDMO(I,I) AND TDMO(A,A)
! SCRATCH SPACES: A1(),A2(), SIZE NEEDED IS NVMAX**2
!                  ALSO F(), SIZE NEEDED IS NVSQ
call CSCALE(INDX,INTSYM,C1,SQ2)
call CSCALE(INDX,INTSYM,C2,SQ2)
NCLIM = 4
if (IFIRST /= 0) NCLIM = 2
! MOVE TRANSITION DENSITY MATRIX TO F IN SYMMETRY BLOCKS
call IPO(IPOF,NVIR,MUL,NSYM,1,-1)
do IASYM=1,NSYM
  IAB = IPOF(IASYM)
  NA1 = NVIRP(IASYM)+1
  NA2 = NVIRP(IASYM)+NVIR(IASYM)
  do NA=NA1,NA2
    do NB=NA1,NA2
      IAB = IAB+1
      F(IAB) = TDMO(LN+NA,LN+NB)
    end do
  end do
end do
II1 = 0
ITAIL = IRC(NCLIM)
do INDA=1,ITAIL
  do I=1,LN
    II1 = II1+1
    IOC(I) = (1+ICUNP(ICSPCK,II1))/2
  end do
  if (INDA <= IRC(1)) then
    TSUM = C1(INDA)*C2(INDA)
  else
    MYSYM = JSUNP(INTSYM,INDA)
    MYL = MUL(MYSYM,LSYM)
    INMY = INDX(INDA)+1
    if (INDA <= IRC(2)) then
      ! DOUBLET-DOUBLET INTERACTIONS
      if (NVIR(MYL) == 0) cycle
      IPF = IPOF(MYL)+1
      NVIRA = NVIR(MYL)
      call DGER(NVIRA,NVIRA,One,C1(INMY),1,C2(INMY),1,F(IPF),NVIRA)
      LNA = LN+NVIRP(MYL)
      TSUM = Zero
      do I=1,NVIRA
        TERM = C1(INMY-1+I)*C2(INMY-1+I)
        IA = LNA+I
        TDMO(IA,IA) = TDMO(IA,IA)+TERM
        TSUM = TSUM+TERM
      end do
    else
      ! TRIPLET-TRIPLET AND SINGLET-SINGLET INTERACTIONS
      IFT = 1
      if (INDA > IRC(3)) IFT = 0
      call IPO(IPOA,NVIR,MUL,NSYM,MYL,IFT)
      TSUM = Zero
      do IASYM=1,NSYM
        IAB = IPOF(IASYM+1)-IPOF(IASYM)
        if (IAB == 0) cycle
        ICSYM = MUL(MYL,IASYM)
        NVIRA = NVIR(IASYM)
        NVIRC = NVIR(ICSYM)
        if (NVIRC == 0) cycle
        if (MYL /= 1) then
          if (IASYM > ICSYM) then
            call MTRANS(C1(INMY+IPOA(IASYM)),A1,NVIRA,NVIRC)
            call MTRANS(C2(INMY+IPOA(IASYM)),A2,NVIRA,NVIRC)
          else
            NAC = NVIRA*NVIRC
            if (IFT == 0) then
              call DCOPY_(NAC,C1(INMY+IPOA(ICSYM)),1,A1,1)
              call DCOPY_(NAC,C2(INMY+IPOA(ICSYM)),1,A2,1)
            else
              call VNEG(NAC,C1(INMY+IPOA(ICSYM)),1,A1,1)
              call VNEG(NAC,C2(INMY+IPOA(ICSYM)),1,A2,1)
            end if
          end if
        else
          if (IFT == 0) then
            call SQUAR(C1(INMY+IPOA(IASYM)),A1,NVIRA)
            call SQUAR(C2(INMY+IPOA(IASYM)),A2,NVIRA)
          else
            call SQUARM(C1(INMY+IPOA(IASYM)),A1,NVIRA)
            call SQUARM(C2(INMY+IPOA(IASYM)),A2,NVIRA)
          end if
        end if
        IPF = IPOF(IASYM)+1
        call DGEMM_('N','T',NVIRA,NVIRA,NVIRC,One,A1,NVIRA,A2,NVIRA,One,F(IPF),NVIRA)
        INN = 1
        LNC = LN+NVIRP(ICSYM)
        do I=1,NVIRC
          TERM = DDOT_(NVIRA,A1(INN),1,A2(INN),1)
          TSUM = TSUM+TERM
          IC = LNC+I
          TDMO(IC,IC) = TDMO(IC,IC)+TERM
          INN = INN+NVIRA
        end do
      end do
      TSUM = TSUM/2
    end if
  end if
  do I=1,LN
    TDMO(I,I) = TDMO(I,I)+IOC(I)*TSUM
  end do
end do
do IASYM=1,NSYM
  IAB = IPOF(IASYM)
  NA1 = NVIRP(IASYM)+1
  NA2 = NVIRP(IASYM)+NVIR(IASYM)
  do NA=NA1,NA2
    do NB=NA1,NA2
      IAB = IAB+1
      if (NA /= NB) TDMO(LN+NA,LN+NB) = F(IAB)
    end do
  end do
end do
call CSCALE(INDX,INTSYM,C1,SQ2INV)
call CSCALE(INDX,INTSYM,C2,SQ2INV)

return

end subroutine ABTD
