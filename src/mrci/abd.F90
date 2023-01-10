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

subroutine ABD(ICSPCK,INTSYM,INDX,C,DMO,A,B,F)

use mrci_global, only: ENP, IFIRST, IRC, IROW, LN, LSYM, NCSHT, NELEC, NSYM, NVIR, NVIRP, SQ2, SQ2INV
use Symmetry_Info, only: Mul
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ICSPCK(*), INTSYM(*), INDX(*)
real(kind=wp), intent(inout) :: C(*)
real(kind=wp), intent(_OUT_) :: DMO(*), A(*), B(*), F(*)
integer(kind=iwp) :: I, IAB, IASYM, ICSYM, IFT, II1, IIA, IIC, IIN, IJ, INDA, INMY, INN, IOC(55), IPF, IPOA(9), IPOF(9), ITAIL, &
                     LNA, LNC, MYL, MYSYM, NA, NA1, NA2, NAB, NAC, NB, NCLIM, NVIRA, NVIRC
real(kind=wp) :: ENPINV, RSUM, TR, TSUM
integer(kind=iwp), external :: ICUNP, JSUNP
real(kind=wp), external :: DDOT_

! SCRATCH SPACE: A(),B(),F().
call CSCALE(INDX,INTSYM,C,SQ2)
NCLIM = 4
if (IFIRST /= 0) NCLIM = 2
ENPINV = One/ENP
! MOVE DENSITY MATRIX TO F IN SYMMETRY BLOCKS
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
      F(IAB) = DMO(NAB)
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
    TSUM = ENPINV*C(INDA)**2
  else
    MYSYM = JSUNP(INTSYM,INDA)
    MYL = MUL(MYSYM,LSYM)
    INMY = INDX(INDA)+1
    if (INDA <= IRC(2)) then
      ! DOUBLET-DOUBLET INTERACTIONS
      if (NVIR(MYL) == 0) cycle
      call FMUL2(C(INMY),C(INMY),A,NVIR(MYL),NVIR(MYL),1)
      IPF = IPOF(MYL)+1
      IIN = IPOF(MYL+1)-IPOF(MYL)
      F(IPF:IPF+IIN-1) = F(IPF:IPF+IIN-1)+ENPINV*A(1:IIN)
      NVIRA = NVIR(MYL)
      LNA = LN+NVIRP(MYL)
      IIA = IROW(LNA+1)
      TSUM = Zero
      do I=1,NVIRA
        RSUM = ENPINV*C(INMY)**2
        INMY = INMY+1
        TSUM = TSUM+RSUM
        IIA = IIA+LNA+I
        DMO(IIA) = DMO(IIA)+RSUM
      end do
    else
      ! TRIPLET-TRIPLET AND SINGLET-SINGLET INTERACTIONS
      IFT = 1
      if (INDA > IRC(3)) IFT = 0
      call IPO(IPOA,NVIR,MUL,NSYM,MYL,IFT)
      IIN = 0
      TSUM = Zero
      do IASYM=1,NSYM
        IAB = IPOF(IASYM+1)-IPOF(IASYM)
        if (IAB == 0) cycle
        ICSYM = MUL(MYL,IASYM)
        if (NVIR(ICSYM) == 0) cycle
        if (MYL /= 1) then
          if (IASYM > ICSYM) then
            call MTRANS(C(INMY+IPOA(IASYM)),A,NVIR(IASYM),NVIR(ICSYM))
          else
            NAC = NVIR(IASYM)*NVIR(ICSYM)
            if (IFT == 0) call DCOPY_(NAC,C(INMY+IPOA(ICSYM)),1,A,1)
            if (IFT == 1) call VNEG(NAC,C(INMY+IPOA(ICSYM)),1,A,1)
          end if
        else
          if (IFT == 0) call SQUAR(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
          if (IFT == 1) call SQUARM(C(INMY+IPOA(IASYM)),A,NVIR(IASYM))
        end if
        NVIRA = NVIR(IASYM)
        NVIRC = NVIR(ICSYM)
        call FMUL2(A,A,B,NVIR(IASYM),NVIR(IASYM),NVIR(ICSYM))
        IPF = IPOF(IASYM)+1
        F(IPF:IPF+IAB-1) = F(IPF:IPF+IAB-1)+ENPINV*B(1:IAB)
        INN = 1
        LNC = LN+NVIRP(ICSYM)
        IIC = IROW(LNC+1)
        do I=1,NVIRC
          RSUM = ENPINV*DDOT_(NVIRA,A(INN),1,A(INN),1)
          TSUM = TSUM+RSUM
          IIC = IIC+LNC+I
          DMO(IIC) = DMO(IIC)+RSUM
          INN = INN+NVIRA
        end do
      end do
      TSUM = TSUM/2
    end if
  end if
  IJ = 0
  do I=1,LN
    IJ = IJ+I
    DMO(IJ) = DMO(IJ)+IOC(I)*TSUM
  end do
end do
do IASYM=1,NSYM
  IAB = IPOF(IASYM)
  NA1 = NVIRP(IASYM)+1
  NA2 = NVIRP(IASYM)+NVIR(IASYM)
  do NA=NA1,NA2
    do NB=NA1,NA2
      IAB = IAB+1
      if (NA >= NB) cycle
      NAB = IROW(LN+NB)+LN+NA
      DMO(NAB) = F(IAB)
    end do
  end do
end do
TR = Zero
IJ = 0
do I=1,NCSHT
  IJ = IJ+I
  TR = TR+DMO(IJ)
end do
if (abs(TR-real(NELEC,kind=wp)) > 1.0e-8_wp) write(u6,310) TR
call CSCALE(INDX,INTSYM,C,SQ2INV)

return

310 format(/,6X,'TRACE OF DENSITY MATRIX',F16.8)

end subroutine ABD
