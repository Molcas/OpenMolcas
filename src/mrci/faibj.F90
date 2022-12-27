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

subroutine FAIBJ(INTSYM,INDX,C,S,ABIJ,AIBJ,AJBI,A,B,F,FSEC)

use mrci_global, only: IRC, IREST, IROW, ITER, LASTAD, LN, Lu_60, LUSYMB, NBITM3, NSM, NSYM, NVIR, NVIRT, SQ2, SQ2INV
use guga_util_global, only: COP, IAD10, ICOP1, nCOP
use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: INTSYM(*), INDX(*)
real(kind=wp), intent(inout) :: C(*), S(*), ABIJ(*), AIBJ(*), AJBI(*), FSEC(*)
real(kind=wp), intent(_OUT_) :: A(*), B(*), F(*)
integer(kind=iwp) :: IAB, IADD10, IADR, IASYM, IBSYM, ICHK, iCoup, iCoup1, IFAB, IFT, IFTA, IFTB, II, IIN, IJ1, ILIM, INDA, INDB, &
                     INDCOP, INDI, INMY, INNY, INS, IPF, IPF1, IPOA(9), IPOB(9), IPOF(9), ISTAR, ITURN, iTyp, JTURN, LENBUF, &
                     LENCOP, MYL, MYSYM, NI, NJ, NOT2, NOVST, NSIJ, NVIRA, NVIRB, NVIRC, NYL, NYSYM
real(kind=wp) :: COPI, CPL, CPLA, FAC, FACS, TERM
logical(kind=iwp) :: Skip
integer(kind=iwp), allocatable :: iBuf(:)
real(kind=wp), allocatable :: Buf(:)
real(kind=wp), external :: DDOT_

call mma_allocate(Buf,NBITM3,label='Buf')
call mma_allocate(iBuf,NBITM3+2,label='iBuf')

!vv this code is a real compiler killer!

! POW: Unnecessary but warningstopping initializations
iTyp = -1234567
iCoup = -1234567
iCoup1 = -1234567

call CSCALE(INDX,INTSYM,C,SQ2)
call CSCALE(INDX,INTSYM,S,SQ2INV)
ICHK = 0
IFAB = 0
NOVST = LN*NVIRT+1+(NVIRT*(NVIRT+1))/2
NOT2 = IROW(LN+1)

IADD10 = IAD10(6)

! Long loop, reading buffers until end of buffers is signalled
! by length field holding a negative number.
do
  call dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
  call iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
  LENCOP = ICOP1(nCOP+1)
  if (LENCOP < 0) exit

  ! Loop over the elements of this buffer
  do II=1,LENCOP
    INDCOP = ICOP1(II)
    if (ICHK == 0) then
      if (INDCOP == 0) then
        ICHK = 1
      else
        if (IFAB == 1) then
          CPLA = COP(II)
          IFAB = 0
        else
          IFAB = ibits(INDCOP,0,1)
          ITURN = ibits(INDCOP,1,1)
          ITYP = ibits(INDCOP,2,3)
          ICOUP = ibits(INDCOP,5,13)
          ICOUP1 = ibits(INDCOP,18,13)
          CPL = COP(II)
          CPLA = Zero
          if (IFAB /= 0) cycle
          if (ITURN == 0) then
            ! FIRST ORDER INTERACTION
            INDA = ICOUP
            INDB = IRC(ITYP+1)+ICOUP1
            ISTAR = 1
            if (ITYP == 1) ISTAR = INS+1
            if (INS /= 0) then
              COPI = CPL*C(INDA)
              S(INDX(INDB)+1:INDX(INDB)+INS) = S(INDX(INDB)+1:INDX(INDB)+INS)+COPI*FSEC(ISTAR:ISTAR+INS-1)
              TERM = DDOT_(INS,FSEC(ISTAR),1,C(INDX(INDB)+1),1)
              S(INDA) = S(INDA)+CPL*TERM
            end if
            cycle
          end if
        end if

        ! INTERACTIONS BETWEEN DOUBLES AND
        ! INTERACTIONS BETWEEN SINGLES
        if ((ITER /= 1) .or. (IREST /= 0)) then

          call faibj2(IFTA,IFTB,ICOUP1,ICOUP,INDA,INDB,MYSYM,INTSYM,NYSYM,NSIJ,MYL,NYL,FACS,IPOA,IPOB,INMY,INNY,INDX,iTYP)

          if (ITYP == 5) then
            ! DOUBLET-DOUBLET INTERACTIONS
            IIN = IPOF(MYL+1)-IPOF(MYL)
            if (IIN /= 0) then
              IPF = IPOF(MYL)+1
              F(1:IIN) = CPL*AIBJ(IPF:IPF+IIN-1)+CPLA*ABIJ(IPF:IPF+IIN-1)
              if (INDA == INDB) call DCOPY_(NVIR(MYL),[Zero],0,F,NVIR(MYL)+1)
              call DGEMV_('T',NVIR(MYL),NVIR(NYL),FACS,F,NVIR(MYL),C(INMY),1,One,S(INNY),1)
              if (INDA /= INDB) then
                call DGEMV_('N',NVIR(MYL),NVIR(NYL),FACS,F,NVIR(MYL),C(INNY),1,One,S(INMY),1)
              end if
            end if
          else
            ! TRIPLET-SINGLET, SINGLET-TRIPLET,
            ! TRIPLET-TRIPLET AND SINGLET-SINGLET INTERACTIONS

            call loop70(C,S,ABIJ,AIBJ,AJBI,A,B,F,IPOF,IPOA,IPOB,MYL,NYL,INDA,INDB,INMY,INNY,IFTB,IFTA,FACS,IAB,CPL,CPLA,NVIRA, &
                        NVIRC,NVIRB)

          end if
        end if
      end if
    else
      ICHK = 0

      ! Unpack indices NI and NJ from INDCOP
      INDI = INDCOP
      NI = ibits(INDI,0,10)
      NJ = ibits(INDI,10,10)

      NSIJ = MUL(NSM(NI),NSM(NJ))
      call IPO(IPOF,NVIR,MUL,NSYM,NSIJ,-1)
      IJ1 = IROW(NI)+NJ
      ILIM = IPOF(NSYM+1)
      ! Clear matrices ABIJ, AIBJ, and AJBI.
      ABIJ(1:ILIM) = Zero
      AIBJ(1:ILIM) = Zero
      AJBI(1:ILIM) = Zero
      if ((ITER == 1) .and. (IREST == 0)) then
        Skip = .true.
      else
        ! READ (AB/IJ) INTEGRALS

        IADR = LASTAD(NOVST+IJ1)
        JTURN = 0
        Skip = .false.
      end if
      do
        if (Skip) then
          Skip = .false.
        else
          call iDAFILE(Lu_60,2,iBuf,NBITM3+2,IADR)
          call dDAFILE(Lu_60,2,Buf,NBITM3,IADR)
          LENBUF = iBuf(NBITM3+1)
          IADR = iBuf(NBITM3+2)
          call faibj5(LENBUF,JTURN,iBuf,Buf,AIBJ,ABIJ)

          if (IADR /= -1) cycle
          if (JTURN == 1) exit
        end if

        ! READ (AI/BJ) INTEGRALS

        IADR = LASTAD(NOVST+NOT2+IJ1)
        JTURN = 1
      end do

      ! CONSTRUCT FIRST ORDER MATRICES

      FAC = One
      if (NI == NJ) FAC = Half
      IIN = 0

      IFT = 0
      call faibj3(NSIJ,IFT,AIBJ,FSEC,FAC,IIN,INS,IPOA,IPOF)

      if ((ITER /= 1) .or. (IREST /= 0)) then
        do IASYM=1,NSYM
          NVIRA = NVIR(IASYM)
          if (NVIRA == 0) cycle
          IBSYM = MUL(NSIJ,IASYM)
          NVIRB = NVIR(IBSYM)
          if (NVIRB == 0) cycle
          IPF = IPOF(IASYM)+1
          IPF1 = IPOF(IBSYM)+1
          if (IASYM > IBSYM) then
            call MTRANS(AIBJ(IPF1),AJBI(IPF),NVIRA,NVIRB)
          else if (NSIJ /= 1) then
            call MTRANS(ABIJ(IPF1),ABIJ(IPF),NVIRA,NVIRB)
            call MTRANS(AIBJ(IPF1),AJBI(IPF),NVIRA,NVIRB)
          else
            call SQUAR2(ABIJ(IPF),NVIRA)
            if (NI == NJ) call SQUAR2(AIBJ(IPF),NVIRA)
            call MTRANS(AIBJ(IPF),AJBI(IPF),NVIRA,NVIRB)
          end if
        end do
      end if
    end if
  end do
end do

call CSCALE(INDX,INTSYM,C,SQ2INV)
call CSCALE(INDX,INTSYM,S,SQ2)
call mma_deallocate(Buf)
call mma_deallocate(iBuf)

return

end subroutine FAIBJ
