!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!               2021, Ignacio Fdez. Galvan                             *
!***********************************************************************
! 2021: Remove GOTOs

subroutine TAB2(NREF,IOCR,nIOCR,L0,L1,L2,L3,INTNUM,LV,LSYM,ICIALL,IFCORE,ICOR,NONE_,JONE)

use guga_global, only: BL1, BL2, BS1, BS2, BS3, BS4, COUP, COUP1, IA, IAF, IB, IBF, IFIRST, IJ, IJF, ILIM, IPO, IPRINT, IRC, IV0, &
                       IX, IY, JNDX, K0, K1, K2, K3, LN, MXVERT, N, NIORB, S
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nIOCR, INTNUM, LV, LSYM, ICIALL, IFCORE, ICOR(*), NONE_, JONE(*)
integer(kind=iwp), intent(inout) :: NREF, IOCR(nIOCR)
integer(kind=iwp), intent(_OUT_) :: L0(*), L1(*), L2(*), L3(*)
integer(kind=iwp) :: I, I0, I1, IA1, IAC, IAT, IB1, IBMAX, IBS, IBT, IEL, II, IIJ, IIJF, IIM, IIM2, IJD, IJFL, IJFS, IJL, IJR, &
                     IJRL, IJS, IJS1, IL, IN_, INUM, ISTA, ISTOP, ISUM, ITTT, IUT, IUT1, J, J11, J3, J4, JJ, JJ1, JJ2, JL, JMAX, &
                     K, LN1, NAC, NACU, NIJ, NIJ1
real(kind=wp) :: FB, FBB
logical(kind=iwp) :: skip
integer(kind=iwp), allocatable :: IORB(:), K00(:), K11(:), K22(:), K33(:), L00(:), L11(:), L22(:), L33(:)

IEL = 2
if (IFIRST /= 0) IEL = 1

! NUMBER OF ACTIVE ELECTRONS
NAC = N-2*NIORB

! UPPER LIMIT FOR NUMBER OF ELECTRONS IN ACTIVE SPACE
NACU = NAC+IEL

IUT = 0
IB(1) = int(2*S)
IA(1) = int(N-2*S)/2
IJ(LN+1) = 0
IJ(LN) = 1
NIJ = 1
IJR = 1
IJS = 2
IJRL = IJR
call mma_allocate(IORB,MXVERT,label='IORB')
IORB(1) = 0
do II=1,LN
  IIM = LN-II+1-LV
  IAC = IIM-NIORB-1
  IIM2 = (IIM-1)*2

  do

    ! S=0

    INUM = N-2*IA(IJR)-IB(IJR)
    skip = .false.
    if ((INTNUM /= 0) .and. (IIM > 0)) then
      if (IAC < 0) then
        if ((IB(IJR) /= 0) .and. (INUM+IIM2 == N-2)) skip = .true.
      else
        if (IB(IJR) > IAC+3) skip = .true.
      end if
    end if
    if (.not. skip) then
      NIJ = NIJ+1
      IA(NIJ) = IA(IJR)
      IB(NIJ) = IB(IJR)
      IORB(NIJ) = IORB(IJR)+2
      if (IIM > 0) then
        call CHEL(IA(NIJ),IB(NIJ),IIM,IEL,ISTOP)
        if (ISTOP == 1) NIJ = NIJ-1
      end if
    end if

    ! S=1

    if (IB(IJR) /= 0) then
      skip = .false.
      if ((INTNUM /= 0) .and. (IIM > 0)) then
        if (IAC >= 0) then
          if (INUM+1 > NACU) then
            skip = .true.
          else if (INUM+1 == NACU) then
            IBS = IB(IJR)-1
            if ((IBS /= 0) .and. (IBS /= 2)) skip = .true.
          end if
        end if
      end if
      if (.not. skip) then
        NIJ = NIJ+1
        IA(NIJ) = IA(IJR)
        IB(NIJ) = IB(IJR)-1
        IORB(NIJ) = IORB(IJR)+1
        if (IIM > 0) then
          call CHEL(IA(NIJ),IB(NIJ),IIM,IEL,ISTOP)
          if (ISTOP == 1) NIJ = NIJ-1
        end if
      end if
    end if

    ! S=2

    if (IA(IJR) /= 0) then
      skip = .false.
      if ((INTNUM /= 0) .and. (IIM > 0)) then
        if (IAC < 0) then
          if (IB(IJR) >= 2) skip = .true.
        else
          if ((IB(IJR)+1 > IAC+3) .or. (INUM+1 > NACU)) then
            skip = .true.
          else if (INUM+1 == NACU) then
            IBS = IB(IJR)+1
            if ((IBS /= 0) .and. (IBS /= 2)) skip = .true.
          end if
        end if
      end if
      if (.not. skip) then
        NIJ = NIJ+1
        IA(NIJ) = IA(IJR)-1
        IB(NIJ) = IB(IJR)+1
        IORB(NIJ) = IORB(IJR)+1
        if (IIM > 0) then
          call CHEL(IA(NIJ),IB(NIJ),IIM,IEL,ISTOP)
          if (ISTOP == 1) NIJ = NIJ-1
        end if
      end if
    end if

    ! S=3

    if (IA(IJR) /= 0) then
      skip = .false.
      if ((INTNUM /= 0) .and. (IAC >= 0)) then
        if ((IB(IJR) > IAC+3) .or. (INUM+2 > NACU)) then
          skip = .true.
        else if (INUM+2 == NACU) then
          IBS = IB(IJR)
          if ((IBS /= 0) .and. (IBS /= 2)) skip = .true.
        end if
      end if
      if (.not. skip) then
        NIJ = NIJ+1
        IA(NIJ) = IA(IJR)-1
        IB(NIJ) = IB(IJR)
        IORB(NIJ) = IORB(IJR)
        if (IIM > 0) then
          call CHEL(IA(NIJ),IB(NIJ),IIM,IEL,ISTOP)
          if (ISTOP == 1) NIJ = NIJ-1
        end if
      end if
    end if

    if (IJR == IJRL) exit
    IJR = IJR+1
  end do

  ! DELETE VERTICES
  NIJ1 = NIJ-1
  IN_ = IJS
  IUT = IJS
  do IJD=IJS,NIJ1
    JJ1 = NIJ-IJD+IJS-1
    J = JJ1+1
    do K=IJS,JJ1
      if ((IA(J) == IA(K)) .and. (IB(J) == IB(K))) then
        IA(J) = -1
        IB(J) = -1
        exit
      end if
    end do
  end do
  ! PACK VERTICES
  IJS1 = IJS+1
  do J=IJS1,NIJ
    if ((IA(J) == -1) .and. (IB(J) == -1)) then
      IN_ = IN_+1
    else
      IN_ = IN_+1
      IUT = IUT+1
      IA(IUT) = IA(IN_)
      IB(IUT) = IB(IN_)
      IORB(IUT) = IORB(IN_)
    end if
  end do
  ! ORDER VERTICES
  IUT1 = IUT-1
  do J=IJS,IUT1
    J11 = J+1
    do K=J11,IUT
      if (IA(J)-IA(K) <= 0) then
        if ((IA(J)-IA(K) < 0) .or. (IB(J) <= IB(K))) then
          IAT = IA(J)
          IBT = IB(J)
          IA(J) = IA(K)
          IB(J) = IB(K)
          IA(K) = IAT
          IB(K) = IBT
        end if
      end if
    end do
  end do
  if (II /= LN) IJ(LN-II) = IUT
  IJR = IJS
  IJS = IUT+1
  IJRL = IUT
  NIJ = IUT
end do
call mma_deallocate(IORB)
if (N == 2) then
  IUT = IUT+1
  IA(IUT) = 0
  IB(IUT) = 0
  IA(IUT-1) = 0
  IB(IUT-1) = 1
  IA(IUT-2) = 0
  IB(IUT-2) = 2
end if
JJ2 = 0
call mma_allocate(K00,MXVERT,label='K00')
call mma_allocate(K11,MXVERT,label='K11')
call mma_allocate(K22,MXVERT,label='K22')
call mma_allocate(K33,MXVERT,label='K33')
call mma_allocate(L00,MXVERT,label='L00')
call mma_allocate(L11,MXVERT,label='L11')
call mma_allocate(L22,MXVERT,label='L22')
call mma_allocate(L33,MXVERT,label='L33')
do II=1,LN
  I = LN-II+1-LV
  IIM2 = (I-1)*2
  I0 = I+LV
  JJ1 = IJ(I0+1)+1
  JJ2 = IJ(I0)
  J3 = JJ2+1
  if (I0 /= 1) J4 = IJ(I0-1)
  if (I0 == 1) J4 = IUT
  ! DETERMINE CASE DOWN
  do J=JJ1,JJ2
    IA1 = IA(J)
    IB1 = IB(J)
    INUM = 2*IA1+IB1
    K00(J) = 0
    K11(J) = 0
    K22(J) = 0
    K33(J) = 0
    do JJ=J3,J4
      if (IA1 /= IA(JJ)) then
        if ((IA1-IA(JJ)) == 1) then
          if (IB1 /= IB(JJ)) then
            if ((IB(JJ)-IB1) == 1) then
              if ((I > NIORB) .or. (INTNUM == 0) .or. (I <= 0) .or. (IB1 < 2)) K22(J) = JJ
            end if
          else
            K33(J) = JJ
          end if
        end if
      else if (IB1 /= IB(JJ)) then
        if ((IB1-IB(JJ)) == 1) K11(J) = JJ
      else
        if ((I <= NIORB) .and. (INTNUM /= 0) .and. (I > 0) .and. (IB1 /= 0)) then
          if ((INUM-IIM2 == 2) .or. (IA1 < I-1)) cycle
        end if
        K00(J) = JJ
      end if
    end do
  end do
  ! DETERMINE CASE UP
  do J=J3,J4
    IA1 = IA(J)
    IB1 = IB(J)
    INUM = 2*IA1+IB1
    L00(J) = 0
    L11(J) = 0
    L22(J) = 0
    L33(J) = 0
    do JJ=JJ1,JJ2
      if (IA(JJ) /= IA1) then
        if ((IA(JJ)-IA1) == 1) then
          if (IB(JJ) /= IB1) then
            if ((IB1-IB(JJ)) == 1) then
              if ((I > NIORB) .or. (INTNUM == 0) .or. (I <= 0) .or. (IB1 < 3)) L22(J) = JJ
            end if
          else
            L33(J) = JJ
          end if
        end if
      else if (IB(JJ) /= IB1) then
        if ((IB(JJ)-IB1) == 1) L11(J) = JJ
      else
        if ((I <= NIORB) .and. (INTNUM /= 0) .and. (I > 0) .and. (IB1 /= 0)) then
          if ((INUM-IIM2 == 2) .or. (IA1 < I-1)) cycle
        end if
        L00(J) = JJ
      end if
    end do
  end do
end do
IV0 = IUT
K00(IUT) = 0
K11(IUT) = 0
K22(IUT) = 0
K33(IUT) = 0
K00(IUT+1) = 0
K11(IUT+1) = 0
K22(IUT+1) = 0
K33(IUT+1) = 0
if (ICIALL /= 0) call CIALL(LSYM,NREF,IOCR,nIOCR,L00,L11,L22,L33,LV)
call DELTAB(NREF,IOCR,L0,L1,L2,L3,INTNUM,LV,IFCORE,ICOR,NONE_,JONE,K00,K11,K22,K33,L00,L11,L22,L33)
call mma_deallocate(K00)
call mma_deallocate(K11)
call mma_deallocate(K22)
call mma_deallocate(K33)
call mma_deallocate(L00)
call mma_deallocate(L11)
call mma_deallocate(L22)
call mma_deallocate(L33)
do I=1,ILIM
  ISTA = (I-1)*MXVERT
  if (IPRINT >= 5) write(u6,101)
  if (IPRINT >= 5) write(u6,100) (J,IA(J),IB(J),K0(ISTA+J),K1(ISTA+J),K2(ISTA+J),K3(ISTA+J),L0(ISTA+J),L1(ISTA+J),L2(ISTA+J), &
                                  L3(ISTA+J),J=1,IUT)
end do
IBMAX = 0
do J=1,IUT
  if (IB(J) > IBMAX) IBMAX = IB(J)
end do
call mma_allocate(BS1,IBMAX+3,label='BS1')
call mma_allocate(BS2,IBMAX+3,label='BS2')
call mma_allocate(BS3,IBMAX+3,label='BS3')
call mma_allocate(BS4,IBMAX+3,label='BS4')
call mma_allocate(BL1,IBMAX+3,label='BL1')
call mma_allocate(BL2,IBMAX+3,label='BL2')
IUT1 = IUT-1
do IL=1,ILIM
  ISTA = (IL-1)*MXVERT
  IX(ISTA+1) = 1
  do II=1,LN
    if (II /= LN) then
      I = LN-II
      IJL = IJ(I)
      IJS = IJ(I+1)+1
    else
      IJL = IUT
      IJS = IUT-3
      if (IFIRST /= 0) IJS = IUT-1
    end if
    do J=IJS,IJL
      ISUM = 0
      if (L0(ISTA+J) /= 0) ISUM = ISUM+IX(ISTA+L0(ISTA+J))
      if (L1(ISTA+J) /= 0) then
        IY(ISTA+L1(ISTA+J),1) = ISUM
        ISUM = ISUM+IX(ISTA+L1(ISTA+J))
      end if
      if (L2(ISTA+J) /= 0) then
        IY(ISTA+L2(ISTA+J),2) = ISUM
        ISUM = ISUM+IX(ISTA+L2(ISTA+J))
      end if
      if (L3(ISTA+J) /= 0) then
        IY(ISTA+L3(ISTA+J),3) = ISUM
        ISUM = ISUM+IX(ISTA+L3(ISTA+J))
      end if
      IX(ISTA+J) = ISUM
    end do
  end do
end do
do I=1,ILIM
  ISTA = (I-1)*MXVERT
  if (IPRINT >= 5) write(u6,102)
  if (IPRINT >= 5) write(u6,200) (J,IY(ISTA+J,1),IY(ISTA+J,2),IY(ISTA+J,3),IX(ISTA+J),J=1,IUT)
end do
if (IPRINT >= 2) write(u6,210) IUT
write(u6,214)

if (IFIRST == 0) then
! This is a normal calculation, both singles and doubles included.
  write(u6,215) (IX(IUT+1-ITTT+MXVERT*(ITTT-1)),ITTT=1,4)
else
  ! "FIRST" keyword has been given. Then this is just a so-called
  ! first-order CI, i.e., singles only.
  write(u6,215) (IX(IUT+1-ITTT+MXVERT*(ITTT-1)),ITTT=1,2)
end if

IRC(1) = IX(IUT)
do I=2,ILIM
  ISTA = (I-1)*MXVERT
  IRC(I) = IX(ISTA+IUT+1-I)+IRC(I-1)
end do
ISUM = IRC(ILIM)
call mma_allocate(JNDX,ISUM,label='JNDX')
do I=1,IBMAX+3
  FBB = real(I-1,kind=wp)
  FB = FBB/(FBB+1)
  BS1(I) = sqrt(FB)
  FB = (FBB+2)/(FBB+1)
  BS2(I) = sqrt(FB)
  if (I > 1) BS3(I) = One/BS1(I)
  BS4(I) = One/BS2(I)
  FB = FBB*FBB-1
  if (I > 1) BL1(I) = sqrt(FB)/FBB
  FB = (FBB+2)**2-1
  BL2(I) = sqrt(FB)/(FBB+2)
end do
! PUT ZEROS IN VECTORS
COUP(1:LN) = Zero
COUP1(1:LN) = Zero
IN_ = 0
do I=1,LN
  II = LN-I+1
  IJS = IJ(II+1)+1
  IJL = IJ(II)
  IJFS = IJF(II+1)+1
  IJFL = IJF(II)
  do IIJ=IJS,IJL
    do IIJF=IJFS,IJFL
      if ((IA(IIJ) == IAF(IIJF)) .and. (IB(IIJ) == IBF(IIJF))) then
        IPO(IIJ) = IIJF
        exit
      end if
    end do
  end do
end do
IPO(IUT) = IJF(1)+1
if (IPRINT >= 10) write(u6,350) (IPO(J),J=1,IUT)
JMAX = 0
LN1 = LN+1
do I=2,LN1
  I1 = I-1
  JL = IJ(I1)-IJ(I)
  if (JL >= JMAX) JMAX = JL
end do

if (IPRINT >= 2) then
  write(u6,411) JMAX
  write(u6,412) IBMAX
  if (JMAX > 31) call Abend()
end if

return

100 format(6X,I3,5X,2I4,5X,8I4)
101 format(///,6X,'TAB2',//,8X,'J',8X,'A',3X,'B',7X,'K0',2X,'K1',2X,'K2',2X,'K3',2X,'L0',2X,'L1',2X,'L2',2X,'L3',/)
102 format(///,6X,'INDEX TABLE',//,8X,'J',8X,'Y1',3X,'Y2',3X,'Y3',9X,'X',/)
200 format(6X,I3,5X,3I5,5X,I5)
210 format(/,6X,'NUMBER OF VERTICES',I10)
214 format(///,6X,'INTERNAL CONFIGURATIONS (FORMAL)')
215 format(/,6X,'NUMBER OF VALENCE STATES',I16,/,6X,'NUMBER OF DOUBLET COUPLED SINGLES',I7,/,6X, &
           'NUMBER OF TRIPLET COUPLED DOUBLES',I7,/,6X,'NUMBER OF SINGLET COUPLED DOUBLES',I7)
!216 format(/,6X,'NUMBER OF VALENCE STATES',I16,/,6X,'NUMBER OF DOUBLET COUPLED SINGLES',I7)
350 format(/,6X,5I5)
411 format(/,6X,'NUMBER OF VERTICES IN ONE ROW',I6,6X,'(PRESENT LIMIT 31)')
412 format(6X,'MAXIMUM B VALUE',I20,6X)

end subroutine TAB2
