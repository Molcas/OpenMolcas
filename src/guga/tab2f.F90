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

subroutine TAB2F(IVER,LV)

use guga_global, only: IAF, IBF, IFIRST, IJF, IPRINT, IVF0, K0F, K1F, K2F, K3F, LN, MXVERT, N, S
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IVER, LV
integer(kind=iwp) :: I, IA1, IAT, IB1, IBT, IEL, II, IIM, IJD, IJR, IJRL, IJS, IJS1, IN_, ISTOP, IUT, IUT1, J, J11, J3, J4, JJ, &
                     JJ1, JJ2, K, NIJ, NIJ1, nijj
integer(kind=iwp), allocatable :: IORB(:)

nijj = 0
IEL = 2
if (IFIRST /= 0) IEL = 1
IEL = IEL+1
IUT = 0
IBF(1) = int(2*S)
IAF(1) = int(N-2*S)/2
IJF(LN+1) = 0
IJF(LN) = 1
NIJ = 1
IJR = 1
IJS = 2
IJRL = IJR
call mma_allocate(IORB,MXVERT,label='IORB')
IORB(1) = 0
do II=1,LN
  IIM = LN-II+1-LV
  do
    ! S=0
    NIJ = NIJ+1
    !if (NIJ > IVER)
    nijj = max(nij,nijj)
    IAF(NIJ) = IAF(IJR)
    IBF(NIJ) = IBF(IJR)
    IORB(NIJ) = IORB(IJR)+2
    if (IIM > 0) then
      call CHEL(IAF(NIJ),IBF(NIJ),IIM,IEL,ISTOP)
      if (ISTOP == 1) NIJ = NIJ-1
    end if
    if (IBF(IJR) /= 0) then
      ! S=1
      NIJ = NIJ+1
      !if (NIJ > IVER)
      nijj = max(nij,nijj)
      IAF(NIJ) = IAF(IJR)
      IBF(NIJ) = IBF(IJR)-1
      IORB(NIJ) = IORB(IJR)+1
      if (IIM > 0) then
        call CHEL(IAF(NIJ),IBF(NIJ),IIM,IEL,ISTOP)
        if (ISTOP == 1) NIJ = NIJ-1
      end if
    end if
    if (IAF(IJR) /= 0) then
      ! S=2
      NIJ = NIJ+1
      !if (NIJ > IVER)
      nijj = max(nij,nijj)
      IAF(NIJ) = IAF(IJR)-1
      IBF(NIJ) = IBF(IJR)+1
      IORB(NIJ) = IORB(IJR)+1
      if (IIM > 0) then
        call CHEL(IAF(NIJ),IBF(NIJ),IIM,IEL,ISTOP)
        if (ISTOP == 1) NIJ = NIJ-1
      end if
    end if
    if (IAF(IJR) /= 0) then
      ! S=3
      NIJ = NIJ+1
      !if (NIJ > IVER)
      nijj = max(nij,nijj)
      IAF(NIJ) = IAF(IJR)-1
      IBF(NIJ) = IBF(IJR)
      IORB(NIJ) = IORB(IJR)
      if (IIM > 0) then
        call CHEL(IAF(NIJ),IBF(NIJ),IIM,IEL,ISTOP)
        if (ISTOP == 1) NIJ = NIJ-1
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
      if ((IAF(J) == IAF(K)) .and. (IBF(J) == IBF(K))) then
        IAF(J) = -1
        IBF(J) = -1
        exit
      end if
    end do
  end do
  ! PACK VERTICES
  IJS1 = IJS+1
  do J=IJS1,NIJ
    if ((IAF(J) == -1) .and. (IBF(J) == -1)) then
      IN_ = IN_+1
    else
      IN_ = IN_+1
      IUT = IUT+1
      IAF(IUT) = IAF(IN_)
      IBF(IUT) = IBF(IN_)
      IORB(IUT) = IORB(IN_)
    end if
  end do
  ! ORDER VERTICES
  IUT1 = IUT-1
  do J=IJS,IUT1
    J11 = J+1
    do K=J11,IUT
      if (IAF(J)-IAF(K) <= 0) then
        if ((IAF(J)-IAF(K) < 0) .or. (IBF(J) <= IBF(K))) then
          IAT = IAF(J)
          IBT = IBF(J)
          IAF(J) = IAF(K)
          IBF(J) = IBF(K)
          IAF(K) = IAT
          IBF(K) = IBT
        end if
      end if
    end do
  end do
  if (II /= LN) IJF(LN-II) = IUT
  IJR = IJS
  IJS = IUT+1
  IJRL = IUT
  NIJ = IUT
end do
call mma_deallocate(IORB)
JJ2 = 0
do II=1,LN
  I = LN-II+1
  JJ1 = IJF(I+1)+1
  JJ2 = IJF(I)
  J3 = JJ2+1
  if (I /= 1) J4 = IJF(I-1)
  if (I == 1) J4 = IUT
  ! DETERMINE CASE DOWN
  do J=JJ1,JJ2
    IA1 = IAF(J)
    IB1 = IBF(J)
    K0F(J) = 0
    K1F(J) = 0
    K2F(J) = 0
    K3F(J) = 0
    do JJ=J3,J4
      if (IA1 /= IAF(JJ)) then
        if ((IA1-IAF(JJ)) == 1) then
          if (IB1 /= IBF(JJ)) then
            if ((IBF(JJ)-IB1) == 1) K2F(J) = JJ
          else
            K3F(J) = JJ
          end if
        end if
      else if (IB1 /= IBF(JJ)) then
        if ((IB1-IBF(JJ)) == 1) K1F(J) = JJ
      else
        K0F(J) = JJ
      end if
    end do
  end do
end do
IVF0 = IUT
K0F(IUT) = 0
K1F(IUT) = 0
K2F(IUT) = 0
K3F(IUT) = 0
K0F(IUT+1) = 0
K1F(IUT+1) = 0
K2F(IUT+1) = 0
K3F(IUT+1) = 0
if (IPRINT >= 5) write(u6,101)
if (IPRINT >= 5) write(u6,100) (J,IAF(J),IBF(J),K0F(J),K1F(J),K2F(J),K3F(J),J=1,IUT)
write(u6,*) ' Number of vertices',nijj,iut
if (nijj > iver) then
  write(u6,310) IVER
  call Abend()
end if

return

100 format(6X,I3,5X,2I4,5X,4I4)
101 format(///,6X,'TAB2F',//,8X,'J',7X,'AF',2X,'BF',6X,'K0F',1X,'K1F',1X,'K2F',1X,'K3F',/)
310 format(/,6X,'NUMBER OF VERTICES EXCEEDS',I7)

end subroutine TAB2F
