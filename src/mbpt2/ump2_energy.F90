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
! Copyright (C) 2026, Okan Koeksal                                     *
!***********************************************************************

subroutine UMP2_Energy(La,Lb,Eoa,Eva,Eob,Evb, &
                       nOA,nVA,nOB,nVB,nCho,DenTol, &
                       Eaa,Ebb,Eab,ierr)

use Definitions, only: wp, iwp
use Constants, only: Zero
use Para_Info, only: nProcs, MyRank, Is_Real_Par

implicit none

integer(kind=iwp), intent(in) :: nOA,nVA,nOB,nVB,nCho

real(kind=wp), intent(in) :: La(nVA,nOA,nCho)
real(kind=wp), intent(in) :: Lb(nVB,nOB,nCho)
real(kind=wp), intent(in) :: Eoa(nOA),Eva(nVA)
real(kind=wp), intent(in) :: Eob(nOB),Evb(nVB)
real(kind=wp), intent(in) :: DenTol

real(kind=wp), intent(out) :: Eaa,Ebb,Eab
integer(kind=iwp), intent(out) :: ierr

integer(kind=iwp) :: i,j,a,b,PairIndex,LocalErr
real(kind=wp) :: Coul,Den,WorkE(3)
logical :: Parallel

Eaa = Zero
Ebb = Zero
Eab = Zero
ierr = 0
LocalErr = 0
Parallel = (nProcs > 1) .and. Is_Real_Par()

! A positive denominator tolerance must be supplied.
if (DenTol <= Zero) LocalErr = 1

if (LocalErr == 0) then
  call SameSpin(La,Eoa,Eva,nOA,nVA,Eaa,LocalErr)
end if

if (LocalErr == 0) then
  call SameSpin(Lb,Eob,Evb,nOB,nVB,Ebb,LocalErr)
end if

! Opposite-spin contribution: no exchange term.
!
! MPI-v2 requires La/Lb to contain the COMPLETE globally assembled
! Cholesky-vector dimension on every rank.  The complete Cholesky-vector
! dot product is then evaluated on the rank that owns a particular
! occupied alpha/beta pair.  This preserves the exact serial Cholesky
! contraction while distributing the expensive excitation loops.
if (LocalErr == 0) then
AB_Occupied: do i=1,nOA
    do j=1,nOB

      PairIndex = (i-1)*nOB + (j-1)
      if (Parallel) then
        if (mod(PairIndex,nProcs) /= MyRank) cycle
      end if

      do a=1,nVA
        do b=1,nVB

          Den = Eoa(i)+Eob(j)-Eva(a)-Evb(b)

          if (abs(Den) <= DenTol) then
            LocalErr = 2
            exit AB_Occupied
          end if

          Coul = dot_product(La(a,i,:),Lb(b,j,:))
          Eab = Eab+Coul*Coul/Den

        end do
      end do
    end do
  end do AB_Occupied
end if

! A denominator failure on any rank invalidates the complete energy.
! All ranks must participate in this reduction before any rank returns.
ierr = LocalErr
if (Parallel) call gaIgOP_SCAL(ierr,'max')

if (ierr /= 0) then
  call ClearEnergies()
  return
end if

! Eaa/Ebb/Eab currently contain rank-local partial sums.  Reduce them
! to the complete UMP2 spin components on every rank.
if (Parallel) then
  WorkE(1) = Eaa
  WorkE(2) = Ebb
  WorkE(3) = Eab
  call GADGOp(WorkE(1),3,'+')
  Eaa = WorkE(1)
  Ebb = WorkE(2)
  Eab = WorkE(3)
end if

contains

subroutine SameSpin(L,Eo,Ev,nO,nV,Ess,Err)

  integer(kind=iwp), intent(in) :: nO,nV
  real(kind=wp), intent(in) :: L(nV,nO,nCho)
  real(kind=wp), intent(in) :: Eo(nO),Ev(nV)
  real(kind=wp), intent(out) :: Ess
  integer(kind=iwp), intent(inout) :: Err

  integer(kind=iwp) :: ii,jj,aa,bb,Pair
  real(kind=wp) :: Jint,Kint,Delta,Anti

  Ess = Zero

  ! No same-spin double excitations in these cases.
  if ((nO < 2) .or. (nV < 2)) return

SS_Occupied: do ii=1,nO
    do jj=1,nO
      if (ii == jj) cycle

      Pair = (ii-1)*nO + (jj-1)
      if (Parallel) then
        if (mod(Pair,nProcs) /= MyRank) cycle
      end if

      do aa=1,nV
        do bb=1,nV
          if (aa == bb) cycle

          Delta = Eo(ii)+Eo(jj)-Ev(aa)-Ev(bb)

          if (abs(Delta) <= DenTol) then
            Err = 2
            exit SS_Occupied
          end if

          Jint = dot_product(L(aa,ii,:),L(bb,jj,:))
          Kint = dot_product(L(bb,ii,:),L(aa,jj,:))
          Anti = Jint-Kint

          ! Full ordered sums require the factor 1/4.
          Ess = Ess+0.25_wp*Anti*Anti/Delta

        end do
      end do
    end do
  end do SS_Occupied

end subroutine SameSpin

subroutine ClearEnergies()

  Eaa = Zero
  Ebb = Zero
  Eab = Zero

end subroutine ClearEnergies

end subroutine UMP2_Energy
