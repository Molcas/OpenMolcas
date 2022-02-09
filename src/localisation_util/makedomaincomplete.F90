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
! Copyright (C) 2006, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine MakeDomainComplete(iDomain,f,S,T,Threshold,nBas_per_Atom,nBas_Start,nBas,nAtom)
! Thomas Bondo Pedersen, January 2006.
!
! Purpose: Boughton-Pulay completeness check of an orbital domain.
!          Extend domain as needed to obtain completeness.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtom, nBas_per_Atom(nAtom), nBas_Start(nAtom), nBas
integer(kind=iwp), intent(inout) :: iDomain(0:nAtom)
real(kind=wp), intent(out) :: f
real(kind=wp), intent(in) :: S(nBas,nBas), T(nBas), Threshold
integer(kind=iwp) :: iA, iB, iCol, irc, iRow, kTi, lnu, mu1, nA, nmu, nnu, nSize, nu, nu1
character(len=80) :: Txt
logical(kind=iwp) :: Complete
real(kind=wp), allocatable :: Scr(:), Si(:,:), Sl(:,:), Ti(:)
real(kind=wp), external :: ddot_
character(len=*), parameter :: SecNam = 'MakeDomainComplete'

nA = iDomain(0)
f = Zero
Complete = nA == nAtom
do while ((nA < nAtom) .and. (.not. Complete))

  ! Allocate S[i], T[i], and Scr[i].
  ! --------------------------------

  nSize = nBas_per_Atom(iDomain(1))
  do iA=2,nA
    nSize = nSize+nBas_per_Atom(iDomain(iA))
  end do
  call mma_allocate(Si,nSize,nSize,label='MkDmC_Si')
  call mma_allocate(Sl,nSize,nSize,label='MkDmC_Sl')
  call mma_allocate(Ti,nSize,label='MkDmC_Ti')
  call mma_allocate(Scr,nSize,label='MkDmC_Scr')

  ! Get S[i] and T[i].
  ! ------------------

  kTi = 1
  iCol = 0
  do iB=1,nA
    nu1 = nBas_Start(iDomain(iB))
    nnu = nBas_per_Atom(iDomain(iB))
    do lnu=0,nnu-1
      nu = nu1+lnu
      iCol = iCol+1
      iRow = 0
      do iA=1,nA
        mu1 = nBas_Start(iDomain(iA))
        nmu = nBas_per_Atom(iDomain(iA))
        Si(:,iCol) = S(mu1:mu1+nmu-1,nu)
        iRow = iRow+nmu
      end do
    end do
    Ti(kTi:kTi+nnu-1) = T(nu1:nu1+nnu-1)
    kTi = kTi+nnu
  end do

  ! Solve S[i]Y=T[i] (Y stored in T[i] on exit).
  ! Use a scratch array for S[i], as it is ruined on exit from
  ! LinEqSolv.
  ! ----------------------------------------------------------

  irc = 0
  Sl(:,:) = Si
  call LinEqSolv(irc,'N',Sl,nSize,Ti,nSize,nSize,1)
  if (irc /= 0) then
    write(Txt,'(A,I9)') 'LinEqSolv returned',irc
    if (irc < 0) then
      call SysAbendMsg(SecNam,Txt,'LinEqSolv input error!')
    else
      call SysAbendMsg(SecNam,Txt,'Singular domain overlap matrix!')
    end if
  end if

  ! Compute f=1-Y(T)S[i]Y.
  ! ----------------------

  call dGeMV_('N',nSize,nSize,One,Si,nSize,Ti,1,Zero,Scr,1)
  f = One-dDot_(nSize,Ti,1,Scr,1)

  ! Deallocation.
  ! -------------

  call mma_deallocate(Si)
  call mma_deallocate(Sl)
  call mma_deallocate(Ti)
  call mma_deallocate(Scr)

  ! Check completeness (f<=Threshold).
  ! If not complete, add next atom to domain.
  ! (If complete, we break the while loop.)
  ! -----------------------------------------

  Complete = f <= Threshold
  if (.not. Complete) nA = nA+1

end do

! Set new domain size.
! --------------------

iDomain(0) = nA

end subroutine MakeDomainComplete
