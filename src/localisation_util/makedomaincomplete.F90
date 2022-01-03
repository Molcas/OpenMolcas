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

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nAtom, iDomain(0:nAtom), nBas_per_Atom(nAtom), nBas_Start(nAtom), nBas
real(kind=wp) :: f, S(nBas,nBas), T(nBas), Threshold
#include "WrkSpc.fh"
integer(kind=iwp) :: iA, iB, iCol, ip_Scr, ip_Si, ip_Sl, ip_Ti, irc, iRow, kSi, kTi, l_Scr, l_Si, l_Sl, l_Ti, lnu, mu1, nA, nmu, &
                     nnu, nSize, nu, nu1
character(len=80) :: Txt
logical(kind=iwp) :: Complete
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
  l_Si = nSize*nSize
  l_Sl = l_Si
  l_Ti = nSize
  l_Scr = nSize
  call GetMem('MkDmC_Si','Allo','Real',ip_Si,l_Si)
  call GetMem('MkDmC_Sl','Allo','Real',ip_Sl,l_Sl)
  call GetMem('MkDmC_Ti','Allo','Real',ip_Ti,l_Ti)
  call GetMem('MkDmC_Scr','Allo','Real',ip_Scr,l_Scr)

  ! Get S[i] and T[i].
  ! ------------------

  kTi = ip_Ti
  iCol = 0
  do iB=1,nA
    nu1 = nBas_Start(iDomain(iB))
    nnu = nBas_per_Atom(iDomain(iB))
    do lnu=0,nnu-1
      nu = nu1+lnu
      iCol = iCol+1
      iRow = 0
      kSi = ip_Si+nSize*(iCol-1)
      do iA=1,nA
        mu1 = nBas_Start(iDomain(iA))
        nmu = nBas_per_Atom(iDomain(iA))
        call dCopy_(nmu,S(mu1,nu),1,Work(kSi+iRow),1)
        iRow = iRow+nmu
      end do
    end do
    call dCopy_(nnu,T(nu1),1,Work(kTi),1)
    kTi = kTi+nnu
  end do

  ! Solve S[i]Y=T[i] (Y stored in T[i] on exit).
  ! Use a scratch array for S[i], as it is ruined on exit from
  ! LinEqSolv.
  ! ----------------------------------------------------------

  irc = 0
  call dCopy_(l_Si,Work(ip_Si),1,Work(ip_Sl),1)
  call LinEqSolv(irc,'N',Work(ip_Sl),nSize,Work(ip_Ti),nSize,nSize,1)
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

  call dGeMV_('N',nSize,nSize,One,Work(ip_Si),nSize,Work(ip_Ti),1,Zero,Work(ip_Scr),1)
  f = One-dDot_(nSize,Work(ip_Ti),1,Work(ip_Scr),1)

  ! Deallocation.
  ! -------------

  call GetMem('MkDmC_Scr','Free','Real',ip_Scr,l_Scr)
  call GetMem('MkDmC_Ti','Free','Real',ip_Ti,l_Ti)
  call GetMem('MkDmC_Sl','Free','Real',ip_Sl,l_Sl)
  call GetMem('MkDmC_Si','Free','Real',ip_Si,l_Si)

  ! Check completeness (f<=Threshold).
  ! If not complete, add next atom to domain.
  ! (If complete, we break the while loop.)
  ! -----------------------------------------

  Complete = f <= Threshold
  if (.not. Complete) then
    nA = nA+1
  end if

end do

! Set new domain size.
! --------------------

iDomain(0) = nA

end subroutine MakeDomainComplete
