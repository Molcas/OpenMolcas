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
!  DefineDomain
!
!> @brief
!>   Define orbital domains
!> @author Thomas Bondo Pedersen
!>
!> @details
!> Orbital domains are defined by summing up gross atomic Mulliken
!> charges for each orbital. Once the sum is greater than or equal
!> to the threshold \p ThrDomain(1), the domain is defined. The
!> second threshold \p ThrDomain(2) is then used to check
!> completeness (Pulay-style) of the definition and if needed,
!> more atoms are added to the domain. To avoid the second step
!> (i.e. the completeness check), simply put \p ThrDomain(2) > ``1.0`` in
!> which case array \p f is undefined on exit.
!>
!> On exit, the contents of iDomain array are:
!>
!> - \p iDomain(0,i): number of atoms in domain \c i.
!> - \p iDomain(n,i): id of atom \c n (\c n = ``1``, ``2``, ..., \p iDomain(0,i)) in domain \c i.
!>
!> Return codes:
!>
!> - \p irc = ``-1``: input error(s) detected.
!> - \p irc =  ``0``: all OK.
!> - \p irc =  ``1``: this indicates a bug in the charge setup.
!> - \p irc =  ``2``: can only be set if debug is turned on; in that case, this code means
!>                    that the computed charges do not sum up to the number of occupied orbitals, \p nOcc.
!> - \p irc =  ``3``: can only be set if debug is turned on; in that case, this code means
!>                    that at least one domain is defined as empty or having at least one atom index out of bounds.
!>
!> @param[out] irc           Return code
!> @param[out] iDomain       Domain definition
!> @param[out] QD            Final total charges for each domain
!> @param[out] f             Function values from completeness check
!> @param[in]  C             MO coefficients
!> @param[in]  ThrDomain     Thresholds for domain definition
!> @param[in]  nBas_per_Atom Number of basis functions on each atom
!> @param[in]  nBas_Start    Start index for basis functions on each atom
!> @param[in]  nAtom         Number of atoms
!> @param[in]  nBas          Number of basis functions
!> @param[in]  nOcc          Number of occupied orbitals
!***********************************************************************

subroutine DefineDomain(irc,iDomain,QD,f,C,ThrDomain,nBas_per_Atom,nBas_Start,nAtom,nBas,nOcc)

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6, r8

implicit none
integer(kind=iwp) :: irc, nAtom, nOcc, iDomain(0:nAtom,nOcc), nBas_per_Atom(nAtom), nBas_Start(nAtom), nBas
real(kind=wp) :: QD(nOcc), f(nOcc), C(nBas,nOcc), ThrDomain(2)
#include "WrkSpc.fh"
integer(kind=iwp) :: i, iAt, iAtom, iCount, iOff, iOff0, ip_absQ, ip_iPivot, ip_Q, ip_S, ip_T, jOff, jOff0, kOff, kOff0, kOffT, &
                     l_absQ, l_iPivot, l_Q, l_S, l_T, nB(1), nErr, nSrt
real(kind=wp) :: Charge, Chrg, Diff
#if defined (_DEBUGPRINT_)
#define DBG .true.
#else
#define DBG .false.
#endif
logical(kind=iwp), parameter :: LocDbg = DBG
real(kind=r8), external :: ddot_

! Check input.
! ------------

irc = 0
if ((nAtom < 1) .or. (nBas < 1) .or. (nOcc < 1)) return

! Allocate and read overlap matrix stored as full square.
! -------------------------------------------------------

l_S = nBas*nBas
call GetMem('DfDm_S','Allo','Real',ip_S,l_S)

nB(1) = nBas
call GetOvlp_Localisation(Work(ip_S),'Sqr',nB,1)

! Allocations.
! ------------

l_T = nBas*nOcc
l_Q = nAtom*nOcc
call GetMem('DfDm_T','Allo','Real',ip_T,l_T)
call GetMem('DfDm_Q','Allo','Real',ip_Q,l_Q)

! Compute T=SC.
! -------------

call DGEMM_('N','N',nBas,nOcc,nBas,One,Work(ip_S),nBas,C,nBas,Zero,Work(ip_T),nBas)

! Compute atomic contributions to Mulliken charges.
! -------------------------------------------------

call dCopy_(l_Q,[Zero],0,Work(ip_Q),1)
iOff0 = ip_T-1
jOff0 = ip_Q-1
do i=1,nOcc
  iOff = iOff0+nBas*(i-1)
  jOff = jOff0+nAtom*(i-1)
  do iAtom=1,nAtom
    Work(jOff+iAtom) = Work(jOff+iAtom)+dDot_(nBas_per_Atom(iAtom),C(nBas_Start(iAtom),i),1,Work(iOff+nBas_Start(iAtom)),1)
  end do
end do

! Debug: check charges.
! The sum of them should equal the trace of the occupied MO overlap
! matrix, i.e. nOcc.
! -----------------------------------------------------------------

if (LocDbg) then
  write(u6,*)
  write(u6,*) 'DefineDomain: checking charge calculation:'
  Charge = Zero
  kOff0 = ip_Q-1
  do i=1,nOcc
    Chrg = Zero
    kOff = kOff0+nAtom*(i-1)
    do iAtom=1,nAtom
      Chrg = Chrg+Work(kOff+iAtom)
    end do
    Charge = Charge+Chrg
    Diff = Chrg-One
    if (abs(Diff) > 1.0e-10_wp) then
      write(u6,*)
      write(u6,*) '  Orbital ',i,':'
      write(u6,*) '  Charge    : ',Chrg
      write(u6,*) '  Expected  : ',One
      write(u6,*) '  Difference: ',Diff
    end if
  end do
  Diff = Charge-real(nOcc,kind=wp)
  write(u6,*)
  write(u6,*) '  Total charge: ',Charge
  write(u6,*) '  Expected    : ',real(nOcc,kind=wp)
  write(u6,*) '  Difference  : ',Diff
  if (abs(Diff) > 1.0e-10_wp) then
    irc = 2
    Go To 1 ! return after deallocation
  end if
end if

! For each orbital, create an index array in descending order of
! absolute contributions to charges. Array iDomain will then contain
! the list of atoms ordered according to contributions.
! ------------------------------------------------------------------

l_iPivot = nAtom
l_absQ = nAtom
call GetMem('DfDm_iPivot','Allo','Inte',ip_iPivot,l_iPivot)
call GetMem('DfDm_absQ','Allo','Real',ip_absQ,l_absQ)
do i=1,nOcc
  iOff = ip_Q+nAtom*(i-1)
  nSrt = nAtom
  do iAtom=0,nAtom-1
    Work(ip_absQ+iAtom) = abs(Work(iOff+iAtom))
  end do
  call CD_DiaMax(Work(ip_absQ),nAtom,iWork(ip_iPivot),iDomain(1,i),nSrt,-One)
  if (nSrt /= nAtom) then
    call GetMem('DfDm_iPivot','Free','Inte',ip_iPivot,l_iPivot)
    irc = 1 ! ooops: something is fishy here...
    Go To 1 ! return after deallocations
  end if
end do
call GetMem('DfDm_absQ','Free','Real',ip_absQ,l_absQ)
call GetMem('DfDm_iPivot','Free','Inte',ip_iPivot,l_iPivot)

! For each orbital, define domain according to charge threshold.
! --------------------------------------------------------------

iOff0 = ip_Q-1
do i=1,nOcc
  iOff = iOff0+nAtom*(i-1)
  iCount = 1
  iAtom = iDomain(iCount,i)
  Charge = Work(iOff+iAtom)
  do while ((iCount < nAtom) .and. (Charge < ThrDomain(1)))
    iCount = iCount+1
    iAtom = iDomain(iCount,i)
    Charge = Charge+Work(iOff+iAtom)
  end do
  iDomain(0,i) = iCount
end do

! Debug.
! ------

if (LocDbg) then
  nErr = 0
  write(u6,*)
  write(u6,*) 'DefineDomain: domains and charges after step 1:'
  write(u6,*) 'Threshold: ',ThrDomain(1)
  do i=1,nOcc
    write(u6,*)
    write(u6,*) 'Domain ',i,': ',iDomain(0,i),' atoms:'
    if (iDomain(0,i) < 1) then
      write(u6,*) 'No atoms in domain !?!?!'
      nErr = nErr+1
    else if (iDomain(0,i) > nAtom) then
      write(u6,*) 'Number of atoms > nAtom in domain !?!?!'
      nErr = nErr+1
    else
      Charge = Zero
      kOff0 = ip_Q-1+nAtom*(i-1)
      do iAt=1,iDomain(0,i)
        iAtom = iDomain(iAt,i)
        if ((iAtom < 1) .or. (iAtom > nAtom)) then
          write(u6,*) '  Atom: ',iAtom,' !?!?!'
          nErr = nErr+1
        else
          kOff = kOff0+iAtom
          write(u6,*) '  Atom: ',iAtom,'  Charge: ',Work(kOff)
          Charge = Charge+Work(kOff)
        end if
      end do
      write(u6,*) '  Total charge: ',Charge
    end if
  end do
  if (nErr /= 0) then
    irc = 3
    Go To 1 ! return after deallocation
  end if
end if

! For each orbital, check completeness and add atoms as needed to
! meet the requirement f<=threshold.
! ---------------------------------------------------------------

if (ThrDomain(2) < One) then
  do i=1,nOcc
    kOffT = ip_T+nBas*(i-1)
    call MakeDomainComplete(iDomain(0,i),f(i),Work(ip_S),Work(kOffT),ThrDomain(2),nBas_per_Atom,nBas_Start,nBas,nAtom)
  end do
end if

! Compute total charges for each domain.
! --------------------------------------

do i=1,nOcc
  kOff = ip_Q-1+nAtom*(i-1)
  QD(i) = Zero
  do iAt=1,iDomain(0,i)
    iAtom = iDomain(iAt,i)
    QD(i) = QD(i)+Work(kOff+iAtom)
  end do
end do

! Debug.
! ------

if (LocDbg) then
  nErr = 0
  write(u6,*)
  write(u6,*) 'DefineDomain: domains and charges after step 2:'
  write(u6,*) 'Threshold: ',ThrDomain(2)
  if (ThrDomain(2) < One) then
    do i=1,nOcc
      write(u6,*)
      write(u6,*) 'Domain ',i,': ',iDomain(0,i),' atoms:'
      if (iDomain(0,i) < 1) then
        write(u6,*) 'No atoms in domain !?!?!'
        nErr = nErr+1
      else if (iDomain(0,i) > nAtom) then
        write(u6,*) 'Number of atoms > nAtom in domain !?!?!'
        nErr = nErr+1
      else
        Charge = Zero
        kOff0 = ip_Q-1+nAtom*(i-1)
        do iAt=1,iDomain(0,i)
          iAtom = iDomain(iAt,i)
          if ((iAtom < 1) .or. (iAtom > nAtom)) then
            write(u6,*) '  Atom: ',iAtom,' !?!?!'
            nErr = nErr+1
          else
            kOff = kOff0+iAtom
            write(u6,*) '  Atom: ',iAtom,'  Charge: ',Work(kOff)
            Charge = Charge+Work(kOff)
          end if
        end do
        write(u6,*) '  Total charge: ',Charge
        if (abs(Charge-QD(i)) > 1.0e-12_wp) then
          write(u6,*) 'Total charge is inconsistent with QD !'
          nErr = nErr+1
        end if
      end if
    end do
  else
    write(u6,*) 'Threshold >= 1.0: step 2 was skipped.'
    write(u6,*) 'Domains are unchanged from step 1.'
  end if
  if (nErr /= 0) then
    irc = 3
    Go To 1 ! return after deallocation
  end if
end if

! Deallocations.
! --------------

1 call GetMem('DfDm_Q','Free','Real',ip_Q,l_Q)
call GetMem('DfDm_T','Free','Real',ip_T,l_T)
call GetMem('DfDm_S','Free','Real',ip_S,l_S)

end subroutine DefineDomain
