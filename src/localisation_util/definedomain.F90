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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: nAtom, nOcc, nBas_per_Atom(nAtom), nBas_Start(nAtom), nBas
integer(kind=iwp), intent(inout) :: iDomain(0:nAtom,nOcc)
real(kind=wp), intent(out) :: QD(nOcc), f(nOcc)
real(kind=wp), intent(in) :: C(nBas,nOcc), ThrDomain(2)
integer(kind=iwp) :: i, iAt, iAtom, iCount, nB(1), nErr, nSrt
real(kind=wp) :: Charge, Chrg, Diff
integer(kind=iwp), allocatable :: iPivot(:)
real(kind=wp), allocatable :: absQ(:), Q(:,:), S(:,:), T(:,:)
#ifdef _DEBUGPRINT_
#define DBG .true.
#else
#define DBG .false.
#endif
logical(kind=iwp), parameter :: LocDbg = DBG
real(kind=wp), external :: ddot_

! Check input.
! ------------

irc = 0
if ((nAtom < 1) .or. (nBas < 1) .or. (nOcc < 1)) return

! Allocate and read overlap matrix stored as full square.
! -------------------------------------------------------

call mma_allocate(S,nBas,nBas,label='DfDm_S')

nB(1) = nBas
call GetOvlp_Localisation(S,'Sqr',nB,1)

! Allocations.
! ------------

call mma_allocate(T,nBas,nOcc,label='DfDm_T')
call mma_allocate(Q,nAtom,nOcc,label='DfDm_Q')

! Compute T=SC.
! -------------

call DGEMM_('N','N',nBas,nOcc,nBas,One,S,nBas,C,nBas,Zero,T,nBas)

! Compute atomic contributions to Mulliken charges.
! -------------------------------------------------

Q(:,:) = Zero
do i=1,nOcc
  do iAtom=1,nAtom
    Q(iAtom,i) = Q(iAtom,i)+dDot_(nBas_per_Atom(iAtom),C(nBas_Start(iAtom),i),1,T(nBas_Start(iAtom),i),1)
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
  do i=1,nOcc
    Chrg = Zero
    do iAtom=1,nAtom
      Chrg = Chrg+Q(iAtom,i)
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
    call FreeMem()
    return
  end if
end if

! For each orbital, create an index array in descending order of
! absolute contributions to charges. Array iDomain will then contain
! the list of atoms ordered according to contributions.
! ------------------------------------------------------------------

call mma_allocate(iPivot,nAtom,label='DfDm_iPivot')
call mma_allocate(absQ,nAtom,label='DfDm_absQ')
do i=1,nOcc
  nSrt = nAtom
  do iAtom=1,nAtom
    absQ(iAtom) = abs(Q(iAtom,i))
  end do
  call CD_DiaMax(absQ,nAtom,iPivot,iDomain(1,i),nSrt,-One)
  if (nSrt /= nAtom) then
    call mma_deallocate(iPivot)
    irc = 1 ! ooops: something is fishy here...
    call FreeMem()
    return
  end if
end do
call mma_deallocate(absQ)
call mma_deallocate(iPivot)

! For each orbital, define domain according to charge threshold.
! --------------------------------------------------------------

do i=1,nOcc
  iCount = 1
  iAtom = iDomain(iCount,i)
  Charge = Q(iAtom,i)
  do while ((iCount < nAtom) .and. (Charge < ThrDomain(1)))
    iCount = iCount+1
    iAtom = iDomain(iCount,i)
    Charge = Charge+Q(iAtom,i)
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
      do iAt=1,iDomain(0,i)
        iAtom = iDomain(iAt,i)
        if ((iAtom < 1) .or. (iAtom > nAtom)) then
          write(u6,*) '  Atom: ',iAtom,' !?!?!'
          nErr = nErr+1
        else
          write(u6,*) '  Atom: ',iAtom,'  Charge: ',Q(iAtom,i)
          Charge = Charge+Q(iAtom,i)
        end if
      end do
      write(u6,*) '  Total charge: ',Charge
    end if
  end do
  if (nErr /= 0) then
    irc = 3
    call FreeMem()
    return
  end if
end if

! For each orbital, check completeness and add atoms as needed to
! meet the requirement f<=threshold.
! ---------------------------------------------------------------

if (ThrDomain(2) < One) then
  do i=1,nOcc
    call MakeDomainComplete(iDomain(0,i),f(i),S,T(:,i),ThrDomain(2),nBas_per_Atom,nBas_Start,nBas,nAtom)
  end do
end if

! Compute total charges for each domain.
! --------------------------------------

do i=1,nOcc
  QD(i) = Zero
  do iAt=1,iDomain(0,i)
    iAtom = iDomain(iAt,i)
    QD(i) = QD(i)+Q(iAtom,i)
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
        do iAt=1,iDomain(0,i)
          iAtom = iDomain(iAt,i)
          if ((iAtom < 1) .or. (iAtom > nAtom)) then
            write(u6,*) '  Atom: ',iAtom,' !?!?!'
            nErr = nErr+1
          else
            write(u6,*) '  Atom: ',iAtom,'  Charge: ',Q(iAtom,i)
            Charge = Charge+Q(iAtom,i)
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
    call FreeMem()
    return
  end if
end if

call FreeMem()

return

contains

subroutine FreeMem()

  ! Deallocations.
  ! --------------

  call mma_deallocate(S)
  call mma_deallocate(T)
  call mma_deallocate(Q)

end subroutine FreeMem

end subroutine DefineDomain
