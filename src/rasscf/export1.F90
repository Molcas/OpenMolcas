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
! Copyright (C) 1994, Markus P. Fuelscher                              *
!***********************************************************************

subroutine Export1(iFinal,CMO,DA,PA,DAO,Focc)
!***********************************************************************
!                                                                      *
!     purpose: Save all information relevant to geometry               *
!              optimizations                                           *
!                                                                      *
!     calling arguments:                                               *
!     iFinal  : Switch including routing information                   *
!     CMO     : MO coefficients in last CI                             *
!     DA      : 1-density of active orbitals                           *
!     PA      : 2-density of active orbitals                           *
!                                                                      *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1994                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************
!
use gas_data, only: iDoGAS
use rasscf_global, only: DoDMRG, iRLXRoot, KSDFT, NAC, NACPAR, NACPR2, nRoots, ThrSX, ThrTE, Weight
use general_data, only: NACTEL, NSYM, NHOLE1, NELEC3, NASH, NDEL, NFRO, NISH, NTOT1, NTOT2
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iFinal
real(kind=wp) :: CMO(*), DA(*), PA(*), DAO(*), Focc(*)
integer(kind=iwp) :: i, iR, iRLXRoot1, iRLXRoot2, iS, iSA, nTemp(8), nW
real(kind=wp) :: Dum(1), Tmp
logical(kind=iwp) :: Found, SCF
character(len=16) :: mstate
character(len=8) :: Method, RlxLbl

!----------------------------------------------------------------------*
!     Save information pertinent to the gradient calculation           *
!----------------------------------------------------------------------*
!...  Add elementary information ......................................*
SCF = .false.
if ((nac == 0) .or. (2*nac == nactel)) SCF = .true.

if (SCF) then
  do iS=1,nSym
    nTemp(is) = nAsh(is)+nIsh(is)
  end do
  call Put_iArray('nIsh',nTemp,nSym)
  do iS=1,nSym
    nTemp(is) = 0
  end do
  call Put_iArray('nAsh',nTemp,nSym)
else
  call Put_iArray('nIsh',nIsh,nSym)
  call Put_iArray('nAsh',nAsh,nSym)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Find the correct method label

Method = 'CASSCF  '
if (KSDFT /= 'SCF') Method = 'CASDFT  '

! Set flags for State-Average or Single root

if (nRoots /= 1) then

  ! iSA=-1 non-equivalent multi state SA-CASSCF
  ! iSA=0  equivalent multi state SA-CASSCF
  ! iSA=2  single root SA-CASSCF

  Method = 'CASSCFSA'

  ! Check if equal weight SA-CASSCF

  iSA = 0
  do iR=2,nRoots
    if (Weight(1) /= Weight(iR)) iSA = -1
  end do
  if (iSA /= 0) then

    ! Check if SA-CASSCF is optimized for just one root.

    nW = 0
    do iR=1,nRoots
      if (Weight(iR) /= Zero) nW = nW+1
    end do
    if (nW == 1) iSA = 2
  end if
  call Put_iScalar('SA ready',iSA)
  if ((iSA == 0) .or. (iSA == -1)) then
    mstate = '****************'
    call Put_cArray('MCLR Root',mstate,16)
  end if
end if

! Check if it is a RASSCF function and not a CASSCF
if ((nHole1 /= 0) .or. (nElec3 /= 0)) Method(1:1) = 'R'
! Check if it is a GASSCF function
if (iDoGAS) Method(1:1) = 'G'
! Check if it is a DMRGSCF function
if (doDMRG) then
  Method = 'DMRGSCF '
  if (nroots /= 1) Method = 'DMRGSCFS'
end if

call Put_cArray('Relax Method',Method,8)
!                                                                      *
!***********************************************************************
!                                                                      *
call Get_iScalar('nSym',i)
call Put_iArray('nFro',nFro,i)
call Put_iArray('nDel',nDel,i)
!...  Add MO-coefficients .............................................*
call Put_dArray('Last orbitals',CMO,NTOT2)
!...  Add one body density matrix in AO/SO basis ......................*
call Put_dArray('D1ao',DAO,NTOT1)
!...  Remove the variational density if it exists .....................*
call Put_dArray('D1aoVar',DAO,0)
!...  Add one body density matrix in MO, active orbitals only .........*
call Put_dArray('D1mo',DA,NACPAR)
!...  Add two body density matrix in MO basis, active orbitals only ...*
if (.not. SCF) call Put_dArray('P2mo',PA,NACPR2)
!...  Next version of MOLCAS add the state to relax file ..............*
call Qpg_iScalar('Relax Original root',Found)
if (Found) then
  call Get_iScalar('Relax Original root',irlxroot1)
  call Get_iScalar('Relax CASSCF root',irlxroot2)
  if (irlxroot1 == irlxroot2) call Put_iScalar('Relax Original root',irlxroot)
else
  call Put_iScalar('Relax Original root',irlxroot)
end if
call Put_iScalar('Relax CASSCF root',irlxroot)
!...  Remove overlaps (computed by rassi) .............................*
call Put_darray('State Overlaps',Dum,0)
call Put_lscalar('Track Done',.false.)
!...  Add generalized Fock matrix .....................................*
if (ifinal >= 1) then
  call Put_dArray('FockOcc',Focc,ntot1)
  RlxLbl = 'Thrs    '
  tmp = max(thrte,thrsx)
  call put_dscalar(RlxLbl,tmp)
end if
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

end subroutine Export1
