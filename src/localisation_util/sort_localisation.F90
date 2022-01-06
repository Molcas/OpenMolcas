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
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Sort_Localisation(CMO,nBas,nOcc,nFro,nSym)
! Thomas Bondo Pedersen, November 2005.
!
! Purpose: sort CMOs according to Cholesky orbital ordering.

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nOcc(nSym), nFro(nSym)
real(kind=wp), intent(inout) :: CMO(*)
#include "WrkSpc.fh"
integer(kind=iwp) :: iComp, iOpt, ipDen, ipOaux, ipOvlp, ipScr, ipU, ipX, irc, iSyLbl, iSym, k1, kC, kS, kSqr, kTri, kX, lDen, &
                     lOAux, lOvlp, lScr, lU, lX
real(kind=wp) :: ThrCho, xNrm
character(len=80) :: Txt
character(len=8) :: Label
character(len=*), parameter :: SecNam = 'Sort_Localisation'

! Static setting of decomposition threshold.
! ------------------------------------------

ThrCho = 1.0e-12_wp

! Get a copy of the occupied orbitals: X=CMO.
! -------------------------------------------

lX = nBas(1)*nOcc(1)
do iSym=2,nSym
  lX = lX+nBas(iSym)*nOcc(iSym)
end do
call GetMem('XCho','Allo','Real',ipX,lX)
k1 = 1
kX = ipX
do iSym=1,nSym
  kC = k1+nBas(iSym)*nFro(iSym)
  call dCopy_(nBas(1)*nOcc(1),CMO(kC),1,Work(kX),1)
  k1 = k1+nBas(iSym)**2
  kX = kX+nBas(iSym)*nOcc(iSym)
end do

! Get the overlap matrix.
! -----------------------

lOAux = nBas(1)*(nBas(1)+1)/2
lOvlp = nBas(1)*nBas(1)
do iSym=1,nSym
  lOaux = lOaux+nBas(iSym)*(nBas(iSym)+1)/2
  lOvlp = lOvlp+nBas(iSym)*nBas(iSym)
end do
lOaux = lOaux+4
call GetMem('Ovlp','Allo','Real',ipOvlp,lOvlp)
call GetMem('AuxOvlp','Allo','Real',ipOaux,lOaux)

irc = -1
iOpt = 2
iComp = 1
iSyLbl = 1
Label = 'Mltpl  0'
call RdOne(irc,iOpt,Label,iComp,Work(ipOaux),iSyLbl)
if (irc /= 0) then
  write(u6,*) SecNam,': RdOne returned ',irc
  write(u6,*) 'Label = ',Label,'  iSyLbl = ',iSyLbl
  call SysAbendMsg(SecNam,'I/O error in RdOne',' ')
end if

kTri = ipOaux
kSqr = ipOvlp
do iSym=1,nSym
  call Tri2Rec(Work(kTri),Work(kSqr),nBas(iSym),.false.)
  kTri = kTri+nBas(iSym)*(nBas(iSym)+1)/2
  kSqr = kSqr+nBas(iSym)*nBas(iSym)
end do
call GetMem('AuxOvlp','Free','Real',ipOaux,lOaux)

! Sort each symmetry block.
! -------------------------

kX = ipX
kC = 1
kS = ipOvlp
do iSym=1,nSym

  ! Cycle loop for empty symmetry blocks.
  ! -------------------------------------

  if ((nBas(iSym) < 1) .or. (nOcc(iSym) < 1)) Go To 100

  ! Allocations.
  ! ------------

  lDen = nBas(iSym)*nBas(iSym)
  lU = nOcc(iSym)*nOcc(iSym)
  lScr = nBas(iSym)*nOcc(iSym)
  call GetMem('SrtDen','Allo','Real',ipDen,lDen)
  call GetMem('SrtU','Allo','Real',ipU,lU)
  call GetMem('SrtScr','Allo','Real',ipScr,lScr)

  ! Cholesky decompose D=C^TC and thus define the ordering.
  ! At this stage, X contains the original MOs (CMO).
  ! After the decomposition, X contains the Cholesky MOs.
  ! -------------------------------------------------------

  call GetDens_Localisation(Work(ipDen),Work(kX),nBas(iSym),nOcc(iSym))
  irc = -1
  call ChoLoc(irc,Work(ipDen),Work(kX),ThrCho,xNrm,nBas(iSym),nOcc(iSym))
  if (irc /= 0) then
    write(u6,*) SecNam,': ChoLoc returned ',irc
    write(u6,*) 'Symmetry block: ',iSym
    write(u6,*) 'Unable to continue...'
    write(Txt,'(A,I6)') 'ChoLoc return code:',irc
    call SysAbendMsg(SecNam,'Density Cholesky decomposition failed!',Txt)
  end if

  ! Compute address in CMO, skipping frozen orbitals.
  ! -------------------------------------------------

  k1 = kC+nBas(iSym)*nFro(iSym)

  ! Compute U = X^TSC.
  ! ------------------

  call GetUmat_Localisation(Work(ipU),Work(kX),Work(kS),CMO(k1),Work(ipScr),lScr,nBas(iSym),nOcc(iSym))

  ! Sort.
  ! -----

  call Sort_Localisation_1(CMO(k1),Work(ipU),nBas(iSym),nOcc(iSym))

  ! Update counters.
  ! ----------------

  kX = kX+nBas(iSym)*nOcc(iSym)
  kC = kC+nBas(iSym)**2
  kS = kS+nBas(iSym)**2

  ! De-allocations.
  ! ---------------

  call GetMem('SrtScr','Free','Real',ipScr,lScr)
  call GetMem('SrtU','Free','Real',ipU,lU)
  call GetMem('SrtDen','Free','Real',ipDen,lDen)

  ! Loop cycling (empty symmetries jump here).
  ! ------------------------------------------

100 continue

end do

! De-allocations.
! ---------------

call GetMem('XCho','Free','Real',ipX,lX)
call GetMem('Ovlp','Free','Real',ipOvlp,lOvlp)

end subroutine Sort_Localisation
