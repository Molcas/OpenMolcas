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
! Copyright (C) 1991, Manuela Merchan                                  *
!***********************************************************************

subroutine AutoCut()
!***********************************************************************
!                                                                      *
! Purpose:                                                             *
! This routine counts the number of orbitals that have an              *
! occupation number smaller than a given threshold and                 *
! marks them as deleted orbitals.                                      *
!                                                                      *
!**** M. Merchan, University of Valencia, Spain, 1991 ******************

use motra_global, only: CutThrs, nBas, nDel, nFro, nOrb, nOrbt, nOrbtt, nSym, Occ
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: iBas, iDel, ipBas, iSym

! FIXME: Occ was never read or allocated
call Untested('AutoCut')

!----------------------------------------------------------------------*
! Start procedure                                                      *
!----------------------------------------------------------------------*
ipBas = 0
do iSym=1,nSym
  iDel = 0
  do iBas=1,nBas(iSym)
    if (Occ(ipBas+iBas) <= abs(CutThrs(iSym))) iDel = iDel+1
  end do
  ipBas = ipBas+nBas(iSym)
  if (nDel(iSym) < iDel) nDel(iSym) = iDel
  if ((nDel(iSym)+nFro(iSym)) > nBas(iSym)) then
    write(u6,*) 'AutoCut:nDel(iSym)+nFro(iSym)) > nBas(iSym)'
    write(u6,*) 'iSym=',iSym
    write(u6,*) 'nDel(iSym)=',nDel(iSym)
    write(u6,*) 'nFro(iSym)=',nFro(iSym)
    write(u6,*) 'nBas(iSym)=',nBas(iSym)
    call Abend()
  end if
end do
!----------------------------------------------------------------------*
! Normal termination                                                   *
!----------------------------------------------------------------------*
nOrbt = 0
nOrbtt = 0
do iSym=1,nSym
  nOrb(iSym) = nBas(iSym)-nFro(iSym)-nDel(iSym)
  nOrbt = nOrbt+nOrb(iSym)
  nOrbtt = nOrbtt+nOrb(iSym)*(nOrb(iSym)+1)/2
end do

return

end subroutine AutoCut
