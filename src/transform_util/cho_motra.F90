!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Cho_MOtra(CMO,nCMOs,Do_int,ihdf5)

implicit none
integer nCMOs, ihdf5
real*8 CMO(nCMOs)
logical Do_int
character*6 BName
integer iSym
integer nSym
integer nBas(8), nOrb(8)
integer nFro(8), nIsh(8), nAsh(8)
integer nSsh(8), nDel(8)
logical InitChoEnv

call Get_iScalar('nSym',nSym)
call Get_iArray('nBas',nBas,nSym)
call Get_iArray('nOrb',nOrb,nSym)
call Get_iArray('nFro',nFro,nSym)
call Get_iArray('nIsh',nIsh,nSym)
call Get_iArray('nAsh',nAsh,nSym)
call Get_iArray('nDel',nDel,nSym)
do iSym=1,nSym
  nSsh(iSym) = nBas(iSym)-nDel(iSym)-nAsh(iSym)-nIsh(iSym)-nFro(iSym)
end do
BName = '_CHMOT'
InitChoEnv = .true.
call Cho_MOTra_Inner(CMO,nCMOs,nSym,nBas,nOrb,nFro,nIsh,nAsh,nSsh,nDel,BName,Do_Int,ihdf5,InitChoEnv)

end subroutine Cho_MOTra
