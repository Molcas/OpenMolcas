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

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nCMOs, ihdf5
real(kind=wp), intent(in) :: CMO(nCMOs)
logical(kind=iwp), intent(in) :: Do_int
integer(kind=iwp) :: nAsh(8), nBas(8), nDel(8), nFro(8), nIsh(8), nSsh(8), nSym
character(len=6), parameter :: BName = '_CHMOT'
logical(kind=iwp), parameter :: InitChoEnv = .true.

call Get_iScalar('nSym',nSym)
call Get_iArray('nBas',nBas,nSym)
call Get_iArray('nFro',nFro,nSym)
call Get_iArray('nIsh',nIsh,nSym)
call Get_iArray('nAsh',nAsh,nSym)
call Get_iArray('nDel',nDel,nSym)
nSsh(1:nSym) = nBas(1:nSym)-nDel(1:nSym)-nAsh(1:nSym)-nIsh(1:nSym)-nFro(1:nSym)
call Cho_MOTra_Inner(CMO,nCMOs,nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,BName,Do_Int,ihdf5,InitChoEnv)

end subroutine Cho_MOTra
