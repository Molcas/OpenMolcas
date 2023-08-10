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

subroutine Cho_PrintPointers(irc,Lunit)
!
! Purpose: print entries in the Cholesky module

use Cholesky, only: iAtomShl, iBasSh, iiBstRSh, IndRed, IndRSh, InfRed, InfVec, IntMap, iQuAB, iRS2F, iScr, iShlSO, iSOShl, iSP2F, &
                    nBasSh, nBstSh, nDimRS, nnBstRSh
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: irc, Lunit

write(Lunit,*) '*** Contents of Cholesky:'
write(Lunit,*) '    (dimension)'
write(Lunit,*)
call XFlush(Lunit)
write(Lunit,*) 'InfRed  : ',size(InfRed)
write(Lunit,*) 'InfVec  : ',size(InfVec)
write(Lunit,*) 'IndRed  : ',size(IndRed)
write(Lunit,*) 'IndRSh  : ',size(IndRsh)
write(Lunit,*) 'iScr    : ',size(iScr)
write(Lunit,*) 'iiBstRSh: ',size(iiBstRSh)
write(Lunit,*) 'nnBstRSh: ',size(nnBstRSh)
write(Lunit,*) 'IntMap  : ',size(IntMap)
write(Lunit,*) 'nDimRS  : ',size(nDimRS)
write(Lunit,*) 'iRS2F   : ',size(iRS2F)
write(Lunit,*) 'iSOShl  : ',size(iSOShl)
write(Lunit,*) 'iShlSO  : ',size(iShlSO)
write(Lunit,*) 'iQuab   : ',size(iQuab)
write(Lunit,*) 'iBasSh  : ',size(iBasSh)
write(Lunit,*) 'nBasSh  : ',size(nBasSh)
write(Lunit,*) 'nBstSh  : ',size(nBstSh)
write(Lunit,*) 'iAtomShl: ',size(iAtomShl)
write(Lunit,*) 'iSP2F   : ',size(iSP2F)
write(Lunit,*)
call XFlush(Lunit)

irc = 0

end subroutine Cho_PrintPointers
