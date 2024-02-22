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
Subroutine CXClose(SGS,CIS,EXS)
use stdalloc, only: mma_deallocate
use Struct, only: SGStruct, CIStruct, EXStruct
Type (SGStruct) SGS
Type (CIStruct) CIS
Type (ExStruct) ExS
! Unpack structure SGS:
nSym   =SGS%nSym
! Unpack structure CIS:
nMidV =CIS%nMidV
nIpWlk=CIS%nIpWlk
nWalk =CIS%nWalk
! Unpack structure EXS:
MxEO  =EXS%MxEO
lNOCP =EXS%lNOCP
lIOCP =EXS%lIOCP
nICoup=EXS%nICoup
lICoup=EXS%lICoup
nVTab =EXS%nVTab
lVTab =EXS%lVTab
lMVL  =EXS%lMVL
lMVR  =EXS%lMVR
nNOW=2*nMidV*nSym
Call GetMem('MVR','Free','Inte',lMVR,2*nMidV)
Call GetMem('MVL','Free','Inte',lMVL,2*nMidV)
Call mma_deallocate(CIS%NOCSF)
Call mma_deallocate(CIS%IOCSF)
Call mma_deallocate(CIS%NOW)
Call mma_deallocate(CIS%IOW)
nNOCP=MxEO*nMidV*nSym
nIOCP=nNOCP
Call GetMem('NOCP','Free','Inte',lNOCP,nNOCP)
Call GetMem('IOCP','Free','Inte',lIOCP,nIOCP)
Call mma_deallocate(CIS%NCSF)
Call mma_deallocate(CIS%ICase)
nnICoup=(3*nICoup+1)/2
Call GetMem('ICoup','Free','Inte',lICoup,nnICoup)
Call GetMem('VTab','Free','Real',lVtab,nVtab)

end subroutine CXClose
