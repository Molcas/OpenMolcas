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
Subroutine CXClose(iSGStruct,iCIStruct,iXStruct)
use Struct, only: nSGSize, nCISize, nXSize
Dimension iSGStruct(nSGSize)
Dimension iCIStruct(nCISize)
Dimension iXStruct (nXSize)
! Unpack structure iSGStruct:
nSym   =iSGStruct(1)
!UNUSED      nLev   =iSGStruct(2)
!UNUSED      lISm   =iSGStruct(3)
!UNUSED      nVert  =iSGStruct(4)
!UNUSED      lDRT   =iSGStruct(5)
!UNUSED      lDown  =iSGStruct(6)
!UNUSED      lUp    =iSGStruct(7)
!UNUSED      MidLev =iSGStruct(8)
!UNUSED      MVSta  =iSGStruct(9)
!UNUSED      MVEnd  =iSGStruct(10)
!UNUSED      lMAW   =iSGStruct(11)
!UNUSED      lLTV   =iSGStruct(12)
! Unpack structure iCIStruct:
nMidV =iCIStruct(1)
nIpWlk=iCIStruct(2)
lNOW  =iCIStruct(3)
lIOW  =iCIStruct(4)
lNCSF =iCIStruct(5)
lNOCSF=iCIStruct(6)
lIOCSF=iCIStruct(7)
nWalk =iCIStruct(8)
lICase=iCIStruct(9)
! Unpack structure iXStruct:
MxEO  =iXStruct(1)
lNOCP =iXStruct(2)
lIOCP =iXStruct(3)
nICoup=iXStruct(4)
lICoup=iXStruct(5)
nVTab =iXStruct(6)
lVTab =iXStruct(7)
lMVL  =iXStruct(8)
lMVR  =iXStruct(9)
!UNUSED      NT1MX =IXSTRUCT(10)
!UNUSED      NT2MX =IXSTRUCT(11)
!UNUSED      NT3MX =IXSTRUCT(12)
!UNUSED      NT4MX =IXSTRUCT(13)
!UNUSED      NT5MX =IXSTRUCT(14)
nNOW=2*nMidV*nSym
!UNUSED      nIOW=nNOW
Call GetMem('MVR','Free','Inte',lMVR,2*nMidV)
Call GetMem('MVL','Free','Inte',lMVL,2*nMidV)
nNOCSF=nMidV*(nSym**2)
nIOCSF=nNOCSF
Call GetMem('NOCSF','Free','Inte',lNOCSF,nNOCSF)
Call GetMem('IOCSF','Free','Inte',lIOCSF,nIOCSF)
Call GetMem('NOW','Free','Inte',lNOW,nNOW)
Call GetMem('IOW','Free','Inte',lIOW,nNOW)
nNOCP=MxEO*nMidV*nSym
nIOCP=nNOCP
Call GetMem('NOCP','Free','Inte',lNOCP,nNOCP)
Call GetMem('IOCP','Free','Inte',lIOCP,nIOCP)
Call GetMem('NCSF','Free','Inte',lNCSF,nSym)
nICase=nWalk*nIpWlk
Call GetMem('ICASE','Free','Inte',lICase,nICase)
nnICoup=(3*nICoup+1)/2
Call GetMem('ICoup','Free','Inte',lICoup,nnICoup)
Call GetMem('VTab','Free','Real',lVtab,nVtab)

end subroutine CXClose
