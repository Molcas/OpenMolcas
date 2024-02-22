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
SUBROUTINE CXSVAL(CIS,IXS,NMIDV,NIPWLK,&
                  NWALK,MXEO,LNOCP,LIOCP,NICOUP,LICOUP,NVTAB,&
                  LVTAB,LMVL,LMVR,NT1MX,NT2MX,NT3MX,NT4MX,NT5MX)
!
! Purpose: Dereference the CI structure and the excitation
! structure arrays and return values and pointers.
!
use Struct, only: nXSize, CIStruct
IMPLICIT REAL*8 (A-H,O-Z)
Type (CIStruct) CIS
Dimension IXS(nXSize)

! CI structure, sizes, addresses...
nMidV  =CIS%nMidV
nIpWlk =CIS%nIpWlk
nWalk  =CIS%nWalk
! Excitation operators, coupling coefficients,...
MxEO =IXS(1)
lNOCP =IXS(2)
lIOCP =IXS(3)
nICoup =IXS(4)
lICoup =IXS(5)
nVTab =IXS(6)
lVTab =IXS(7)
lMVL =IXS(8)
lMVR =IXS(9)
NT1MX =IXS(10)
NT2MX =IXS(11)
NT3MX =IXS(12)
NT4MX =IXS(13)
NT5MX =IXS(14)

END SUBROUTINE CXSVAL
