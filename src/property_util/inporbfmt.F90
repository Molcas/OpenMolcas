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
! INPORB formatting parameters
! to add a new version, increase version number and
! add appropriate new Magic string and format descriptors

module InpOrbFmt

use Definitions, only: iwp

implicit none
private

integer(kind=iwp), parameter :: mxVer = 5, &
                                iVer10 = 1, &
                                iVer11 = 2, &
                                iVer20 = 3, &
                                iVer21 = 4, &
                                iVer22 = 5, &
                                nDivOrb(mxVer) = [4,4,5,5,5], &
                                nDivOcc(mxVer) = [4,4,10,4,5], &
                                nDivOccHR(mxVer) = [0,0,0,10,10], &
                                nDivEne(mxVer) = [4,4,10,10,10], &
                                nDivInd(mxVer) = [4,10,10,10,10], &
                                nSkpInd(mxVer) = [0,1,1,1,1]
character(len=*), parameter :: Magic(mxVer) = ['#INPORB 1.0', &
                                               '#INPORB 1.1', &
                                               '#INPORB 2.0', &
                                               '#INPORB 2.1', &
                                               '#INPORB 2.2' &
                                              ], &
                               FmtOrb(mxVer) = [character(len=40) :: &
                                                '(4E18.12)', &
                                                '(4E18.12)', &
                                                '(5(1X,ES21.14))', &
                                                '(5(1X,ES21.14))', &
                                                '(5(1X,ES21.14))' &
                                               ], &
                               FmtOcc(mxVer) = [character(len=40) :: &
                                                '(4E18.12)', &
                                                '(4E18.12)', &
                                                '(10(1X,F7.4))', &
                                                '(4E18.12)', &
                                                '(5(1X,ES21.14))' &
                                               ], &
                               FmtOccHR(mxVer) = [character(len=40) :: &
                                                  '', &
                                                  '', &
                                                  '', &
                                                  '(10(1X,F7.4))', &
                                                  '(10(1X,F7.4))' &
                                                 ], &
                               FmtEne(mxVer) = [character(len=40) :: &
                                                '(4E18.12)', &
                                                '(4E18.12)', &
                                                '(10(1X,ES11.4))', &
                                                '(10(1X,ES11.4))', &
                                                '(10(1X,ES11.4))' &
                                               ], &
                               FmtInd(mxVer) = [character(len=40) :: &
                                                '(A4)', &
                                                '(2X,A10)', &
                                                '(2X,A10)', &
                                                '(2X,A10)', &
                                                '(2X,A10)' &
                                               ]

public :: FmtEne, FmtInd, FmtOcc, FmtOccHR, FmtOrb, Magic, iVer10, iVer11, iVer20, iVer21, iVer22, mxVer, nDivEne, nDivInd, &
          nDivOcc, nDivOccHR, nDivOrb, nSkpInd

end module InpOrbFmt
