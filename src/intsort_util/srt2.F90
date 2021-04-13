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
!----------------------------------------------------------------------*
!                                                                      *
!     Information pertinent to the bin sort algorithm                  *
!                                                                      *
!     *--------------------------------------------------------*       *
!     *                                                        *       *
!     *  Note:                                                 *       *
!     *  These definitions depend on the record structure of   *       *
!     *       2el integral file which are given in the         *       *
!     *         common /TWODEF/ of the RDORD utility           *       *
!     *                                                        *       *
!     *--------------------------------------------------------*       *
!                                                                      *
!     Parameter definitions:                                           *
!     mxBin  : maximum number of bins allowed. The number should       *
!              not be smaller than:                                    *
!              If nSyOp=1 mxBin=  1 and If Square=.true. mxBin=  1     *
!                 nSyOp=2 mxBin=  4                      mxBin=  5     *
!                 nSyOp=4 mxBin= 19                      mxBin= 28     *
!                 nSyOp=8 mxBin=106                      mxBin=176     *
!     lBin   : buffer length of a bin                                  *
!              Note: if lBin is not a multiple of 4*lDaRec             *
!              the buffers used for intermediate storing of            *
!              labels on LuTmp are used very inefficiently             *
!     mxSrtA : maximum sorting area to be used                         *
!                                                                      *
!     Entries to common SRT2:                                          *
!     lwIBin : array used to store index Bins                          *
!     lwVBin : array used to store value Bins                          *
!     nRec   : number of records per slice                             *
!     iDIBin : disk addresses of index bins                            *
!     iDVBin : disk addresses of value bins                            *
!     n_Int  : number of integrals in a bin                            *
!     mInt   : number of integrals and bytes in a slice                *
!     nOff1  : memory allocation offset for index bins                 *
!     nOff2  : memory allocation offset for value bins                 *
!     LuTwo  : logical unit number of ordered 2el file                 *
!     LuTmp  : logical unit number of temporary file                   *
!     iDaTw0 : first disk address after header of ordered 2el file     *
!     iDaTwo : current disk position of LuTwo                          *
!     iDaTmp : current disk position of LuTmp                          *
!     mDaTwo : highest accessed disk position of LuTwo                 *
!     mDaTwo : highest accessed disk position of LuTmp                 *
!     MxOrd  : total number of records on LuTwo                        *
!                                                                      *
!----------------------------------------------------------------------*

module Srt2

use Definitions, only: wp, iwp

implicit none
private

#include "TwoDef.fh"
integer(kind=iwp), parameter :: mxBin = 2048

integer(kind=iwp) :: iDIBin(3,mxBin), iDVBin(4,mxBin), nRec(mxBin), n_Int(mxBin), mInt(3,mxBin)
integer(kind=iwp) :: LuTwo, LuTmp, iDaTw0, iDaTwo, iDaTmp, mDaTwo, mDaTmp, MxOrd, lBin

integer(kind=iwp), allocatable :: IndBin(:), lIndx(:), lInts(:), lwIBin(:,:)
real(kind=wp), allocatable :: lwVBin(:,:), ValBin(:)

public :: &
iDaTmp, &
iDaTw0, &
iDaTwo, &
iDIBin, &
iDVBin, &
IndBin, &
lBin, &
lDaRec, &
lIndx, &
lInts, &
lStRec, &
lTop, &
LuTmp, &
LuTwo, &
lwIBin, &
lwVBin, &
mDaTmp, &
mDaTwo, &
mInt, &
mxBin, & ! *
MxOrd, &
n_Int, &
nRec, &
nSect, &
ValBin

end module Srt2
