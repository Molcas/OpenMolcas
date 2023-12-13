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
!     Information pertinent to the generation of integral adresses     *
!                                                                      *
!     Square : flag which indicates the desired ordering scheme        *
!     nSyOp  : number of irreducible representations                   *
!     mxSyP  : max number of pairs of symmetry indices                 *
!     nBs    : number of atomic orbitals per irred. rep.               *
!     nSkip  : flag to exclude symmetry combinations                   *
!     DimSyB : dimension of a symmetry block                           *
!     TriSyB : symmetry block number of pairs of symmetry indices      *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     Information pertinent to translation from symmetry block         *
!     to bin numbers.                                                  *
!                                                                      *
!     Parameter definitions:                                           *
!     mSyBlk : highest possible symmetry block number                  *
!                                                                      *
!     Entries to common SRT1:                                          *
!     nBin   : total number of Bins                                    *
!     nSln   : number of submatrices (slices) per symmetry block       *
!     lSll   : length of each slice                                    *
!     IstBin : offset for bin number per symmetry block                *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     Information pertinent to the bin sort algorithm                  *
!                                                                      *
!     *--------------------------------------------------------*       *
!     *                                                        *       *
!     *  Note:                                                 *       *
!     *  These definitions depend on the record structure of   *       *
!     *       2el integral file which are given in the         *       *
!     *         TwoDat module                                  *       *
!     *                                                        *       *
!     *--------------------------------------------------------*       *
!                                                                      *
!     Parameter definitions:                                           *
!     lBin   : buffer length of a bin                                  *
!              Note: if lBin is not a multiple of 4*lDaRec             *
!              the buffers used for intermediate storing of            *
!              labels on LuTmp are used very inefficiently             *
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

module Sort_data

use Definitions, only: wp, iwp

implicit none
private

!integer(kind=iwp), parameter :: mSyBlk = 36*36 ! where 36 = 8*(8+1)/2
integer(kind=iwp) :: mSyBlk

integer(kind=iwp) :: DimSyB(8,8), iDaTmp, iDaTw0, iDaTwo,                 lBin,               LuTmp, LuTwo, mDaTmp, mDaTwo, MxOrd, &
                     mxSyP, nBin, nBs(8), nSkip(8),               nSyOp, TriSyB(8,8)
logical(kind=iwp) :: Square

integer(kind=iwp), allocatable :: iStBin(:), lSll(:), nSln(:)

integer(kind=iwp), allocatable :: iDIBin(:,:), iDVBin(:,:), IndBin(:), lIndx(:), lInts(:), lwIBin(:,:), mInt(:,:), n_Int(:), nRec(:)
real(kind=wp), allocatable :: lwVBin(:,:), ValBin(:)

public :: DimSyB, iDaTmp, iDaTw0, iDaTwo, iDIBin, iDVBin, IndBin, iStBin, lBin, lIndx, lInts, lSll, LuTmp, LuTwo, lwIBin, lwVBin, &
          mDaTmp, mDaTwo, mInt, mSyBlk, MxOrd, mxSyP, n_Int, nBin, nBs, nRec, nSkip, nSln, nSyOp, Square, TriSyB, ValBin

end module Sort_data
