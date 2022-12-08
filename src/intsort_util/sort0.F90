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
! Copyright (C) 1991,1996, Markus P. Fuelscher                         *
!               1991, Per Ake Malmqvist                                *
!***********************************************************************

subroutine SORT0()
!***********************************************************************
!                                                                      *
!     Purpose: Prepare all information needed to sort 2el integral     *
!              into canonical order. This is the interface which       *
!              matches the sorting algorithm to SEWARD.                *
!     Method:  The method of choice is a modification of the bin       *
!              sorting algorithm. The sorting proceeds in three        *
!              steps: 1). The integrals are  distributed into          *
!                         the bins where the buffers to save integral  *
!                         values and labels are kept separated.        *
!                         If a bin is filled the 'value buffer' is     *
!                         written to LuOrd but the 'label buffer' is   *
!                         stored on LuTmp. Prior to writing all data   *
!                         are compressed.                              *
!                         (cf. subroutines SORT1A and SORT1B)          *
!                     2). All integrals belonging to the same slice    *
!                         are picked from the disk and sorted in       *
!                         memory. When sorting is completed they       *
!                         overwrite the old records and, if necessary, *
!                         new records are appended to LuOrd. Forward   *
!                         chaining indices are constructed and added to*
!                         the buffers. Hereafter LuTmp is not needed   *
!                         anymore.                                     *
!                         (cf. subroutines SORT2, SORT2A and SORT2B)   *
!                     3). Finally, the records are brought into        *
!                         sequential order and the file header is      *
!                         completed                                    *
!                         (cf. subroutines SORT3 and SORT4)            *
!     References: M. Yoshimine in 'Numerical Algorithms in Chemistry:  *
!                 Algebraic Methods', NRCC wokshop 1978                *
!                                                                      *
!     Called from: Seward_main                                         *
!                                                                      *
!     Calls to : MkSrt0,MkSrt1,MkSrt2,OpnOrd,ErrOrd,IniPkR8            *
!                                                                      *
!     Calling parameters: none                                         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M. P. Fuelscher and P.-Aa. Malmqvist                             *
!     University of Lund, Sweden, 1991                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history:                                                         *
!     - build in option to use multiple file units                     *
!       M. P. Fuelscher, University of Lund, Sweden, 1996              *
!     - build in option to use a virtual disk                          *
!       M. P. Fuelscher, University of Lund, Sweden, 1996              *
!                                                                      *
!***********************************************************************

use Basis_Info, only: nBas
use TwoDat, only: isForm, lDaRec, RAMD, sNew, TocTwo
use sort_data, only: iDaTmp, iDaTw0, iDaTwo, IndBin, lBin, lIndx, lInts, LuTmp, LuTwo, lwIBin, lwVBin, mDaTmp, mDaTwo, nBin, ValBin
use Symmetry_Info, only: nIrrep, iSkip
use Gateway_global, only: iPack
use Gateway_Info, only: PkAcc
use stdalloc, only: mma_allocate
use Definitions, only: iwp, u6

implicit none
#include "print.fh"
integer(kind=iwp) :: iDisk, iOpt, iPrint, iRc, iRout, Kase
logical(kind=iwp) :: PkMode
integer(kind=iwp), external :: isfreeunit

!----------------------------------------------------------------------*
!     pick up the print level                                          *
!----------------------------------------------------------------------*
iRout = 80
iPrint = nPrint(iRout)
if (iPrint > 10) write(u6,*) ' >>> Enter SORT0 <<<'
!----------------------------------------------------------------------*
!     assume there is no virtual disk available                        *
!----------------------------------------------------------------------*
RAMD%act = .false.
!----------------------------------------------------------------------*
!     open the two electron integral file                              *
!----------------------------------------------------------------------*

LuTwo = isfreeunit(40)
iOpt = ibset(0,sNew)
iRc = 0
call OPNORD(iRc,iOpt,'ORDINT',LuTwo)
if (iRc /= 0) then
  write(u6,*) 'SORT0: Error opening ORDINT'
  call Abend()
end if

iOpt = TocTwo(IsForm)
Kase = iand(iOpt,15)
if (Kase == 0) then
  lBin = 4*lDaRec  ! tce
else
  lBin = 32*lDaRec ! rle
end if

!----------------------------------------------------------------------*
!     Initialize the common /SRT0/                                     *
!     (Get total no. of irred. representations, basis dimensions ect.) *
!----------------------------------------------------------------------*

call MKSRT0(0,nIrrep,nBas,iSkip)

!----------------------------------------------------------------------*
!     Initialize the common /SRT1/, i.e.,                              *
!     determine the partitioning of the 2el integrals into             *
!     submatrices (slices) where the size of the slices is             *
!     determined by the maximum number of bins (mxBin)                 *
!----------------------------------------------------------------------*

call MKSRT1()

!----------------------------------------------------------------------*
!     allocate the space required in phase1 of the bin sort algoritm   *
!----------------------------------------------------------------------*

call mma_allocate(lwVBin,lBin,nBin,Label='lwVBin')
call mma_allocate(lwIBin,lBin,nBin,Label='lwIBin')

call mma_allocate(lIndx,lBin,Label='lIndx')
call mma_allocate(lInts,lBin,Label='lInts')
call mma_allocate(ValBin,lBin,Label='ValBin')
call mma_allocate(IndBin,lBin,Label='IndBin')

!----------------------------------------------------------------------*
!     compute various offsets for each Bin and also                    *
!     initialize various pointers, counters and disk adresses          *
!----------------------------------------------------------------------*

call MKSRT2()

!----------------------------------------------------------------------*
!     set packing mode and initialize packing table                    *
!----------------------------------------------------------------------*

PkMode = .true.
if (iPack /= 0) PkMode = .false.
call INIPKR8(PkAcc,PkMode)

!----------------------------------------------------------------------*
!     generate table of contents on 2el integral file                  *
!----------------------------------------------------------------------*

call MKORD(iDisk)

!----------------------------------------------------------------------*
!     Because the bin sorting algorithm heavily relies on random       *
!     access to the disks from here on the disk access control is      *
!     transfered to the sorting algorithm directly.                    *
!----------------------------------------------------------------------*

iDaTw0 = iDisk
iDaTwo = iDaTw0
mDaTwo = iDaTw0

!----------------------------------------------------------------------*
!     open the scratch file                                            *
!----------------------------------------------------------------------*

LuTmp = isfreeunit(50)
call DANAME_MF(LuTmp,'TEMP01')
iDaTmp = 0
mDaTmp = 0

return

end subroutine SORT0
