************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1991,1996, Markus P. Fuelscher                         *
*               1991, Per Ake Malmqvist                                *
************************************************************************
      Subroutine SORT0
************************************************************************
*                                                                      *
*     Purpose: Prepare all information needed to sort 2el integral     *
*              into canonical order. This is the interface which       *
*              matches the sorting algorithm to SEWARD.                *
*     Method:  The method of choice is a modification of the bin       *
*              sorting algorithm. The sorting processes in three       *
*              steps: 1). The integrals are  distributed into          *
*                         the bins where the buffers to save integral  *
*                         values and labels are kept separated.        *
*                         If a bin is filled the 'value buffer' is     *
*                         written to LuOrd but the 'label buffer' is   *
*                         stored on LuTmp. Prior to writing all data   *
*                         are compressed.                              *
*                         (c.f. subroutines SORT1A and SORT1B)         *
*                     2). All integrals belonging to the same slice    *
*                         are picked from the disk and sorted in       *
*                         memory. When sorting is completed the        *
*                         overwrite the old records and, if necessary, *
*                         new records are appended to LuOrd. Forward   *
*                         chaining indices are constructed an added to *
*                         the buffers. Hereafter LuTmp is not needed   *
*                         anymore.                                     *
*                         (c.f. subroutines SORT2, SORT2A and SORT2B)  *
*                     3). Finally, the records are brought into        *
*                         sequential order and the file header is      *
*                         completed                                    *
*                         (c.f. subroutines SORT3 and SORT4)           *
*     References: M. Yoshimine in 'Numerical Algorithms in Chemistry:  *
*                 Algebraic Methods', NRCC wokshop 1978                *
*                                                                      *
*     Called from: Seward_main                                         *
*                                                                      *
*     Calls to : MkSrt0,MkSrt1,MkSrt2,OpnOrd,ErrOrd,IniPkR8            *
*                                                                      *
*     Calling parameters: none                                         *
*                                                                      *
*     Global data declarations (Include files) :                       *
*     TwoDef  : definitions of the record structure                    *
*     Srt0    : common block containing information pertinent to       *
*               the calculation of 2el integral sequence numbers       *
*     Srt1    : common block containing information the number of      *
*               bins and partitioning of symmetry blocks               *
*     Srt2    : common block containing information pertinent to       *
*               the bin sorting algorithm                              *
*     PkCtl   : packing table                                          *
*     iTMax   : SEWARD's definition of highest angular momentum        *
*     Info    : SEWARD's tables of definitions                         *
*                                                                      *
*     Global data declarations (Include files) :                       *
*                                                                      *
*     Local data declarations: none                                    *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. P. Fuelscher and P.-Aa. Malmqvist                             *
*     University of Lund, Sweden, 1991                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history:                                                         *
*     - build in option to use multiple file units                     *
*       M. P. Fuelscher, University of Lund, Sweden, 1996              *
*     - build in option to use a virtual disk                          *
*       M. P. Fuelscher, University of Lund, Sweden, 1996              *
*                                                                      *
************************************************************************
*
      use Basis_Info, only: nBas
      use srt2
      use Symmetry_Info, only: iSkip
      use Integral_parameters, only: iPack
      Implicit Integer (A-Z)
*
#include "itmax.fh"
#include "info.fh"
#include "TwoDat.fh"
#include "srt0.fh"
#include "srt1.fh"
#include "stdalloc.fh"
#include "print.fh"
#include "warnings.fh"
      Logical PkMode
*----------------------------------------------------------------------*
*     pick up the print level                                          *
*----------------------------------------------------------------------*
      iRout = 80
      iPrint = nPrint(iRout)
      If ( iPrint.gt.10) Write (6,*) ' >>> Enter SORT0 <<<'
*----------------------------------------------------------------------*
*     start timer                                                      *
*----------------------------------------------------------------------*
      Call qEnter('Sort0')
*----------------------------------------------------------------------*
*     assume there is no vortual diak available                        *
*----------------------------------------------------------------------*
      RAMD = .false.
*----------------------------------------------------------------------*
*     open the two electron integral file                              *
*----------------------------------------------------------------------*
*
      LuTwo=isfreeunit(40)
      iOpt=1
      iRc=0
      Call OPNORD(iRc,iOpt,'ORDINT',LuTwo)
      If ( iRc.ne.0 ) Then
         Write (6,*) 'SORT0: Error opening ORDINT'
         Call Abend()
      End If
*
      iOpt=Toctwo(IsForm)
      Kase=iAnd(iOpt,15)
      If(Kase.eq.0) Then
         lBin=lBin_tce
      Else
         lBin=lBin_rle
      End If
*
*----------------------------------------------------------------------*
*     Initialize the common /SRT0/                                     *
*     (Get total no. of irred. representations, basis dimensions ect.) *
*----------------------------------------------------------------------*
*
      Call MKSRT0(0,nIrrep,nBas,iSkip)
*
*----------------------------------------------------------------------*
*     Initialize the common /SRT1/, i.e.,                              *
*     determine the partitioning of the 2el integrals into             *
*     submatrices (slices) where the size of the slices is             *
*     determined by the maximum number of bins (mxBin)                 *
*----------------------------------------------------------------------*
*
      Call MKSRT1()
*
*----------------------------------------------------------------------*
*     allocate the space required in phase1 of the bin sort algoritm   *
*----------------------------------------------------------------------*
*
      Call mma_allocate(lwVBin,lBin,nBin,Label='lwVBin')
      Call mma_allocate(lwIBin,lBin,nBin,Label='lwIBin')
*
      Call mma_allocate(lIndx,lBin,Label='lIndx')
      Call mma_allocate(lInts,lBin,Label='lInts')
      Call mma_allocate(ValBin,lBin,Label='ValBin')
      Call mma_allocate(IndBin,lBin,Label='IndBin')
*
*----------------------------------------------------------------------*
*     compute various offsets for each Bin and also                    *
*     initialize various pointers, counters and disk adresses          *
*----------------------------------------------------------------------*
*
      Call MKSRT2
*
*----------------------------------------------------------------------*
*     set packing mode and initialize packing table                    *
*----------------------------------------------------------------------*
*
      PkMode=.True.
      If ( iPack.ne.0 ) PkMode=.false.
      Call INIPKR8(PkAcc,PkMode)
*
*----------------------------------------------------------------------*
*     generate table of contents on 2el integral file                  *
*----------------------------------------------------------------------*
*
      Call MKORD(iDisk)
*
*----------------------------------------------------------------------*
*     Because the bin sorting algorithm heavily relies on random       *
*     access to the disks from here on the disk access control is      *
*     transfered to the sorting algorithm directly.                    *
*----------------------------------------------------------------------*
*
      iDaTw0=iDisk
      iDaTwo=iDaTw0
      mDaTwo=iDaTw0
*
*----------------------------------------------------------------------*
*     open the scratch file                                            *
*----------------------------------------------------------------------*
*
      LuTmp=isfreeunit(50)
      Call DANAME_MF(LuTmp,'TEMP01')
      iDaTmp=0
      mDaTmp=0

      Call qExit('Sort0')
      Return
      End
