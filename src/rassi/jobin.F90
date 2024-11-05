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
Module JobIn
# include "rasdim.fh"
! The paramaters defined in Molcas should be private
Private MaxBfn,MaxBfn_Aux, MxAO, mxAtom, mxroot, mxNemoAtom, Mxdbsc, lCache, mxact, mxina, mxbas, mxOrb, &
        mxSym, mxGAS, LENIN, LENIN1, LENIN2, LENIN3, LENIN4, LENIN5, LENIN6, LENIN8
Private mxRef,mxIter,mxCiIt,mxSxIt,mxTit

! TEMPORARY DATA FROM JOBIPHS
REAL*8 POTNU1
Integer NACTE1,MPLET1,NSYM1,NFRO1(8),NISH1(8),NASH1(8),NDEL1(8),NBAS1(8),NRS11(8),NRS21(8), &
        NRS31(8),LROT1,NROOT1,IROOT1(mxRoot),NHOL11,NELE31
CHARACTER(LEN=LENIN8) NAME(mxOrb)
CHARACTER(LEN=2) HEAD1(72)
CHARACTER(LEN=4) TITLE1(18,mxTit)
End Module JobIn
