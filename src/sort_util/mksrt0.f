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
* Copyright (C) 1991, Markus P. Fuelscher                              *
*               1991, Per Ake Malmqvist                                *
************************************************************************
      Subroutine MkSrt0(iSquar,nIrrep,nBas,iSkip)
************************************************************************
*                                                                      *
*     Purpose: Set up all information needed to compute 2el integral   *
*              adresses in the subroutine PLF2 and INDSFT2.            *
*                                                                      *
*     Called from: Sort0                                               *
*                                                                      *
*     Calls to : none                                                  *
*                                                                      *
*     Calling parameters:                                              *
*     iSquar  : ordering flag                                          *
*     nIrrep  : number of irreducible representations                  *
*     nBas    : number of basis functions per irred. rep.              *
*     nSkip   : flag to exlude symmetry combinations                   *
*                                                                      *
*     Global data declarations (Include files) :                       *
*     Srt0    : common block containing information pertinent to       *
*               the calculation of 2el integral sequence numbers       *
*                                                                      *
*     local data declarations: none                                    *
*                                                                      *
**** M. Fuelscher and P.-Aa. Malmqvist, Univ. of Lund, Sweden, 1991 ****
*
      Implicit Integer (A-Z)
*
#include "srt0.fh"
#include "print.fh"
#include "SysCtl.fh"
      Dimension nBas(*),iSkip(*)
*
      iRout = 80
      iPrint = nPrint(iRout)
      if ( iPrint.gt.10) Write(6,*) ' >>> Enter MKSRT0 <<<'
      Call qEnter('MkSrt0')
*----------------------------------------------------------------------*
*     Gather information on desired ordering scheme                    *
*----------------------------------------------------------------------*
*
      Square=.true.
      If( iSquar.eq.0 )  Square=.false.
*
*----------------------------------------------------------------------*
*     Gather data on the number of symmetry operations                 *
*----------------------------------------------------------------------*
*
      nSyOp=nIrrep
      mxSyP=nSyOp*(nSyOp+1)/2
*
*----------------------------------------------------------------------*
*     Gather data on the number of basis functions                     *
*----------------------------------------------------------------------*
*
      Do iSymi=1,nSyOp
        nBs(iSymi)=nBas(iSymi)
      End Do
*
*----------------------------------------------------------------------*
*     Put flags to exclude symmetry combinations into common block     *
*----------------------------------------------------------------------*
*
      Do 20 iSymi=1,nSyOp
        nSkip(iSymi)=iSkip(iSymi)
20    Continue
*
*----------------------------------------------------------------------*
*     Precompute the dimension of the symmetry blocks                  *
*     and symmetry block numbers for pairs of symmtry indices          *
*----------------------------------------------------------------------*
*
      Do 30 iSymi=1,nSyOp
        iBsi=nBs(iSymi)
        DimSyB(iSymi,iSymi)=iBsi*(iBsi+1)/2
        TriSyB(iSymi,iSymi)=iSymi*(iSymi+1)/2
        Do 40 jSymj=1,iSymi-1
          jBsj=nBs(jSymj)
          DimSyB(iSymi,jSymj)=iBsi*jBsj
          DimSyB(jSymj,iSymi)=jBsj*iBsi
          TriSyB(iSymi,jSymj)=jSymj+iSymi*(iSymi-1)/2
          TriSyB(jSymj,iSymi)=jSymj+iSymi*(iSymi-1)/2
40      Continue
30    Continue
*
      Call qExit('MkSrt0')
      Return
      End
