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
* Copyright (C) 1996, Markus P. Fuelscher                              *
************************************************************************
      Subroutine SORT1C(nUt,vInt,nSqNum,nSyBlk)
************************************************************************
*                                                                      *
*     Purpose: Provided SEWARD was allowed to use a virtual disk       *
*              scatter the 2el integrals into the approriate place.    *
*                                                                      *
*     Called from: SORT1A                                              *
*                                                                      *
*     Calling Parameters:                                              *
*     nUt    : number of 2el integrals in the buffers vInt,nSqNum      *
*              and nSyBlk                                              *
*     vInt   : Buffer of 2el integral values                           *
*     nSqNum : sequence number of the integral relative to             *
*              the first adress of the symmetry block                  *
*     nSyBlk : symmetry block number of an integral                    *
*                                                                      *
*     Global data declarations (Include files) :                       *
*     TwoDat : definitions of sorting flags and address tables         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. P. Fuelscher                                                  *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
*
      Implicit Real*8 (A-H,O-Z)
*
#include "Molcas.fh"
#include "TwoDat.fh"
*
      Real*8 vInt(nUt),nSqNum(nUt),nSyBlk(nUt)
*
*----------------------------------------------------------------------*
*     Turn timing ON                                                   *
*----------------------------------------------------------------------*
*
      Call QEnter('Sort1C')
*
*----------------------------------------------------------------------*
*     scatter the 2el integrals to the appropriate position            *
*     of the virtual disk                                              *
*----------------------------------------------------------------------*
*
*     Write (*,'(2X,5(I4,I8,F12.8))')
*    &      (nSyBlk(iUt),
*    &       RAMD_adr(nBatch(nSyBlk(iUt)))+nSqNum(iUt)-1,
*    &       vInt(iUt),
*    &       iUt=1,nUt)
      Do iUt=1,nUt
        iSyBlk=INT(nSyBlk(iUt))
        iBatch=nBatch(iSyBlk)
        iOff=RAMD_adr(iBatch)+INT(nSqNum(iUt))
        RAMD_ints(iOff) = vInt(iUt)
      End Do
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Call QExit('Sort1C')
      Return
      End
