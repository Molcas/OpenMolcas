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
      Subroutine MkOrd(iDisk)
*
************************************************************************
*                                                                      *
*    Purpose: Create the table of content of the OrdInt file           *
*                                                                      *
*    Called from: Sort0 ans Sort3                                      *
*                                                                      *
*    Calls to : DaFile                                                 *
*                                                                      *
*    Calling parameters:                                               *
*    iDisk  : First free entry on disk after table of contents         *
*                                                                      *
*    Global data declarations (Include files) :                        *
*    TwoDat  : table of contents and auxiliary information             *
*              on the ordered 2el file                                 *
*    TowID   : Table of file identifiers                               *
*    TwoDef  : definitions of the record structure                     *
*    Srt0    : common block containing information pertinent to        *
*              the calculation of 2el integral sequence numbers        *
*    Srt1    : common block containing information the number of       *
*              bins and partitioning of symmetry blocks                *
*    Srt2    : common block containing information pertinent to        *
*              the bin sorting algorithm                               *
*    PkCtl   : packing table                                           *
*                                                                      *
*    Local data declarations: none                                     *
*                                                                      *
**** M. Fuelscher and P.-Aa. Malmqvist, Univ. of Lund, Sweden, 1991 ****
*
      Implicit Integer (A-Z)
*

#include "itmax.fh"
#include "info.fh"
#include "FileIDs.fh"
#include "TwoDat.fh"
#include "TwoDef.fh"
#include "srt0.fh"
#include "srt1.fh"
#include "srt2.fh"
#include "PkCtl.fh"
*
C     Call qEnter('MkOrd')
*---------------------------------------------------------------------*
*     Initialize table of content                                     *
*---------------------------------------------------------------------*
*
      Do iToc=1,lTocTwo
         TocTwo(iToc)=iNoNum
      End Do
*
*---------------------------------------------------------------------*
*     Write file identifier                                           *
*---------------------------------------------------------------------*
*
      TocTwo(isId)=IDtwo
      TocTwo(isVer)=VNtwo
      TocTwo(isForm)=iOrdFm
*
*---------------------------------------------------------------------*
*     Write ordring mode                                              *
*---------------------------------------------------------------------*
*
      TocTwo(isOrd)=iSquar
*
*---------------------------------------------------------------------*
*     Write symmetry and basis set information                        *
*---------------------------------------------------------------------*
*
      TocTwo(isSym)=nSyOp
      Do iSymi=0,nSyOp-1
        TocTwo(isSkip+iSymi)=nSkip(iSymi+1)
      End Do
      nPairs=nSyOp*(nSyOp+1)/2
      iBatch=0
      Do iSymi=1,nSyOp
        Do jSymj=1,iSymi
          Do kSymk=1,nSyOp
            Do lSyml=1,kSymk
              If ( ieor(iSymi-1,jSymj-1).eq.ieor(kSymk-1,lSyml-1) ) then
                ijPair=jSymj+iSymi*(iSymi-1)/2
                klPair=lSyml+kSymk*(kSymk-1)/2
                iSyBlk=(ijPair-1)*nPairs+klPair
                iBatch=iBatch+1
                nBatch(iSyBlk)=iBatch
              End If
            End Do
          End Do
        End Do
      End Do
*
*---------------------------------------------------------------------*
*     Write the skip parameters                                       *
*---------------------------------------------------------------------*
*
      TocTwo(isSym)=nSyOP
      Do iSymi=0,nSyOp-1
        TocTwo(isBas+iSymi)=nBs(iSymi+1)
      End Do
*
*---------------------------------------------------------------------*
*     Write disk access table                                         *
*---------------------------------------------------------------------*
*
      Do iBatch=0,175
         TocTwo(isDAdr+iBatch)=0
      End Do
*
      iBin=1
      Do iSymi=1,nSyOp
         iFlit=nSkip(iSymi)
         Do jSymj=1,iSymi
            jFlit=nSkip(jSymj)
            iSymj=1+ieor(iSymi-1,jSymj-1)
            iSyblj=jSymj+iSymi*(iSymi-1)/2
            kSymMx=iSymi
            If ( Square ) kSymMx=nSyOp
            Do kSymk=1,kSymMx
               kFlit=nSkip(kSymk)
               lSymMx=jSymj
               If( kSymk.ne.iSymi .or. Square ) lSymMx=kSymk
               Do lSyml=1,lSymMx
                  lFlit=nSkip(lSyml)
                  kSyml=1+ieor(kSymk-1,lSyml-1)
*
                  kSybll=lSyml+kSymk*(kSymk-1)/2
                  If ( ieor(iSymj-1,kSyml-1).eq.0 ) then
                     If ( (iFlit+jFlit+kFlit+lFlit).eq.0 ) then
                        iSyBlk=kSybll+mxSyP*(iSyblj-1)
                        iBatch=nBatch(iSyBlk)
                        TocTwo(isDAdr+iBatch-1)=iDVBin(2,iBin)
                        iBin=iBin+nSln(iSyBlk)
                     End If
                  End If
               End Do
            End Do
         End Do
      End Do
*
      TocTwo(isMxDa)=mDaTwo
*
*---------------------------------------------------------------------*
*     Write packing information                                       *
*---------------------------------------------------------------------*
*
      Call Real2Int(PkThrs,TocTwo(isPkTh))
      Call Real2Int(PkCutof,TocTwo(isPkCt))
      Call Real2Int(PkScal,TocTwo(isPkSc))
      iPack=1
      If ( Pack ) iPack=0
      TocTwo(isPkPa)=iPack
      iAssm=1
      If ( Assm ) iAssm=0
      TocTwo(isPkAs)=iAssm
      Do iExp=0,4095
        TocTwo(isPkTb+iExp)=PkTab(iExp)
      End Do
*
*---------------------------------------------------------------------*
*     Transfer table of content to disk                               *
*---------------------------------------------------------------------*
*
      iDisk=0
      iOpt=1
      LuTwo=AuxTwo(isUnit)
      Call iDAFILE(LuTwo,iOpt,TocTwo,lTocTwo,iDisk)
      AuxTwo(isDaDa)=iDisk
*
C     Call qExit('MkOrd')
      Return
      End
