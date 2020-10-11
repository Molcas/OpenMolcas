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
* Copyright (C) 1993, Markus P. Fuelscher                              *
*               1993, Per-Olof Widmark                                 *
************************************************************************
      Subroutine GetOrd(rc,Square,nSym,nBas,nSkip)
************************************************************************
*                                                                      *
*    Purpose: Read the table of content from the OrdInt file           *
*                                                                      *
*    Calling parameters:                                               *
*    nSym   : contains on return the number of irred. rep.             *
*    nBas   : contains on return the number of basis functions per     *
*             irred. rep.                                              *
*    nSkip  : contains on return the skipping flag per irred. rep.     *
*    rc     : return code                                              *
*                                                                      *
*    Global data declarations (Include files) :                        *
*    TwoDat : table of contents and auxiliary information              *
*    TowRc  : Table of return code                                     *
*    PkCtl  : packing table                                            *
*    TowID  : Table of file identifiers                                *
*                                                                      *
*    Local data declarations: none                                     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. P. Fuelscher and P.O. Widmark                                 *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Integer (A-Z)
*
#include "TwoDat.fh"
#include "TwoRc.fh"
#include "PkCtl.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "Molcas.fh"
#include "warnings.fh"
*
      Dimension nBas(0:7),nSkip(0:7)
*
      Logical DoCholesky
      Logical Square
      Character TheName*16
      Data TheName/'GetOrd'/
*----------------------------------------------------------------------*
*     Start the procedure                                              *
*----------------------------------------------------------------------*
      rc=rc0000

*----------------------------------------------------------------------*
*     Cholesky cheating: set output variables regardless of OrdInt file*
*     Note that all elements of nSkip are zero (no skipping)           *
*     Square is set false (arbitrarily)                                *
*     Skip data relating to this file (Toc etc.) entirely....          *
*----------------------------------------------------------------------*
      Call decideoncholesky(DoCholesky)
      If (DoCholesky) Then
         Call Get_iScalar('nSym',nSym)
         Call Get_iArray('nBas',nBas,nSym)
         Do iSym = 0,7
            nSkip(iSym) = 0
         End Do
         Square = .false.
         Return
      End If
*----------------------------------------------------------------------*
*     Pick up the file definitions                                     *
*----------------------------------------------------------------------*
      LuTwo=AuxTwo(isUnit)
      Open=AuxTwo(isStat)
*----------------------------------------------------------------------*
*     Check the file status                                            *
*----------------------------------------------------------------------*
      If ( Open.ne.1 ) Then
        rc=rcTC01
            Call SysAbendMsg(TheName,
     * 'The ORDINT file has not been opened',' ')
      End If
*---------------------------------------------------------------------*
*     Check the ordering parameter                                    *
*---------------------------------------------------------------------*
      If ( TocTwo(isOrd).lt.0 .or. TocTwo(isOrd).gt.1 ) then
        rc=rcTC03
            Call SysWarnMsg(TheName,
     * 'The file carries an invalid ordering parameter',' ')
       Call SysValueMsg ('TocTwo(isOrd)', TocTwo(isOrd))

      End If
      Square=TocTwo(isOrd).eq.1
*---------------------------------------------------------------------*
*     Check the number of symmetry operations                         *
*---------------------------------------------------------------------*
      nSym=TocTwo(isSym)
      If ( nSym.ne.1 .and. nSym.ne.2 .and.
     &     nSym.ne.4 .and. nSym.ne.8 ) then
        rc=rcTC04
            Call SysWarnMsg(TheName,
     * 'The file carries an invalid number '//
     * 'of irreducible representations',' ')
       Call SysValueMsg ('nSym', nSym)

      End If
      nPairs=nSym*(nSym+1)/2
      iBatch=0
      Do iSym=1,nSym
        Do jSym=1,iSym
          Do kSym=1,nSym
            Do lSym=1,kSym
              If ( ieor(iSym-1,jSym-1).eq.ieor(kSym-1,lSym-1) ) then
                ijPair=jSym+iSym*(iSym-1)/2
                klPair=lSym+kSym*(kSym-1)/2
                iSyBlk=(ijPair-1)*nPairs+klPair
                iBatch=iBatch+1
                nBatch(iSyBlk)=iBatch
              End If
            End Do
          End Do
        End Do
      End Do
*---------------------------------------------------------------------*
*     Check number of basis function                                  *
*---------------------------------------------------------------------*
      ntBas=0
      Do iSym=0,(nSym-1)
         nBas(iSym)=TocTwo(isBas+iSym)
         ntBas=ntBas+nBas(iSym)
         If ( nBas(iSym).lt.0) then
            Call SysWarnMsg(TheName,
     *                      'Invalid number of basis functions',' ')
            Call SysValueWarnMsg ('iSym', iSym)
            Call SysCondMsg('nBas(iSym).ge.0',nBas(iSym),'<',0)
         End If
*
         if( nBas(iSym).gt.mxBas ) then
            Call SysWarnMsg(TheName,
     *                      'Invalid number of basis functions',' ')
            Call SysValueWarnMsg ('iSym', iSym)
         Call SysCondMsg('nBas(iSym).lt.mxBas',nBas(iSym),'>',mxBas)
        End If
      End Do
      If ( ntBas.le.0) then
            Call SysWarnMsg(TheName,
     * 'Invalid number of basis functions',' ')
       Call SysCondMsg('ntBas.gt.0',ntBas,'<=',0)
      End If
      If ( ntBas.gt.mxOrb ) then
            Call SysWarnMsg(TheName,
     * 'Invalid number of basis functions',' ')
       Call SysCondMsg('ntBas.lt.mxOrb',ntBas,'>',mxOrb)
      End If
*---------------------------------------------------------------------*
*     Check the skip parameters                                       *
*---------------------------------------------------------------------*
      Do iSym=0,(nSym-1)
        nSkip(iSym)=TocTwo(isSkip+iSym)
        If ( nSkip(iSym).lt.0  ) then
            Call SysAbendMsg(TheName,
     * 'The table of skiping parameters is spoiled',' ')
        End If
      End Do
*---------------------------------------------------------------------*
*     Check table of disk adresses                                    *
*---------------------------------------------------------------------*
      mxDAdr=TocTwo(isMxDa)
      If ( mxDAdr.lt.0 ) then
            Call SysWarnMsg(TheName,
     * 'The file carries an invalid disk address',' ')
       Call SysCondMsg('mxDAdr.ge.0',mxDAdr,'<',0)

      End If
      Do iTab=0,175
        If ( TocTwo(isDAdr+iTab).lt.0 .or.
     &       TocTwo(isDAdr+iTab).gt. mxDAdr ) then
            Call SysWarnMsg(TheName,
     * 'The table of disk adresses is spoiled',' ')
       Call SysValueWarnMsg ('iTab', iTab)
       Call SysCondMsg('TocTwo(isDAdr+iTab).lt.mxDAdr',
     *   TocTwo(isDAdr+iTab),'>',mxDAdr)

        End If
      End Do
*---------------------------------------------------------------------*
*     Generate and check packing table                                *
*---------------------------------------------------------------------*
      Call Int2Real(TocTwo(isPkTh),PkThrs)
      Call Int2Real(TocTwo(isPkCt),PkCutof)
      If ( PkThrs.lt.0 ) then
            Call SysAbendMsg(TheName,
     * 'The accuracy threshold for unpacking is spoiled',' ')
      End If
      Call Int2Real(TocTwo(isPkSc),PkScal)
      If ( PkScal.ne.1.0D0 .and. PkScal.ne.2.0D0 .and.
     &     PkScal.ne.4.0D0 .and. PkScal.ne.8.0D0 ) then
            Call SysAbendMsg(TheName,
     * 'The scaling constant for unpacking is spoiled',' ')
      End If
      iPack=TocTwo(isPkPa)
      If ( iPack.lt.0 .or. iPack.gt.1 ) then
            Call SysWarnMsg(TheName,
     * 'The packing flag is spoiled',' ')
       Call SysValueMsg ('iPack', iPack)
      End If
      If( TocTwo(isPkPa).eq.0 ) Pack=.true.
      If( TocTwo(isPkPa).eq.1 ) Pack=.false.
      iAssm=TocTwo(isPkAs)
      If ( iAssm.lt.0 .or. iAssm.gt.1 ) then
            Call SysWarnMsg(TheName,
     * 'The assembler flag is spoiled',' ')
       Call SysValueMsg ('iAssm', iAssm)
      End If
      If( TocTwo(isPkAs).eq.0 ) Assm=.true.
      If( TocTwo(isPkAs).eq.1 ) Assm=.false.
      Do iExp=0,4095
        PkTab(iExp)=TocTwo(isPkTb+iExp)
        If ( PkTab(iExp).lt.1 ) then
            Call SysWarnMsg(TheName,
     * 'The packing table is spoiled',' ')
       Call SysValueWarnMsg ('iExp', iExp)
       Call SysCondMsg('PkTab(iExp).gt.0 ',
     *   PkTab(iExp),'<',1)
        End If
      End Do
*---------------------------------------------------------------------*
*     exit                                                            *
*---------------------------------------------------------------------*
      Return
      End
