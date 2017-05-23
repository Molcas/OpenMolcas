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
* Copyright (C) 1993, Per-Olof Widmark                                 *
*               1993, Markus P. Fuelscher                              *
*               1995, Anders Bernhardsson                              *
************************************************************************
      Subroutine WrMCK(rc,Option,InLab,iComp,Data,iSymLab)
************************************************************************
*                                                                      *
*     Purpose: write data to one-electron integral file                *
*                                                                      *
*     Calling parameters:                                              *
*     rc      : Return code                                            *
*     Option  : Switch to set options                                  *
*     InLab   : Identifier for the data to write                       *
*     Comp    : Composite identifier to select components              *
*     Data    : contains on input the data to store on disk            *
*     SymLab  : symmetry label of the provided data                    *
*                                                                      *
*     Global data declarations (Include file):                         *
*     Parm    : Table of paramaters                                    *
*     rcParm  : Table of return codes                                  *
*     Switch  : Table of options                                       *
*     Common  : Common block containing TocOne                         *
*     Data    : Data definitions                                       *
*                                                                      *
*     Local data declarations:                                         *
*     Label   : character*8, used to covert incoming names             *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark and M.P. Fuelscher                                  *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*     Modified by:                                                     *
*     Anders Bernhardsson                                              *
*     University of Lund, Sweden, 1995                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Integer (A-Z)
*
#include "MckRc.fh"
#include "MckFlags.fh"
#include "SysDef.fh"
#include "MckDat.fh"
      Character*(*) InLab
      Dimension Data(*)
C     Real*8 DDot_, Check
C     external ddot_
*
      Logical Debug
C     Character*8 TmpLab,Label, Label_Add*11
      Character*8 TmpLab,Label
      Dimension LabTmp(2)
*     Equivalence (TmpLab,LabTmp)
      Character*16 TheName
      Data TheName/'WrMck'/
      Data Debug /.False./
*----------------------------------------------------------------------*
*     Start procedure:                                                 *
*     Define inline function (symmetry multiplication)                 *
*----------------------------------------------------------------------*
      MulTab(i,j)=iEor(i-1,j-1)+1
*----------------------------------------------------------------------*
*     Pick up the file definitions                                     *
*----------------------------------------------------------------------*
*     Call qEnter('WrMCK')
      SymLab=iSymLab
      icpi=itob
      Comp=iComp
      Length=rc
      rc    = rc0000
      LuMCK = AuxMCK(pLu  )
      Open  = AuxMCK(pOpen)
*----------------------------------------------------------------------*
*     Check the file status                                            *
*----------------------------------------------------------------------*
      If(Open.ne.1) Then
          Call SysFileMsg(TheName,'MSG: open',LuMck,' ')
      End If
*----------------------------------------------------------------------*
*     Truncate the label to 8 characters and convert it to upper case  *
*----------------------------------------------------------------------*
*     Call StdFmt(InLab,Label)
      Label=InLab
      Call UpCase(Label)
      TmpLab=Label
      Call ByteCopy(TmpLab,LabTmp,8)
*----------------------------------------------------------------------*
*     Print debugging information                                      *
*----------------------------------------------------------------------*
      If(iAnd(option,1024).ne.0) debug=.true.
      If(Debug) Then
         Write(6,*) '<<< Entering WrMck >>>'
         Write(6,'(a,z8)') ' rc on entry:     ',rc
         Write(6,'(a,a)')  ' Label on entry:  ',Label
         Write(6,'(a,z8)') ' Comp on entry:   ',Comp
         Write(6,'(a,z8)') ' SymLab on entry: ',SymLab
         Write(6,'(a,z8)') ' Option on entry: ',Option
         Write(6,*) ' Contents of the Toc'
         Write(6,*) ' ==================='
         Write (6,'(i6,z8)') 'pFID,TocOne(pFID)=',pFID,TocOne(pFID)
         Write (6,'(i6,z8)') pVersN,TocOne(pVersN)
         Write (6,'(i6,z8)') pTitle,TocOne(pTitle)
         Write (6,'(i6,z8)') pOp,TocOne(pOp)
         Write (6,'(i6,z8)') pSym,TocOne(pSym)
         Write (6,'(i6,z8)') pSymOp,TocOne(pSymOp)
         Write (6,'(i6,z8)') pBas,TocOne(pBas)
         Write (6,'(i6,z8)') pNext,TocOne(pNext)
         Write (6,'(i6,z8)') pEnd,TocOne(pEnd)
      End If
*----------------------------------------------------------------------*
*     Store data in TocOne                                             *
*----------------------------------------------------------------------*
*
      If(Label.eq.'TITLE') Then
*
         TocOne(pTitle)=NotNaN
         Call iCopy(nTitle,Data(1),1,TocOne(pTitle+1),1)
*----------------------------------------------------------------------*
*
      Else If(Label.eq.'NSYM') Then
*
         If(Data(1).gt.MxSym .or. Data(1).lt.1) Then
         Call SysWarnMsg(TheName,'Label=',Label)
         Call SysValueMsg('Data(1)=',Data(1))
         End If
         TocOne(pSym)=Data(1)
*----------------------------------------------------------------------*
*
      Else If(Label.eq.'NBAS') Then
*
         If(TocOne(pSym).eq.NaN) Then
         Call SysAbendMsg(TheName,'Undefined Label:',Label)
         End If
         Len=TocOne(pSym)
         Call iCOPY(Len,Data,1,TocOne(pbas),1)
*----------------------------------------------------------------------*
*
      Else If (label.eq.'NISH') THEN
*
         If(TocOne(pSym).eq.NaN) Then
         Call SysAbendMsg(TheName,'Undefined Label:',Label)
         End If
         Len=TocOne(pSym)
         Call iCOPY(Len,Data,1,TocOne(pISH),1)
*----------------------------------------------------------------------*
*
      Else If (label.eq.'NASH') THEN
*
         If(TocOne(pSym).eq.NaN) Then
         Call SysAbendMsg(TheName,'Undefined Label:',Label)
         End If
         Len=TocOne(pSym)
         Call iCOPY(Len,Data,1,TocOne(pASH),1)
*----------------------------------------------------------------------*
*
      Else If (label.eq.'LDISP') THEN
*
         If(TocOne(pSym).eq.NaN) Then
         Call SysAbendMsg(TheName,'Undefined Label:',Label)
         End If
         Len=TocOne(pSym)
         Call iCOPY(Len,Data,1,TocOne(pldisp),1)
*----------------------------------------------------------------------*
*
      Else If (label.eq.'TDISP') THEN
*
         If(TocOne(pndisp).eq.NaN) Then
         Call SysAbendMsg(TheName,'Undefined Label:',Label)
         End If
         Len=TocOne(pndisp)
         Call iCOPY(Len,Data,1,TocOne(ptdisp),1)
*----------------------------------------------------------------------*
*
      Else If (label.eq.'NDISP') THEN
*
         ToCOne(pnDisp)=Data(1)
*----------------------------------------------------------------------*
*
      Else If (label.eq.'CHDISP') THEN
*
         If(TocOne(pndisp).eq.NaN) Then
         Call SysAbendMsg(TheName,'Undefined Label:',Label)
         End If
         Len=TocOne(pndisp)*30/icpi+1
         Call iCOPY(Len,Data,1,TocOne(pchdisp),1)
*----------------------------------------------------------------------*
*
      Else If (label.eq.'NRDISP') THEN
*
         If(TocOne(pndisp).eq.NaN) Then
         Call SysAbendMsg(TheName,'Undefined Label:',Label)
         End If
         Len=TocOne(pndisp)
         Call iCOPY(Len,Data,1,TocOne(pnrdisp),1)
*----------------------------------------------------------------------*
*
      Else If (label.eq.'DEGDISP ') THEN
*
         If(TocOne(pndisp).eq.NaN) Then
         Call SysAbendMsg(TheName,'Undefined Label:',Label)
         End If
         Len=TocOne(pndisp)
         Call iCOPY(Len,Data,1,TocOne(pdegdisp),1)
*----------------------------------------------------------------------*
*
      Else If(Label.eq.'SYMOP') Then
*
         If (TocOne(pSym).eq.NaN) Then
         Call SysAbendMsg(TheName,'Undefined Label:',Label)
         End If
         Len=(3*TocOne(pSym)+ItoB-1)/ItoB
         Call iCopy(Len,Data,1,TocOne(pSymOp),1)
*----------------------------------------------------------------------*
*
      Else If (label.eq.'PERT') Then
*
         Call icopy(16/icpi,Data,1,TocOne(pPert),1)
*----------------------------------------------------------------------*
      Else
*----------------------------------------------------------------------*
*     Store operators as individual records on disk                    *
*     Note: If the incoming operator has already been stored           *
*     previously (label, component and symmetry labels are identical)  *
*     it will replace the existing one.                                *
*----------------------------------------------------------------------*
         If ((label.eq.'STATHESS').or.
     &       (label.eq.'RESPHESS').or.
     &       (label.eq.'CONNHESS').or.
     &       (label.eq.'HESS').or.
     &       (label.eq.'NUCGRAD').or.
     &       (label.eq.'TWOGRAD')) Then
           Comp=1
           SymLab=1
         End If
*
         If(TocOne(pBas).eq.NaN) Then
            Call SysAbendMsg(TheName,'Undefined Label:',Label)
         End If
*
         k=0
         Do 500 i=MxOp,1,-1
#ifdef _I8_
            If(TocOne(pOp+LenOp*(i-1)+oLabel  ).eq.LabTmp(1) .and.
     &         TocOne(pOp+LenOp*(i-1)+oComp   ).eq.Comp  ) k=i
*     &    .and. TocOne(pOp+LenOp*(i-1)+oSymLb  ).eq.SymLab ) k=i
#else
            If(TocOne(pOp+LenOp*(i-1)+oLabel  ).eq.LabTmp(1) .and.
     &         TocOne(pOp+LenOp*(i-1)+oLabel+1).eq.LabTmp(2) .and.
     &         TocOne(pOp+LenOp*(i-1)+oComp   ).eq.Comp  ) k=i
*     &    .and.TocOne(pOp+LenOp*(i-1)+oSymLb  ).eq.SymLab) k=i
#endif
500      Continue
*
         iDisk=TocOne(pOp+LenOp*(k-1)+oAddr   )
         If (k.eq.0) Then
            Do 505 i=MxOp,1,-1
               If(TocOne(pOp+LenOp*(i-1)+oLabel).eq.NaN) k=i
505         Continue
            iDisk=TocOne(pNext)
            If (Debug) Then
               Write(6,*) ' This is a new field!'
               Write(6,*) ' iDisk=',iDisk
               Write(6,*) ' FldNo=',k
               Write(6,*) ' pNext=',pNext
            End If
         Else
            If(Debug) Then
               Write (6,*) ' This is an old field!'
               Write(6,*) ' iDisk=',iDisk
               Write(6,*) ' FldNo=',k
               Write(6,*) ' pNext=',pNext
            End If
         End If
         If(k.eq.0) Call SysAbendMsg(TheName,'Undefined Label:',Label)
C     Write(*,*) isymlab,label
*                                                                      *
************************************************************************
*                                                                      *
         If (Label.eq.'MOPERT') Then

            nA=0
            Do i=0,TocOne(psym)-1
               nA=TocOne(pAsh+i)+nA
            End Do
            Len=nA*(na+1)/2
            Len=len*(len+1)/2
C           write(*,*) len
*                                                                      *
************************************************************************
*                                                                      *
         Else If (label.eq.'NUCGRAD') Then
             Comp=1
             SymLab=1
             If(TocOne(pldisp).eq.NaN) Then
         Call SysAbendMsg(TheName,'Undefined Label:',Label)
             End If
             Len=TocOne(pldisp)
*                                                                      *
************************************************************************
*                                                                      *
         Else If (label.eq.'TWOGRAD') Then
             Comp=1
             SymLab=1
             If(TocOne(pldisp).eq.NaN) Then
         Call SysAbendMsg(TheName,'Undefined Label:',Label)
             End If
             Len=TocOne(pldisp)
*                                                                      *
************************************************************************
*                                                                      *
         Else If ((label.eq.'STATHESS').or.
     &            (label.eq.'RESPHESS').or.
     &            (label.eq.'CONNHESS').or.
     &            (label.eq.'HESS')) Then
             Comp=1
             SymLab=1
             If(TocOne(pndisp).eq.NaN) Then
         Call SysAbendMsg(TheName,'Undefined Label:',Label)
             End If
             Len=0
             Do iSym=0,TocOne(pSym)-1
               Len=Len+TocOne(pldisp+isym)*
     &                (TocOne(pldisp+isym)+1)/2
             End Do
*                                                                      *
************************************************************************
*                                                                      *
         Else If (label.eq.'INACTIVE') THEN
             Len=0
             Do iS=1,TocOne(pSym)
              Do jS=1,TocOne(pSym)
               ijS=MulTab(iS,jS)
               If(iAnd(2**(ijS-1),iSymLab).ne.0) Then
                 jBas=TocOne(pbas+jS-1)
                 iBas=TocOne(pbas+iS-1)
           If (jBas.eq.NaN)
     * Call SysAbendMsg(TheName,'jBas.eq.NaN at label',Label)
           If (iBas.eq.NaN)
     * Call SysAbendMsg(TheName,'iBas.eq.NaN at label',Label)
                 Len=len+iBas*jBas
               End If
              End Do
             End Do
*                                                                      *
************************************************************************
*                                                                      *
         Else If (label.eq.'TOTAL') THEN
             Len=0
             Do iS=1,TocOne(pSym)
              Do jS=1,TocOne(pSym)
               ijS=MulTab(iS,jS)
               If(iAnd(2**(ijS-1),iSymLab).ne.0) Then
                 jBas=TocOne(pbas+jS-1)
                 iBas=TocOne(pbas+iS-1)
           If (jBas.eq.NaN)
     * Call SysAbendMsg(TheName,'jBas.eq.NaN at label',Label)
           If (iBas.eq.NaN)
     * Call SysAbendMsg(TheName,'iBas.eq.NaN at label',Label)
                 Len=len+iBas*jBas
               End If
              End Do
             End Do
*                                                                      *
************************************************************************
*                                                                      *
         Else
             Len=0
             Do 510 i=1,TocOne(pSym)
              Do 510 j=1,i
               ij=MulTab(i,j)-1
               If(iAnd(2**ij,iSymLab).ne.0) Then
                If(i.eq.j) Then
                  Len=Len+TocOne(pBas-1+i)*(TocOne(pBas-1+i)+1)/2
                Else
                  Len=Len+TocOne(pBas-1+i)*TocOne(pBas-1+j)
                End If
               End If
510          Continue
             If (iAnd(Option,slength).ne.0) Len=Length
         End If
*                                                                      *
************************************************************************
*                                                                      *
         Len=rtoi*Len
         TocOne(pOp+LenOp*(k-1)+oLabel  )=LabTmp(1)
#ifndef _I8_
         TocOne(pOp+LenOp*(k-1)+oLabel+1)=LabTmp(2)
#endif
         TocOne(pOp+LenOp*(k-1)+oComp   )=Comp
         TocOne(pOp+LenOp*(k-1)+oSymLb  )=iSymLab
         TocOne(pOp+LenOp*(k-1)+oAddr   )=iDisk
C           write(*,*) leN,idisk,nauxdt
         Call iDAFile(LuMCK,1,Data,Len+nAuxDt,iDisk)
C         If (Label.eq.'TOTAL'    .or.
C     &       Label.eq.'INACTIVE' .or.
C     &       Label.eq.'MOPERT'       ) Then
CC            Write (*,*) 'iComp=',iComp
CC            Call RecPrt(Label,' ',Data,1,Len/RtoI)
C             Check=DDot_(Len/RtoI,Data,1,Data,1)
C             Label_Add=' '
C             Label_Add=Label
C             Write(Label_Add(9:11),'(I3.3)') iComp
C*            Write (*,*) Label_Add,Label,iComp,Check(1)
C             If (Label.eq.'TOTAL') Then
C                Call Add_Info(Label_Add,Check,1,6)
C             Else If (Label.eq.'INACTIVE') Then
C                Call Add_Info(Label_Add,Check,1,7)
C             Else
C                Call Add_Info(Label_Add,Check,1,1)
C             End If
C         End If
*        TocOne(pNext)=iDisk
         TocOne(pNext)=Max(TocOne(pNext),iDisk)
      End If
*----------------------------------------------------------------------*
*     Finally copy the TocOne back to disk                             *
*----------------------------------------------------------------------*
      iDisk=0
      Call iDaFile(LuMCK,1,TocOne,lToc,iDisk)
*----------------------------------------------------------------------*
*     Print debugging information                                      *
*----------------------------------------------------------------------*
      If(Debug) Then
         Write(6,*) '<<< Exiting WrMck >>>'
         Write(6,'(a,z8)') ' rc on exit:     ',rc
         Write(6,'(a,a)')  ' Label on exit:  ',Label
         Write(6,'(a,z8)') ' Comp on exit:   ',Comp
         Write(6,'(a,z8)') ' SymLab on exit: ',SymLab
         Write(6,'(a,z8)') ' Option on exit: ',Option
      End If
*----------------------------------------------------------------------*
*     Terminate procedure                                              *
*----------------------------------------------------------------------*
*     Call qExit('WrMck')
      Return
      End
