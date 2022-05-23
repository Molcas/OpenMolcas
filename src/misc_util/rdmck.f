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
      Subroutine RdMCK(rc,Option,InLab,iComp,Data,iSymLab)
************************************************************************
*                                                                      *
*     Purpose: Read data from one-electron integral file               *
*                                                                      *
*     Calling parameters:                                              *
*     rc      : Return code                                            *
*     Option  : Switch to set options                                  *
*     InLab   : Identifier for the data to read                        *
*     Comp    : Composite identifier to select components              *
*     Data    : contains on output the requested data                  *
*     SymLab  : symmetry label of the requested data                   *
*                                                                      *
*     Global data declarations (Include file):                         *
*     Parm    : Table of paramaters                                    *
*     rcParm  : Table of return codes                                  *
*     Switch  : Table of options                                       *
*     Common  : Common block containing ToC                            *
*     Data    : Data definitions                                       *
*                                                                      *
*     Local data declarations:                                         *
*     Label   : character*8, used to covert incoming names             *
*     TmpBuf  : I/O buffer                                             *
*     HldBuf  : I/O buffer                                             *
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
*
      Logical :: Debug=.False.
      Character(LEN=8) TmpLab,Label
*     Dimension LabTmp(2)
*     Equivalence (TmpLab,LabTmp)
      Dimension TmpBuf(nBuf),HldBuf(1)
      Character(LEN=16) :: TheName= 'RdMck'
      Integer :: CurrOp=1
*----------------------------------------------------------------------*
*     Start procedure:                                                 *
*     Define inline function (symmetry multiplication)                 *
*----------------------------------------------------------------------*
      MulTab(i,j)=iEor(i-1,j-1)+1
*----------------------------------------------------------------------*
*     Pick up the file definitions                                     *
*----------------------------------------------------------------------*
      icpi=itob
      Length=rc
      rc    = rc0000
      LuMck = AuxMck(pLu  )
      Open  = AuxMck(pOpen)
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
      Comp=icomp
      Symlab=isymlab
      If ((label.eq.'STATHESS').or.
     &    (label.eq.'RESPHESS').or.
     &    (label.eq.'CONNHESS').or.
     &    (label.eq.'HESS').or.
     &    (label.eq.'NUCGRAD').or.
     &    (label.eq.'TWOGRAD')) Then
       Comp=1
       Symlab=1
      End If
*----------------------------------------------------------------------*
*     Print debugging information                                      *
*----------------------------------------------------------------------*
      If(iAnd(option,1024).ne.0) debug=.true.
      If(Debug) Then
         Write(6,*) '<<< Entering RdMck  >>>'
         Write(6,'(a,z8)') ' rc on entry:     ',rc
         Write(6,'(a,a)')  ' Label on entry:  ',Label
         Write(6,'(a,z8)') ' Comp on entry:   ',Comp
         Write(6,'(a,z8)') ' SymLab on entry: ',SymLab
         Write(6,'(a,z8)') ' Option on entry: ',Option
      End If
*----------------------------------------------------------------------*
*     Check reading mode                                               *
*----------------------------------------------------------------------*
      If(iAnd(iAnd(option,sRdFst),sRdNxt).ne.0) Then
         Call SysWarnMsg(TheName,'Invalid value',' ')
         Call SysCondMsg('iAnd(iAnd(option,sRdFst),sRdNxt).eq.0',
     &           iAnd(iAnd(option,sRdFst),sRdNxt),'<>',0)
      Else If(iAnd(iAnd(option,sRdFst),sRdCur).ne.0) Then
         Call SysWarnMsg(TheName,'Invalid value',' ')
         Call SysCondMsg('iAnd(iAnd(option,sRdFst),sRdCur).eq.0',
     &           iAnd(iAnd(option,sRdFst),sRdCur),'<>',0)
      Else If(iAnd(iAnd(option,sRdNxt),sRdCur).ne.0) Then
         Call SysWarnMsg(TheName,'Invalid value',' ')
         Call SysCondMsg('iAnd(iAnd(option,sRdNxt),sRdCur).eq.0',
     &           iAnd(iAnd(option,sRdNxt),sRdCur),'<>',0)
      End If
*----------------------------------------------------------------------*
*     Read data from ToC                                               *
*----------------------------------------------------------------------*
      NoGo=sRdFst+sRdNxt+sRdCur
      If(Label.eq.'TITLE' .and. iAnd(option,NoGo).eq.0) Then
         If(TocOne(pTitle).eq.NaN) Then
         Call SysAbendMsg(TheName,'Undefined Label:',Label)

         End If
         If(iAnd(option,sOpSiz).eq.0) Then
            Call iCopy(nTitle,TocOne(pTitle+1),1,Data(1),1)
            If(debug) Then
               Write(6,'(a,z8)') ' Reading Title:'
               Write(6,'(8(1x,z8))') (Data(k),k=1,nTitle)
            End If
         Else
            Data(1)=nTitle
            If(debug) Then
               Write(6,'(a,z8)') ' Reading Title:',Data(1)
            End If
         End If
      else If(Label.eq.'CHDISP' .and. iAnd(option,NoGo).eq.0) Then
         If(TocOne(pCHDISP).eq.NaN) Then
         Call SysAbendMsg(TheName,'Undefined Label:',Label)
         End If
         If(iAnd(option,sOpSiz).eq.0) Then
            Len=TocOne(pnDisp)*30/icpi+1
            Call iCopy(Len,TocOne(pchdisp),1,Data(1),1)
            If(debug) Then
               Write(6,'(a,z8)') ' Reading perturbations:'
               Write(6,'(8(1x,z8))') (Data(k),k=1,nTitle)
            End If
         Else
            Data(1)=TocOne(pnDisp)*30/icpi+1
            If(debug) Then
               Write(6,'(a,z8)') ' Reading perturbations:',Data(1)
            End If
         End If
      Else If(Label.eq.'NDISP' .and. iAnd(option,NoGo).eq.0) Then
         If(TocOne(pndisp).eq.NaN) Then
         Call SysAbendMsg(TheName,'Undefined Label:',Label)
         End If
         If(iAnd(option,sOpSiz).eq.0) Then
            Data(1)=TocOne(pndisp)
         Else
            Data(1)=1
         End If
         If(debug) Then
            Write(6,'(a,z8)') ' Reading nSym: ',Data(1)
         End If
      Else If(Label.eq.'NSYM' .and. iAnd(option,NoGo).eq.0) Then
         If(TocOne(pSym).eq.NaN) Then
         Call SysAbendMsg(TheName,'Undefined Label:',Label)
         End If
         If(iAnd(option,sOpSiz).eq.0) Then
            Data(1)=TocOne(pSym)
         Else
            Data(1)=1
         End If
         If(debug) Then
            Write(6,'(a,z8)') ' Reading nSym: ',Data(1)
         End If
      Else If(Label.eq.'NBAS' .and. iAnd(option,NoGo).eq.0) Then
         If(TocOne(pBas).eq.NaN) Then
         Call SysAbendMsg(TheName,'Undefined Label:',Label)
         End If
         If(iAnd(option,sOpSiz).eq.0) Then
            Len=TocOne(pSym)
            Call iCopy(Len,TocOne(pBas),1,Data,1)
            If(debug) Then
               Write(6,'(a,8z8)') ' Reading nBas: ',(Data(k),k=1,Len)
            End If
         Else
            Data(1)=TocOne(pSym)
            If(debug) Then
               Write(6,'(a,z8)') ' Reading nBas: ',Data(1)
            End If
         End If
      Else If(Label.eq.'LDISP' .and. iAnd(option,NoGo).eq.0) Then
         If(TocOne(pldisp).eq.NaN) Then
         Call SysAbendMsg(TheName,'Undefined Label:',Label)
         End If
         If(iAnd(option,sOpSiz).eq.0) Then
            Len=TocOne(psym)
            Call iCopy(Len,TocOne(pldisp),1,Data,1)
            If(debug) Then
               Write(6,'(a,8z8)') ' Reading ldisp: ',(Data(k),k=1,Len)
            End If
         Else
            Data(1)=TocOne(pSym)
            If(debug) Then
               Write(6,'(a,z8)') ' Reading ldisp: ',Data(1)
            End If
         End If
         Else If(Label.eq.'TDISP' .and. iAnd(option,NoGo).eq.0) Then
          If(TocOne(ptdisp).eq.NaN) Then
         Call SysAbendMsg(TheName,'Undefined Label:',Label)
          End If
          If(iAnd(option,sOpSiz).eq.0) Then
             Len=TocOne(pndisp)
             Call iCopy(Len,TocOne(ptdisp),1,Data,1)
             If(debug) Then
                Write(6,'(a,8z8)') ' Reading nBas: ',(Data(k),k=1,Len)
             End If
          Else
             Data(1)=TocOne(pSym)
             If(debug) Then
                Write(6,'(a,z8)') ' Reading nBas: ',Data(1)
             End If
          End If
      Else If(Label.eq.'NASH' .and. iAnd(option,NoGo).eq.0) Then
         If(TocOne(pASH).eq.NaN) Then
         Call SysAbendMsg(TheName,'Undefined Label:',Label)
         End If
         If(iAnd(option,sOpSiz).eq.0) Then
            Len=TocOne(pSYM)
            Call iCopy(Len,TocOne(pASH),1,Data,1)
            If(debug) Then
               Write(6,'(a,8z8)') ' Reading nASH: ',(Data(k),k=1,Len)
            End If
         Else
            Data(1)=TocOne(pSym)
            If(debug) Then
               Write(6,'(a,z8)') ' Reading nASH: ',Data(1)
            End If
         End If
      Else If (label.eq.'PERT') Then
         Call icopy(16/icpi,TocOne(pPert),1,Data,1)
      Else If (label.eq.'NRCTDISP') THEN
*
         If(TocOne(pndisp).eq.NaN) Then
         Call SysAbendMsg(TheName,'Undefined Label:',Label)
         End If
         Len=TocOne(pndisp)
         Call iCOPY(Len,TocOne(pnrdisp),1,Data,1)
*
      Else If(Label.eq.'SYMOP' .and. iAnd(option,NoGo).eq.0) Then
         If (TocOne(pSym).eq.NaN) Then
         Call SysAbendMsg(TheName,'Undefined Label:',Label)
         End If
         If(iAnd(option,sOpSiz).eq.0) Then
*           Len=(3*TocOne(pSym)-1)/icpi+1
            Len=(3*TocOne(pSym)+ItoB-1)/ItoB
            Call iCopy(Len,TocOne(pSymOp),1,Data,1)
            If(debug) Then
               Write(6,'(a)') ' Reading symmetry operators:'
               Write(6,'(8(1x,z8))') (Data(k),k=1,Len)
            End If
         Else
            Data(1)=TocOne(pSym)
            If(debug) Then
               Write(6,'(a,z8)') ' Reading symmetry operators:',Data(1)
            End If
         End If
      Else If (label.eq.'DEGDISP ') THEN
*
         If(TocOne(pndisp).eq.NaN) Then
         Call SysAbendMsg(TheName,'Undefined Label:',Label)
         End If
         Len=TocOne(pndisp)
         Call iCOPY(Len,TocOne(pdegdisp),1,Data,1)
*
      Else
*----------------------------------------------------------------------*
*     Read operators from integral records                             *
*----------------------------------------------------------------------*
         If(iAnd(option,sRdNxt).ne.0) Then
            If(debug) Then
               Write(6,'(a)') ' Reading next item'
            End If
            CurrOp=CurrOp+1
            If(CurrOp.gt.MxOp) Then
               CurrOp=0
            Else If(TocOne(pOp+LenOp*(CurrOp-1)+oLabel).eq.Nan) Then
               CurrOp=0
            Else
               i=CurrOp
*               LabTmp(1)=TocOne(pOp+LenOp*(i-1)+oLabel  )
*#ifndef _I8_
*               LabTmp(2)=TocOne(pOp+LenOp*(i-1)+oLabel+1)
*#endif
               Call ByteCopy(TocOne(pOp+LenOp*(i-1)+oLabel),TmpLab,8)
               Label=TmpLab
               InLab=Label
               SymLab=TocOne(pOp+LenOp*(i-1)+oSymLb)
               Comp=TocOne(pOp+LenOp*(i-1)+oComp   )
            End If
         Else If(iAnd(option,sRdFst).ne.0) Then
            If(debug) Then
               Write(6,'(a)') ' Reading first item'
            End If
            CurrOp=1
            If(TocOne(pOp+LenOp*(CurrOp-1)+oLabel).eq.Nan) Then
               CurrOp=0
            Else
               i=CurrOp
*               LabTmp(1)=TocOne(pOp+LenOp*(i-1)+oLabel  )
*#ifndef _I8_
*               LabTmp(2)=TocOne(pOp+LenOp*(i-1)+oLabel+1)
*#endif
               Call ByteCopy(TocOne(pOp+LenOp*(i-1)+oLabel),TmpLab,8)
               Label=TmpLab
               InLab=Label
               SymLab=TocOne(pOp+LenOp*(i-1)+oSymLb)
               Comp=TocOne(pOp+LenOp*(i-1)+oComp   )
            End If
         Else If(iAnd(option,sRdCur).ne.0) Then
            If(debug) Then
               Write(6,'(a)') ' Reading current item'
            End If
            If(CurrOp.lt.1 .or. CurrOp.gt.MxOp) Then
               CurrOp=0
            Else If(TocOne(pOp+LenOp*(CurrOp-1)+oLabel).eq.Nan) Then
               CurrOp=0
            Else
               i=CurrOp
*               LabTmp(1)=TocOne(pOp+LenOp*(i-1)+oLabel  )
*#ifndef _I8_
*               LabTmp(2)=TocOne(pOp+LenOp*(i-1)+oLabel+1)
*#endif
               Call ByteCopy(TocOne(pOp+LenOp*(i-1)+oLabel),TmpLab,8)
               Label=TmpLab
               InLab=Label
               SymLab=TocOne(pOp+LenOp*(i-1)+oSymLb)
               Comp=TocOne(pOp+LenOp*(i-1)+oComp   )
            End If
         Else
            CurrOp=0
            Do 500 i=MxOp,1,-1
*               LabTmp(1)=TocOne(pOp+LenOp*(i-1)+oLabel  )
*#ifndef _I8_
*               LabTmp(2)=TocOne(pOp+LenOp*(i-1)+oLabel+1)
*#endif
               Call ByteCopy(TocOne(pOp+LenOp*(i-1)+oLabel),TmpLab,8)
               CmpTmp=TocOne(pOp+LenOp*(i-1)+oComp   )
               TmpCmp=Comp
               If(TmpLab.eq.Label .and. CmpTmp.eq.TmpCmp) CurrOp=i
500         Continue
         End If
         If(CurrOp.eq.0) Then
         Call SysAbendMsg(TheName,'Current Operation .eq. 0',' ')
         End If
         SymLab=TocOne(pOp+LenOp*(CurrOp-1)+oSymLb)
*
         If (Label.eq.'MOPERT') Then
            iTmp=SymLab
            Do iIrr=1,TocOne(pSym)
             iTmp=iTmp/2
            End Do
            na=0
            Do i=0,TocOne(psym)-1
             na=TocOne(pAsh+i)+na
            End Do
            Len=na*(na+1)/2
            Len=len*(len+1)/2
         Else If (label.eq.'NUCGRAD') Then
             Comp=1
             SymLab=1
             If(TocOne(pldisp).eq.NaN) Then
         Call SysAbendMsg(TheName,'Undefined Label:',Label)
             End If
             Len=TocOne(pldisp)
         Else If (label.eq.'TWOGRAD') Then
             Comp=1
             SymLab=1
             If(TocOne(pldisp).eq.NaN) Then
         Call SysAbendMsg(TheName,'Undefined Label:',Label)
             End If
             Len=TocOne(pldisp)
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
*
         Else If (label.eq.'INACTIVE') THEN
             Len=0
             Do iS=1,TocOne(pSym)
              Do jS=1,TocOne(pSym)
               ijS=MulTab(iS,jS)
               If(iAnd(2**(ijS-1),SymLab).ne.0) Then
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

         Else If (label.eq.'TOTAL') THEN
             Len=0
             Do iS=1,TocOne(pSym)
              Do jS=1,TocOne(pSym)
               ijS=MulTab(iS,jS)
               If(iAnd(2**(ijS-1),SymLab).ne.0) Then
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
          Else
*
            Len=0
            Do 510 i=1,TocOne(pSym)
            Do 511 j=1,i
              ij=MulTab(i,j)-1
              If(iAnd(2**ij,SymLab).ne.0) Then
                 If(i.eq.j) Then
                    Len=Len+TocOne(pBas-1+i)*(TocOne(pBas-1+i)+1)/2
                 Else
                    Len=Len+TocOne(pBas-1+i)*TocOne(pBas-1+j)
                 End If
              End If
511         Continue
510         Continue
            If (iAnd(Option,slength).ne.0) Len=Length
         End If

         Data(1)=Len
         If ( IAND(option,sOpSiz).eq.0 ) Then
           Len=rtoi*Len
           IndDta=1
           IndHld=1
           iDisk=TocOne(pOp+LenOp*(CurrOp-1)+oAddr   )
           Do 550 k=0,Len+nAuxDt-1,nBuf
              iBuf=Max(0,Min(nBuf,Len+nAuxDt-k))
              tBuf=Max(0,Min(nBuf,Len-k))
              eBuf=iBuf-tBuf
              Call iDaFile(LuMCK,2,TmpBuf,iBuf,iDisk)
              IndTmp=1
              Call iCopy(tBuf,TmpBuf(IndTmp),1,Data(IndDta),1)
              If(debug) Then
                 Write(6,'(a,z8)') ' Reading buffer to: ',IndDta
                 Write(6,'(8(1x,z8))') (Data(IndDta+m),m=0,tBuf-1)
              End If
              IndTmp=IndTmp+tBuf
              IndDta=IndDta+tBuf
              If(tBuf.lt.iBuf) Then
                 Call iCopy(eBuf,TmpBuf(IndTmp),1,HldBuf(IndHld),1)
                 IndTmp=IndTmp+eBuf
                 IndHld=IndHld+eBuf
              End If
550        Continue
*          If(iAnd(sNoOri,option).eq.0) Then
*             Call iCopy(6,HldBuf(1),1,Data(IndDta),1)
*             If(debug) Then
*                Write(*,'(a,z8)') ' Reading buffer to: ',IndDta
*                Write(*,'(8(1x,z8))') (Data(IndDta+m),m=0,5)
*             End If
*          End If
           IndDta=IndDta+6
*          If(iAnd(sNoNuc,option).eq.0) Then
*             Call iCopy(2,HldBuf(7),1,Data(IndDta),1)
*             If(debug) Then
*                Write(*,'(a,z8)') ' Reading buffer to: ',IndDta
*                Write(*,'(8(1x,z8))') (Data(IndDta+m),m=0,1)
*             End If
*          End If
*          IndDta=IndDta+2
        End If
      End If
*----------------------------------------------------------------------*
*     Print debugging information                                      *
*----------------------------------------------------------------------*
      If(Debug) Then
         Write(6,*) '<<< Exiting RdMck >>>'
         Write(6,'(a,z8)') ' rc on exit:     ',rc
         Write(6,'(a,a)')  ' Label on exit:  ',Label
         Write(6,'(a,z8)') ' Comp on exit:    ',Comp
         Write(6,'(a,z8)') ' SymLab on exit: ',SymLab
         Write(6,'(a,z8)') ' Option on exit: ',Option
      End If
*----------------------------------------------------------------------*
*     Terminate procedure                                              *
*----------------------------------------------------------------------*
      Return
      End
