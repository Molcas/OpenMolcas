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
************************************************************************
      Subroutine RdOne(rc,Option,InLab,Comp,Data,SymLab)
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
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Integer (A-Z)
*
#include "OneRc.fh"
#include "OneFlags.fh"

#include "OneDat.fh"
*
      Character*(*) InLab
      Dimension Data(*)
*
      Character*8 TmpLab,Label
*     Dimension LabTmp(2)
*     Equivalence (TmpLab,LabTmp)
*
      Parameter (lBuf=1024)
      Real*8    TmpBuf(lBuf),AuxBuf(4)
      Logical debug, Close
      Data CurrOp/1/
      Save CurrOp
*----------------------------------------------------------------------*
*     Start procedure:                                                 *
*     Define inline function (symmetry multiplication)                 *
*----------------------------------------------------------------------*
      MulTab(i,j)=iEor(i-1,j-1)+1
*----------------------------------------------------------------------*
*     Pick up the file definitions                                     *
*----------------------------------------------------------------------*
*     Call qEnter('RdOne')
      rc    = rc0000
      LuOne = AuxOne(pLu  )
      Open  = AuxOne(pOpen)
*----------------------------------------------------------------------*
*     Check the file status                                            *
*----------------------------------------------------------------------*
      Close=.False.
      If ( Open.ne.1 ) Then
*        Write (*,*) ' I will open the file for you!'
*
*------- Well, I'll open and close it for you under the default name
*
         LuOne=77
         LuOne=isFreeUnit(LuOne)
         Label='ONEINT  '
*        Write(*,*) 'RdOne: opening OneInt'
         iRC=-1
         iOpt=0
         Call OpnOne(iRC,iOpt,Label,LuOne)
         If (iRC.ne.0) Then
            Write (6,*) 'RdOne: Error opening file'
            Call Abend
         End If
         Close=.True.
      End If
*----------------------------------------------------------------------*
*     Truncate the label to 8 characters and convert it to upper case  *
*----------------------------------------------------------------------*
      Label=InLab
      Call UpCase(Label)
      TmpLab=Label
*----------------------------------------------------------------------*
*     Print debugging information                                      *
*----------------------------------------------------------------------*
      debug=.false.
      If(iAnd(option,1024).ne.0) debug=.true.
      If(debug) Then
         Write(6,*) '<<< Entering RdOne >>>'
         Write(6,'(a,z8)') ' rc on entry:     ',rc
         Write(6,'(a,a)')  ' Label on entry:  ',Label
         Write(6,'(a,z8)') ' Comp on entry:   ',Comp
         Write(6,'(a,z8)') ' SymLab on entry: ',SymLab
         Write(6,'(a,z8)') ' Option on entry: ',Option
      End If
*----------------------------------------------------------------------*
*     Check reading mode                                               *
*----------------------------------------------------------------------*
      If((iAnd(iAnd(option,sRdFst),sRdNxt)).ne.0) then
         Write (6,*) 'RdOne: Invalid option(s)'
         Write (6,*) 'option=',option
         Call QTrace()
         Call Abend()
      Else If((iAnd(iAnd(option,sRdFst),sRdCur)).ne.0) then
         Write (6,*) 'RdOne: Invalid option(s)'
         Write (6,*) 'option=',option
         Call QTrace()
         Call Abend()
      Else If((iAnd(iAnd(option,sRdNxt),sRdCur)).ne.0) then
         Write (6,*) 'RdOne: Invalid option(s)'
         Write (6,*) 'option=',option
         Call QTrace()
         Call Abend()
      End If
*----------------------------------------------------------------------*
*     Load back TocOne                                                 *
*----------------------------------------------------------------------*
      iDisk=0
      Call iDaFile(LuOne,2,TocOne,lToc,iDisk)
*----------------------------------------------------------------------*
*     Read data from ToC                                               *
*----------------------------------------------------------------------*
      NoGo=sRdFst+sRdNxt+sRdCur
*----------------------------------------------------------------------*
*     Read operators from integral records                             *
*----------------------------------------------------------------------*
      If (iAnd(option,sRdNxt).ne.0) Then
         CurrOp=CurrOp+1
         If (CurrOp.gt.MxOp) Then
            CurrOp=0
         Else If(TocOne(pOp+LenOp*(CurrOp-1)+oLabel).eq.Nan) Then
            CurrOp=0
         Else
            i=CurrOp
*            LabTmp(1)=TocOne(pOp+LenOp*(i-1)+oLabel  )
*#ifndef _I8_
*            LabTmp(2)=TocOne(pOp+LenOp*(i-1)+oLabel+1)
*#endif
            Call ByteCopy(TocOne(pOp+LenOp*(i-1)+oLabel),TmpLab,8)
            Label=TmpLab
            InLab=Label
            SymLab=TocOne(pOp+LenOp*(i-1)+oSymLb)
            Comp=TocOne(pOp+LenOp*(i-1)+oComp   )
         End If
      Else If(iAnd(option,sRdFst).ne.0) Then
         CurrOp=1
         If(TocOne(pOp+LenOp*(CurrOp-1)+oLabel).eq.Nan) Then
            CurrOp=0
         Else
            i=CurrOp
*            LabTmp(1)=TocOne(pOp+LenOp*(i-1)+oLabel  )
*#ifndef _I8_
*            LabTmp(2)=TocOne(pOp+LenOp*(i-1)+oLabel+1)
*#endif
            Call ByteCopy(TocOne(pOp+LenOp*(i-1)+oLabel),TmpLab,8)
            Label=TmpLab
            InLab=Label
            SymLab=TocOne(pOp+LenOp*(i-1)+oSymLb)
            Comp=TocOne(pOp+LenOp*(i-1)+oComp   )
         End If
      Else If(iAnd(option,sRdCur).ne.0) Then
         If(CurrOp.lt.1 .or. CurrOp.gt.MxOp) Then
            CurrOp=0
         Else If(TocOne(pOp+LenOp*(CurrOp-1)+oLabel).eq.Nan) Then
            CurrOp=0
         Else
            i=CurrOp
*            LabTmp(1)=TocOne(pOp+LenOp*(i-1)+oLabel  )
*#ifndef _I8_
*            LabTmp(2)=TocOne(pOp+LenOp*(i-1)+oLabel+1)
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
*            LabTmp(1)=TocOne(pOp+LenOp*(i-1)+oLabel  )
*#ifndef _I8_
*            LabTmp(2)=TocOne(pOp+LenOp*(i-1)+oLabel+1)
*#endif
            Call ByteCopy(TocOne(pOp+LenOp*(i-1)+oLabel),TmpLab,8)
            CmpTmp=TocOne(pOp+LenOp*(i-1)+oComp   )
            TmpCmp=Comp
            If(TmpLab.eq.Label .and. CmpTmp.eq.TmpCmp) CurrOp=i
500      Continue
      End If
      If(CurrOp.eq.0) Then
         rc=rcRD03
*        Write (*,*) 'RdOne: Information not available'
*        Write (*,*) 'Option=',Option
*        Write (*,*) 'Comp=',Comp
*        Write (*,*) 'SymLab=',SymLab
*        Write (*,*) 'Label=',Label
         Go To 999
      End If
      SymLab=TocOne(pOp+LenOp*(CurrOp-1)+oSymLb)
      Len=0
      Do 510 i=1,nSym
      Do 510 j=1,i
         ij=MulTab(i,j)-1
         If(iAnd(2**ij,SymLab).ne.0) Then
            If(i.eq.j) Then
               Len=Len+nBas(i)*(nBas(i)+1)/2
            Else
               Len=Len+nBas(i)*nBas(j)
            End If
         End If
510   Continue
      Data(1)=Len
      If ( IAND(option,sOpSiz).eq.0 ) Then
         IndAux = 0
         IndDta = 0
         iDisk=TocOne(pOp+LenOp*(CurrOp-1)+oAddr)
         Do i = 0,Len+3,lBuf
           nCopy  = MAX(0,MIN(lBuf,Len+4-i))
           nSave  = MAX(0,MIN(lBuf,Len-i))
           Call dDaFile(LuOne,2,TmpBuf,nCopy,iDisk)
           Call dCopy_(nSave,TmpBuf,1,Data(IndDta+1),1)
           IndDta = IndDta+RtoI*nSave
           Do j = nSave+1,nCopy
             IndAux = IndAux+1
*            AuxBuf(IndAux) = TmpBuf(nSave+IndAux)
             AuxBuf(IndAux) = TmpBuf(j)
           End Do
         End Do
         If(iAnd(sNoOri,option).eq.0) Then
            Call dCopy_(3,AuxBuf,1,Data(IndDta+1),1)
         End If
         If(iAnd(sNoNuc,option).eq.0) Then
            Call dCopy_(1,AuxBuf(4),1,Data(IndDta+RtoI*3+1),1)
         End If
      End If
*#define _DMET_
#ifdef _DMET_
      Call PrMtrx("ONEINT inside Rdone",1,Comp,1,Data)
#endif
*
 999  Continue
      If (Close) Then
*        Write (*,*) ' I will close the file for you!'
         iRC=-1
         iOpt=0
         Call ClsOne(iRC,iOpt)
         If (iRC.ne.0) Then
            Write (6,*) 'RdOne: Error closing file'
            Call Abend
         End If
      End If
*----------------------------------------------------------------------*
*     Terminate procedure                                              *
*----------------------------------------------------------------------*
*     Call qExit('RdOne')
      Return
      End

      Subroutine iRdOne(rc,Option,InLab,Comp,iData,SymLab)
      Implicit Integer (A-Z)
*
      Character*(*) InLab
      Dimension iData(*)
*
      Call RdOne(rc,Option,InLab,Comp,iData,SymLab)
      return
      end
