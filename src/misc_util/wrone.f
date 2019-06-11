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
      Subroutine iWrOne(rc,Option,InLab,Comp,Data,SymLab)
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
      Dimension LabTmp(2)
*     Equivalence (TmpLab,LabTmp)
      Logical debug, Close
*----------------------------------------------------------------------*
*     Start procedure:                                                 *
*     Define inline function (symmetry multiplication)                 *
*----------------------------------------------------------------------*
      MulTab(i,j)=iEor(i-1,j-1)+1
*----------------------------------------------------------------------*
*     Pick up the file definitions                                     *
*----------------------------------------------------------------------*
*     Call qEnter('WrOne')
      rc    = rc0000
      LuOne = AuxOne(pLu  )
      Open  = AuxOne(pOpen)
*----------------------------------------------------------------------*
*     Check the file status                                            *
*----------------------------------------------------------------------*
      Close=.False.
      If ( Open.ne.1 ) Then
*
*------- Well, I'll open and close it for you under the default name
*
         LuOne=77
         LuOne=isFreeUnit(LuOne)
         Label='ONEINT  '
*        Write(*,*) 'WrOne: opening OneInt'
         iRC=-1
         iOpt=0
         Call OpnOne(iRC,iOpt,Label,LuOne)
         If (iRC.ne.0) Then
            Write (6,*) 'WrOne: Error opening file'
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
      Call ByteCopy(TmpLab,LabTmp,8)
*----------------------------------------------------------------------*
*     Print debugging information                                      *
*----------------------------------------------------------------------*
      debug=.false.
      If(iAnd(option,1024).ne.0) debug=.true.
      If(debug) Then
         Call DmpOne
         Write(6,*) '<<< Entering WrOne >>>'
         Write(6,'(a,z8)') ' rc on entry:     ',rc
         Write(6,'(a,a)')  ' Label on entry:  ',Label
         Write(6,'(a,z8)') ' Comp on entry:   ',Comp
         Write(6,'(a,z8)') ' SymLab on entry: ',SymLab
         Write(6,'(a,z8)') ' Option on entry: ',Option
      End If
*----------------------------------------------------------------------*
*     Store operators as individual records on disk                    *
*     Note: If the incoming operator has already been stored           *
*     previously (label, component and symmetry labels are identical)  *
*     it will replace the existing one.                                *
*----------------------------------------------------------------------*
         k=0
         Do 500 i=MxOp,1,-1
#ifdef _I8_
            If(TocOne(pOp+LenOp*(i-1)+oLabel  ).eq.LabTmp(1) .and.
     &         TocOne(pOp+LenOp*(i-1)+oComp   ).eq.Comp      .and.
     &         TocOne(pOp+LenOp*(i-1)+oSymLb  ).eq.SymLab          ) k=i
#else
            If(TocOne(pOp+LenOp*(i-1)+oLabel  ).eq.LabTmp(1) .and.
     &         TocOne(pOp+LenOp*(i-1)+oLabel+1).eq.LabTmp(2) .and.
     &         TocOne(pOp+LenOp*(i-1)+oComp   ).eq.Comp      .and.
     &         TocOne(pOp+LenOp*(i-1)+oSymLb  ).eq.SymLab          ) k=i
#endif
500      Continue
         iDisk=TocOne(pOp+LenOp*(k-1)+oAddr   )
         If(k.eq.0) Then
           Do 505 i=MxOp,1,-1
              If(TocOne(pOp+LenOp*(i-1)+oLabel).eq.NaN) k=i
505        Continue
           iDisk=TocOne(pNext)
         End If
         If(k.eq.0) Then
            rc=rcWR11
            Write (6,*) 'WrOne: The total number of operators',
     &                  ' exceeds the limit'
            Write (6,*) 'k.eq.0'
            Call QTrace()
            Call Abend()
         End If
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
510      Continue
         Len=RtoI*(Len+nAuxDt)
         TocOne(pOp+LenOp*(k-1)+oLabel  )=LabTmp(1)
#ifndef _I8_
         TocOne(pOp+LenOp*(k-1)+oLabel+1)=LabTmp(2)
#endif
         TocOne(pOp+LenOp*(k-1)+oComp   )=Comp
         TocOne(pOp+LenOp*(k-1)+oSymLb  )=SymLab
         TocOne(pOp+LenOp*(k-1)+oAddr   )=iDisk
         Call iDaFile(LuOne,1,Data,Len,iDisk)
         TocOne(pNext)=Max(TocOne(pNext),iDisk)
*----------------------------------------------------------------------*
*     Finally copy the TocOne back to disk                             *
*----------------------------------------------------------------------*
      iDisk=0
      Call iDaFile(LuOne,1,TocOne,lToc,iDisk)
*
      If (Close) Then
         iRC=-1
         iOpt=0
         Call ClsOne(iRC,iOpt)
         If (iRC.ne.0) Then
            Write (6,*) 'WrOne: Error closing file'
            Call Abend
         End If
      End If
*----------------------------------------------------------------------*
*     Terminate procedure                                              *
*----------------------------------------------------------------------*
*     Call qExit('WrOne')
      Return
      End

      Subroutine WrOne(rc,Option,InLab,Comp,Data,SymLab)
      Implicit Integer (A-Z)
*
      Character*(*) InLab
      Real*8 Data(*)
*
      Call WrOne_Internal(Data)
*
*     This is to allow type punning without an explicit interface
      Contains
      Subroutine WrOne_Internal(Data)
      Use Iso_C_Binding
      Real*8, Target :: Data(*)
      Integer, Pointer :: iData(:)
      Call C_F_Pointer(C_Loc(Data(1)),iData,[1])
      Call iWrOne(rc,Option,InLab,Comp,iData,SymLab)
      Nullify(iData)
      return
      End Subroutine WrOne_Internal
*
      end
