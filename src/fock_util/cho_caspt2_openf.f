************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SubRoutine Cho_CASPT2_OpenF(iOpt,iTyp,iSym,nBatch)
C
C
C     Purpose: open (iOpt=1), close and keep (iOpt=2), or close and
C              delete (iOpt=3) Cholesky vector files for CASPT2 program
C              (full vectors).
C              For iOpt=0, the units are initialized (to -1).
C              iTyp=1: transformed Cholesky vectors (Inactive).
C              iTyp=2: transformed Cholesky vectors (Active).
C              iTyp=3: vectors from (pi|qj) decomposition (Not
C                      implemented yet!)
C
#include "implicit.fh"
#include "WrkSpc.fh"
      Character*16 SecNam
      Parameter (SecNam = 'Cho_CASPT2_OpenF')
      Character*3 BaseNm
      Character*7 FullNm
#include "chocaspt2.fh"
      DIMENSION NUMCHO(8)

      Save NCALLS
      Integer NCALLS
      Data NCALLS /0/

*******************************************************************
      lUnit_F(j,k,l) = iWork(ipUnit_F(j)+(l-1)*nIsplit(j)+k-1)
*******************************************************************

      If (nBatch .gt. 999) Then
         Call Cho_x_Quit(SecNam,' nBatch limited to 999 !!!',' ')
      End If
      Call Get_iScalar('nSym',nSym)
      Call Get_iArray('NumCho',NumCho,nSym)

      If (NCALLS.EQ.0) THEN
        Do iB=1,nBatch
          iWork(ipUnit_F(iSym)+(iTyp-1)*nIsplit(iSym)+iB-1) = -1
        End Do
      End If

C     Initialize units and return for iOpt=0.
C     ---------------------------------------

      If (iOpt .eq. 0) Then
         Do iB=1,nBatch
            iWork(ipUnit_F(iSym)+(iTyp-1)*nIsplit(iSym)+iB-1) = -1
         End Do
         Return
      End If

C     Open or close files.
C     --------------------
      If (iTyp.lt.1 .or. iTyp.gt.2) Then
         Call Cho_x_Quit(SecNam,'iTyp error',' ')
      End If

      If (iOpt .eq. 1) Then
         If (NumCho(iSym).gt.0) Then
            Do iB=1,nBatch
               If(lUnit_F(iSym,iB,iTyp) .lt. 1) Then
                 Call Cho_caspt2_GetBaseNm(BaseNm,iTyp)
                 Write(FullNm,'(A3,I1,I3)') BaseNm,iSym,iB
                 LuV=7 ! initial guess
                 Call daName_MF_WA(LuV,FullNm)!handle inquire/free unit
                 iWork(ipUnit_F(iSym)+(iTyp-1)*nIsplit(iSym)+iB-1) = LuV
      write(6,*)' Opened file ''',FullNm,''' as unit nr LuV=',LuV
      iaddr=ipUnit_F(iSym)+(iTyp-1)*nIsplit(iSym)+iB-1
      write(6,*)' Unit number LuV is stored at address ',iaddr
               End If
            End Do
         Else
            Do iB=1,nBatch
               iWork(ipUnit_F(iSym)+(iTyp-1)*nIsplit(iSym)+iB-1) = -1
            End Do
         End If
      Else If (iOpt .eq. 2) Then
         Do iB=1,nBatch
            If (lUnit_F(iSym,iB,iTyp) .gt. 0) Then
      write(6,*)' Closing lUnit_F(iSym,iB,iTyp)=',lUnit_F(iSym,iB,iTyp)
               Call daClos(lUnit_F(iSym,iB,iTyp))
               iWork(ipUnit_F(iSym)+(iTyp-1)*nIsplit(iSym)+iB-1) = -1
            End If
         End Do
      Else If (iOpt .eq. 3) Then
         Do iB=1,nBatch
            If (lUnit_F(iSym,iB,iTyp) .gt. 0) Then
      write(6,*)' Erasing lUnit_F(iSym,iB,iTyp)=',lUnit_F(iSym,iB,iTyp)
               Call daEras(lUnit_F(iSym,iB,iTyp))
               iWork(ipUnit_F(iSym)+(iTyp-1)*nIsplit(iSym)+iB-1) = -1
            End If
         End Do
      Else
         Call Cho_x_Quit(SecNam,'iOpt out of bounds',' ')
      End If

      End


      SubRoutine Cho_caspt2_GetBaseNm(BaseNm,iTyp)
      Implicit None
      Character*3 BaseNm
      Integer     iTyp

      If (iTyp .eq. 1) Then
         BaseNm = '_PI'
      Else If (iTyp .eq. 2) Then
         BaseNm = '_PW'
      Else If (iTyp .eq. 3) Then
         BaseNm = '_CD'
      Else
         BaseNm = '_un'
      End If

      End
