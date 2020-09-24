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
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
*               1998, Roland Lindh                                     *
*               2003, Valera Veryazov                                  *
************************************************************************
      SubRoutine SOrb(LuOrb,SIntTh,iTerm)
      use SCF_Arrays
      Implicit Real*8 (a-h,o-z)
#include "mxdm.fh"
#include "infscf.fh"
*
      nD = iUHF + 1
      Call SOrb_(LuOrb,SIntTh,iTerm,CMO,TrM,nBB,nD,OneHam,Fock,Ovrlp,
     &           nBT,Eorb,OccNo,nnB)
*
      Return
      End
      SubRoutine SOrb_(LuOrb,SIntTh,iTerm,CMO,TrM,mBB,nD,OneHam,Fock,
     &                 Ovrlp,mBT,EOrb,OccNo,mmB)
************************************************************************
*                                                                      *
*     purpose: Get starting orbitals from:                             *
*              -1) default choice                                      *
*               0) diagonalizaton of the core                          *
*               1) via intermediate calculation of HF AOs              *
*               2) input orbitals                                      *
*               3) input density matrix                                *
*                                                                      *
*     called from: SCF                                                 *
*                                                                      *
*     calls to: Start0, Start2, Start3, SorbChk                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history:                                                         *
*                                                                      *
*     Modified by R. Lindh, May  98, Tokyo, Japan                      *
*     UHF - V.Veryazov, 2003                                           *
*                                                                      *
************************************************************************
*
      Implicit Real*8 (a-h,o-z)
*
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "infso.fh"
#include "file.fh"
#ifdef _HDF5_
#  include "mh5.fh"
#endif
      Real*8 CMO(mBB,nD), TrM(mBB,nD), OneHam(mBT), Fock(mBT,nD),
     &       Ovrlp(mBT), EOrb(mmB,nD), OccNo(mmB,nD)
      Character FName*512, KSDFT_save*16
      Logical FstItr
      Logical found
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
#ifdef _DEBUGPRINT_
      Call qEnter('SOrb')
#endif
*
      CALL DecideonCholesky(DoCholesky)
*-------- Cholesky and NDDO are incompatible
      IF (DoCholesky.and.InVec.eq.1) THEN
         call WarningMessage(1,
     &    ' In SORB: Cholesky and NDDO not implemented !!!; '//
     &    ' NDDO option ignored')
         InVec=-1
      ENDIF
*
*---- Is default clause choosen?
*
      If(InVec.eq.-1) Then
         Call qpg_darray('SCF orbitals',found,ndata)
         If(found .and. nData.eq.nBB) Then
            Call qpg_darray('OrbE',found,ndata)
            If(found) Then
               InVec=8
            End If
         End If
      End If
      If(InVec.eq.-1) Then
         Call qpg_darray('Guessorb',found,ndata)
         If(found .and. nData.eq.nBB) Then
            Call qpg_darray('Guessorb energies',found,ndata)
            If(found) Then
               InVec=9
            End If
         End If
      End If
      If(InVec.eq.-1) Then
         InVec=0
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Reset some parameters depending on which set of starting orbitals
*     we are using. If we have a good starting set of vectors we want
*     the Quasi Newton Raphson to kick in as soon as possible. Further,
*     we would like to use a single numerical grid in the calculation
*     and we would like to a single threshold for the direct SCF.
*
* No no no, not good for difficult cases!!!
*
*     if(iUHF.eq.0) Then
*     If (InVec.ne.0) Then
*        QNRTh =1.0D0
*        DiisTh=1.0D0
*        One_Grid=.True.
*        Two_Thresholds=.False.
*     End If
*     End If
*
*                                                                      *
************************************************************************
*                                                                      *

*
*---- Has the user selected a method?
*
100   Continue
      If (InVec.eq.0) Then
*-------- Diagonalize core
          Call Start0(CMO,TrM,mBB,nD,OneHam,Ovrlp,mBT,EOrb,mmB)
      Else If (InVec.eq.6) Then
         write(6,*)
         write(6,*) '     Constrained SCF calculation '
         write(6,*)
         StVec='Constrained orbitals'
         One_Grid=.True.
         FName=SCF_FileOrb
         Call Chk_Vec_UHF(FName,LuOrb,isUHF)
         If (IsUHF.eq.1) Then
            InVec=2
            Go To 100
         EndIf
         Call Start6(FName,LuOrb,CMO,mBB,nD,EOrb,OccNo,mmB)
      Else If (InVec.eq.8) Then
         StVec='Detected old SCF orbitals'
         One_Grid=.True.
         Call start0y(CMO,mBB,nD,EOrb,mmB)
      Else If (InVec.eq.9) Then
         StVec='Detected guessorb starting orbitals'
         One_Grid=.True.
         Call start0x(CMO,mBB,nD,EOrb,mmB)
      Else If (InVec.eq.2) Then
*-------- Read INPORB
         One_Grid=.True.
         FName=SCF_FileOrb
         Call Start2(FName,LuOrb,CMO,mBB,nD,Ovrlp,mBT,
     &               EOrb,OccNo,mmB)
      Else If (InVec.eq.3) Then
*-------- Read COMOLD
         One_Grid=.True.
         Call Start3(CMO,TrM,mBB,nD,OneHam,Ovrlp,mBT)
*-------- Only if not Cholesky do NDDO
      Else If (InVec.eq.1) Then
*
*-------- HF AO orbitals as intermediate step...
*
*------- NDDO, always none-DFT
*
         Call SwiOpt(.False.,OneHam,Ovrlp,mBT,CMO,mBB,nD)
         Call Start0(CMO,TrM,mBB,nD,OneHam,Ovrlp,mBT,EOrb,mmB)
         InVec=0
         Call SOrbCHk(OneHam,Ovrlp,Fock,mBT,nD,CMO,mBB)
         KSDFT_save=KSDFT
         KSDFT='SCF'
         Call WrInp_SCF(SIntTh)
         FstItr=.True.
         Call WfCtl_SCF(iTerm,'NDDO      ',FstItr,SIntTh)
         KSDFT=KSDFT_save
         Call Free_TLists
         If (iTerm.ne.0) Call Quit(iTerm)
         Write(6,*)
         Write(6,'(A)') 'Generation of NDDO vectors completed!'
         Write(6,*)
         Write(6,*) '2nd step: optimizing HF MOs...'
         Write(6,*) '------------------------------'
         Call SwiOpt(.TRUE.,OneHam,Ovrlp,mBT,CMO,mBB,nD)
*------- Reset to to start from the current MO set
         Call Init_SCF()
         InVec=5
!IFG: I presume the arguments after LuOut in these two calls are correct,
!     they were missing!
         If(iUHF.eq.0) Then
            FName='SCFORB'
            Call Start2(FName,LuOut,CMO,mBB,nD,Ovrlp,mBT,
     &               EOrb,OccNo,mmB)
         Else
            FName='UHFORB'
            Call Start2(FName,LuOut,CMO,mBB,nD,Ovrlp,mBT,
     &               EOrb,OccNo,mmB)
         End If
      End If
      If (Scrmbl) Then
         Do iD = 1, nD
            Call Scram(CMO(1,iD),nSym,nBas,nOrb,ScrFac)
         End Do
      End If
      Call SOrbCHk(OneHam,Ovrlp,Fock,mBT,nD,CMO,mBB)
#ifdef _HDF5_
      If (isHDF5) Call mh5_close_file(fileorb_id)
#endif
#ifdef _DEBUGPRINT_
      Call qExit('SOrb')
#endif
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
      Return
      End
