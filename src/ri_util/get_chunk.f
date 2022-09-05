!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      Subroutine Get_Chunk(LenVec,NumVec_,iChoVec,iSym,iVec_Global)
      use Chunk_Mod
#ifdef _MOLCAS_MPP_
      Use Para_Info, Only: MyRank, Is_Real_Par
#endif
      Implicit Real*8 (A-H,O-Z)
#ifdef _MOLCAS_MPP_
#include "mafdecls.fh"
#endif
!
#ifdef _MOLCAS_MPP_
      If (Is_Real_Par()) Then
         Call GA_Sync()
         Call Cho_RI_SwapVecUnit(iSym)
!
!        Get the subrange and check that we are within range.
!
         J_s=iMap(1+MyRank)
         If (J_s.gt.NumVec_) Go To 999
         J_e=iMap(1+MyRank+1)-1
         J_e=Min(J_e,NumVec_)
         nJ=J_e-J_s+1
         If (nJ.gt.0) Then
            MuNu_s=1
            MuNu_e=LenVec
            Call GA_Access(ip_Chunk,MuNu_s,MuNu_e,J_s,J_e,Index,ld)
!           Call RecPrt('Dbl_Mb(Index)',' ',Dbl_Mb(Index),LenVec,nJ)
            Call Cho_PutVec(DBL_MB(Index),LenVec,nJ,iChoVec+1,iSym)
            Call Cho_RI_SetInfVec_5(iVec_Global,iChoVec+1,J_s,J_e,iSym)
            Call GA_Release(ip_Chunk,MuNu_s,MuNu_e,J_s,J_e)
            iChoVec = iChoVec + nJ
         End If
  999    Continue
         Call Cho_RI_SwapVecUnit(iSym)
      Else
         Call Cho_PutVec(Chunk,LenVec,NumVec_,iChoVec+1,iSym)
         iChoVec = iChoVec + NumVec_
      End If
!
#else
      Call Cho_PutVec(Chunk,LenVec,NumVec_,iChoVec+1,iSym)
      iChoVec = iChoVec + NumVec_
! Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(iVec_Global)
      End If
#endif
!
      Return
      End
#if defined (_MOLCAS_MPP_)
      SubRoutine Cho_RI_SetInfVec_5(iVec_Global,iVec_Local,J_s,J_e,iSym)
!
!     Set mapping from local to global vector index (needed in parallel
!     RI gradient code).
!
      use ChoSwp, only: InfVec
      Implicit None
      Integer iVec_Global, iVec_Local, J_s, J_e, iSym
#include "cholesky.fh"

      Integer iOff, nVec, iVec

      iOff = iVec_Global + J_s - 2
      nVec = J_e - J_s + 1
      Do iVec = 1,nVec
         InfVec(iVec_Local-1+iVec,5,iSym) = iOff + iVec
      End Do

      End
#endif
