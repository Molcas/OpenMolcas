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
      Subroutine Put_Chunk(ip_ChoVec,MuNu_s,MuNu_e,j_s,j_e,Rv,nMuNu,
     &                     LenVec)
#ifdef _MOLCAS_MPP_
      Use Para_Info, Only: Is_Real_Par
#endif
      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"
      Real*8 Rv(nMuNu,(j_e-j_s+1))
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _MOLCAS_MPP_
*
      NumVec_ = j_e - j_s + 1
      If (NumVec_ .gt. 0) Then
         If (Is_Real_Par()) Then
            Call GA_Put(ip_ChoVec,MuNu_s,MuNu_e,j_s,j_e,Rv,nMuNu)
         Else
            mMuNu=MuNu_s-1
            jp_ChoVec=ip_ChoVec+mMuNu
            Do jVec = 1, NumVec_
               call dcopy_(nMuNu,Rv(1,jVec),1,Work(jp_ChoVec),1)
               jp_ChoVec = jp_ChoVec + LenVec
            End Do
         End If
      End If
*
#else
*
      mMuNu=MuNu_s-1
      NumVec_ = j_e - j_s + 1
*
      jp_ChoVec=ip_ChoVec+mMuNu
      Do jVec = 1, NumVec_
         call dcopy_(nMuNu,Rv(1,jVec),1,Work(jp_ChoVec),1)
         jp_ChoVec = jp_ChoVec + LenVec
      End Do
*
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(MuNu_e)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
