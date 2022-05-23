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
* PAM2009 Who wrote this? What purpose?
      Subroutine Check_InVec(InVec)
#ifdef _MOLCAS_MPP_
      Use Para_Info, Only: nProcs, Is_Real_Par
      If (.Not.Is_Real_Par()) Return
      InVec_Tot = InVec
      Call GAIGOP_SCAL(InVec_Tot,'+')
      If (InVec_Tot.ne.0) Then
         If (InVec.eq.0) Then
            mProcs=-1
         Else
            mProcs = InVec_Tot/InVec
         End If
         If (mProcs.ne.nProcs) Then
            Write (6,*) 'Check_InVec: different orbital options on '
     &                //'different nodes'
            Write (6,*) 'Sets InVec to 0'
            InVec=0
         End If
      End If
#else
c Avoid unused argument warnings
#ifdef _WARNING_WORKAROUND_
      If (.False.) Call Unused_integer(InVec)
#endif
#endif
      Return
      End
