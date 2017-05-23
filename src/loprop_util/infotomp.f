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
      Subroutine InfoToMp(nSym,nBas,Energy_Without_FFPT,ip_Ene_Occ,
     &                   nOcOb,UserDen,Restart)
      Implicit real*8(a-h,o-z)

#include "WrkSpc.fh"
      Integer nBas(8)
      Character*8 Method
      Character*40 VTitle
      Character*6 FName
      Logical UserDen,Restart
*                                                                      *
************************************************************************
*                                                                      *
* Set some variables needed for creating the MpProp file.              *
*                                                                      *
************************************************************************
*                                                                      *
      NOCOB=0
      If(.not.UserDen) then
         nVec=0
         nOcc=0
         Do iSym = 1, nSym
            nVec = nVec + nBas(iSym)**2
            nOcc = nOcc + nBas(iSym)
         End Do
         Call Allocate_Work(ip_Ene_Occ,nOcc)
         If (Restart) Then
            Call Get_dScalar('MpProp Energy',Energy_Without_FFPT)
            Call Get_dArray('MpProp Orb Ener',Work(ip_Ene_occ),nOcc)
            Call Get_iScalar('MpProp nOcOb',nOcOb)
         Else
            Call Get_dScalar('Last energy',Energy_Without_FFPT)
            Call Put_dScalar('MpProp Energy',Energy_Without_FFPT)

            Call Allocate_Work(ip_Vec,nVec)
            Call Allocate_Work(ip_Occ,nOcc)

            Lu_   = 11
            FName ='INPORB'
            iDum  = 0
            iWarn = 2
            Call RdVec(FName,Lu_,'COE',nSym,nBas,nBas,Work(ip_Vec),
     &           Work(ip_Occ),Work(ip_Ene_Occ),iDum,VTitle,iWarn,iErr)
            Close(Lu_)

            Do i=0,nOcc-1
               If (Work(ip_Occ+i).ne.0.0D0) Then
               nOcOb = nOcOb + 1
               End If
            End Do
            Call Put_dArray('MpProp Orb Ener',Work(ip_Ene_occ),nOcc)
            Call Put_iScalar('MpProp nOcOb',nOcOb)

            Call Free_Work(ip_Vec)
            Call Free_Work(ip_Occ)
         End If
      Else  !Here we go if user give density as input. Need to put some
            !dummy values for later; also we give a value to the relax
            !method, which for now is called 'external'.
        Energy_Without_FFPT=0.0D0
        nOcc=0
        Do iSym = 1, nSym
           nOcc = nOcc + nBas(iSym)
        End Do
        Call Allocate_Work(ip_Ene_Occ,nOcc)
        Do i=0,nOcc-1
          Work(ip_Ene_Occ+i)=0.0D0
        Enddo
        Write(Method,'(A)')'External'
        Call Put_cArray('Relax Method',Method,8)
      Endif

      Return
      End
