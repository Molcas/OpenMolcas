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
      Subroutine espf_write(MltOrd,iRMax,DeltaR,iGrdTyp,nGrdPt,DoTinker,
     &                      DoGromacs,lMorok,ipMltp,nMult,ipIsMM,natom,
     &                      Show_espf,Forces,DoDirect)
      Implicit Real*8 (A-H,O-Z)
#include "espf.fh"
#include "stdalloc.fh"
      Real*8, Allocatable:: Grad(:,:)
*
      Logical DoTinker,DoGromacs,lMorok,Show_espf,Forces,DoDirect,Exist
*
* Espf data are saved
*
      iPL = iPL_espf()
*
* Save data in the ESPF.DATA file
*
      IPotFl=12
      IPotFl=IsFreeUnit(IPotFl)
      Call Molcas_Open(IPotFl,'ESPF.DATA')
      Write(IPotFl,'(A10,I10)')    'MLTORD    ',MltOrd/4
      Write(IPotFl,'(A10,I10)')    'IRMAX     ',iRMax
      Write(IPotFl,'(A10,F12.9)')  'DELTAR    ',DeltaR
      Write(IPotFl,'(A10,I10)')    'GRIDTYPE  ',iGrdTyp
      Write(IPotFl,'(A10,I10)')    'GRID      ',nGrdPt
      If (DoTinker)  Write(IPotFl,'(A10)') 'TINKER    '
      If (DoGromacs) Write(IPotFl,'(A10)') 'GROMACS   '
      If (lMorok)    Write(IPotFl,'(A10)') 'LA_MOROK  '
      If (DoDirect)  Write(IPotFl,'(A10)') 'DIRECT    '
      If (ipMltp.ne.ip_Dummy) Then
         Write(IPotFl,'(A10,I10)') 'MULTIPOLE ',nMult
         iMlt = 0
         If (MltOrd.eq.1) Then
            Do iAt = 0, natom-1
               If (iWork(ipIsMM+iAt).eq.0) Then
                  Write(IPotFl,'(I6,4F15.8)') iAt+1,Work(ipMltp+iMlt),
     &                                        Zero,Zero,Zero
                  iMlt = iMlt+1
               End If
            End Do
         Else
            Do iAt = 0, natom-1
               If (iWork(ipIsMM+iAt).eq.0) Then
                  Write(IPotFl,'(I6,4F15.8)') iAt+1,
     &                                   (Work(ipMltp+iMlt+j),j=0,3)
                  iMlt = iMlt+4
               End If
            End Do
         End If
      End If
      Write(IPotFl,'(A10)')        'ENDOFESPF '
      Close (IPotFl)
*
      If (Show_espf .or. iPL.ge.4) Then
         Write(6,'(/,A,/)') ' Informations found in the ESPF data file:'
         Write(6,'(A10,I10)')    ' MLTORD   ',MltOrd/4
         Write(6,'(A10,I10)')    ' IRMAX    ',iRMax
         Write(6,'(A10,F12.9)')  ' DELTAR   ',DeltaR
         Write(6,'(A10,I10)')    ' GRIDTYPE ',iGrdTyp
         Write(6,'(A10,I10)')    ' GRID     ',nGrdPt
         If (DoTinker)  Write(6,'(A10)') ' TINKER   '
         If (DoGromacs) Write(6,'(A10)') ' GROMACS  '
         If (lMorok)    Write(6,'(A10)') ' LA_MOROK '
         If (DoDirect)  Write(6,'(A10)') ' DIRECT   '
         If (ipMltp.ne.ip_Dummy) Then
            Write(6,'(A10,I10)') ' MULTIPOLE ',nMult
            iMlt = 0
            If (MltOrd.eq.1) Then
               Do iAt = 0, natom-1
                  If (iWork(ipIsMM+iAt).eq.0) Then
                     Write(6,'(I6,4F15.8)') iAt+1,
     &                                           Work(ipMltp+iMlt),
     &                                           Zero,Zero,Zero
                     iMlt = iMlt+1
                  End If
               End Do
            Else
               Do iAt = 0, natom-1
                  If (iWork(ipIsMM+iAt).eq.0) Then
                     Write(6,'(I6,4F15.8)') iAt+1,
     &                                      (Work(ipMltp+iMlt+j),j=0,3)
                     iMlt = iMlt+4
                  End If
               End Do
            End If
         End If
         Write(6,'(A10)')        ' ENDOFESPF'
      End If
*
* Special case: Tinker is the driver of the QM/MM calculation.
* QM energy + gradient + ESPF multipoles are stored into the QMMM file
*
      Call F_Inquire('QMMM',Exist)
      If (Exist .and. Forces .and. .not. DoTinker) Then
         ITkQMMM=IsFreeUnit(15)
         Call Molcas_Open(ITkQMMM,'QMMM')
         Call Get_dScalar('Last energy',EQMMM)
         Write(ITkQMMM,'(F12.7,I5)') EQMMM,MltOrd/4
         Call mma_allocate(Grad,3,nAtom,Label='Grad')
         Call Get_Grad(Grad,3*nAtom)
         Do iAt = 1, natom
            iBlaQ = ipMltp + MltOrd*(iAt-1)
            Write(ITkQMMM,'(7F12.7)') Grad(1:3,iAt),
     &                    (Work(iBlaQ+J),J=0,MltOrd-1)
         End Do
         Close(ITkQMMM)
         Call mma_deallocate(Grad)
         Close(ITkQMMM)
      End If
*
      Return
      End
