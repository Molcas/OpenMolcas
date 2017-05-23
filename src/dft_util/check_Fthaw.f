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
      Subroutine check_Fthaw(iRC)

      Implicit Real*8 (a-h,o-z)
#include "warnings.fh"
      COMMON  / OFembed_T / ThrFThaw
      Character*16 NamRfil
      Logical ok
      Real*8 Ene(1000,4)
*
      If (ThrFThaw.le.0.0d0) Return
*
      Call f_inquire('AUXRFIL',ok)
      If ( .not.ok ) Return
*
      Call Get_NameRun(NamRfil)
      Call NameRun('AUXRFIL')
      Call Get_dScalar('Last energy',EneB)
      Call NameRun(NamRfil)
      Call Get_dScalar('Last energy',EneA)
*
      iSeed=7
      Lu=IsFreeUnit(iSeed)
      Call f_inquire('FRETHAW',ok)
      If ( .not.ok ) Then
         open(Lu,file='FRETHAW')
         write(Lu,'(I4,2F18.10)') 1, EneA, EneB
         Go To 99
      Else
         open(Lu,file='FRETHAW',status='old')
      EndIf
*
      read(Lu,'(I4,2F18.10)') iter0, Ene(1,1), Ene(1,3)
      If (iter0.eq.1000) Then
         write(6,*) ' Error! check_Fthaw: maxIter reached! '
         Call Abend
      EndIf
      Do i=2,iter0
         read(Lu,'(I4,4F18.10)') kiter, Ene(i,1), Ene(i,2), Ene(i,3),
     &                                  Ene(i,4)
      End Do
*
      iter=iter0+1
      DEneA=EneA-Ene(iter0,1)
      DEneB=EneB-Ene(iter0,3)
*
      Rewind Lu
      write(Lu,'(I4,2F18.10)') iter, Ene(1,1), Ene(1,3)
      Do i=2,iter0
         write(Lu,'(I4,4F18.10)') iter, Ene(i,1), Ene(i,2), Ene(i,3),
     &                                  Ene(i,4)
      End Do
      write(Lu,'(I4,4F18.10)') iter, EneA, DEneA, EneB, DEneB
*
      write(6,*)
      write(6,*) '**************************************************'//
     &           '*****************************'
      write(6,*) '*************** Energy Statistics for Freeze-n-Thaw'//
     &           ' ***************************'
      write(6,*) '**************************************************'//
     &           '*****************************'
      write(6,*) '         Energy_A       Delta(Energy_A)      '//
     &           'Energy_B       Delta(Energy_B)'
      write(6,'(I3,1X,F18.10,18X,F18.10)') 1, Ene(1,1), Ene(1,3)
      Do i=2,iter0
         write(6,'(I3,1X,4F18.10)') i, Ene(i,1), Ene(i,2), Ene(i,3),
     &                                Ene(i,4)
      End Do
      write(6,'(I3,1X,4F18.10)') iter, EneA, DEneA, EneB, DEneB
      write(6,*) '**************************************************'//
     &           '*****************************'
*
      If ( abs(DEneA).lt.ThrFThaw .and.
     &     abs(DEneB).lt.ThrFThaw ) Then
         write(6,'(A,E9.2,A)')' Convergence reached ! (Thr = ',ThrFThaw,
     &                        ')'
         write(6,*)
         iRC=_RC_ALL_IS_WELL_
         Close(Lu,status='delete')
         Return
      Else
         write(6,'(A,E9.2,A)')' Convergence NOT reached yet ! (Thr = ',
     &                        ThrFThaw,')'
         write(6,*)
      EndIf
*
99    Continue
      Close(Lu,status='keep')
*
      Return
      End
