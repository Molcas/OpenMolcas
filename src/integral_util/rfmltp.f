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
* Copyright (C) 1992, Roland Lindh                                     *
*               1994, Markus P. Fuelscher                              *
************************************************************************
      SubRoutine RFmltp()
      use PCM_arrays, only: MM
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "rctfld.fh"
#include "WrkSpc.fh"
*
      If (.Not.lRF) Return
      nComp = (lMax+1)*(lMax+2)*(lMax+3)/6
      Call GetMem('VTot','Allo','Real',ipVTot,nComp)
      Call GetMem('QTot','Allo','Real',ipQTot,nComp)
*
      Call RFmltp_(MM,Work(ipVTot),Work(ipQTot),nComp)
*
      Call GetMem('QTot','Free','Real',ipQTot,nComp)
      Call GetMem('VTot','Free','Real',ipVTot,nComp)
*
      Return
      End
      Subroutine RFmltp_(Qs,QTot,VTot,nComp)
************************************************************************
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*                                                                      *
*     modified by M. P. Fuelscher, 94/04/28                            *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Real*8 Qs(nComp,2), QTot(nComp), VTot(nComp)
#include "print.fh"
#include "real.fh"
#include "rctfld.fh"
*
      iRout = 4
      iPrint = nPrint(iRout)
      If ( lRF .and. .Not.PCM .and. lRFCav) then
         call dcopy_(nComp,Qs(1,1),1,QTot,1)
         Call DaXpY_(nComp,One,Qs(1,2),1,QTot,1)
         If (iPrint.ge.99) Call RecPrt('Total Multipole Moments',' ',
     &                                 QTot,1,nComp)
         call dcopy_(nComp,QTot,1,VTot,1)
*--------Compute the electric field due to the total charge
*        distribution.
         Call AppFld(VTot,rds,Eps,lMax,EpsInf,NonEq_ref)
         If (iPrint.ge.99) Call RecPrt('Total Electric Field',
     &                                 ' ',VTot,1,nComp)
*                                                                      *
************************************************************************
*                                                                      *
         Write (6,*)
     &   '     Multipole analysis of the contributions to the '//
     &          'dielectric solvation energy'
         Write (6,*)
         Write (6,*) '     --------------------------------------'
         Write (6,*) '        l             dE '
         Write (6,*) '     --------------------------------------'
         Esolv=Zero
         iOff = 1
         Do l = 0, lMax
            nElem = (l+1)*(l+2)/2
            dEsolv= -Half*DDot_(nElem,QTot(iOff),1,VTot(iOff),1)
            Write (6,'(8X,I2,10X,F13.10)') l,dEsolv
            iOff = iOff + nElem
            Esolv = Esolv + dEsolv
         End Do
         Write (6,*) '     --------------------------------------'
         Write (6,*)
         Write (6,*)
         Write (6,*) '     Total Multipole Moments (cartesian)'
         Write (6,*) '     -----------------------------------'
         iM = 1
         Do l = 0, lMax
            nElem = (l+1)*(l+2)/2
            jM = iM
            Do iElem = 1, nElem, 7
               nM=Min(7,nElem-iElem+1)
               Write (6,'(8X,7E14.5)') (QTot(i),i=jM,jM+nM-1)
               jM = jM + nM
            End Do
            iM = iM + nElem
         End Do
         Write (6,*) '     -----------------------------------'
         Write (6,*)
         Write (6,*)
         Write (6,*) '     Total Electric Field (cartesian)'
         Write (6,*) '     --------------------------------'
         iM = 1
         Do l = 0, lMax
            nElem = (l+1)*(l+2)/2
            jM = iM
            Do iElem = 1, nElem, 7
               nM=Min(7,nElem-iElem+1)
               Write (6,'(8X,7E14.5)') (VTot(i),i=jM,jM+nM-1)
               jM = jM + nM
            End Do
            iM = iM + nElem
         End Do
         Write (6,*) '     -----------------------------------'
         Write (6,*)
      End If
*
      Return
      End
