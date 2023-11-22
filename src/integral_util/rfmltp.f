!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1992, Roland Lindh                                     *
!               1994, Markus P. Fuelscher                              *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine RFmltp()
      use PCM_arrays, only: MM
      use stdalloc, only: mma_allocate, mma_deallocate
      use rctfld_module
      Implicit None

      Real*8, Allocatable:: VTot(:), QTot(:)
      Integer nComp
!
      If (.Not.lRF) Return
      nComp = (lMax+1)*(lMax+2)*(lMax+3)/6
      Call mma_allocate(VTot,nComp,Label='VTot')
      Call mma_allocate(QTot,nComp,Label='QTot')
!
      Call RFmltp_Internal(MM,nComp)
!
      Call mma_deallocate(VTot)
      Call mma_deallocate(QTot)
!
      Contains

      Subroutine RFmltp_Internal(Qs,nComp)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!                                                                      *
!     modified by M. P. Fuelscher, 94/04/28                            *
!***********************************************************************
      use Constants, only: Zero, One, Half
      use rctfld_module
      Implicit None
      Integer nComp
      Real*8 Qs(nComp,2)

      Integer l, nElem, iM, jM, iElem, iOff, i, nM
      Real*8 ESolv, dESolv
      Real*8, External:: DDot_
!
      If ( lRF .and. .Not.PCM .and. lRFCav) then
         call dcopy_(nComp,Qs(1,1),1,QTot,1)
         Call DaXpY_(nComp,One,Qs(1,2),1,QTot,1)
#ifdef _DEBUGPRINT_
         Call RecPrt('Total Multipole Moments',' ',QTot,1,nComp)
#endif
         call dcopy_(nComp,QTot,1,VTot,1)
!--------Compute the electric field due to the total charge
!        distribution.
         Call AppFld(VTot,rds,Eps,lMax,EpsInf,NonEq_ref)
#ifdef _DEBUGPRINT_
         Call RecPrt('Total Electric Field',' ',VTot,1,nComp)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
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
               Write (6,'(8X,7ES14.5)') (QTot(i),i=jM,jM+nM-1)
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
               Write (6,'(8X,7ES14.5)') (VTot(i),i=jM,jM+nM-1)
               jM = jM + nM
            End Do
            iM = iM + nElem
         End Do
         Write (6,*) '     -----------------------------------'
         Write (6,*)
      End If
!
      End SubRoutine RFmltp_Internal

      End SubRoutine RFmltp
