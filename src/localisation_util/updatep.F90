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
! Copyright (C) Yannick Carissan                                       *
!               2005, Thomas Bondo Pedersen                            *
!***********************************************************************
      Subroutine UpdateP(PACol,Name,nBas_Start,                         &
     &                   nOrb2Loc,nAtoms,PA,gamma_rot,                  &
     &                   iMO_s,iMO_t,Debug)
!
!     Author: Yannick Carissan.
!
!     Modifications:
!        - October 6, 2005 (Thomas Bondo Pedersen):
!          Reduce operation count and use BLAS.
!
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "Molcas.fh"
      Integer nBas_Start(*)
      Real*8 PACol(nOrb2Loc,2)
      Real*8 PA(nOrb2Loc,nOrb2Loc,nAtoms)
      Character*(LENIN8) Name(*),PALbl
      Logical Debug
!
      cosg   = cos(gamma_rot)
      sing   = sin(gamma_rot)
      cos2g  = cosg*cosg
      sin2g  = sing*sing
      cosing = cosg*sing
!
      Do iAt=1,nAtoms
!       Call RecPrt('PA(1,1,iAt)',' ',PA(1,1,iAt),nOrb2Loc,nOrb2Loc)
!
!------ Copy out the PAss, PAtt, and PAst elements.
!
        PAss = PA(iMO_s,iMO_s,iAt)
        PAst = PA(iMO_s,iMO_t,iAt)
        PAtt = PA(iMO_t,iMO_t,iAt)
!       Write (6,*) 'updateP:',PAss,PAst,PAtt
#if defined (_DEBUGPRINT_)
        PAts = PA(iMO_t,iMO_s,iAt)
        Tst  = PAst - PAts
        If (abs(Tst) .gt. 1.0d-14) Then
           Write(6,*) 'Broken symmetry in UpdateP!!'
           Write(6,*) 'MOs s and t: ',iMO_s,iMO_t
           Write(6,*) 'PAst = ',PAst
           Write(6,*) 'PAts = ',PAts
           Write(6,*) 'Diff = ',Tst
           Call SysAbendMsg('UpdateP','Broken symmetry!',' ')
        End If
#endif
!
!------ Copy out columns s and t of PA.
!
        Call dCopy_(nOrb2Loc,PA(1,iMO_s,iAt),1,PACol(1,1),1)
        Call dCopy_(nOrb2Loc,PA(1,iMO_t,iAt),1,PACol(1,2),1)
!
!------ Compute transformed columns.
!
        Call dScal_(nOrb2Loc,cosg,PA(1,iMO_s,iAt),1)
        Call dAXPY_(nOrb2Loc, sing,PACol(1,2),1,PA(1,iMO_s,iAt),1)
        Call dScal_(nOrb2Loc,cosg,PA(1,iMO_t,iAt),1)
        Call dAXPY_(nOrb2Loc,-sing,PACol(1,1),1,PA(1,iMO_t,iAt),1)
!
!------ Compute PAss, PAtt, PAst, and PAts (= PAst).
!
        PA(iMO_s,iMO_s,iAt)= PAss*cos2g + PAtt*sin2g                    &
     &                     + Two*PAst*cosing
        PA(iMO_t,iMO_s,iAt)= (PAtt-PAss)*cosing + PAst*(cos2g-sin2g)
        PA(iMO_s,iMO_t,iAt)= PA(iMO_t,iMO_s,iAt)
        PA(iMO_t,iMO_t,iAt)= PAtt*cos2g + PAss*sin2g                    &
     &                     - Two*PAst*cosing
!
!------ Copy columns to rows.
!
        Call dCopy_(nOrb2Loc,PA(1,iMO_s,iAt),1,PA(iMO_s,1,iAt),nOrb2Loc)
        Call dCopy_(nOrb2Loc,PA(1,iMO_t,iAt),1,PA(iMO_t,1,iAt),nOrb2Loc)
!
      End Do
!
      If (Debug) Then
        Write(6,*) 'In UpdateP'
        Write(6,*) '----------'
        Do iAt=1,nAtoms
          PALbl='PA__'//Name(nBas_Start(iAt))(1:LENIN)
          Call RecPrt(PALbl,' ',PA(1,1,iAt),nOrb2Loc,nOrb2Loc)
        End Do
      End If
!
      Return
      End
