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
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      SubRoutine GeoNew(Print)
************************************************************************
*                                                                      *
* Object: to pick up the geometry from a special file. This will only  *
*         make any difference of there exist a file otherwise SEWARD   *
*         will use the geometry as specified by the standard input     *
*         file.                                                        *
*                                                                      *
* Called from: Input                                                   *
*                                                                      *
* Calling    : qEnter                                                  *
*              qExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
*             University of Lund, SWEDEN                               *
*             March 1991                                               *
************************************************************************
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
      Logical Print
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "SysDef.fh"
      Logical Exist
      Real*8, Dimension (:,:), Allocatable :: CN
      Interface
        Subroutine Get_Coord_New(CN,lBuf)
        Real*8, Dimension (:,:), Allocatable :: CN
        Integer lBuf
        End Subroutine
      End Interface

*
*     Prologue
*
*     Check if there is a data field called 'NewGeom'
*
      Call Get_Coord_New(CN,lBuf)
*
*     Quit if the datadfield 'NewGeom' is not available. However,
*     if the field is available on RUNOLD pick it up there.
*
      If ( lBuf.eq.0 ) Then
*
*        Check RUNOLD
*
         Call f_Inquire('RUNOLD',Exist)
         If (Exist) Then
            Call NameRun('RUNOLD')
            Call Get_Coord_New(CN,lBuf)
            If (lBuf.eq.0) Then
               Call qExit('GeoNew')
               nNuc=0
               Call NameRun('RUNFILE')
               Return
            Else
               Call Get_iScalar('Unique atoms',nNuc)
               Call NameRun('RUNFILE')
               If (Print) Then
                  Write (6,*)
                  Write (6,'(A)') '    Geometry read from RUNOLD'
                  Write (6,*)
               End If
            End If
         Else
            Call qExit('GeoNew')
            nNuc=0
            Return
         End If
      Else
         Call Get_iScalar('Unique atoms',nNuc)
         If (Print) Then
            Write (6,*)
            Write (6,'(A)') '    Geometry read from RUNFILE'
            Write (6,*)
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Replace coodinates read in subroutine input
*
      iDC = 1
      iNuc = 0
      Call RecPrt('CN',' ',CN,3,nNuc)
      Do iCnttp = 1, nCnttp
         If (.Not.pChrg(iCnttp).and..Not.FragCnttp(iCnttp) .and.
     &       .Not.AuxCnttp(iCnttp)) Then
            Do iCnt = 1, dbsc(iCnttp)%nCntr
               dbsc(iCnttp)%Coor(1:3,iCnt)=CN(1:3,iDC)
               iDC = iDC + 1
               iNuc = iNuc + 1
               If (iNuc.eq.nNuc) Go To 999
            End Do
         End If
      End Do
 999  Continue
*
*     Now put in the updated coordinates for the RI/CD basis too, note
*     they are pick up from the parent valence basis set.
*
      Do kCnttp = 1, nCnttp
         If (AuxCnttp(kCnttp)) Then
            iCnttp=Parent_iCnttp(kCnttp)
            If (iCnttp.ne.0)
     &         dbsc(kCnttp)%Coor(:,:)=dbsc(iCnttp)%Coor(:,:)
         End If
      End Do
*     Do iCnttp = 1, nCnttp
*        Call RecPrt('dbsc(iCnttp)%Coor',' ',dbsc(iCnttp)%Coor(1,1),
*    &               3,dbsc(iCnttp)%nCntr)
*     End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Epilogue, end
*
      Call mma_deallocate(CN)
      Return
      End
