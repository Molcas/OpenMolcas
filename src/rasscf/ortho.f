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
* Copyright (C) 1992, Per Ake Malmqvist                                *
*               1992, Markus P. Fuelscher                              *
************************************************************************
      Subroutine Ortho_RASSCF(SMAT,SCRATCH,CMO,TEMP)
************************************************************************
*                                                                      *
*     purpose: Orthogonalize MOs (one symmetry block at a time)        *
*                                                                      *
*     calling arguments:                                               *
*     Smat    : overlap matrix                                         *
*     CMO     : MO-coefficients                                        *
*     Temp    : temporary work space                                   *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.-AA. Malmqvist and M.P. Fuelscher                              *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
*
      Implicit Real*8 (A-H,O-Z)
*
#include "rasdim.fh"
#include "warnings.fh"
#include "rasscf.fh"
#include "general.fh"
#include "output_ras.fh"
#include "orthonormalize.fh"
      Parameter (ROUTINE='ORTHO   ')
*
      Dimension Smat(*),SCRATCH(*),CMO(*),Temp(*)
*                                                                      *
************************************************************************
*                                                                      *
      Call qEnter('ORTHO')
*                                                                      *
************************************************************************
*                                                                      *
*     Read overlap matrix SMAT:
*
      i_Rc=0
      i_Opt=6
      i_Component=1
      i_SymLbl=1
      Call RdOne(i_Rc,i_Opt,'Mltpl  0',i_Component,Smat,i_SymLbl)
      If ( i_Rc.ne.0 ) Then
        Write(LF,*)' ORTHO could not read overlaps from ONEINT.'
        Write(LF,*)' RASSCF is trying to orthonormalize orbitals but'
        Write(LF,*)' could not read overlaps from ONEINT. Something'
        Write(LF,*)' is wrong with the file, or possibly with the'
        Write(LF,*)' program. Please check.'
        Call Quit(_RC_IO_ERROR_READ_)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Orthonormalize symmetry blocks:
*
      ip_Smat=1
      ip_CMO=1
      Do  iSym=1,nSym
         iBas=nBas(iSym)
         iOcc=iBas-nDel(iSym)
         If ( iBas.gt.0 ) then
*
            Call SQUARE(SMAT(ip_Smat),Temp,1,iBas,iBas)
*
C           Call RecPrt('S',' ',Temp,iBas,iBas)
C           Call RecPrt('CMO',' ',CMO(ip_CMO),iBas,iBas)

            If (Lowdin_ON) Then
*
* --- compute C^T*S*C = W  (overlap in MO basis)
*
               Call DGEMM_('T','N',
     &                     iOcc,iBas,iBas,
     &                     1.0d0,CMO(ip_CMO),iBas,
     &                     Temp,iBas,
     &                     0.0d0,SCRATCH,iOcc)
               Call DGEMM_('N','N',
     &                     iOcc,iOcc,iBas,
     &                     1.0d0,SCRATCH,iOcc,
     &                     CMO(ip_CMO),iBas,
     &                     0.0d0,Temp,iOcc)
*
* --- compute W^-1/2
*
               Call Lowdin(Temp,SCRATCH,iOcc)
*
* --- compute C' = C*W^-1/2
*
               Call DGEMM_('N','N',
     &                     iBas,iOcc,iOcc,
     &                     1.0d0,CMO(ip_CMO),iBas,
     &                     SCRATCH,iOcc,
     &                     0.0d0,Temp,iBas)

* PAM March 2016: Probable bugfix needed (Thanks, Liviu!)
* not affecting any tests (!)
* by adding the following line:
              CALL DCOPY_(iBas*iOcc,Temp,1,CMO(ip_CMO),1)
            Else

               Call ORTHO1(Temp,CMO(ip_CMO),SCRATCH,iBas,iOcc)
            End If
*
C           Call RecPrt('CMO',' ',CMO(ip_CMO),iBas,iBas)
*
            ip_Smat=ip_Smat+(iBas*iBas+iBas)/2
            ip_CMO=ip_CMO+iBas*iBas
         End If
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Call qExit('ORTHO')
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
