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
      Subroutine Print_OpInfo()
#ifdef _EFP_
      Use EFP_Module
      Use EFP
#endif
      use External_Centers
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "rmat.fh"
      Character*72 tempStr
      Character*14 Format_XF
      Logical PrintOperators
      Real*8 A(3)
      Integer iStb(0:7), jCoSet(8,8)
*                                                                      *
************************************************************************
*                                                                      *
*---- Statement Function
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*                                                                      *
************************************************************************
*                                                                      *
      iRout=2
      iPrint = nPrint(iRout)
      If (iPrint.eq.0) Return
      Call qEnter('Print_OpInfo')
      LuWr=6
*                                                                      *
************************************************************************
*                                                                      *
      PrintOperators=.False.
      PrintOperators=PrintOperators.or.(nEF.ne.0)
      PrintOperators=PrintOperators.or.(nDMS.ne.0)
      PrintOperators=PrintOperators.or.(nWel.ne.0)
      PrintOperators=PrintOperators.or.lXF
      PrintOperators=PrintOperators.or.RMat_On
      If (PrintOperators) Then
        Write (LuWr,*)
        Call CollapseOutput(1,'   Operator info:')
        Write (LuWr,'(3X,A)') '   --------------'
        Write (LuWr,*)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (nEF.ne.0) Then
        If (nOrdEF.eq.0) Then
          Write(LuWr,'(2X,A,1X,I8)')
     &      'Centers for electric potential option:',nEF
        Else If (nOrdEF.eq.1) Then
          Write(LuWr,'(2X,A,1X,I8)')
     &      'Centers for electric field option:',nEF
        Else If (nOrdEF.eq.2) Then
          Write(LuWr,'(2X,A,1X,I8)')
     &      'Centers for electric field gradient and contact option:',
     &      nEF
        End If
        Do i=1,nEF
          Write(LuWr,'(4X,I8,3(1X,F14.8))') i,(EF_Centers(j,i),j=1,3)
        End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (nDMS.ne.0) Then
         Call RecPrt(' Gauge Origin for diamagnetic shielding',' ',
     &               Dxyz,1,3)
         Call RecPrt(' Centers for diamagnetic shielding',
     &               ' ',DMS_Centers,3,nDMS)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (nWel.ne.0) Then
         Write (LuWr,*)
         Write (LuWr,*) ' Spherical well specification in au'
         Write (LuWr,*) ' =================================='
         Write (LuWr,*) '   Coeff.      Exp.        R0      '
         Do iWel = 1, nWel
            Write (LuWr,'(3(F10.6,2x))') Wel_Info(3,iWel),
     &                                   Wel_Info(2,iWel),
     &                                   Wel_Info(1,iWel)
         End Do
         Write (LuWr,*)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (lXF) Then
*
         If (nPrint(2).lt.6) Go To 666
         If (iXPolType.gt.0) Then
            tempStr = '       a(xx)       a(xy)    '
     &                //'   a(xz)       a(yy)       a(yz)    '
     &                //'   a(zz)'
         Else
            tempStr = ' '
         EndIf
         Write (LuWr,*)
         Write (LuWr,*) ' External field specification in au'
         Write (LuWr,*) ' =================================='
         If (nOrd_XF.eq.0) Then
            Write (LuWr,*) '     x           y           z      '
     &                //'     Z' // tempStr
         ElseIf (nOrd_XF.eq.1) Then
            Write (LuWr,*) '     x           y           z      '
     &                //'     Z         my(x)       my(y)    '
     &                //'   my(z)' // tempStr
         ElseIf (nOrd_XF.eq.2) Then
            Write (LuWr,*) '     x           y           z      '
     &                //'     Z         my(x)       my(y)    '
     &                //'   my(z)       Q(xx)       Q(xy)    '
     &                //'   Q(xz)       Q(yy)       Q(yz)    '
     &                //'   Q(zz)' // tempStr
         ElseIf (nOrd_XF.eq.-1) Then
            Write (LuWr,*) '     x           y           z ' // tempstr
         Else
            Call WarningMessage(2,'Option not implemented yet!')
            Call Abend
         End If
*
666      Continue
*
         Write(Format_XF,'(A,I2.2,A)') '(',nData_XF,'(F10.6,2x))'
         XnetCharg=0.0
         Do iXF = 1, nXF
            A(1:3)=XF(1:3,iXF)
            Charge_iXF=XF(4,iXF)
            iChxyz=iChAtm(A,iChBas(2))
            iDum=0
            Call Stblz(iChxyz,nStab_iXF,iStb,iDum,jCoSet)
            If (nPrint(2).ge.6)
     &         Write(LuWr,Format_XF) (XF(i,iXF),i=1,nData_XF)
            XnetCharg=XnetCharg+DBLE(nIrrep/nStab_iXF)*Charge_iXF
         End do
         Write (LuWr,*)
         Write (LuWr,*) ' Net charge from external field: ',XnetCharg
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (RMat_On) Then
         Write (LuWr,*)
     &          ' Parameters for radial integration (R-matrix option)'
         Write (LuWr,*)
     &          ' ==================================================='
         Write (LuWr,'(A,G12.5)') '   rmatr     :', RmatR
         Write (LuWr,'(A,G12.5)') '   epsabs    :', Epsabs
         Write (LuWr,'(A,G12.5)') '   epsrel    :', Epsrel
         Write (LuWr,'(A,G12.5)') '   qcoul     :', qCoul
         Write (LuWr,'(A,G12.5)') '   dipol(1)  :', dipol(1)
         Write (LuWr,'(A,G12.5)') '   dipol(2)  :', dipol(2)
         Write (LuWr,'(A,G12.5)') '   dipol(3)  :', dipol(3)
         Write (LuWr,'(A,G12.5)') '   epsq      :', epsq
         Write (LuWr,'(A,G12.5)') '   bparm     :', bParm
      End If
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _EFP_
      If (nEFP_fragments.ne.0) Then
         Call EFP_PRINT_BANNER()
         Write (LuWr,*)
         Write (LuWr,*)
     &          ' Specification of Effective Fragment Potentials'
         Write (LuWr,*)
         If     (Coor_Type.eq.XYZABC_type) Then
            Write (LuWr,*) 'In XYZABC format'
         ElseIf (Coor_Type.eq.POINTS_type) Then
            Write (LuWr,*) 'In Points format'
         ElseIf (Coor_Type.eq.ROTMAT_type) Then
            Write (LuWr,*) 'In RotMat format'
         Else
            Write (LuWr,*) 'Illegal Coor_type:',Coor_Type
            Call Abend()
         End If
         Do i = 1, nEFP_Fragments
            Write (LuWr,*)
            Write (LuWr,*) 'Fragment:',FRAG_TYPE(i)
            If     (Coor_Type.eq.XYZABC_type) Then
            ElseIf (Coor_Type.eq.POINTS_type) Then
               Do j = 1, 3
                  Write (LuWr,'(A10,3F20.10)')
     &               ABC(j,i)(1:10),
     &               (EFP_Coors((j-1)*3+k,i),k=1,3)
               End Do
            ElseIf (Coor_Type.eq.ROTMAT_type) Then
            End If
         End Do
      End If
#endif
*                                                                      *
************************************************************************
*                                                                      *
      If (PrintOperators) Then
        Call CollapseOutput(0,'   Operator info:')
        Write(LuWr,*)
      End If
      Call qExit('Print_OpInfo')
      Return
      End
