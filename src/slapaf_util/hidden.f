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
      Subroutine Hidden(mTtAtm,ipCoor,ipAN,nHidden,rHidden,nMDstep)
      Implicit Real*8 (a-h,o-z)
*
*  Add to the Grand atom list some hidden atoms, coming e.g.
*  from the MM part of a QM/MM system
*
#include "Molcas.fh"
#include "angstr.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "periodic_table.fh"
      Logical Do_ESPF,Exist,Exist2
      Character(LEN=180) Line
      Character(LEN=180), External::  Get_Ln
      Real*8 XYZ(3)
      Character(LEN=2) Symbol
      Character(LEN=LENIN), Allocatable:: LabMMO(:)
      Real*8, Allocatable:: h_xyz(:,:)
      Integer, Allocatable:: h_AN(:)
*
      iPL = iPrintLevel(-1)
*     iPL=99
*
*#define _DEBUGPRINT_
*
#ifdef _DEBUGPRINT_
      iPL = 4
#endif
      nMDstep = 0
      iHidden = 0
      Do_ESPF = .False.
*
*  Is there a ESPF/QMMM file ?
*
      Call DecideOnESPF(Do_ESPF)
      If (Do_ESPF) Then
*
*        Try MM atoms from Tinker QM/MM interface file
*
         Call F_Inquire('QMMM',Exist)
         If (Exist) Then
            ITkQMMM=IsFreeUnit(25)
            Call Molcas_Open(ITkQMMM,'QMMM')
            Line = ' '
            Do While (Index(Line,'TheEnd ') .eq. 0)
               Line=Get_Ln(ITkQMMM)
*
*  Read the maximum number of MM atoms + some temporary arrays allocation
*
               If (Index(Line,'NMM').ne.0) Then
                  Call Get_I1(2,nHidden)
                  If (iPL.gt.3) Write(6,'(A,I5,A)')'Found ',nHidden,
     &                                    ' hidden atoms.'
                  If(nHidden.gt.0) Then
                     Call mma_allocate(h_xyz,3,nHidden,Label='h_xyz')
                     Call mma_allocate(h_AN,nHidden,Label='h_AN')
                     Do iHidden = 1, nHidden
                        Line=Get_Ln(ITkQMMM)
                        If (Index(Line,'MMCoord').eq.0) Then
                           Write(6,*) 'Error in hidden.',
     &                      ' Last line does not start with MMCoord:'
                           Write(6,*) Line
                           Call Quit_onUserError()
                        End If
                        Call Get_I1(2,iAtNum)
                        h_AN(iHidden) = -iAtNum
                        Call Get_F(3,XYZ,3)
                        h_xyz(1,iHidden) = XYZ(1)/Angstr
                        h_xyz(2,iHidden) = XYZ(2)/Angstr
                        h_xyz(3,iHidden) = XYZ(3)/Angstr
                     End Do
                  End If
               Else If (Index(Line,'MD ').ne.0) Then
                  Call Get_I1(2,nMDstep)
               End If
            End Do
            Close (ITkQMMM)
         End If
*
*        Try outer MM atoms stored on runfile
*
         If (.Not.Exist) Then
            Call Qpg_dArray('MMO Coords',Exist2,nHidden)
         Else
            Exist2 = .False.
         End If
         If (Exist2) Then
            nHidden=nHidden/3
            Call mma_allocate(h_xyz,3,nHidden,Label='h_xyz')
            Call mma_allocate(h_AN,nHidden,Label='h_AN')
            Call mma_allocate(LabMMO,nHidden,Label='LabMMO')
            Call Get_dArray('MMO Coords',h_xyz,nHidden*3)
            Call Get_cArray('MMO Labels',LabMMO,LENIN*nHidden)
            Do iHidden = 1, nHidden
               Symbol(1:1) = LabMMO(iHidden)(1:1)
               Symbol(2:2) = LabMMO(iHidden)(2:2)
               If (Symbol(2:2).Eq.'_') Symbol = ' '//Symbol(1:1)
               Do i = 0, Num_Elem
                  If (Ptab(i) == Symbol) Then
                     h_AN(iHidden) = -i
                     Exit
                  End If
               End Do
            End Do
            Call mma_deallocate(LabMMO)
         End If
      End If
      If (iPL.gt.3) Call RecPrt('Hidden coord:',' ',h_xyz,3,nHidden)
*
*  Select the hidden atoms to be kept.
*
      nKept = 0
      If (nHidden .gt. 0)
     &   Call Select_hidden(mTtAtm,nHidden,Work(ipCoor),h_xyz,h_AN,
     &                      nKept,rHidden,iPL)
*
*  Copy all the arrays needed by box and nlm
*
      If (nKept .gt. 0) Then
         If (iPL .gt. 3) Then
            Write(6,'(A8,I5,A)') 'Hidden: ',nKept,' atoms are kept.'
            If (nMDstep.gt.0) Write(6,'(8X,I5,A)') nMDstep,' mean Hess'
         End If
         mTot = mTtAtm + nKept
         Call Allocate_Work(ipCoor_h,3*mTot)
         Call Allocate_iWork(ipAN_h,mTot)
*
         Call dCopy_(3*mTtAtm,Work(ipCoor),1,Work(ipCoor_h),1)
         Call iCopy(mTtAtm,iWork(ipAN),1,iWork(ipAN_h),1)
*
*  Copy the kept hidden atom coordinates, atom numbers and masses
*
         iKept = 0
         Do iHidden = 0, (nHidden-1)
            If (h_AN(1+iHidden) .gt. 0) Then
               Call dCopy_(3,h_xyz(:,1+iHidden),1,
     &                     Work(ipCoor_h+3*(mTtAtm+iKept)),1)
               iAN = h_AN(1+iHidden)
               iWork(ipAN_h+mTtAtm+iKept) = iAN
               iKept = iKept + 1
            End If
         End Do
         If (iKept .ne. nKept) Then
            Write(6,'(A)') ' Hidden: wrong number of kept hidden atoms.'
            Call Quit_OnUserError()
         End If
         Call mma_deallocate(h_AN)
         Call mma_deallocate(h_xyz)
         Call GetMem('Coor','Free','Real',ipCoor,3*mTtAtm)
         Call GetMem('AN','Free','Inte',ipAN,mTtAtm)
*
*  The end
*
         ipCoor  = ipCoor_h
         ipAN    = ipAN_h
*
         If (iPL .gt. 3) Then
            Call RecPrt('Hidden: Coor',' ',Work(ipCoor),3,mTot)
         End If
      End If
      nHidden = nKept
*
      Return
      End
