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
* Copyright (C) 2010, Mickael G. Delcey                                *
*               2017, Ignacio Fdez. Galvan                             *
************************************************************************
      subroutine processRP(KeepGroup,SymThr,DInf,nDInf)
#ifndef _HAVE_EXTRA_
      Use XYZ
#endif
      Implicit Real*8 (a-h,o-z)
      Character KeepGroup*180,KWord*180,KG*180
      Real*8 DInf(nDInf)
#ifdef _HAVE_EXTRA_
      Character*180 minGroup,Key
      Character*180, External :: Get_Ln
#else
      Real*8, Dimension(:,:,:), Allocatable :: DumRot
      Real*8, Dimension(:,:), Allocatable :: DumTrans
#endif
#include "itmax.fh"
#include "info.fh"
************************************************************************
*                                                                      *
*    A silly routine to try to handle symmetry in RP-Coord section     *
*    verify that R, P and the standard structures have one common      *
*    symmetry. The user should preferably use a defined symmetry if    *
*    the TS symmetry is lower than R and P                             *
*                                                                      *
*    M.G. Delcey     June 2010                                         *
*    Lund University                                                   *
*                                                                      *
*    Adaptation to OpenMolcas                                          *
*    I. Fdez. Galvan   July 2017                                       *
*                                                                      *
************************************************************************
*
**   If C1, all is already done
*
      KG=KeepGroup
      Call UpCase(KG)
      if (KG(1:1).eq.'E'.or.KG(1:2).eq.'C1') then
           KG='NOSYM'
      endif
      If (KG(1:5).eq.'NOSYM') Go To 30
*
#ifdef _HAVE_EXTRA_
      LuRP=IsFreeUnit(10)
      Call molcas_open(LuRP,'findsym.out')
      KWord = Get_Ln(LuRP)
      minGroup = Get_Ln(LuRP)
      close(LuRP)
      If (index(KWord,'C1').ne.0) Go To 30
*
**    Else findsym will find the group
*
      Call findsymf('findsym.RP1',KG,SymThr,ierr)
      If (ierr.ne.0) Go To 20
      Call molcas_open(LuRP,'findsym.out')
      Key = Get_Ln(LuRP)
      If (index(Key,'C1').ne.0) Then
        Close(LuRP)
        Go To 30
      End If
      Key = Get_Ln(LuRP)
      KWord = Get_Ln(LuRP)
      Call Get_I(1,nRP,1)
*  Only read the second structure in findsym.out
      Do i=1,nRP+3
        Read(LuRP,*)
      End Do
      Do i=1,nRP
         Read(LuRP,*) KWord,j,(DInf(ipRP1+3*(i-1)+j),j=0,2)
      EndDo
      nRP=nRP*3
      close(LuRP)
*
**    Compare with the one of the current structure
**    no need if user defined
*
      Write(6,*) Key,minGroup
      If (KG(1:4).eq.'FULL') Then
        If (Key.ne.minGroup) Go To 21
      EndIf
*
**
*
      Call findsymf('findsym.RP2',KG,SymThr,ierr)
      If (ierr.ne.0) Go To 20
      Call molcas_open(LuRP,'findsym.out')
      Read(LuRP,*)
      KWord = Get_Ln(LuRP)
*
      If (KG(1:4).eq.'FULL') Then
        If (KWord.ne.minGroup) Then
          Close(LuRP)
          Go To 21
        End If
      EndIf
*
      KWord = Get_Ln(LuRP)
      Call Get_I(1,i,1)
      If (3*i.ne.nRP) Go To 20
      Do i=1,nRP/3+3
        Read(LuRP,*)
      End Do
      Do i=1,nRP/3
         Read(LuRP,*) KWord,j,(DInf(ipRP1+nRP+3*(i-1)+j),j=0,2)
      EndDo
#else
      ! If manual symmetry, trust it
      ! Only check if automatic symmetry requested
      If (KG(1:4) .ne. 'FULL') Go To 30
      ! The detected symmetry for the two structures must match
      LuRP = 10
      LuRP = IsFreeUnit(LuRP)
      Call Molcas_Open(LuRP,'findsym.RP1')
      Call Read_XYZ(LuRP,DumRot,DumTrans)
      Close(LuRP)
      Call Parse_Group(KeepGroup,SymThr)
      nRP = Out_Raw(DInf(ipRP1))
      Call Clear_XYZ()
      !
      KG = Trim(Symmetry)
      LuRP = IsFreeUnit(LuRP)
      Call Molcas_Open(LuRP,'findsym.RP2')
      Call Read_XYZ(LuRP,DumRot,DumTrans)
      Close(LuRP)
      Call Parse_Group(KeepGroup,SymThr)
      i = Out_Raw(DInf(ipRP1+nRP))
      If (i .ne. nRP) Go To 20
      Call Clear_XYZ()
      If (Symmetry .ne. KG) Go To 21
#endif
 30   Continue
      return
************************************************************************
 20   Continue
      Call WarningMessage(2,
     &       'Error in RP-Coord section, check symmetry')
      Call Quit_OnUserError()
 21   Continue
      KWord='Error in RP-Coord section, structures do not have the '//
     &      'same symmetry. Please define manually the symmetry group.'
      Call WarningMessage(2,KWord)
      Call Quit_OnUserError()
      end
