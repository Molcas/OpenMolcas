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
* Copyright (C) Per-Olof Widmark                                       *
************************************************************************
************************************************************************
*                                                                      *
* This routine perform some initialization after reading the first     *
* one electron file.                                                   *
* These parameters are checked when reading subsequent one electron    *
* files.                                                               *
*                                                                      *
*======================================================================*
*                                                                      *
* Author: Per-Olof Widmark                                             *
*         IBM Sweden                                                   *
*                                                                      *
************************************************************************
      Subroutine Init_GenANO
      Implicit Real*8 (a-h,o-z)
#include "parm.fh"
#include "common.fh"
#include "symlab.fh"
      Character*(LENIN4) LblCnt(MxAtom)
      Logical Found

      Call Get_nAtoms_all(nAtoms)
      Call Get_cArray('LP_L',LblCnt,(LENIN4)*nAtoms)
      Found=.FALSE.
      Do i=1,nAtoms
        If (LblCnt(i)(1:LENIN).eq.Center) Found=.TRUE.
      End Do
      If (.not.Found) Then
        Call WarningMessage(2,'Center '//Center//' not found')
        Call Quit_OnUserError()
      End IF
      ind=0
      Do 100 iSym=1,nSym
*        Write(*,'(a,i1,a)') ' <<< Symmetry ',iSym,' >>>'
         Do 110 iBas=1,nBas(iSym)
            ind=ind+1
*           Write(*,*) Name(ind)
            Do 111 l=0,MxLqn
               If(Name(ind)(1:LENIN).eq.Center) Then
                  If(Name(ind)(LENIN1:LENIN4).eq.type(l*(l+1)+1)) Then
                     nPrim(l)=nPrim(l)+1
                  End If
               End If
111         Continue
110      Continue
100   Continue
      Write(6,*)
      Write(6,'(a,8i5)') 'Number of primitives per shell:',nPrim
      nDsym=0
      Do 200 l=0,MxLqn
*        Write(*,'(1x,a,i5)') type(l*(l+1)+1),nPrim(l)
         nDsym=nDsym+(2*l+1)*nPrim(l)*(nPrim(l)+1)/2
200   Continue
      Do 300 i=1,nDsym
         pDsym(i)=0.0d0
         tDsym(i)=0.0d0
300   Continue
*----------------------------------------------------------------------*
      ind=0
      len=0
      next=1
      Do 400 iShell=0,MxLqn
         len=nPrim(iShell)*(nPrim(iShell)+1)/2
         Do 410 iComp=-iShell,iShell
            ind=ind+1
            iSymBk(ind)=next
            next=next+len
410      Continue
400   Continue
*     Write(*,*) 'In Init: density block pointers'
*     Write(*,'(1x,10i5)') (iSymBk(i),i=1,(MxLqn+1)*(MxLqn+1))
      Return
      End
