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
* Copyright (C) Ben Swerts                                             *
*               2016, Liviu Ungur                                      *
************************************************************************
      SubRoutine GetFragment(lUnit,nCnttp)
************************************************************************
*                                                                      *
*    Objective: To read frozen fragment information.                   *
*               This means that we read (from input stream) the unique *
*               centers with their basis sets, the coordinates of all  *
*               atoms of the fragment, the orbital energies, the       *
*               orbital coefficients and the Mulliken charges.         *
*                                                                      *
* Called from: Input                                                   *
*                                                                      *
*     Author: Ben Swerts                                               *
*   Modified: Liviu Ungur                                              *
************************************************************************
      Use Basis_Info
      Implicit None
#include "stdalloc.fh"
#include "angstr.fh"
      Integer lUnit,nCnttp
      Character*180 Line, Get_Ln
      integer storageSize,LineWords
      parameter(storageSize = 200, LineWords=storageSize/8)
      Integer nFragType,nFragCoor,nFragEner,nFragDens
!     LineWords = 25
      Character*(storageSize) sBasis
      Real*8 eqBasis(LineWords)
      Equivalence(sBasis,eqBasis)
      Integer iPrint,i,j,iBasis,ierr
*                                                                      *
************************************************************************
*                                                                      *
      Interface
      subroutine Read_v(lunit,Work,istrt,iend,inc,ierr)
      Integer lUnit, iStrt, iEnd, Inc, iErr
      Real*8 Work(iend)
      End subroutine Read_v
      End Interface
*                                                                      *
************************************************************************
*                                                                      *
!#define _DEBUG_
#ifdef _DEBUG_
      iPrint=99
#else
      iPrint=5
#endif
*                                                                      *
************************************************************************
*                                                                      *
      nFrag_LineWords=LineWords ! Needed in Module Basis_Info
*                                                                      *
************************************************************************
*                                                                      *
*     Keyword: LBASIS
*                                                                      *
************************************************************************
*                                                                      *
*     Local basis sets for unique centers
*
      if(iPrint.ge.99) write(6,*) 'Reading LBASIS'
      Line = Get_Ln(lUnit)
      If (Index(Line,'LBASIS').eq.0) Then
        write (6,*) 'ERROR: Keyword LBASIS expected, offending line:'
        write (6,*) Line
        Call Quit_OnUserError()
      Endif
*
*     Fragment types
*
      Line = Get_Ln(lUnit)
      Call Get_i1(1,nFragType)
      dbsc(nCnttp)%nFragType=nFragType
      Call mma_allocate(dbsc(nCnttp)%FragType,LineWords,nFragType,
     &                  Label='FragType')
      if(iPrint.ge.99) write(6,*) 'number of LBASIS = ',nFragType
*
*     read the basis sets labels
*
      do i = 1,nFragType
          sBasis=Get_Ln(lUnit)
          do j = 1,LineWords
             dbsc(nCnttp)%FragType(j,i) = eqBasis(j)
          enddo
          if(iPrint.ge.49) write(6,*) 'GetFragment: basis set ', sBasis
      enddo
*                                                                      *
************************************************************************
*                                                                      *
*     Keyword: RELCOORDS
*                                                                      *
************************************************************************
*                                                                      *
*     All atoms: index of the associated basis set and coordinates
*
      If(iPrint.ge.99) write(6,*) 'Reading RELCOORDS'

      Line = Get_Ln(lUnit)
      If (Index(Line,'RELCOORDS').eq.0) Then
        write (6,*) 'ERROR: Keyword RELCOORDS expected, offending line:'
        write (6,*) Line
        Call Quit_OnUserError()
      Endif

      Line = Get_Ln(lUnit)
      Call Get_i1(1,nFragCoor)
      dbsc(nCnttp)%nFragCoor=nFragCoor
      if(iPrint.ge.99) write(6,*) 'number of RELCOORDS = ',nFragCoor
*
*     read all centers, but reserve space for the Mulliken charges
*
*     FragCoor(1,i): Index of the FragType
*     FragCoor(2,i): x coordinate
*     FragCoor(3,i): y coordinate
*     FragCoor(4,i): z coordinate
*     FragCoor(5,i): Mulliken charge
*
      Call mma_allocate(dbsc(nCnttp)%FragCoor,5,nFragCoor,
     &                  Label='FragCoor')
      do i = 1,nFragCoor
        Line = Get_Ln(lUnit)
        Call Get_i1(1,iBasis)
        dbsc(nCnttp)%FragCoor(1,i) = dble(iBasis)
        Call Get_f(2,dbsc(nCnttp)%FragCoor(2,i),3)
        If (Index(Line,'ANGSTROM').ne.0) Then
           If (iPrint.ge.49)
     &        write(6,*) 'Reading the relcoords in Angstrom'
              dbsc(nCnttp)%FragCoor(2:4,i)=
     &            dbsc(nCnttp)%FragCoor(2:4,i)/Angstr
        End If
      enddo
*                                                                      *
************************************************************************
*                                                                      *
*     keyword: ENERGIES
*                                                                      *
************************************************************************
*                                                                      *
*     Orbital energies (taken from the ONE ELECTRON ENERGIES in ScfOrb)
*
      if(iPrint.ge.99) write(6,*) 'Reading ENERGIES'
      Line = Get_Ln(lUnit)
      If (Index(Line,'ENERGIES').eq.0) Then
        write (6,*) 'ERROR: Keyword ENERGIES expected, offending line:'
        write (6,*) Line
        Call Quit_OnUserError()
      Endif

      Line = Get_Ln(lUnit)
      Call Get_i1(1,nFragEner)
      dbsc(nCnttp)%nFragEner=nFragEner
      Call mma_Allocate(dbsc(nCnttp)%FragEner,nFragEner,
     &                   Label='FragEner')
*
      Call Read_v(lUnit,dbsc(nCnttp)%FragEner,1,nFragEner,1,ierr)
      if(iPrint.ge.99) Call RecPrt('Fragment MO energies',' ',
     &                              dbsc(nCnttp)%FragEner,1,nFragEner)
      If (ierr.ne.0) Then
        write(6,*) 'ERROR: number of energy values is not correct'
        write(6,*) ierr
        Call Quit_OnUserError()
      Endif
*                                                                      *
************************************************************************
*                                                                      *
*     keyword: MOCOEFF
*                                                                      *
************************************************************************
*                                                                      *
*     MO coefficients (taken from the ORBITALs in ScfOrb)
*
      if(iPrint.ge.99) write(6,*) 'Reading MOCOEFF'
      Line = Get_Ln(lUnit)
      If(Index(Line,'MOCOEFF').eq.0) Then
        write(6,*) 'ERROR: Keyword MOCOEFF expected, offending line:'
        write(6,*) Line
        Call Quit_OnUserError()
      Endif

      Line=Get_Ln(lUnit)
      Call Get_i1(1,nFragDens)
      dbsc(nCnttp)%nFragDens=nFragDens
      Call mma_allocate(dbsc(nCnttp)%FragCoef,nFragDens,nFragEner,
     &                  Label='FragCoef')
*
      Call Read_v(lUnit,dbsc(nCnttp)%FragCoef,1,nFragDens*nFragEner,1,
     &                       iErr)
      If(iPrint.ge.99) Call RecPrt('Fragment MO coefficients',' ',
     &                      dbsc(nCnttp)%FragCoef,nFragDens,nFragEner)
      If(ierr.ne.0) Then
        write(6,*) 'ERROR: number of coefficients is not correct'
        Call Quit_OnUserError()
      Endif
*                                                                      *
************************************************************************
*                                                                      *
*     keyword: MULLIKEN
*                                                                      *
************************************************************************
*                                                                      *
*     Mulliken charges
*
      if(iPrint.ge.99) write(6,*) 'Reading MULLIKEN'
      Line = Get_Ln(lUnit)
      If (Index(Line,'MULLIKEN').eq.0) Then
        write(6,*) 'ERROR: Keyword MULLIKEN expected, offending line:'
        write(6,*) Line
        Call Quit_OnUserError()
      Endif
      Call Read_v(lUnit,dbsc(nCnttp)%FragCoor,5,5*nFragCoor,5,ierr)
      If(iPrint.ge.99) Call RecPrt('Fragment Mulliken charges',' ',
     &                      dbsc(nCnttp)%FragCoor,5,nFragCoor)
      If(ierr.ne.0) Then
        write(6,*) 'ERROR: number of Mulliken charges is not correct'
        Call Quit_OnUserError()
      Endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
