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
      SubRoutine GetFragment(lUnit,ipExp,MxShll,iShll,nFragType,
     &                       nFragCoor,nFragEner,nFragDens,ipFragType,
     &                       ipFragCoor,ipFragEner,ipFragCoef,
     &                       DInf,nDInf)
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
      Implicit None
#include "print.fh"
#include "angstr.fh"
      Integer nDInf,lUnit,MxShll,iShll
      Integer nFragType,nFragCoor,nFragEner,nFragDens
      Integer ipFragType,ipFragCoor,ipFragEner,ipFragCoef
      Real*8 DInf(nDInf)
      Character*180 Line, Get_Ln
      Integer ipExp(MxShll)
      integer storageSize,LineWords
      parameter(storageSize = 200, LineWords=storageSize/8)
!     LineWords = 25
      Character*(storageSize) sBasis
      Real*8 eqBasis(LineWords)
      Equivalence(sBasis,eqBasis)
      Integer iRout,iPrint,iStart,i,j,iBasis,ierr
*                                                                      *
************************************************************************
*                                                                      *
      iRout=30
      iPrint = nPrint(iRout)
      Call QEnter('GetFragment')
*
      iStart = ipExp(iShll+1)
*                                                                      *
************************************************************************
*                                                                      *
* Local basis sets for unique centers
*
      if(iPrint.ge.99) write(6,*) 'Reading LBASIS'
      Line = Get_Ln(lUnit)
      If (Index(Line,'LBASIS').eq.0) Then
        write (6,*) 'ERROR: Keyword LBASIS expected, offending line:'
        write (6,*) Line
        Call Quit_OnUserError()
      Endif

      Line = Get_Ln(lUnit)
      Call Get_i(1,nFragType,1)
      ipFragType = iStart
      if(iPrint.ge.99) write(6,*) 'number of LBASIS = ',nFragType

* read the basis sets
      do i = 1,nFragType
          sBasis=Get_Ln(lUnit)
          do j = 1,LineWords
             DInf(iStart+j-1) = eqBasis(j)
          enddo
          if(iPrint.ge.49) write(6,*) 'GetFragment: basis set ', sBasis
             iStart = iStart + LineWords
      enddo
*
* All atoms: index of the associated basis set and coordinates
*
      if(iPrint.ge.99) write(6,*) 'Reading RELCOORDS'

      Line = Get_Ln(lUnit)
      If (Index(Line,'RELCOORDS').eq.0) Then
        write (6,*) 'ERROR: Keyword RELCOORDS expected, offending line:'
        write (6,*) Line
        Call Quit_OnUserError()
      Endif

      Line = Get_Ln(lUnit)
      Call Get_i(1,nFragCoor,1)
      ipFragCoor = iStart
      if(iPrint.ge.99) write(6,*) 'number of RELCOORDS = ',nFragCoor
*
* read all centers, but reserve space for the Mulliken charges
*
      do i = 1,nFragCoor
        Line = Get_Ln(lUnit)
        Call Get_i(1,iBasis,1)
        DInf(iStart) = dble(iBasis)
        Call Get_f(2,DInf(iStart+1),3)
        If (Index(Line,'ANGSTROM').ne.0) Then
         If(iPrint.ge.49) write(6,*) 'Reading the relcoords in Angstrom'
          Do j = 1, 3
            DInf(iStart+j) = DInf(iStart+j)/Angstr
          End Do
        End If
        iStart = iStart + 5
      enddo
*
* Orbital energies (taken from the ONE ELECTRON ENERGIES in ScfOrb)
*
      if(iPrint.ge.99) write(6,*) 'Reading ENERGIES'
      Line = Get_Ln(lUnit)
      If (Index(Line,'ENERGIES').eq.0) Then
        write (6,*) 'ERROR: Keyword ENERGIES expected, offending line:'
        write (6,*) Line
        Call Quit_OnUserError()
      Endif

      Line = Get_Ln(lUnit)
      Call Get_i(1,nFragEner,1)
      ipFragEner = iStart

      Call Read_v(lUnit,DInf,iStart,iStart+nFragEner-1,1,ierr)
      if(iPrint.ge.99) Call RecPrt('Fragment MO energies',' ',
     &  DInf(ipFragEner),1,nFragEner)
      If (ierr.ne.0) Then
        write(6,*) 'ERROR: number of energy values is not correct'
        write(6,*) ierr
        Call Quit_OnUserError()
      Endif
      iStart = iStart + nFragEner
*
* MO coefficients (taken from the ORBITALs in ScfOrb)
*
      if(iPrint.ge.99) write(6,*) 'Reading MOCOEFF'
      Line = Get_Ln(lUnit)
      If(Index(Line,'MOCOEFF').eq.0) Then
        write(6,*) 'ERROR: Keyword MOCOEFF expected, offending line:'
        write(6,*) Line
        Call Quit_OnUserError()
      Endif

      Line=Get_Ln(lUnit)
      Call Get_i(1,nFragDens,1)
      ipFragCoef = iStart
      Call Read_v(lUnit,DInf,ipFragCoef,
     &            ipFragCoef+nFragDens*nFragEner-1,1,ierr)
      If(iPrint.ge.99) Call RecPrt('Fragment MO coefficients',' ',
     &  DInf(ipFragCoef),nFragDens,nFragEner)
      If(ierr.ne.0) Then
        write(6,*) 'ERROR: number of coefficients is not correct'
        Call Quit_OnUserError()
      Endif
      iStart = iStart + nFragDens*nFragEner
*
* Mulliken charges
*
      if(iPrint.ge.99) write(6,*) 'Reading MULLIKEN'
      Line = Get_Ln(lUnit)
      If (Index(Line,'MULLIKEN').eq.0) Then
        write(6,*) 'ERROR: Keyword MULLIKEN expected, offending line:'
        write(6,*) Line
        Call Quit_OnUserError()
      Endif
* Temporarily use the DInf array to store the charges, but they are moved
* together with the coordinates afterwards
      Call Read_v(lUnit,DInf,iStart,
     &            iStart+nFragCoor-1,1,ierr)
      If(iPrint.ge.99) Call RecPrt('Fragment Mulliken charges',' ',
     &  DInf(iStart),nFragCoor,1)
      If(ierr.ne.0) Then
        write(6,*) 'ERROR: number of Mulliken charges is not correct'
        Call Quit_OnUserError()
      Endif
      Do i = 1, nFragCoor
        DInf(ipFragCoor + 4 + (i-1)*5) = DInf(iStart + i - 1)
      End Do

      ipExp(iShll+1) = iStart
      Call QExit('GetFragment')

      Return
      End
