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
* Copyright (C) 2006, Roland Lindh                                     *
************************************************************************
      SubRoutine Print_Geometry(iOpt)
************************************************************************
*                                                                      *
* Object: to print the molecular coordinates, bonds, angles and        *
*         torsional angles.                                            *
*                                                                      *
*     Author: Roland Lindh, Dept Chem. Phys., Lund University, Sweden  *
*             September 2006                                           *
************************************************************************
      use Basis_Info
      use Center_Info
      use Period
      use Temporary_Parameters, only: Expert
      use Sizes_of_Seward, only: S
      use Real_Info, only: Rtrnc
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "print.fh"
      Character(LEN=1) help_c
      Character(LEN=16) FMT
      Character(LEN=LENIN), Allocatable:: Lblxxx(:)
      Real*8, Dimension (:,:), Allocatable :: Centr
#include "angstr.fh"
*                                                                      *
************************************************************************
*                                                                      *
      iRout=2
      iPrint = nPrint(iRout)
      If (iPrint.eq.0) Return
      LuWr=6
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_allocate(Centr,3,S%mCentr,Label='Centr')
      Call mma_allocate(Lblxxx,S%mCentr,Label='Lblxxx')
*                                                                      *
************************************************************************
*                                                                      *

      Write (LuWr,*)
      Call CollapseOutput(1,'   Molecular structure info:')
      Write (LuWr,'(3X,A)') '   -------------------------'
      Write (LuWr,*)
      if(iOpt.eq.0) FMT='(19X,A)'
      if(iOpt.eq.1) FMT='(11X,A)'
*
      Write (LuWr,FMT)
     &       ' ************************************************ '
      if(iOpt.eq.0) then
      Write (LuWr,FMT)
     &       ' **** Cartesian Coordinates / Bohr, Angstrom **** '
      else
      Write (LuWr,FMT)
     &       ' **** Cartesian Coordinates / Angstrom       **** '
            endif
      Write (LuWr,FMT)
     &       ' ************************************************ '
      Write (LuWr,*)
      if(iOpt.eq.0) then
      Write (LuWr,'(A,A,A)')  '     Center  Label ',
     &             '               x              y              z',
     &  '                     x              y              z'
      else
      Write (LuWr,'(A,A,A)')  '     Center  Label ',
     &             '               x              y              z'
      endif
*                                                                      *
************************************************************************
*                                                                      *
      ndc = 0
      nc = 1
      Do jCnttp = 1, nCnttp
         mCnt = dbsc(jCnttp)%nCntr
         If (dbsc(jCnttp)%Aux.or.dbsc(jCnttp)%Frag)Then
            ndc = ndc + mCnt
            Go To 32
         End If
         Do jCnt = 1, mCnt
            ndc = ndc + 1
            Do i = 0, nIrrep/dc(ndc)%nStab - 1
               Call OA(dc(ndc)%iCoSet(i,0),dbsc(jCnttp)%Coor(1:3,jCnt),
     &                 Centr(1:3,nc))
               If (Show) Then
                  help_c = ' '
                  If(Cell_l) Then
                     Do j=1,lthCell
                         If(AdCell(j).EQ.ndc) help_c = '*'
                     End Do
                  End If
                  if(iOpt.eq.0) then
                     Write (LuWr,
     &                   '(6X,I3,A1,5X,A,3F15.6,7X,3F15.6)')
     &                    nc, help_c, dc(ndc)%LblCnt,
     &                    Centr(1:3,nc),
     &                    Centr(1:3,nc)*angstr
                  else
                     Write (LuWr,
     &                   '(6X,I3,A1,5X,A,3F15.6)')
     &                    nc, help_c, dc(ndc)%LblCnt,
     &                    Centr(1:3,nc)*angstr
                  endif
               End If
               if (nc.gt.8*MxAtom) Then
                  Call WarningMessage(2,'lblxxx too small')
                  Call Abend()
               End If
               lblxxx(nc)=dc(ndc)%LblCnt(1:LENIN)
               nc = nc + 1
            End Do
         End Do
32       Continue
      End Do
      nc=nc-1
*                                                                      *
************************************************************************
*                                                                      *
*     Compute distances
*
      If (S%mCentr.le.2) Go To 55
      Call Dstncs(lblxxx,Centr,nc,angstr,S%Max_Center,6)
      If (.Not.Expert) Call DstChk(Centr,lblxxx,nc)
*
*     Compute valence bond angels
*
      If (iPrint.lt.5.or.S%mCentr.lt.3.or.iOpt.eq.1) Go To 55
      Call Angles(lblxxx,Centr,nc,rtrnc,S%Max_Center)
*
*     Compute dihedral angles
*
      If (iPrint.lt.5.or.S%mCentr.lt.4) Go To 55
      Call Dihedr(lblxxx,Centr,nc,rtrnc,S%Max_Center)
*                                                                      *
************************************************************************
*                                                                      *
 55   Continue
*
      Call mma_deallocate(Lblxxx)
      Call mma_deallocate(Centr)
*                                                                      *
************************************************************************
*                                                                      *
      Call CollapseOutput(0,'   Molecular structure info:')
      Write (LuWr,*)
      Return
      End
