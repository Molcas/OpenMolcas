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
      Subroutine PlaceIt(Coord,iQ_Atoms,iCNum)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"

      Dimension Coord(MxAt*3)
      Dimension AvstPart(MxPut),IndexSet(MxPut)
      Dimension CordstTemp(MxPut*MxCen,3)
      Character Head*200
      Logical Changed


      Do 11, i=1,nPart !For each solvent particle, compute the
        Sbig=1D+20   !smallest distance to any QM-atom from the
        Do 12, j=1,iQ_Atoms  !oxygen on water.
          S=0
          Do 13, k=1,3
            S=S+(Coord((j-1)*3+k)-Cordst(nCent*(i-1)+1,k))**2
13        Continue
          If(S.le.Sbig) then
            Sbig=S
            AvstPart(i)=S
          Endif
12      Continue
11    Continue

      Do 21, i=1,MxPut
        IndexSet(i)=i
21    Continue

31    Continue  !Order the indeces suchwise that smallest distance
        Changed=.false.  !goes first. The sorting routine is blunt
        Do 32, i=1,nPart-1 !but at this stage of the execution time
          If(AvstPart(i+1).lt.AvstPart(i)) then !is not a problem.
            Atemp=AvstPart(i)
            AvstPart(i)=AvstPart(i+1)
            AvstPart(i+1)=Atemp
            iTemp=IndexSet(i)
            IndexSet(i)=IndexSet(i+1)
            IndexSet(i+1)=iTemp
            Changed=.true.
          Endif
32      Continue
      If(Changed) GoTo 31

      Do 41, i=1,nPart  !Put coordinates of solvent suchwise that
        Do 42, j=1,nCent !smallest distances goes first.
          CordstTemp((i-1)*nCent+j,1)=Cordst((i-1)*nCent+j,1)
          CordstTemp((i-1)*nCent+j,2)=Cordst((i-1)*nCent+j,2)
          CordstTemp((i-1)*nCent+j,3)=Cordst((i-1)*nCent+j,3)
42      Continue
41    Continue
      Do 43, i=1,nPart
        ind=IndexSet(i)
        Do 44, j=1,nCent
          Cordst((i-1)*nCent+j,1)=CordstTemp((ind-1)*nCent+j,1)
          Cordst((i-1)*nCent+j,2)=CordstTemp((ind-1)*nCent+j,2)
          Cordst((i-1)*nCent+j,3)=CordstTemp((ind-1)*nCent+j,3)
44      Continue
43    Continue

      Do 51,iz=1,iQ_Atoms  !Substitute the first coordinate slots with
        Cordst(iz,1)=Coord((iz-1)*3+1) !QM-molecule, or since we have
        Cordst(iz,2)=Coord((iz-1)*3+2) !ordered above, this is
        Cordst(iz,3)=Coord((iz-1)*3+3) !equivalent with removing closest
51    Continue                         !solvents and there put QM-mol.
      Do 52, iextr=iQ_Atoms+1,iCnum*nCent  !Just dummy-coordinates
        Cordst(iextr,1)=Coord(1)        !added to empty slots.
        Cordst(iextr,2)=Coord(2)
        Cordst(iextr,3)=Coord(3)
52    Continue

      If(iPrint.ge.10) then  !Optional printing.
        Write(Head,*)'Coordinates of the system after substitution and'
     &//' reordening of solvent molecules.'
        Call Cooout(Head,Cordst,nPart,nCent)
      Endif

      Return
      End

*----------------------------------------------------------------------*
* With this function we wish to place the QM-molecule properly when we *
* run with solvetn configurations from the sampfile. This we do by     *
* making the center-of-masses to coincide.                             *
*----------------------------------------------------------------------*
      Subroutine PlaceIt9(Coord,Cordst,info_atom,iQ_Atoms)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"

      Dimension Coord(MxAt*3),Cordst(MxCen*MxPut,3)
      Dimension info_atom(MxAt)

      CMSewx=0
      CMSewy=0
      CMSewz=0
      CMSamx=0
      CMSamy=0
      CMSamz=0
      Wtot=0
      Do 1, i=1,iQ_Atoms
        CMSewx=CMSewx+Coord((i-1)*3+1)*info_atom(i)
        CMSewy=CMSewy+Coord((i-1)*3+2)*info_atom(i)
        CMSewz=CMSewz+Coord((i-1)*3+3)*info_atom(i)
        CMSamx=CMSamx+Cordst(i,1)*info_atom(i)
        CMSamy=CMSamy+Cordst(i,2)*info_atom(i)
        CMSamz=CMSamz+Cordst(i,3)*info_atom(i)
        Wtot=Wtot+dble(info_atom(i))
1     Continue
      CMSewx=CMSewx/Wtot
      CMSewy=CMSewy/Wtot
      CMSewz=CMSewz/Wtot
      CMSamx=CMSamx/Wtot
      CMSamy=CMSamy/Wtot
      CMSamz=CMSamz/Wtot
      Tx=CMSewx-CMSamx
      Ty=CMSewy-CMSamy
      Tz=CMSewz-CMSamz
      Do 2, i=1,iQ_Atoms
        Cordst(i,1)=Coord((i-1)*3+1)-Tx
        Cordst(i,2)=Coord((i-1)*3+2)-Ty
        Cordst(i,3)=Coord((i-1)*3+3)-Tz
2     Continue

      Return
      End
