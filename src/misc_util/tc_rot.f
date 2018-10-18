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
      Subroutine Rotation(TotalM,TRotA,TRotB,TRotC,
     &                       nsRot,nFAtoms,lSlapaf)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "Molcas.fh"
#include "constants2.fh"
      Real*8 dEV, dVec(3)
#include "WrkSpc.fh"
      Real*8  CM(3)                         ! Center of masses
      Real*8  FCoor(3,mxAtom)               ! Full Coords
      Real*8  CCoor(3,mxAtom)               ! Mass-centered Coords
      Real*8  SOCoor(3,mxAtom)              ! Symmetry-Oriented Coords
      Real*8  Inrt(3,3), RotE(3), Vec(3,3)  ! Inertia components
      Real*8  Mass(mxAtom)                  ! Masses
      Character*(LENIN) FAtLbl(mxAtom)      ! Atomic labes
      Integer nAtomicN(mxAtom)              ! Atomic numbers
      Logical lSlapaf

*
* --- Get Atomic Full Labels, Coordinates & Mass - FAtLbl, FCoor, Mass
*
      Call GetFullCoord(FCoor,Mass,FAtLbl,nFAtoms,mxAtm,lSlapaf)
*
* --- Define the Center of Masses - CM()
*
      CM(1)  = 0.0d0
      CM(2)  = 0.0d0
      CM(3)  = 0.0d0
      TotalM = 0.0d0
      Do i=1, nFAtoms
        TotalM = TotalM + Mass(i)
        CM(1)  = CM(1) + Mass(i)*FCoor(1,i)
        CM(2)  = CM(2) + Mass(i)*FCoor(2,i)
        CM(3)  = CM(3) + Mass(i)*FCoor(3,i)
      EndDo
      CM(1)  = CM(1)/TotalM
      CM(2)  = CM(2)/TotalM
      CM(3)  = CM(3)/TotalM
*
* --- Shift coordinates in the center of masses - CCoord()
*
      Do i=1, nFAtoms
        CCoor(1,i) = FCoor(1,i) - CM(1)
        CCoor(2,i) = FCoor(2,i) - CM(2)
        CCoor(3,i) = FCoor(3,i) - CM(3)
      EndDo
*
* --- Compute the Inertia-matrix - Inrt(3,3)
*
      Do i = 1, 3
        Do j = 1, 3
          Inrt(i,j) = 0.0d0
        EndDo
      EndDo
      Do i = 1, nFAtoms
        dX = CCoor(1,i)
        dY = CCoor(2,i)
        dZ = CCoor(3,i)
        Inrt(1,1) = Inrt(1,1) + Mass(i) * (dY*dY+dZ*dZ) ! YY ZZ
        Inrt(2,1) = Inrt(2,1) - Mass(i) * dX * dY       ! XY
        Inrt(2,2) = Inrt(2,2) + Mass(i) * (dX*dX+dZ*dZ) ! XX ZZ
        Inrt(3,1) = Inrt(3,1) - Mass(i) * dX * dZ       ! XZ
        Inrt(3,2) = Inrt(3,2) - Mass(i) * dY * dZ       ! YZ
        Inrt(3,3) = Inrt(3,3) + Mass(i) * (dX*dX+dY*dY) ! XX YY
      EndDo
      Inrt(1,2) = Inrt(2,1)
      Inrt(1,3) = Inrt(3,1)
      Inrt(2,3) = Inrt(3,2)
*
* --- and diagonalize it - Inrt(3,3)
*
      Call GetMem('EVal','Allo','Real',ipEVal,3*(3+1)/2)
      Call GetMem('EVec','Allo','Real',ipEVec,3*3)
      Do i = 1, 3
        Do j = 1, 3
          ij=i*(i-1)/2 + j + ipEval -1
          Work(ij)=Inrt(i,j)
        End Do
      End Do
      call dcopy_(3*3,Zero,0,Work(ipEVec),1)
      call dcopy_(3,One,0,Work(ipEVec),3+1)
      Call Jacob (Work(ipEVal),Work(ipEVec),3,3)
      Call Jacord(Work(ipEVal),Work(ipEVec),3,3)
      Do i = 1, 3
        RotE(i)=Work(i*(i+1)/2+ipEVal-1)
        Do j = 1, 3
           Vec(i,j) = Work(ipEVec-1+i+(j-1)*3)
        EndDO
      End Do
*
*     Sort the principal axis such that z' is the one with the lowest
*     eigenvalue.
*
      Do i = 1, 2
        Do j = i+1, 3
          If (RotE(i).LT.RotE(j)) Then
            dEV     = RotE(i)
            dVec(1) = Vec(1,i)
            dVec(2) = Vec(2,i)
            dVec(3) = Vec(3,i)
              RotE(i)  = RotE(j)
              Vec(1,i) = Vec(1,j)
              Vec(2,i) = Vec(2,j)
              Vec(3,i) = Vec(3,j)
            RotE(j)  = dEV
            Vec(1,j) = dVec(1)
            Vec(2,j) = dVec(2)
            Vec(3,j) = dVec(3)
          End If
        EndDo
      EndDo
*
*     Rotate coords to Symmetry-Oriented
*
      Do iAtom=1, nFAtoms
        Do i = 1, 3
          dSum = 0.0d0
          Do j = 1, 3
            dSum = dSum + CCoor(j,iAtom) * Vec(j,i)
          EndDo
          SOCoor(i,iAtom) = dSum
        EndDo
      EndDo
*
* --- Rotational Symmetry factor - nsRot
*
      If (nsRot.EQ.0) nsRot = 1
      If (nFAtoms.EQ.2) then
        If (Mass(1).EQ.Mass(2)) nsRot = 2
      EndIf
*
      TRotA = 8.661377d01/(RotE(3)+1.0d-99)
      TRotB = 8.661377d01/(RotE(2)+1.0d-99)
      TRotC = 8.661377d01/(RotE(1)+1.0d-99)
*
**    Check if linear molecule
*
      nrot=3
      if (TRotA.gt.1.0d99) nrot=nrot-1
      if (TRotB.gt.1.0d99) nrot=nrot-1
      if (TRotC.gt.1.0d99) nrot=nrot-1
*
* --- Print results
*
      Write(6,'(A)') ' Mass-centered Coordinates (Angstrom):'
      Write(6,'(1X,60A)')
     &    '***********************************************************'
      Write(6,'(1X,60A)')
     &    'Label   N         X           Y           Z          Mass  '
      Write(6,'(1X,60A)')
     &    '-----------------------------------------------------------'
      Do i=1,nFAtoms
         Write(6,'(1X,A,1X,I4,1X,3F12.6,1x,F12.5)')
     & FAtLbl(i),nAtomicN(i),(Angstrom*SOCoor(j,i),j=1,3),Mass(i)
      EndDo
      Write(6,'(1X,60A)')
     &    '-----------------------------------------------------------'
      Write(6,'(A,F12.6)') ' Molecular mass:',TotalM
      Write(6,'(A,3F10.4)') ' Rotational Constants (cm-1):',
     &              (auTocm*Half/(uToau*RotE(i)),i=1,nrot)
      Write(6,'(A,3F10.4)') ' Rotational Constants (GHz) :',
     &              (1.0D-9*auToHz*Half/(uToau*RotE(i)),i=1,nrot)
      Write(6,'(A,3F10.4)') ' Rotational temperatures (K):',
     &                          (8.661377d01/RotE(i),i=1,nrot)
      Write(6,'(A,I2)') ' Rotational Symmetry factor: ',nsRot
*
      Call GetMem('EVec','Free','Real',ipEVec,3*3)
      Call GetMem('EVal','Free','Real',ipEVal,3*(3+1)/2)
*
      Return
      End
