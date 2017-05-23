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
      Subroutine Wr_Cord(nAtoms)

      Implicit real*8 (a-h,o-z)

#include "WrkSpc.fh"
#include "MpParam.fh"
#include "Address.fh"
#include "MolProp.fh"

      iStdOut = 6

      Write (iStdOut,*)
      Write (iStdOut,'(10X,A)')
     &      ' ************************************************ '
      Write (iStdOut,'(10X,A)')
     &      ' **** Cartesian Coordinates / Bohr, Angstrom **** '
      Write (iStdOut,'(10X,A)')
     &      ' ************************************************ '
      Write (iStdOut,*)
      Write (iStdOut,'(A11,A7,A)')  'Center',  'Label',
     &      '              X              Y              Z',
     &      '                     X              Y              Z'
      Do i=1,nAtoms
*Jose    Write (iStdOut,'(6X,I3,5X,A4,3(5X,F10.6),7X,3(5X,F10.6))')
         Write (iStdOut,'(6X,I3,5X,A6,3(5X,F10.6),7X,3(5X,F10.6))')
     &   i,Labe(i),Cor(1,i,i),Cor(2,i,i),Cor(3,i,i),
     &   Cor(1,i,i)*Bohr_to_Ang,
     &   Cor(2,i,i)*Bohr_to_Ang,Cor(3,i,i)*Bohr_to_Ang
      EndDo
      Write (iStdOut,*)
      Write (iStdOut,*)
      Write (iStdOut,'(A11,A7,A13,A25,A20)') 'Center','Label',
     &'Atomtype','Effective core charge','Atomic Parameter'
      Write (iStdOut,*)
      Do i=1,nAtoms
*Jose Write (iStdOut,'(6X,I3,5X,A4,10X,I3,19X,F6.1,10X,I10)')i,Labe(i),
      Write (iStdOut,'(6X,I3,5X,A6,10X,I3,19X,F6.1,10X,I10)')i,Labe(i),
     &iAtomType(i),Work(iQnuc+i-1), iAtomPar(i)
      EndDo
      Write (iStdOut,*)
      Return
      End
