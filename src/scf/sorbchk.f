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
      SubRoutine SOrbChk(OneHam,Ovrlp,Fock,mBT,nD,CMO,mBB)
      Implicit Real*8 (a-h,o-z)
*
      Real*8 OneHam(mBT), Ovrlp(mBT),Fock(mBT,nD), CMO(mBB,nD)
*
      Do iD = 1, nD
*----    Check orthonormality of start orbitals
         Call ChkOrt(CMO(1,iD),mBB,Ovrlp,nBT,Whatever)
*
*----    Form the first Fock matrix
         Call DCopy_(mBT,OneHam,1,Fock(1,iD),1)
      End Do
*
      Return
      End
