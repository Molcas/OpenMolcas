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
      Subroutine PrDiOp(Text,nSym,nBas,XInt)
************************************************************************
*                                                                      *
*     Object: Print a diagonal block matrix                            *
*                                                                      *
************************************************************************
*
      Implicit Real*8 ( A-H,O-Z )
*
      Character*(*) Text
      Dimension XInt(*)
      Integer nBas(*)
*
*----------------------------------------------------------------------*
*                                                                      *
*     Start procedure                                                  *
*     Loop over pairs of symmetry labels. Skip all offdiagonal         *
*     blocks.                                                          *
*                                                                      *
*----------------------------------------------------------------------*
*
      lText=Min(120,LEN(Text))
      Write(6,'(6X,A)')Text(1:lText)
      iOff=0
      Do iSym=1,nSym
        nBs=nBas(iSym)
        If ( nBs.ne.0 ) Then
          Write(6,'(6X,A,I2)') 'Symmetry species',iSym
          Call TriPrt(' ',' ',XInt(iOff+1),nBs)
        End If
      End Do
*
*----------------------------------------------------------------------*
*     Terminate procedure                                              *
*----------------------------------------------------------------------*
*
      Return
      End
