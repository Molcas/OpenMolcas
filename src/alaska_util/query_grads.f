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
* Copyright (C) 2015, Ignacio Fdez. Galvan                             *
************************************************************************
*  Query_Grads
*
*> @brief Query sizes from a gradients file
*> @author Ignacio Fdez. Galv&aacute;n
*>
*> @details
*> Reads the number of roots and length of the gradients from the
*> gradients file (GRADS).
*>
*> @param[out] Exists Whether or not the gradients file exists
*> @param[out] nRoots Number of roots allowed in the gradients file
*> @param[out] nGrad  Length of the vectors in the gradients file
************************************************************************
      Subroutine Query_Grads(Exists,nRoots,nGrad)
      Implicit None
      Logical :: Exists
      Integer :: nRoots,nGrad
      Integer, Dimension(5) :: TOC
      Integer, Dimension(1) :: iDum
      Integer :: LuGrad,iAd
      Logical :: Found
      Character(Len=5) :: Filename
*
      Filename='GRADS'
      Call f_Inquire(Filename,Found)
      If (.Not.Found) Then
        Exists=.False.
        nRoots=0
        nGrad=0
        Return
      End If
*
      LuGrad=20
      Call DaName(LuGrad,Filename)
      iAd=0
      Call iDaFile(LuGrad,2,TOC,Size(TOC),iAd)
      Call iDaFile(LuGrad,2,iDum,1,iAd)
      nRoots=iDum(1)
      Call iDaFile(LuGrad,2,iDum,1,iAd)
      nGrad=iDum(1)
      Call DaClos(LuGrad)
*
      End Subroutine Query_Grads
