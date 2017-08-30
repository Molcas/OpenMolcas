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
* Copyright (C) Thomas Bondo Pedersen                                  *
************************************************************************
*  ChoLoc
*
*> @brief
*>   Cholesky decompose density to obtain Cholesky localized orbitals
*> @author T. B. Pedersen
*>
*> @details
*> Localize orbitals by Cholesky decomposition of the density
*> matrix. The threshold should be set such that the
*> decomposition is essentially exact (e.g. ``1.0d-12``). If not, you
*> might risk that the localization fails (non-zero return code),
*> since the number of Cholesky vectors will be less than the
*> number of occupied (\p nOcc) molecular orbitals. On sucessful
*> completion, \p irc = ``0`` is returned.
*>
*> @note
*> The density matrix is destroyed during the localization!
*>
*> @param[out]    irc  Return code
*> @param[in,out] Dens Density matrix
*> @param[out]    CMO  Cholesky MO coefficients
*> @param[in]     Thrs Threshold for decomposition
*> @param[out]    xNrm Total norm of Cholesky vectors
*> @param[in]     nBas Number of basis functions
*> @param[in]     nOcc Number of occupied orbitals
************************************************************************
      SubRoutine ChoLoc(irc,Dens,CMO,Thrs,xNrm,nBas,nOcc)
      Implicit None
      Integer irc, nBas, nOcc
      Real*8  Thrs, xNrm
      Real*8  Dens(nBas,nBas), CMO(nBas,nOcc)

      Character*6 SecNam
      Parameter (SecNam = 'ChoLoc')

      Integer  nVec
      real*8 ddot_
      external ddot_

      irc  = 0
      xNrm = -9.9d9

      nVec = 0
      Call CD_InCore(Dens,nBas,CMO,nOcc,nVec,Thrs,irc)
      If (irc .ne. 0) Then
         Write(6,*) SecNam,': CD_InCore returned ',irc
         Return
      Else If (nVec .ne. nOcc) Then
         Write(6,*) SecNam,': nVec.NE.nOcc'
         Write(6,*) '   nVec,nOcc = ',nVec,nOcc
         irc = 1
         Return
      End If

      xNrm = sqrt(dDot_(nBas*nOcc,CMO,1,CMO,1))

      End
