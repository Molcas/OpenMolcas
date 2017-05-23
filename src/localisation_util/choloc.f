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
      SubRoutine ChoLoc(irc,Dens,CMO,Thrs,xNrm,nBas,nOcc)
************************************************************
*
*   <DOC>
*     <Name>ChoLoc</Name>
*     <Syntax>Call ChoLoc(irc,Dens,CMO,Thrs,xNrm,nBas,nOcc)</Syntax>
*     <Arguments>
*       \Argument{irc}{Return code}{Integer}{out}
*       \Argument{Dens}{Density matrix}{Real*8}{inout}
*       \Argument{CMO}{Cholesky MO coefficients}{Real*8}{out}
*       \Argument{Thrs}{Threshold for decomposition}{Real*8}{in}
*       \Argument{xNrm}{Total norm of Cholesky vectors}{Real*8}{out}
*       \Argument{nBas}{Number of basis functions}{Integer}{in}
*       \Argument{nOcc}{Number of occupied orbitals}{Integer}{in}
*     </Arguments>
*     <Purpose>Cholesky decompose density to obtain Cholesky localized
*              orbitals
*     </Purpose>
*     <Dependencies></Dependencies>
*     <Author>T.B. Pedersen</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*         Localize orbitals by Cholesky decomposition of the density
*         matrix. The threshold should be set such that the
*         decomposition is essentially exact (e.g. 1.0d-12). If not, you
*         might risk that the localization fails (non-zero return code),
*         since the number of Cholesky vectors will be less than the
*         number of occupied (nOcc) molecular orbitals. On sucessful
*         completion, irc=0 is returned. NOTE: the density matrix is
*         destroyed during the localization!
*     </Description>
*    </DOC>
*
************************************************************

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
