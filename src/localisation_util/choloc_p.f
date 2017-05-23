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
* Copyright (C) Francesco Aquilante                                    *
************************************************************************
      SubRoutine ChoLoc_p(irc,Dens,CMO,Thrs,xNrm,nBas,nOcc,iD)
******************************************************************
*
*   <DOC>
*     <Name>ChoLoc\_p</Name>
*     <Syntax>Call ChoLoc\_p(irc,Dens,CMO,Thrs,xNrm,nBas,nOcc)</Syntax>
*     <Arguments>
*       \Argument{irc}{Return code}{Integer}{out}
*       \Argument{Dens}{Density matrix}{Array Real*8}{inout}
*       \Argument{CMO}{Cholesky MO coefficients}{Real*8}{out}
*       \Argument{Thrs}{Threshold for decomposition}{Real*8}{in}
*       \Argument{xNrm}{Total norm of Cholesky vectors}{Real*8}{out}
*       \Argument{nBas}{Number of basis functions}{Integer}{in}
*       \Argument{nOcc}{Number of occupied orbitals}{Integer}{in}
*       \Argument{iD}{Index array of parent diagonals}{Integer array}{inout}
*     </Arguments>
*     <Purpose>Cholesky decompose density to obtain Cholesky localized
*              orbitals
*     </Purpose>
*     <Dependencies></Dependencies>
*     <Author>F. Aquilante</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*         Localize orbitals by Cholesky decomposition of the density
*         matrix and returns an index array of the parent diagonals.
*         The threshold should be set such that the
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
      Integer irc, nBas, nOcc, iD(nBas)
      Real*8  Thrs, xNrm
      Real*8  Dens(nBas,nBas), CMO(nBas,nOcc)

      Character*8 SecNam
      Parameter (SecNam = 'ChoLoc_p')

      Integer  nVec
      real*8 ddot_
      external ddot_

      irc  = 0
      xNrm = -9.9d9

      nVec = 0
      Call CD_InCore_p(Dens,nBas,CMO,nOcc,iD,nVec,Thrs,irc)
      If (irc .ne. 0) Then
         Write(6,*) SecNam,': CD_InCore_p returned ',irc
         Return
      Else If (nVec .ne. nOcc) Then
         Write(6,*) SecNam,': nVec.NE.nOcc'
         Write(6,*) '   nVec,nOcc = ',nVec,nOcc
         irc = 1
         Return
      End If

      xNrm = sqrt(dDot_(nBas*nOcc,CMO,1,CMO,1))

      End
*
*
      SubRoutine ChoLoc_xp(irc,Dens,CMO,Thrs,xNrm,nBas,nOcc,iD)
******************************************************************
*
*  Same as ChoLoc_p but handles differently the irc=102
*
************************************************************

      Implicit None
      Integer irc, nBas, nOcc, iD(nBas)
      Real*8  Thrs, xNrm
      Real*8  Dens(nBas,nBas), CMO(nBas,nOcc)

      Character*9 SecNam
      Parameter (SecNam = 'ChoLoc_xp')

      Integer  nVec
      real*8 ddot_
      external ddot_

      irc  = 0
      xNrm = -9.9d9

      nVec = 0
      Call CD_InCore_p(Dens,nBas,CMO,nOcc,iD,nVec,Thrs,irc)
      If (irc.ne.0 .and. irc.ne.102) Then
         Write(6,*) SecNam,': CD_InCore_p returned ',irc
         Return
      Else If (irc.eq.102) Then
         irc=0  ! reset because it is most likely a numerical noise
      Else If (nVec .ne. nOcc) Then
         Write(6,*) SecNam,': nVec.NE.nOcc'
         Write(6,*) '   nVec,nOcc = ',nVec,nOcc
         irc = 1
         Return
      End If

      xNrm = sqrt(dDot_(nBas*nOcc,CMO,1,CMO,1))

      End
