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
      SubRoutine ERLoc(irc,CMO,Thr,ThrGrad,ThrRot,MxIter,nBas,nOcc,nFro,
     &                 nSym,Silent)
************************************************************
*
*   <DOC>
*     <Name>ERLoc</Name>
*     <Syntax>Call ERLoc(irc,CMO,Thr,ThrGrad,ThrRot,MxIter,nBas,nOcc,
*                        nFro,nSym,Silent)
*     </Syntax>
*     <Arguments>
*       \Argument{irc}{Return code}{Integer}{out}
*       \Argument{CMO}{Molecular orbital coefficients}{Real*8}{inout}
*       \Argument{Thr}{Threshold for functional}{Real*8}{in}
*       \Argument{ThrGrad}{Threshold for gradient}{Real*8}{in}
*       \Argument{ThrRot}{Screening threshold}{Real*8}{in}
*       \Argument{MxIter}{Max. number of iterations}{Integer}{in}
*       \Argument{nBas}{Number of basis functions per irrep}{Integer}
*                {in}
*       \Argument{nOcc}{Number of occupied orbitals to localize}
*                {Integer}{in}
*       \Argument{nFro}{Number of frozen occupied orbitals}
*                {Integer}{in}
*       \Argument{nSym}{Number of irreps}{Integer}{in}
*       \Argument{Silent}{Flag to avoid printing}{Logical}{in}
*     </Arguments>
*     <Purpose>Edmiston-Ruedenberg localization of occupied orbitals
*     </Purpose>
*     <Dependencies>Call Cho\_X\_Init</Dependencies>
*     <Author>Thomas Bondo Pedersen</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*        Performs iterative (eta sweeps) Edmiston-Ruedenberg
*        localization of
*        molecular orbitals. CMO is assumed to be the full set of
*        orbitals, but only the first nOcc() will be localized. On exit,
*        the localized MOs are returned in CMO.
*        Note that symmetry is not allowed.
*        If successful, irc=0 is returned. If irc=-1
*        is returned, an error was found in the input and nothing has
*        been done. If irc=1 is returned, the iterative procedure did
*        not converge to within the threshold (Thr,ThrGrad) in the
*        requested max. number of iterations (MxIter). If the user
*        specifies negative thresholds, Thr=1.0d-6 and ThrGrad=1.0d-3
*        will be used.
*        If the user specifies a negative screening threshold,
*        ThrRot=1.0d-10 is used.
*     </Description>
*    </DOC>
*
************************************************************

      Implicit Real*8 (a-h,o-z)
      Real*8  CMO(*)
      Integer nBas(nSym), nOcc(nSym), nFro(nSym)
      Logical Silent
#include "Molcas.fh"

      Character*5 SecNam
      Parameter (SecNam = 'ERLoc')

      Logical Maximization, Converged, Debug

C     Initialization.
C     ---------------

      irc = 0
      If (MxIter .lt. 1) Return

      nBasT = nBas(1)
      Do iSym = 2,nSym
         nBasT = nBasT + nBas(iSym)
      End Do
      If (nBasT .lt. 1) Return

      nOccT = nOcc(1)
      Do iSym = 2,nSym
         nOccT = nOccT + nOcc(iSym)
      End Do
      If (nOccT .lt. 1) Return

C     Edmiston-Ruedenberg does not work with symmetry.
C     TODO/FIXME: we might consider de-symmetrization.
C     ------------------------------------------------

      If (nSym .ne. 1) Then
         irc = -1
         Return
      End If

C     Localize.
C     ---------

      Functional = -9.9d9
      If (Thr .le. 0.0d0) Then
         ThrLoc = 1.0d-6
      Else
         ThrLoc = Thr
      End If
      If (ThrGrad .le. 0.0d0) Then
         ThrGLoc = 1.0d-3
      Else
         ThrGLoc = ThrGrad
      End If
      If (ThrRot .lt. 0.0d0) Then
         ThrRotLoc = 1.0d-10
      Else
         ThrRotLoc = ThrRot
      End If
      Maximization = .True.
      Converged = .False.
      Debug = .False.
      Call EdmistonRuedenberg(Functional,CMO,ThrLoc,ThrRotLoc,ThrGLoc,
     &                        nBas,nOcc,nFro,
     &                        nSym,MxIter,
     &                        Maximization,Converged,Debug,Silent)

C     Check convergence.
C     ------------------

      If (.not. Converged) irc = 1

      End
