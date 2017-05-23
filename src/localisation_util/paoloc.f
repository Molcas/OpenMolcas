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
      SubRoutine PAOLoc(irc,CMO,PAO,Thr,nBas,nOrb,nOcc,nVir,nSym,Mode)
************************************************************
*
*   <DOC>
*     <Name>PAOLoc</Name>
*     <Syntax>Call PAOLoc(irc,CMO,PAO,Thr,nBas,nOrb,nOcc,nVir,nSym,Mode)
*     </Syntax>
*     <Arguments>
*       \Argument{irc}{return code}{Integer}{out}
*       \Argument{CMO}{MO coefficients}{Real*8}{in}
*       \Argument{PAO}{PAOs}{Real*8}{out}
*       \Argument{Thr}{Cholesky decomposition threshold}{Real*8}{in}
*       \Argument{nBas}{Number of basis functions/irrep}{Integer}{in}
*       \Argument{nOrb}{Number of orbitals/irrep}{Integer}{in}
*       \Argument{nOcc}{Number of occupied orbs/irrep}{Integer}{in}
*       \Argument{nVir}{Number of virtual orbs/irrep}{Integer}{in}
*       \Argument{nSym}{Number of irreps}{Integer}{in}
*       \Argument{Mode}{Mode of calculation}{Character*(*)}{in}
*     </Arguments>
*     <Purpose>Compute projected atomic orbitals (PAOs)</Purpose>
*     <Dependencies>The AO overlap matrix on disk</Dependencies>
*     <Author>Thomas Bondo Pedersen</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*        A set of linearly independent nonorthonormal
*        Projected Atomic Orbitals (PAOs) are computed by projecting
*        the AOs onto the virtual space. Subsequently, the linear
*        dependence may be removed through Cholesky decomposition of a
*        density-type matrix constructed from the linearly dependent
*        PAOs according to D(a,b) = sum(c) PAO(a,c) * PAO(b,c). Finally,
*        orthonormal linearly independent PAOs may be obtained by
*        multiplying with the inverse square root of the PAO overlap
*        matrix.
*        The Mode input string controls which set is returned in array
*        PAO:
*        Mode='RAW': return linearly dependent nonorthonormal PAOs.
*        In this case, array PAO should be dimensioned as nBas*nBas (in
*        symmetry blocks).
*        Mode='CHO': return linearly independent nonorthonormal PAOs
*        obtained by Cholesky decomposition of the density-type matrix.
*        In this case, array PAO should be dimensioned as nBas*nVir (in
*        symmetry blocks).
*        Mode='ORT': return linearly independent orthonormal PAOs.
*        In this case, array PAO should be dimensioned as nBas*nVir (in
*        symmetry blocks).
*        If Mode is anything else, Mode='ORT' is assumed.
*        Note that if nOcc+nVir<nBas the remaining nBas-nOcc-nVir
*        orbitals are considered part of the orthogonal complement of
*        the virtual space (and are projected out of the AOs when
*        obtaining the raw PAOs). Thus, one may obtain PAOs for any
*        subset of orbitals by specifying nOcc and nVir arrays
*        appropriately.
*     </Description>
*    </DOC>
*
************************************************************
      Implicit Real*8 (a-h,o-z)
      Real*8  CMO(*), PAO(*)
      Integer nBas(nSym), nOrb(nSym), nOcc(nSym), nVir(nSym)
      Character*(*) Mode
#include "WrkSpc.fh"

      Character*6 SecNam
      Parameter (SecNam = 'PAOLoc')

      Character*3 myMode, DefMode
      Parameter (DefMode = 'ORT')
      Integer DefLevel
      Real*8  DefThr
      Parameter (DefThr = 1.0d-12)

      Logical TestOrth, Normalize

#if defined (_DEBUG_)
      TestOrth = .True.
#else
      TestOrth = .False.
#endif

C     Set return code.
C     ----------------

      irc = 0

C     Set DefLevel from DefMode.
C     --------------------------

      If (DefMode .eq. 'RAW') Then
         DefLevel = 1
      Else If (DefMode .eq. 'CHO') Then
         DefLevel = 2
      Else If (DefMode .eq. 'ORT') Then
         DefLevel = 3
      Else
         Call SysAbendMsg(SecNam,'DefMode not recognized!',
     &                    'Note: this is a programming error...')
         DefLevel = -999999
      End If

C     Interpret Mode input.
C     ---------------------

      lMode = len(Mode)
      If (lMode .lt. 3) Then
         Level = DefLevel
      Else
         myMode = Mode(1:3)
         Call UpCase(myMode)
         If (myMode .eq. 'RAW') Then
            Level = 1
         Else If (myMode .eq. 'CHO') Then
            Level = 2
         Else If (myMode .eq. 'ORT') Then
            Level = 3
         Else
            Level = DefLevel
         End If
      End If

C     Make a dummy allocation to enable de-allocation by flushing.
C     ------------------------------------------------------------

      l_Dum = 1
      Call GetMem('PAOL_Dummy','Allo','Real',ip_Dum,l_Dum)

C     Compute raw projected AOs.
C     --------------------------

      l_R = nBas(1)**2
      Do iSym = 2,nSym
         l_R = l_R + nBas(iSym)**2
      End Do
      Call GetMem('PAOL_R','Allo','Real',ip_R,l_R)

      Normalize = .True.
      Call GetRawPAOs(Work(ip_R),CMO,nBas,nOrb,nOcc,nVir,nSym,Normalize)

      If (Level .eq. 1) Then
         Call dCopy_(l_R,Work(ip_R),1,PAO,1)
         Go To 1 ! return after de-allocation
      End If

C     Use Cholesky decomposition to compute a linearly independent set
C     of nonorthonormal PAOs.
C     ----------------------------------------------------------------

      l_D = nBas(1)**2
      Do iSym = 2,nSym
         l_D = max(l_D,nBas(iSym)**2)
      End Do
      Call GetMem('PAOL_D','Allo','Real',ip_D,l_D)

      If (Thr .le. 0.0d0) Then
         ThrLoc = DefThr
      Else
         ThrLoc = Thr
      End If

      kOffR = ip_R
      kOffP = 1
      Do iSym = 1,nSym
         If (nVir(iSym) .gt. 0) Then
            Call GetDens_Localisation(Work(ip_D),Work(kOffR),
     &                                nBas(iSym),nBas(iSym))
            Call ChoLoc(irc,Work(ip_D),PAO(kOffP),ThrLoc,xNrm,
     &                  nBas(iSym),nVir(iSym))
            If (irc .ne. 0) Go To 1 ! return after de-allocation
         End If
         kOffR = kOffR + nBas(iSym)**2
         kOffP = kOffP + nBas(iSym)*nVir(iSym)
      End Do

      If (Level .eq. 2) Then
         Go To 1 ! return after de-allocation
      End If

C     Orthonormalize the PAOs.
C     ------------------------

      kOffP = 1
      kOffR = ip_R
      Do iSym = 1,nSym
         kOff1 = kOffR + nBas(iSym)*nOcc(iSym)
         Call dCopy_(nBas(iSym)*nVir(iSym),PAO(kOffP),1,Work(kOff1),1)
         kOffP = kOffP + nBas(iSym)*nVir(iSym)
         kOffR = kOffR + nBas(iSym)**2
      End Do
      nOrthPs = 2 ! orthonormalization passes to ensure num. accuracy
      Call OrthoPAO_Localisation(Work(ip_R),nBas,nOcc,nVir,nSym,nOrthPs,
     &                           TestOrth)
      kOffP = 1
      kOffR = ip_R
      Do iSym = 1,nSym
         kOff1 = kOffR + nBas(iSym)*nOcc(iSym)
         Call dCopy_(nBas(iSym)*nVir(iSym),Work(kOff1),1,PAO(kOffP),1)
         kOffP = kOffP + nBas(iSym)*nVir(iSym)
         kOffR = kOffR + nBas(iSym)**2
      End Do

      If (Level .eq. 3) Then
         Go To 1 ! return after de-allocation
      End If

C     De-allocation by flushing.
C     --------------------------

    1 Call GetMem('PAOL_Dummy','Flus','Real',ip_Dum,l_Dum)
      Call GetMem('PAOL_Dummy','Free','Real',ip_Dum,l_Dum)

      End
