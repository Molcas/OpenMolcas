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
* Copyright (C) 2005,2006, Thomas Bondo Pedersen                       *
************************************************************************
      SubRoutine GetRawPAOs(R,C,nBas,nOrb,nFro,nOrb2Loc,nSym,Normalize)
C
C     Thomas Bondo Pedersen, December 2005.
C     - revised january 2006 (Thomas Bondo Pedersen).
C
C     Purpose: compute projected AOs spanning the primary space from the
C              formula
C
C              R = 1 - Do*S = D*S
C
C              where Do is the "AO density" matrix of the orthogonal
C              complement orbitals and D that of the orbital space to be
C              localised (primary space).
C              S is the AO overlap matrix (which is read from disk).
C
C              Which formula is used depends on the dimensions of the
C              two spaces (such that the most economical is computed).
C
C              If (Normalize): normalize each PAO (recommended).
C                              Note that if the norm of the PAO is
C                              smaller than 1.0d-6, it will not be
C                              normalized (it is left unchanged).
C
C-------------------------------------------------------------
C-TODO/FIXME: it is in most cases faster to use Do*S=C*(C^T*S)
C-------------------------------------------------------------
C
      Implicit Real*8 (a-h,o-z)
      Integer nBas(nSym), nOrb(nSym), nFro(nSym), nOrb2Loc(nSym)
      Real*8  R(*), C(*)
      Logical Normalize
#include "WrkSpc.fh"

      Character*10 SecNam
      Parameter (SecNam = 'GetRawPAOs')

      Character*80 Txt

      external ddot_

C     Read the overlap matrix from disk.
C     ----------------------------------

      lOvlp = nBas(1)**2
      Do iSym = 2,nSym
         lOvlp = lOvlp + nBas(iSym)**2
      End Do
      Call GetMem('Ovlp','Allo','Real',ipOvlp,lOvlp)
      Call GetOvlp_Localisation(Work(ipOvlp),'Sqr',nBas,nSym)

C     Compute R.
C     ----------

      lDo = nBas(1)**2
      Do iSym = 2,nSym
         lDo = max(lDo,nBas(iSym)**2)
      End Do
      Call GetMem('Do','Allo','Real',ipDo,lDo)

      kOff  = 1
      kOffS = ipOvlp
      Do iSym = 1,nSym

         nB = nBas(iSym)
         If (nB .gt. 0) Then

            nF    = nFro(iSym)
            nO2L  = nOrb2Loc(iSym)
            nRest = nOrb(iSym) - nF - nO2L
            nOrth = nF + nRest ! dim. of orthogonal complement

            If (nO2L .lt. 1) Then ! R = 0
               Call fZero(R(kOff),nB**2)
            Else If (nOrth .lt. 0) Then ! error
               Call SysAbendMsg(SecNam,
     &              'Dimension of orthogonal complement space < 0',' ')
            Else If (nOrth .eq. 0) Then ! R = 1
               Call fZero(R(kOff),nB**2)
               Do i = 1,nB
                  R(kOff-1+nB*(i-1)+i) = 1.0d0
               End Do
            Else If (nOrth .lt. nO2L) Then ! R = 1 - Do*S
               If (nRest .gt. 0) Then
                 lOff = kOff + nB*(nF+nO2L)
                 Call GetDens_Localisation(Work(ipDo),C(lOff),nB,nRest)
               Else
                 Call fZero(Work(ipDo),nB**2)
               End If
               If (nF .gt. 0) Then
                  Call GetDens_Localisation(R(kOff),C(kOff),nB,nF)
                  Call dAXPY_(nB**2,1.0d0,R(kOff),1,Work(ipDo),1)
               End If
               Call DGEMM_('N','N',nB,nB,nB,
     &                    -1.0d0,Work(ipDo),nB,Work(kOffS),nB,
     &                     0.0d0,R(kOff),nB)
               Do i = 1,nB
                  R(kOff-1+nB*(i-1)+i) = R(kOff-1+nB*(i-1)+i) + 1.0d0
               End Do
            Else ! R = D*S
               lOff = kOff + nB*nF
               Call GetDens_Localisation(Work(ipDo),C(lOff),nB,nO2L)
               Call DGEMM_('N','N',nB,nB,nB,
     &                    1.0d0,Work(ipDo),nB,Work(kOffS),nB,
     &                    0.0d0,R(kOff),nB)
            End If

            kOff  = kOff  + nB**2
            kOffS = kOffS + nB**2

         End If

      End Do

C     If requested, normalize the PAOs.
C     ---------------------------------

      If (Normalize) Then
         kOff  = 1
         kOffS = ipOvlp
         Do iSym = 1,nSym
            nB = nBas(iSym)
            If (nB .gt. 0) Then
               Call DGEMM_('N','N',nB,nB,nB,
     &                    1.0d0,Work(kOffS),nB,R(kOff),nB,
     &                    0.0d0,Work(ipDo),nB)
               Do mu = 0,nB-1
                  kR  = kOff + nB*mu
                  kSR = ipDo + nB*mu
                  Ovlp = dDot_(nB,R(kR),1,Work(kSR),1)
                  If (Ovlp .gt. 1.0d-6) Then
                     Fac = 1.0d0/sqrt(Ovlp)
                     Call dScal_(nB,Fac,R(kR),1)
                  Else If (Ovlp .lt. 0.0d0) Then
                     Write(Txt,'(A,1P,D15.5)') 'Overlap = ',Ovlp
                     Call SysAbendMsg(SecNam,
     &                                'Negative raw PAO overlap!',Txt)
                  End If
               End Do
               kOff  = kOff + nB**2
               kOffS = kOffS + nB**2
            End If
         End Do
      End If

      Call GetMem('Do','Free','Real',ipDo,lDo)
      Call GetMem('Ovlp','Free','Real',ipOvlp,lOvlp)

      End
