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
* Copyright (C) 2006, Thomas Bondo Pedersen                            *
************************************************************************
*  DefineDomain
*
*> @brief
*>   Define orbital domains
*> @author Thomas Bondo Pedersen
*>
*> @details
*> Orbital domains are defined by summing up gross atomic Mulliken
*> charges for each orbital. Once the sum is greater than or equal
*> to the threshold \p ThrDomain(1), the domain is defined. The
*> second threshold \p ThrDomain(2) is then used to check
*> completeness (Pulay-style) of the definition and if needed,
*> more atoms are added to the domain. To avoid the second step
*> (i.e. the completeness check), simply put \p ThrDomain(2) > ``1.0d0`` in
*> which case array \p f is undefined on exit.
*>
*> On exit, the contents of iDomain array are:
*>
*> - \p iDomain(0,i): number of atoms in domain \c i.
*> - \p iDomain(n,i): id of atom \c n (\c n = ``1``, ``2``, ..., \p iDomain(0,i)) in domain \c i.
*>
*> Return codes:
*>
*> - \p irc = ``-1``: input error(s) detected.
*> - \p irc =  ``0``: all OK.
*> - \p irc =  ``1``: this indicates a bug in the charge setup.
*> - \p irc =  ``2``: can only be set if debug is turned on; in that case, this code means
*>                    that the computed charges do not sum up to the number of occupied orbitals, \p nOcc.
*> - \p irc =  ``3``: can only be set if debug is turned on; in that case, this code means
*>                    that at least one domain is defined as empty or having at least one atom index out of bounds.
*>
*> @param[out] irc           Return code
*> @param[out] iDomain       Domain definition
*> @param[out] QD            Final total charges for each domain
*> @param[out] f             Function values from completeness check
*> @param[in]  C             MO coefficients
*> @param[in]  ThrDomain     Thresholds for domain definition
*> @param[in]  nBas_per_Atom Number of basis functions on each atom
*> @param[in]  nBas_Start    Start index for basis functions on each atom
*> @param[in]  nAtom         Number of atoms
*> @param[in]  nBas          Number of basis functions
*> @param[in]  nOcc          Number of occupied orbitals
************************************************************************
      SubRoutine DefineDomain(irc,iDomain,QD,f,C,ThrDomain,
     &                        nBas_per_Atom,nBas_Start,
     &                        nAtom,nBas,nOcc)
      Implicit Real*8 (a-h,o-z)
      Integer iDomain(0:nAtom,nOcc)
      Integer nBas_per_Atom(nAtom), nBas_Start(nAtom)
      Real*8  QD(nOcc), f(nOcc), C(nBas,nOcc), ThrDomain(2)
#include "WrkSpc.fh"

      external ddot_

      Integer nB(1)

      Logical LocDbg
#if defined (_DEBUGPRINT_)
      Parameter (LocDbg = .True.)
#else
      Parameter (LocDbg = .False.)
#endif

C     Check input.
C     ------------

      irc = 0
      If (nAtom.lt.1 .or. nBas.lt.1 .or. nOcc.lt.1) Return

C     Allocate and read overlap matrix stored as full square.
C     -------------------------------------------------------

      l_S = nBas*nBas
      Call GetMem('DfDm_S','Allo','Real',ip_S,l_S)

      nB(1) = nBas
      Call GetOvlp_Localisation(Work(ip_S),'Sqr',nB,1)

C     Allocations.
C     ------------

      l_T = nBas*nOcc
      l_Q = nAtom*nOcc
      Call GetMem('DfDm_T','Allo','Real',ip_T,l_T)
      Call GetMem('DfDm_Q','Allo','Real',ip_Q,l_Q)

C     Compute T=SC.
C     -------------

      Call DGEMM_('N','N',nBas,nOcc,nBas,
     &           1.0d0,Work(ip_S),nBas,C,nBas,
     &           0.0d0,Work(ip_T),nBas)

C     Compute atomic contributions to Mulliken charges.
C     -------------------------------------------------

      Call dCopy_(l_Q,[0.0d0],0,Work(ip_Q),1)
      iOff0 = ip_T - 1
      jOff0 = ip_Q - 1
      Do i = 1,nOcc
         iOff = iOff0 + nBas*(i-1)
         jOff = jOff0 + nAtom*(i-1)
         Do iAtom = 1,nAtom
            Work(jOff+iAtom) = Work(jOff+iAtom)
     &      + dDot_(nBas_per_Atom(iAtom),C(nBas_Start(iAtom),i),1,
     &                                  Work(iOff+nBas_Start(iAtom)),1)
         End Do
      End Do

C     Debug: check charges.
C     The sum of them should equal the trace of the occupied MO overlap
C     matrix, i.e. nOcc.
C     -----------------------------------------------------------------

      If (LocDbg) Then
         Write(6,*)
         Write(6,*) 'DefineDomain: checking charge calculation:'
         Charge = 0.0d0
         kOff0 = ip_Q - 1
         Do i = 1,nOcc
            Chrg = 0.0d0
            kOff = kOff0 + nAtom*(i-1)
            Do iAtom = 1,nAtom
               Chrg = Chrg + Work(kOff+iAtom)
            End Do
            Charge = Charge + Chrg
            Diff = Chrg - 1.0d0
            If (abs(Diff) .gt. 1.0d-10) Then
               x1 = 1.0d0
               Write(6,*)
               Write(6,*) '  Orbital ',i,':'
               Write(6,*) '  Charge    : ',Chrg
               Write(6,*) '  Expected  : ',x1
               Write(6,*) '  Difference: ',Diff
            End If
         End Do
         Diff = Charge - dble(nOcc)
         Write(6,*)
         Write(6,*) '  Total charge: ',Charge
         Write(6,*) '  Expected    : ',dble(nOcc)
         Write(6,*) '  Difference  : ',Diff
         If (abs(Diff) .gt. 1.0d-10) Then
            irc = 2
            Go To 1 ! return after deallocation
         End If
      End If

C     For each orbital, create an index array in descending order of
C     absolute contributions to charges. Array iDomain will then contain
C     the list of atoms ordered according to contributions.
C     ------------------------------------------------------------------

      l_iPivot = nAtom
      l_absQ = nAtom
      Call GetMem('DfDm_iPivot','Allo','Inte',ip_iPivot,l_iPivot)
      Call GetMem('DfDm_absQ','Allo','Real',ip_absQ,l_absQ)
      Do i = 1,nOcc
         iOff = ip_Q  + nAtom*(i-1)
         nSrt = nAtom
         Do iAtom = 0,nAtom-1
            Work(ip_absQ+iAtom) = abs(Work(iOff+iAtom))
         End Do
         Call CD_DiaMax(Work(ip_absQ),nAtom,iWork(ip_iPivot),
     &                  iDomain(1,i),nSrt,-1.0d0)
         If (nSrt .ne. nAtom) Then
            Call GetMem('DfDm_iPivot','Free','Inte',ip_iPivot,l_iPivot)
            irc = 1 ! ooops: something is fishy here...
            Go To 1 ! return after deallocations
         End If
      End Do
      Call GetMem('DfDm_absQ','Free','Real',ip_absQ,l_absQ)
      Call GetMem('DfDm_iPivot','Free','Inte',ip_iPivot,l_iPivot)

C     For each orbital, define domain according to charge threshold.
C     --------------------------------------------------------------

      iOff0 = ip_Q - 1
      Do i = 1,nOcc
         iOff = iOff0 + nAtom*(i-1)
         iCount = 1
         iAtom  = iDomain(iCount,i)
         Charge = Work(iOff+iAtom)
         Do While (iCount.lt.nAtom .and. Charge.lt.ThrDomain(1))
            iCount = iCount + 1
            iAtom  = iDomain(iCount,i)
            Charge = Charge + Work(iOff+iAtom)
         End Do
         iDomain(0,i) = iCount
      End Do

C     Debug.
C     ------

      If (LocDbg) Then
         nErr = 0
         Write(6,*)
         Write(6,*) 'DefineDomain: domains and charges after step 1:'
         Write(6,*) 'Threshold: ',ThrDomain(1)
         Do i = 1,nOcc
            Write(6,*)
            Write(6,*) 'Domain ',i,': ',iDomain(0,i),' atoms:'
            If (iDomain(0,i) .lt. 1) Then
               Write(6,*) 'No atoms in domain !?!?!'
               nErr = nErr + 1
            Else If (iDomain(0,i) .gt. nAtom) Then
               Write(6,*) 'Number of atoms > nAtom in domain !?!?!'
               nErr = nErr + 1
            Else
               Charge = 0.0d0
               kOff0 = ip_Q - 1 + nAtom*(i-1)
               Do iAt = 1,iDomain(0,i)
                  iAtom = iDomain(iAt,i)
                  If (iAtom.lt.1 .or. iAtom.gt.nAtom) Then
                     Write(6,*) '  Atom: ',iAtom,' !?!?!'
                     nErr = nErr + 1
                  Else
                     kOff = kOff0 + iAtom
                     Write(6,*) '  Atom: ',iAtom,
     &                          '  Charge: ',Work(kOff)
                     Charge = Charge + Work(kOff)
                  End If
               End Do
               Write(6,*) '  Total charge: ',Charge
            End If
         End Do
         If (nErr .ne. 0) Then
            irc = 3
            Go To 1 ! return after deallocation
         End If
      End If

C     For each orbital, check completeness and add atoms as needed to
C     meet the requirement f<=threshold.
C     ---------------------------------------------------------------

      If (ThrDomain(2) .lt. 1.0d0) Then
         Do i = 1,nOcc
            kOffT = ip_T + nBas*(i-1)
            Call MakeDomainComplete(iDomain(0,i),f(i),Work(ip_S),
     &                              Work(kOffT),
     &                              ThrDomain(2),nBas_per_Atom,
     &                              nBas_Start,nBas,nAtom)
         End Do
      End If

C     Compute total charges for each domain.
C     --------------------------------------

      Do i = 1,nOcc
         kOff = ip_Q - 1 + nAtom*(i-1)
         QD(i) = 0.0d0
         Do iAt = 1,iDomain(0,i)
            iAtom = iDomain(iAt,i)
            QD(i) = QD(i) + Work(kOff+iAtom)
         End Do
      End Do

C     Debug.
C     ------

      If (LocDbg) Then
         nErr = 0
         Write(6,*)
         Write(6,*) 'DefineDomain: domains and charges after step 2:'
         Write(6,*) 'Threshold: ',ThrDomain(2)
         If (ThrDomain(2) .lt. 1.0d0) Then
            Do i = 1,nOcc
               Write(6,*)
               Write(6,*) 'Domain ',i,': ',iDomain(0,i),' atoms:'
               If (iDomain(0,i) .lt. 1) Then
                  Write(6,*) 'No atoms in domain !?!?!'
                  nErr = nErr + 1
               Else If (iDomain(0,i) .gt. nAtom) Then
                  Write(6,*) 'Number of atoms > nAtom in domain !?!?!'
                  nErr = nErr + 1
               Else
                  Charge = 0.0d0
                  kOff0 = ip_Q - 1 + nAtom*(i-1)
                  Do iAt = 1,iDomain(0,i)
                     iAtom = iDomain(iAt,i)
                     If (iAtom.lt.1 .or. iAtom.gt.nAtom) Then
                        Write(6,*) '  Atom: ',iAtom,' !?!?!'
                        nErr = nErr + 1
                     Else
                        kOff = kOff0 + iAtom
                        Write(6,*) '  Atom: ',iAtom,
     &                             '  Charge: ',Work(kOff)
                        Charge = Charge + Work(kOff)
                     End If
                  End Do
                  Write(6,*) '  Total charge: ',Charge
                  If (abs(Charge-QD(i)) .gt. 1.0d-12) Then
                     Write(6,*) 'Total charge is inconsistent with QD !'
                     nErr = nErr + 1
                  End If
               End If
            End Do
         Else
            Write(6,*) 'Threshold >= 1.0d0: step 2 was skipped.'
            Write(6,*) 'Domains are unchanged from step 1.'
         End If
         If (nErr .ne. 0) Then
            irc = 3
            Go To 1 ! return after deallocation
         End If
      End If

C     Deallocations.
C     --------------

    1 Call GetMem('DfDm_Q','Free','Real',ip_Q,l_Q)
      Call GetMem('DfDm_T','Free','Real',ip_T,l_T)
      Call GetMem('DfDm_S','Free','Real',ip_S,l_S)

      End
      SubRoutine MakeDomainComplete(iDomain,f,S,T,Threshold,
     &                              nBas_per_Atom,nBas_Start,nBas,nAtom)
C
C     Thomas Bondo Pedersen, January 2006.
C
C     Purpose: Boughton-Pulay completeness check of an orbital domain.
C              Extend domain as needed to obtain completeness.
C
      Implicit Real*8 (a-h,o-z)
      Integer iDomain(0:nAtom)
      Real*8  S(nBas,nBas), T(nBas)
      Integer nBas_per_Atom(nAtom), nBas_Start(nAtom)
#include "WrkSpc.fh"

      external ddot_

      Character*18 SecNam
      Parameter (SecNam = 'MakeDomainComplete')

      Character*80 Txt
      Logical Complete

      nA = iDomain(0)
      f = 0.0d0
      Complete = nA.eq.nAtom
      Do While (nA.lt.nAtom .and. .not.Complete)

C        Allocate S[i], T[i], and Scr[i].
C        --------------------------------

         nSize = nBas_per_Atom(iDomain(1))
         Do iA = 2,nA
            nSize = nSize + nBas_per_Atom(iDomain(iA))
         End Do
         l_Si = nSize*nSize
         l_Sl = l_Si
         l_Ti = nSize
         l_Scr = nSize
         Call GetMem('MkDmC_Si','Allo','Real',ip_Si,l_Si)
         Call GetMem('MkDmC_Sl','Allo','Real',ip_Sl,l_Sl)
         Call GetMem('MkDmC_Ti','Allo','Real',ip_Ti,l_Ti)
         Call GetMem('MkDmC_Scr','Allo','Real',ip_Scr,l_Scr)

C        Get S[i] and T[i].
C        ------------------

         kTi = ip_Ti
         iCol = 0
         Do iB = 1,nA
            nu1 = nBas_Start(iDomain(iB))
            nnu = nBas_per_Atom(iDomain(iB))
            Do lnu = 0,nnu-1
               nu = nu1 + lnu
               iCol = iCol + 1
               iRow = 0
               kSi = ip_Si + nSize*(iCol-1)
               Do iA = 1,nA
                  mu1 = nBas_Start(iDomain(iA))
                  nmu = nBas_per_Atom(iDomain(iA))
                  Call dCopy_(nmu,S(mu1,nu),1,Work(kSi+iRow),1)
                  iRow = iRow + nmu
               End Do
            End Do
            Call dCopy_(nnu,T(nu1),1,Work(kTi),1)
            kTi = kTi + nnu
         End Do

C        Solve S[i]Y=T[i] (Y stored in T[i] on exit).
C        Use a scratch array for S[i], as it is ruined on exit from
C        LinEqSolv.
C        ----------------------------------------------------------

         irc = 0
         Call dCopy_(l_Si,Work(ip_Si),1,Work(ip_Sl),1)
         Call LinEqSolv(irc,'N',Work(ip_Sl),nSize,Work(ip_Ti),nSize,
     &                  nSize,1)
         If (irc .ne. 0) Then
            Write(Txt,'(A,I9)') 'LinEqSolv returned',irc
            If (irc .lt. 0) Then
               Call SysAbendMsg(SecNam,Txt,'LinEqSolv input error!')
            Else
               Call SysAbendMsg(SecNam,Txt,
     &                          'Singular domain overlap matrix!')
            End If
         End If

C        Compute f=1-Y(T)S[i]Y.
C        ----------------------

         Call dGeMV_('N',nSize,nSize,1.0d0,Work(ip_Si),nSize,
     &              Work(ip_Ti),1,0.0d0,Work(ip_Scr),1)
         f = 1.0d0 - dDot_(nSize,Work(ip_Ti),1,Work(ip_Scr),1)

C        Deallocation.
C        -------------

         Call GetMem('MkDmC_Scr','Free','Real',ip_Scr,l_Scr)
         Call GetMem('MkDmC_Ti','Free','Real',ip_Ti,l_Ti)
         Call GetMem('MkDmC_Sl','Free','Real',ip_Sl,l_Sl)
         Call GetMem('MkDmC_Si','Free','Real',ip_Si,l_Si)

C        Check completeness (f<=Threshold).
C        If not complete, add next atom to domain.
C        (If complete, we break the while loop.)
C        -----------------------------------------

         Complete = f.le.Threshold
         If (.not.Complete) Then
            nA = nA + 1
         End If

      End Do

C     Set new domain size.
C     --------------------

      iDomain(0) = nA

      End
