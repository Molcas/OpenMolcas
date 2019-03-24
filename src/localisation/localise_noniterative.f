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
      SubRoutine Localise_Noniterative(irc,Model,xNrm)
C
C     Author: T.B. Pedersen
C
C     Purpose: Non-iterative localisation of orbitals.
C              Models implemented:
C                Cholesky [MODEL='CHOL']
C                PAO      [MODEL='PAO ']
C
      Implicit Real*8 (a-h,o-z)
      Character*4 Model
#include "Molcas.fh"
#include "inflocal.fh"
#include "WrkSpc.fh"

      Character*21 SecNam
      Parameter (SecNam = 'Localise_Noniterative')

      Character*4  myModel
      Character*6  Namefile
      Character*80 Txt

      Logical Test_OrthoPAO, Normalize
      Parameter (Test_OrthoPAO = .False.)

      Dimension dum(1),idum(1)

      irc = 0
      xNrm = 0.0d0

      myModel = Model
      Call UpCase(myModel)
      If (myModel .eq. 'CHOL') Then
*        If (.not.Silent) Then
            Write(6,'(/,1X,A)') 'Cholesky localisation'
            Write(6,'(1X,A,1X,D12.4,A)')
     &      'Convergence threshold:',Thrs,' (decomposition)'
            Write(6,'(1X,A,8(1X,I6))')
     &      'Frozen orbitals      :',(nFro(iSym),iSym=1,nSym)
            Write(6,'(1X,A,8(1X,I6))')
     &      'Orbitals to localise :',(nOrb2Loc(iSym),iSym=1,nSym)
*        End If
         l_Dens = nBas(1)**2
         Do iSym = 2,nSym
            l_Dens = max(l_Dens,nBas(iSym)**2)
         End Do
         Call GetMem('Density','Allo','Real',ip_Dens,l_Dens)
         kOffC = ipCMO
         Do iSym = 1,nSym
            If (nOrb2Loc(iSym) .gt. 0) Then
               kOff1 = kOffC + nBas(iSym)*nFro(iSym)
               Call GetDens_Localisation(Work(ip_Dens),Work(kOff1),
     &                                   nBas(iSym),nOrb2Loc(iSym))
               Call ChoLoc(irc,Work(ip_Dens),Work(kOff1),Thrs,yNrm,
     &                     nBas(iSym),nOrb2Loc(iSym))
               xNrm = xNrm + yNrm*yNrm
               If (irc .ne. 0) Then
                  Call GetMem('Density','Free','Real',ip_Dens,l_Dens)
                  irc  = 1
                  xNrm = -9.9d9
                  Return
               End If
            End If
            kOffC = kOffC + nBas(iSym)**2
         End Do
         xNrm = sqrt(xNrm)
         Call GetMem('Density','Free','Real',ip_Dens,l_Dens)
      Else If (myModel .eq. 'PAO ') Then
*        If (.not.Silent) Then
            Write(6,'(/,1X,A)') 'PAO Cholesky localisation'
            Write(6,'(1X,A,1X,D12.4,A)')
     &      'Convergence threshold:',Thrs,' (decomposition)'
            Write(6,'(1X,A,8(1X,I6))')
     &      'Frozen orbitals      :',(nFro(iSym),iSym=1,nSym)
            Write(6,'(1X,A,8(1X,I6))')
     &      'Orbitals to localise :',(nOrb2Loc(iSym),iSym=1,nSym)
*        End If
         l_Dv = nBas(1)**2
         l_R = nBas(1)**2
         Do iSym = 2,nSym
            l_R = l_R + nBas(iSym)**2
            l_Dv = max(l_Dv,nBas(iSym)**2)
         End Do
         Call GetMem('R','Allo','Real',ip_R,l_R)
         Call GetMem('Dv','Allo','Real',ip_Dv,l_Dv)
         Normalize = .True.
         Call GetRawPAOs(Work(ip_R),Work(ipCMO),nBas,nOrb,nFro,nOrb2Loc,
     &                   nSym,Normalize)
         kSav = 0
         If (AnaPAO) Then
            l_DvSav = l_R
            Call GetMem('DvSav','Allo','Real',ip_DvSav,l_DvSav)
            kSav = ip_DvSav
         End If
         kOffR = ip_R
         kOffC = ipCMO
         Do iSym = 1,nSym
            If (nOrb2Loc(iSym) .gt. 0) Then
               Call GetDens_Localisation(Work(ip_Dv),Work(kOffR),
     &                                   nBas(iSym),nBas(iSym))
               If (AnaPAO) Then
                  Call dCopy_(nBas(iSym)**2,Work(ip_Dv),1,Work(kSav),1)
                  kSav = kSav + nBas(iSym)**2
               End If
               kOff1 = kOffC + nBas(iSym)*nFro(iSym)
               Call ChoLoc(irc,Work(ip_Dv),Work(kOff1),Thrs,yNrm,
     &                     nBas(iSym),nOrb2Loc(iSym))
               xNrm = xNrm + yNrm*yNrm
               If (irc .ne. 0) Then
                  If (AnaPAO) Then
                     Call GetMem('DvSav','Free','Real',ip_DvSav,l_DvSav)
                  End If
                  Call GetMem('Dv','Free','Real',ip_Dv,l_Dv)
                  Call GetMem('R','Free','Real',ip_R,l_R)
                  irc  = 1
                  xNrm = -9.9d9
                  Return
               End If
            End If
            kOffR = kOffR + nBas(iSym)**2
            kOffC = kOffC + nBas(iSym)**2
         End Do
         xNrm = sqrt(xNrm)
         If (AnaPAO) Then
            Call PAO_Analysis(Work(ip_DvSav),Work(ip_R),Work(ipCMO))
            Call GetMem('DvSav','Free','Real',ip_DvSav,l_DvSav)
         End If
         Write(Namefile,'(A)') 'DPAORB'
         Write(Txt,'(80X)')
         Write(Txt,'(A)') 'Linearly dependent PAOs'
         Lu_=isFreeUnit(11)
         Call WrVec_Localisation(Namefile,Lu_,'CO',nSym,nBas,nBas,
     &                           Work(ip_R),Work(ipOcc),dum,idum,Txt)
*        If (.not.Silent) Then
            Write(6,'(1X,A)') 'The DPAORB file has been written.'
*        End If
         Write(Namefile,'(A)') 'IPAORB'
         Write(Txt,'(80X)')
         Write(Txt,'(A)') 'Linearly independent PAOs'
         Lu_=isFreeUnit(11)
         Call WrVec_Localisation(Namefile,Lu_,'CO',nSym,nBas,nBas,
     &                           Work(ipCMO),Work(ipOcc),dum,idum,Txt)
*        If (.not.Silent) Then
            Write(6,'(1X,A)') 'The IPAORB file has been written.'
*        End If
         Call GetMem('Dv','Free','Real',ip_Dv,l_Dv)
         Call GetMem('R','Free','Real',ip_R,l_R)
         nOrPs = 2 ! use 2 orthonorm. passes for num. accuracy
         Call OrthoPAO_Localisation(Work(ipCMO),nBas,nFro,nOrb2Loc,nSym,
     &                              nOrPs,Test_OrthoPAO)
      Else
         Write(Txt,'(A,A4)') 'Model = ',Model
         Call SysAbendMsg(SecNam,'Unknown model',Txt)
      End If

      End
