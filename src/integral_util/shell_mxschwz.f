!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!***********************************************************************
!                                                                      *
!  Subroutine Shell_MxSchwz:  gets max integral estimates for each     *
!                             shell pair...                            *
!                                                                      *
!***********************************************************************
      SubRoutine Shell_MxSchwz(nSkal,Schwz_Shl)
!----------------------------------------------------------------------
      use iSD_data, only: iSD
      use Basis_Info, only: Shells, DBSC
      use Symmetry_Info, only: nIrrep
      use Constants, only: Zero
      use k2_structure, only: k2Data, IndK2
      use k2_arrays, only: DoHess_
      Implicit None
      Integer nSkal
      Real*8 Schwz_Shl(nSkal,nSkal)

      Integer ixyz, nabSz, ik2, iS, iShll, iShell, iCmp, iAng, iCnttp,
     &                          jS, jShll, jShell, jCmp, jAng, jCnttp,
     &        ijS, nDCRR, nHm, lDCRR
      Real*8 Schwz_Tmp
!
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
!
!     loop over shell pair...
      Schwz_Shl(:,:)=Zero
      Do iS = 1, nSkal
        iShll= iSD( 0,iS)
        If (Shells(iShll)%Aux .and. iS.ne.nSkal) Go To 100
        iShell=iSD(11,iS)
        iCmp=iSD(2,iS)
        iAng=iSD(1,iS)
        iCnttp=iSD(13,iS)
        Do jS = 1, iS
          jShll= iSD( 0,jS)
          If (Shells(iShll)%Aux.and..Not.Shells(jShll)%Aux) Go To 200
          If (Shells(jShll)%Aux .and. jS.eq.nSkal) Go To 200
!         Write (*,*) 'Shell_..:iS,jS=',iS,jS
          jShell=iSD(11,jS)
          jCmp=iSD(2,jS)
          jAng=iSD(1,jS)
          jCnttp=iSD(13,jS)
          If (iShell.ge.jShell) Then
            ijS = iShell*(iShell-1)/2 + jShell
          Else
            ijS = jShell*(jShell-1)/2 + iShell
          End If
          nDCRR = Indk2(2,ijS)
          ik2   = Indk2(3,ijS)
          nHm=iCmp*jCmp*(nabSz(iAng+jAng)-nabSz(Max(iAng,jAng)-1))
          nHm=nHm*nIrrep
          If (DoHess_) nHm=0
!         now loop over  R operator...
          If (dbsc(iCnttp)%fMass.eq.dbsc(jCnttp)%fMass) Then
             Schwz_tmp = k2data(1,ik2)%EstI
             Do lDCRR = 2, nDCRR
                Schwz_tmp=Max(Schwz_tmp,k2data(lDCRR,ik2)%EstI)
             End Do
          Else
             Schwz_tmp=Zero
          End If
          Schwz_Shl(jS,iS)=Schwz_tmp
          Schwz_Shl(iS,jS)=Schwz_tmp
 200      Continue
        End Do
 100    Continue
      End Do
!     Call RecPrt('Schwz_shl',' ',Schwz_Shl,nSkal,nSkal)
      end SubRoutine Shell_MxSchwz
