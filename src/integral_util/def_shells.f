
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Def_Shells(iSD,nSD,mSkal)
      use Basis_Info
      use Center_Info
      use Sizes_of_Seward, only: S
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "Basis_Mode_Parameters.fh"
#include "Basis_Mode.fh"
#include "disp.fh"
*
      Integer iSD(0:nSD,mSkal)
      Logical  TF, TstFnc
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
*
      TF(mdc,iIrrep,iComp) = TstFnc(dc(mdc)%iCoSet,
     &                              iIrrep,iComp,dc(mdc)%nStab)
*                                                                      *
************************************************************************
*                                                                      *
      If (Basis_Mode.ne.Valence_Mode .and.
     &    Basis_Mode.ne.Auxiliary_Mode .and.
     &    Basis_Mode.ne.Fragment_Mode .and.
     &    Basis_Mode.ne.With_Auxiliary_Mode .and.
     &    Basis_Mode.ne.With_Fragment_Mode .and.
     &    Basis_Mode.ne.All_Mode) Then
         Call WarningMessage(2,'Def_Shells: Basis_Mode is not defined')
         Call Abend()
      End If
*                                                                      *
************************************************************************
*                                                                      *
      iIrrep=0
      nSkal=0
      iAOttp=0 ! Number of AO functions proceeding a particular shell
      S%m2Max=0
*
      If (Atomic) Go To 300
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Molecular setup                                                  *
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      iCnttp = 0
      mdc=0
      iShell = 0
      Do jCnttp = 1, nCnttp
*
*        Make sure that we process the dummy shell last
*
         If (jCnttp.eq.iCnttp_Dummy .and. jCnttp.ne.nCnttp) Then
            iCnttp = iCnttp + 2
         Else If (jCnttp.eq.nCnttp .and. iCnttp.eq.jCnttp) Then
            iCnttp = iCnttp_Dummy
         Else
            iCnttp = iCnttp + 1
         End If
*
         nTest = dbsc(iCnttp)%nVal-1
         mdci = dbsc(iCnttp)%mdci
         Do iCnt = 1, dbsc(iCnttp)%nCntr
            mdci = mdci + 1
            mdc  = mdc  + 1
            iShell_Set = iShell
*
            Do iAng=0, nTest
               iShell = iShell + 1
               iShll = dbsc(iCnttp)%iVal + iAng
               nExpi=Shells(iShll)%nExp
               nBasisi=Shells(iShll)%nBasis
               If (Shells(iShll)%Prjct) Then
                  iCmp = 2*iAng+1
               Else
                  iCmp  = (iAng+1)*(iAng+2)/2
               End If
               If (nExpi.eq.0)   Go To 200
               If (nBasisi.eq.0) Go To 200
               If (Basis_Mode.eq.Valence_Mode .and.
     &             (Shells(iShll)%Aux.or.
     &              Shells(iShll)%Frag)) Go To 200
               If (Basis_Mode.eq.Auxiliary_Mode .and.
     &             .Not.Shells(iShll)%Aux) Go To 200
               If (Basis_Mode.eq.Fragment_Mode .and.
     &             .Not.Shells(iShll)%Frag) Go To 200
               If (Basis_Mode.eq.With_Auxiliary_Mode .and.
     &             Shells(iShll)%Frag) Go To 200
               If (Basis_Mode.eq.With_Fragment_Mode .and.
     &             Shells(iShll)%Aux) Go To 200

               kSh=dbsc(iCnttp)%iVal+iAng
*
               nSkal = nSkal + 1
*
************************************************************************
*
               iSD(0,nSkal)=iShll                    ! Unique shell ind.
               iSD(1,nSkal)=iAng                     ! l value
               iSD(2,nSkal)=iCmp                     ! # of ang. comp.
               iSD(3,nSkal)=nBasisi                  ! # of cont. func.
               iSD(4,nSkal)= -1                      ! Not used
               iSD(5,nSkal)=  nExpi                  ! # of prim.
               iSD(6,nSkal)= -1                      ! Not used
               iSD(7,nSkal)= iAOttp                  !
     &                     + (iCnt-1)*dbsc(iCnttp)%lOffAO
     &                     + Shells(kSh)%kOffAO      !
               iSD(8,nSkal)= -1
               itemp=0                               !
               If (Shells(iShll)%Prjct ) itemp=itemp+1      !
               If (Shells(iShll)%Transf) itemp=itemp+2      !
               iSD(9,nSkal)=itemp                    ! sph., car., cont.
               iSD(10,nSkal)=mdci                    ! Center index
               iSD(11,nSkal)=iShell_Set + iAng + 1
               If (dbsc(iCnttp)%pChrg) Then
                  iSD(12,nSkal)= 1                   ! pseudo charge
               Else
                  iSD(12,nSkal)= 0                   ! pseudo charge
               End If
               iSD(13,nSkal)= iCnttp
               iSD(14,nSkal)= iCnt
*
               nDisp = IndDsp(mdci,iIrrep)
               iTmp=0
               Do iCar = 0, 2
                  iComp = 2**iCar
                  If (TF(mdci,iIrrep,iComp).and.
     &                .Not.dbsc(iCnttp)%pChrg) Then
                     nDisp = nDisp + 1
                     If (Direct(nDisp)) Then
                        iSD(iCar+16,nSkal) = nDisp
                        iTmp=iOr(iTmp,2**iCar)
                     Else
                        iSD(iCar+16,nSkal) = 0
                     End If
                  Else
                     iSD(iCar+16,nSkal) = 0
                  End If
               End Do
               iSD(15,nSkal) = iTmp
*
               S%m2Max=Max(S%m2Max,nExpi**2)
 200           Continue
*                                                                      *
************************************************************************
*                                                                      *
            End Do                       ! iAng
         End Do                          ! iCnt
         iAOttp = iAOttp + dbsc(iCnttp)%lOffAO*dbsc(iCnttp)%nCntr
      End Do                             ! iCnttp
*
      Return
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Atomic setup                                                     *
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
 300  Continue
*
      iCase=1
      mdci = 1
      iCnt = 1
      nFunctions=0

 301  Continue
      If (iCase.eq.1) Then
*
*        Add the auxiliary basis
*
         iCnttp = kCnttp
      Else
*
*        Add the dummy basis
*
         iCnttp = iCnttp_Dummy
      End If
      nTest = dbsc(iCnttp)%nVal-1
*
      Do 400 iAng=0, nTest
         iShll = dbsc(iCnttp)%iVal + iAng
         nExpi=Shells(iShll)%nExp
         If (nExpi.eq.0)   Cycle
         nBasisi=Shells(iShll)%nBasis
         If (nBasisi.eq.0) Go To 400
         If (Shells(iShll)%Frag) Go To 400
         iCmp  = (iAng+1)*(iAng+2)/2
         If (Shells(iShll)%Prjct ) iCmp = 2*iAng+1
         kSh=dbsc(iCnttp)%iVal+iAng
*
         nSkal = nSkal + 1
*
************************************************************************
*
         iSD(0,nSkal)=iShll                    ! Unique shell ind.
         iSD(1,nSkal)=iAng                     ! l value
         iSD(2,nSkal)=iCmp                     ! # of ang. comp.
         iSD(3,nSkal)=nBasisi                  ! # of cont. func.
         iSD(4,nSkal)= -1                      ! Not used
         iSD(5,nSkal)=  nExpi                  ! # of prim.
         iSD(6,nSkal)= -1                      ! Not used
         iSD(7,nSkal)= iAOttp                  !
     &               + Shells(kSh)%kOffAO      !
         iSD(8,nSkal)= -1                      ! Not used
         itemp=0                               !
         If (Shells(iShll)%Prjct ) itemp=itemp+1      !
         If (Shells(iShll)%Transf) itemp=itemp+2      !
         iSD(9,nSkal)=itemp                    ! sph., car., cont.
         iSD(10,nSkal)=mdci                    ! Center index
         iSD(11,nSkal)= iAng + 1               ! Not used
         If (dbsc(iCnttp)%pChrg) Then
            iSD(12,nSkal)= 1                   ! pseudo charge
         Else
            iSD(12,nSkal)= 0                   ! pseudo charge
         End If
         iSD(13,nSkal)= iCnttp
         iSD(14,nSkal)= 1
*
         iSD(15,nSkal) = 0
         iSD(16,nSkal) = 0
         iSD(17,nSkal) = 0
         iSD(18,nSkal) = 0
*
         S%m2Max=Max(S%m2Max,nExpi**2)
*
         If (Shells(iShll)%Prjct ) Then
            nFunctions = nFunctions + nBasisi*(2*iAng+1)
         Else
            nFunctions = nFunctions + nBasisi*(iAng+1)*(iAng+2)/2
         End If
*                                                                      *
************************************************************************
*                                                                      *
 400  Continue                     ! iAng
*
      iCase=iCase+1
      If (iCase.le.2.and.dbsc(iCnttp)%Aux) Go To 301
*
      If (dbsc(iCnttp)%Aux) Then
         nBas(0)=0
      Else
         nBas(0)=nFunctions
      End If
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Write(6,*) 'in Define_Shells...'
      Do i = 1, mSkal
         Write (6,*) 'i=',i
         Write (6,'(10I8,/,8I8)') (iSD(j,i),j=0,nSD)
      End Do
#endif
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Return
      End
