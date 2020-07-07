************************************************************************
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
      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "info.fh"
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
      TF(mdc,iIrrep,iComp) = TstFnc(iOper,nIrrep,iCoSet(0,0,mdc),
     &                       nIrrep/nStab(mdc),iChTbl,iIrrep,iComp,
     &                       nStab(mdc))
      IndSOff(iCnttp,iCnt)=(iCnttp-1)*Max_Cnt+iCnt
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
      iAOttp=0
      m2Max=0
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
*
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
         nTest = nVal_Shells(iCnttp)-1
         mdci = mdciCnttp(iCnttp)
         Do iCnt = 1, dbsc(iCnttp)%nCntr
            mdci = mdci + 1
*
            Do 200 iAng=0, nTest
               iShll = ipVal(iCnttp) + iAng
               If (nExp(iShll).eq.0)   Go To 200
               If (nBasis(iShll).eq.0) Go To 200
               If (Basis_Mode.eq.Valence_Mode .and.
     &             (AuxShell(iShll).or.FragShell(iShll))) Go To 200
               If (Basis_Mode.eq.Auxiliary_Mode .and.
     &             .Not.AuxShell(iShll)) Go To 200
               If (Basis_Mode.eq.Fragment_Mode .and.
     &             .Not.FragShell(iShll)) Go To 200
               If (Basis_Mode.eq.With_Auxiliary_Mode .and.
     &             FragShell(iShll)) Go To 200
               If (Basis_Mode.eq.With_Fragment_Mode .and.
     &             AuxShell(iShll)) Go To 200
               iCmp  = (iAng+1)*(iAng+2)/2
               If (Prjct(iShll)) iCmp = 2*iAng+1
*
               nSkal = nSkal + 1
*
************************************************************************
*
               iSD(0,nSkal)=iShll                    ! Unique shell ind.
               iSD(1,nSkal)=iAng                     ! l value
               iSD(2,nSkal)=iCmp                     ! # of ang. comp.
               iSD(3,nSkal)=nBasis(iShll)            ! # of cont. func.
               iSD(4,nSkal)= ipCff(iShll)            ! pointer to coeff.
               iSD(5,nSkal)=  nExp(iShll)            ! # of prim.
               iSD(6,nSkal)= -1                      ! Not used
               iSD(7,nSkal)= iAOttp                  ! ? magic
     &                     + (iCnt-1)*lOffAO(iCnttp) !
     &                     + kOffAO(iCnttp,iAng)     !
               iSD(8,nSkal)= -1                      ! Not used
               itemp=0                               !
               If ( Prjct(iShll)) itemp=itemp+1      !
               If (Transf(iShll)) itemp=itemp+2      !
               iSD(9,nSkal)=itemp                    ! sph., car., cont.
               iSD(10,nSkal)=mdci                    ! Center index
               iSD(11,nSkal)=Ind_Shell(IndSOff(iCnttp,iCnt)) + iAng + 1
               If (pChrg(iCnttp)) Then
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
     &                .Not.pChrg(iCnttp)) Then
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
               m2Max=Max(m2Max,nExp(iShll)**2)
*                                                                      *
************************************************************************
*                                                                      *
 200        Continue                     ! iAng
         End Do                          ! iCnt
         iAOttp = iAOttp + lOffAO(iCnttp)*dbsc(iCnttp)%nCntr
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
      nTest = nVal_Shells(iCnttp)-1
*
      Do 400 iAng=0, nTest
         iShll = ipVal(iCnttp) + iAng
         If (nExp(iShll).eq.0)   Go To 400
         If (nBasis(iShll).eq.0) Go To 400
         If (FragShell(iShll)) Go To 400
         iCmp  = (iAng+1)*(iAng+2)/2
         If (Prjct(iShll)) iCmp = 2*iAng+1
*
         nSkal = nSkal + 1
*
************************************************************************
*
         iSD(0,nSkal)=iShll                    ! Unique shell ind.
         iSD(1,nSkal)=iAng                     ! l value
         iSD(2,nSkal)=iCmp                     ! # of ang. comp.
         iSD(3,nSkal)=nBasis(iShll)            ! # of cont. func.
         iSD(4,nSkal)= ipCff(iShll)            ! pointer to coeff.
         iSD(5,nSkal)=  nExp(iShll)            ! # of prim.
         iSD(6,nSkal)= -1                      ! Not used
         iSD(7,nSkal)= iAOttp                  ! ? magic
     &               + kOffAO(iCnttp,iAng)     !
         iSD(8,nSkal)= -1                      ! Not used
         itemp=0                               !
         If ( Prjct(iShll)) itemp=itemp+1      !
         If (Transf(iShll)) itemp=itemp+2      !
         iSD(9,nSkal)=itemp                    ! sph., car., cont.
         iSD(10,nSkal)=mdci                    ! Center index
         iSD(11,nSkal)=Ind_Shell(IndSOff(iCnttp,iCnt)) + iAng + 1 !
         If (pChrg(iCnttp)) Then
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
         m2Max=Max(m2Max,nExp(iShll)**2)
*
         If (Prjct(iShll)) Then
            nFunctions = nFunctions + nBasis(iShll)*(2*iAng+1)
         Else
            nFunctions = nFunctions + nBasis(iShll)*(iAng+1)*(iAng+2)/2
         End If
*                                                                      *
************************************************************************
*                                                                      *
 400  Continue                     ! iAng
*
      iCase=iCase+1
      If (iCase.le.2.and.AuxCnttp(iCnttp)) Go To 301
*
      If (AuxCnttp(iCnttp)) Then
         nBas(0)=0
      Else
         nBas(0)=nFunctions
      End If
*define _DEBUG_
#ifdef _DEBUG_
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
      Subroutine Define_Shells_kext(iSD,ikak,nSkal)
      use Basis_Info
      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "info.fh"
*
      Parameter(nSD=14)
      Integer iSD(0:nSD,1024)
*
*---- Statement function
*
      IndSOff(iCnttp,iCnt)=(iCnttp-1)*Max_Cnt+iCnt
*
      nSkal=0
      iAOttp=0
      Do 100 iAng=0, iAngMx
         If (MaxPrm(iAng).eq.0) Go To 100
         iAOttp=0
         Do 200 iCnttp = 1, nCnttp
            mdci = mdciCnttp(iCnttp)
            nTest = nVal_Shells(iCnttp)-1
            If (iAng.gt.nTest) Go To 201
            iShll = ipVal(iCnttp) + iAng
            If (nExp(iShll).eq.0) Go To 201
            If (nBasis(iShll).eq.0) Go To 201
            iCmp  = (iAng+1)*(iAng+2)/2
            If (Prjct(iShll)) iCmp = 2*iAng+1
            Do iCnt = 1, dbsc(iCnttp)%nCntr
               nSkal = nSkal + 1
*                                                                      *
************************************************************************
*                                                                      *
               iSD(0,nSkal)=iShll                    ! Unique shell ind.
               iSD(1,nSkal)=iAng                     ! l value
               iSD(2,nSkal)=iCmp                     ! # of ang. comp.
               iSD(3,nSkal)=nBasis(iShll)            ! # of cont. func.
               iSD(4,nSkal)= ipCff(iShll)            ! pointer to coeff.
               iSD(5,nSkal)=  nExp(iShll)            ! # of prim.
               iSD(6,nSkal)= -1                      ! Not used
               iSD(7,nSkal)= iAOttp                  ! ? magic
     &                     + (iCnt-1)*lOffAO(iCnttp) !
     &                     + kOffAO(iCnttp,iAng)     !
               iSD(8,nSkal)=-1                       ! Not used
               itemp=0                               !
               If ( Prjct(iShll)) itemp=itemp+1      !
               If (Transf(iShll)) itemp=itemp+2      !
               iSD(9,nSkal)=itemp                    ! sph., car., cont.
               iSD(10,nSkal)=mdci+iCnt               ! Center index
*              iSD(11,nSkal)=iSOff(iCnttp,iCnt)      !
               iSD(11,nSkal)=Ind_Shell(IndSOff(iCnttp,iCnt)) + iAng + 1!
               iSD(12,nSkal)= ipVal(iCnttp) + iAng     !
               iSD(13,nSkal)= iCnttp
               iSD(14,nSkal)= iCnt
*                                                                      *
************************************************************************
*                                                                      *
            End Do                       ! iCnt
 201        Continue
            iAOttp = iAOttp + lOffAO(iCnttp)*dbsc(iCnttp)%nCntr
 200     Continue                        ! iCnttp
 100  Continue                           ! iAng
*
*     The order of the shells could be reordered here!
*debugdebug
c     Write(6,*) 'in Define_Shells...'
c     Do i = 1, nSkal
c        Write (*,'(13I8)') (iSD(j,i),j=0,nSD)
c     End Do
*debugdebug
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(ikak)
      End
