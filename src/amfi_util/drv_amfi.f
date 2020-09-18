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
      Subroutine Drv_AMFI(Label,ip,lOper,nComp,rHrmt,iChO, iAtmNr2,
     &                    Charge2)
      use iSD_data
      use Basis_Info
      Implicit Real*8 (a-h,o-z)
      External Rsv_Tsk
#include "angtp.fh"
#include "info.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "nsd.fh"
#include "setup.fh"
#include "para.fh"
      Integer, Allocatable :: iDel(:)
      Real*8, Allocatable :: SOInt(:)
      Real*8 Coor(3)
      Logical EQ, IfTest, Rsv_Tsk
      Character Label*8
      Integer ip(nComp), lOper(nComp), iChO(nComp)
      Integer iAtmNr2(mxdbsc)
      Real*8 Charge2(mxdbsc)
      Data IfTest/.False./
*
!#define _DEBUG_
#ifdef _DEBUG_
      IfTest=.True.
      Write (6,*) ' In OneEl: Label', Label
      Write (6,*) ' In OneEl: nComp'
      Write (6,'(1X,8I5)') nComp
      Write (6,*) ' In OneEl: lOper'
      Write (6,'(1X,8I5)') lOper
      Write (6,*) ' In OneEl: n2Tri'
      Do iComp = 1, nComp
         ip(iComp) = n2Tri(lOper(iComp))
      End Do
      Write (6,'(1X,8I5)') (ip(iComp),iComp=1,nComp)
#endif
*
      Eta_Nuc=Zero
*
*     Allocate memory for symmetry adapted one electron integrals.
*     Will just store the unique elements, i.e. low triangular blocks
*     and lower triangular elements in the diagonal blocks.
*
      ip(:)=-1
      LenTot=0
      Do iComp = 1, nComp
         ip(iComp)=1+LenTot
         LenInt=n2Tri(lOper(iComp))
         LenTot=LenTot+LenInt+4
      End Do
      Call mma_allocate(SOInt,LenTot,label='SOInt')
      SOInt(:)=Zero
*
*---- Generate list of shell information
*
      Call Nr_Shells(nSkal)
*
*---- Check that there are not several instances of the same center.
*
      nCenter=0
      Do iSkal=1,nSkal
         nCenter=Max(nCenter,iSD(10,iSkal))
         If (iSD(1,iSkal).gt.Lmax) Then
            Write (6,*) ' Shells higher than '
     &                 //Angtp(Lmax)//'-functions not allowed in AMFI.'
            Call Quit_OnUserError()
         End If
         If (iSD(1,iSkal).ge.2.and.iAnd(iSD(9,iSkal),1).ne.1) Then
            Write (6,*) ' Only real spherical harmonics allowed'
            Write (6,*) ' for AMFI.'
            Call Quit_OnUserError()
         End If
      End Do
      Coor(:)=Zero
      Do iCenter=1,nCenter
*
*        Identify which dbsc this center belongs.
*
         iCnttp=0
         Do iSkal=1,nSkal
            If (iSD(10,iSkal).eq.iCenter) Then
               iCnttp=iSD(13,iSkal)
               iCnt  =iSD(14,iSkal)
               Coor(1:3)=dbsc(iCnttp)%Coor(1:3,iCnt)
            End If
         End Do
*
*        If iCenter is not a center in the current list of shells that
*        is to be processed then test the next iCenter.
*
         If (iCnttp.eq.0) Cycle
*
         Do iSkal=1,nSkal
            If (iSD(13,iSkal).ne.iCnttp .and.
     &          iSD(10,iSkal).ne.iCenter) Then
               jCnttp=iSD(13,iSkal)
               jCnt  =iSD(14,iSkal)
               If (  EQ(Coor, dbsc(jCnttp)%Coor(1,jCnt)) ) Then
                  Write (6,*) 'Multiple instances of the same center!'
                  Write (6,*) 'This is not allowed with AMFI.'
                  Call Quit_OnUserError()
               End If
            End If
         End Do
      End Do
      If (MolWgh.ne.0 .and. MolWgh.ne.2) Then
         Write (6,*) ' AMFI integrals not implemented for symmetry'
     &                //' adaptation a la MOLECULE'
         Call Quit_OnUserError()
      End If
*
*---- Loop over unique center. Observe that multiple shells of the same
*     angular momentum is not allowed.
*
      Lu_AMFI=21
      LUPROP=22
      call molcas_open(Lu_AMFI,'AMFI_INP')
      call molcas_binaryopen_vanilla(LUPROP,'AMFI_INT')
      nCenter_node=0
*
      Call Init_Tsk(id_Tsk,nCenter)
 10   Continue
      If (.Not.Rsv_Tsk(id_Tsk,iCenter)) Go To 11
*     Do iCenter = 1, nCenter
         nCenter_node=nCenter_node+1
*
         Write (Lu_AMFI,'(A)') ' &AMFI &END'
         If (IfTest) Write (6,'(A)') ' &AMFI &END'
*
*------- Find atom type
*
         mdci=0
         Do iCnttp = 1, nCnttp
            Do iCnt = 1, dbsc(iCnttp)%nCntr
               mdci=mdci+1
               If (mdci.eq.iCenter) Then
               If ((.Not.DKroll).and.(.Not.dbsc(iCnttp)%SODK))
     &            Write (Lu_AMFI,'(A)') 'Breit-Pauli'
                  If (IfTest) Then
                  If (.Not.DKroll.and..Not.dbsc(iCnttp)%SODK)
     &            Write (6,'(A)') 'Breit-Pauli'
               End If
                  If (iAtmNr2(iCnttp).ge.1) Then
                     charge_x=DBLE(iAtmNr2(iCnttp))
                  Else If (iAtmNr2(iCnttp).le.0 .and.
     &                Charge2(iCnttp).eq.Zero) Then
                     charge_x=0.0D0
                  Else
                     Write (6,*) 'Drv_AMFI: Invalid basis!'
                     Write (6,*) 'iAtmNr=',iAtmNr2(iCnttp)
                     Write (6,*) 'Charge2=',Charge2(iCnttp)
                     Call Abend()
                  End If
                  If (Nuclear_Model.eq.Gaussian_Type) Then
                     Eta_Nuc=dbsc(iCnttp)%ExpNuc
                  End If
                  Go To 99
               End If
            End Do
         End Do
 99      Continue
*
         If (Nuclear_Model.eq.Gaussian_Type) Then
            Write (Lu_AMFI,'(A)') 'Finite'
            If (IfTest) Write (6,'(A)') 'Finite'
            Write (Lu_AMFI,*) Eta_Nuc
            If (IfTest) Write (6,*) Eta_Nuc
         End If
*
*------- Generate input for each atom
*
         l_max=-1
         Do iSkal = 1, nSkal
            If (iSD(10,iSkal).eq.iCenter)
     &         l_Max=Max(l_Max,iSD(1,iSkal))
         End Do
         If (l_max.gt.LMax) Then
            Write (6,*) 'AMFI integrals only implemented up to '
     &                 //Angtp(Lmax)//'-functions.'
            Call Quit_OnUserError()
         End If
*
*        Check if there are any core orbitals to be deleted
*
         nCore=0
         Do l = 0, l_max
            Do iSkal = 1, nSkal
               If (iSD(10,iSkal).eq.iCenter .and.
     &             iSD( 1,iSkal).eq.l) Then
*
                   iCnttp = iSD(13,iSkal)
                   If (dbsc(iCnttp)%nSOC.ne.0) Then
                      iShl  = dbsc(iCnttp)%iSOC+l
                      jShll  = dbsc(iCnttp)%iPrj+l
*                     nCore=nCore+Shells(jShll)%nBasis
                      nCore=nCore+dbsc(iCnttp)%kDel(l)
                   End If
               End If
            End Do
         End Do
         If (IfTest) Write (6,*) 'nCore: ', nCore
*
*        Set up delete array
*
         If (nCore.ne.0) Then
            lDel=l_Max+1
            Call mma_allocate(iDel,lDel,label='iDel')
            Do l = 0, l_max
               Do iSkal = 1, nSkal
                  If (iSD(10,iSkal).eq.iCenter .and.
     &                iSD( 1,iSkal).eq.l) Then
*
                      iCnttp = iSD(13,iSkal)
                      If (dbsc(iCnttp)%nSOC.ne.0) Then
                         iShll  = dbsc(iCnttp)%iSOC+l
                         jShll  = dbsc(iCnttp)%iPrj+l
*                        iDel(ip_iDel+l)=Shells(jShll)%nBasis
                         iDel(1+l)=dbsc(iCnttp)%kDel(l)
                      End If
                  End If
               End Do
            End Do
            Write (Lu_AMFI,'(A)') 'AIMP'
            Write (Lu_AMFI,*) lDel-1,(iDel(i),i=1,lDel)
            If (IfTest) Write (6,'(A)') 'AIMP'
            If (IfTest) Write (6,*) lDel, (iDel(i),i=1,lDel)
            Call mma_deallocate(iDel)
         End If
*
         Write (Lu_AMFI,'(A)') '     '
         Write (Lu_AMFI,'(3X,F5.1,I4)') charge_x, l_max
         If (IfTest) Write (6,*) charge_x, l_max
*
         Do l = 0, l_max
            Do iSkal = 1, nSkal
               If (iSD(10,iSkal).eq.iCenter .and.
     &             iSD( 1,iSkal).eq.l) Then
*
                   iCnttp = iSD(13,iSkal)
                   If (dbsc(iCnttp)%nSOC.eq.0) Then
*
*                     Use valence basis
*
                      iShll  = dbsc(iCnttp)%iVal+l
                      iCase  = 2
*
                   Else
*
*                     Use special valence basis in case of a ECP where the
*                     normal valence might not be adequate.
*
                      iShll  = dbsc(iCnttp)%iSOC+l
                      nBas_y = iSD( 3,iSkal)
                      iCase = 1
*
                   End If
                   nBas_x = Shells(iShll)%nBasis
                   nExp_x = Shells(iShll)%nExp
*
                   If (IfTest) Write (6,*)  'iShll=',iShll
                   Write (Lu_AMFI,*) nExp_x, nBas_x
                   If (IfTest) Write (6,*) nExp_x, nBas_x
                   Write (Lu_AMFI,*) (Shells(iShll)%Exp(iExp_x),
     &                                  iExp_x=1,nExp_x)
                   If (IfTest) Write (6,*) (Shells(iShll)%Exp(iExp_x),
     &                                       iExp_x=1,nExp_x)
                   Do iExp_x = 1, nExp_x
                      Write (Lu_AMFI,*)
     &                     (Shells(iShll)%Cff_c(iExp_x,iCff_x,iCase),
     &                              iCff_x=1,nBas_x)
                      If (IfTest) Write (6,*)
     &                     (Shells(iShll)%Cff_c(iExp_x,iCff_x,iCase),
     &                              iCff_x=1,nBas_x)
                   End Do
*
               End If
            End Do
         End Do
         Write (Lu_AMFI,'(A)') 'End of Input'
         If (IfTest) Write (6,'(A)') 'End of Input'
*
*
*------- Now call AMFI
*
         Rewind(Lu_AMFI)
         Call AMFI(Lu_AMFI,LUPROP,iCenter)
*
         Rewind(Lu_AMFI)
*
         Go To 10
 11   Continue
      Call Free_Tsk(id_Tsk)
*     End Do
      Close(Lu_AMFI)
*
*---- Now symmetry adopt.
*
      Rewind(LUPROP)
      Call SymTrafo(LUPROP,ip,lOper,nComp,nBas,nIrrep,Label,MolWgh,
     &              SOInt,LenTot)
*
      Close(LUPROP)
*
      Call mma_deallocate(SOInt)
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real(rHrmt)
         Call Unused_integer_array(iChO)
      End If
      End
