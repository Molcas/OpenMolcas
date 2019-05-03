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
      Subroutine Local_Properties(Coor,nAtoms,ip_sq_mu,nElem,
     &                            Sq_Temp,Origin,iCenter,Ttot_Inv,
     &                            Temp,nij,nPert,ip_D,rMP,lMax,rMPq,
     &                            C_o_C,EC,iANr,Standard,nBas1,nTemp,
     &                            Q_Nuc,Bond_Threshold,Utility,
     &                            Opt_Method,iPlot,iPrint,nSym)
      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 Coor(3,nAtoms), A(3), B(3), Sq_Temp(nTemp),
     &        Origin(3,0:lMax), Ref(3), C_o_C(3), EC(3,nij),
     &       Ttot_Inv(nBas1**2), Temp(nTemp), Q_Nuc(nAtoms),
     &       rMP(nij,0:nElem-1,0:nPert-1), rMPq(0:nElem-1)
      Integer ip_sq_Mu(0:nElem-1), iCenter(nBas1), ip_D(0:6),
     &        iANr(nAtoms)
      Logical Standard,Utility,Reduce_Prt
      External Reduce_Prt
      Character*12  Opt_Method

#include "Molcas.fh"
      CHARACTER*(LENIN) CNAME(MXATOM)
      integer tNuc
      Data Ref/0.0D0,0.0D0,0.0D0/
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
      Call Set_Binom()
*                                                                      *
************************************************************************
*                                                                      *
*---- Loop over perturbations
*
      Do iPert = 0, nPert-1
*
*        Transform the density matrix to the LoProp basis
*
         iOffD=ip_D(iPert)-1
         Do i = 1, nBas1
            Do j = 1, i-1
               ij = i * (i-1)/2 + j
               ij_ = (j-1)*nBas1+i
               ji_ = (i-1)*nBas1+j
               Sq_Temp(ij_)=Half*Work(ij+iOffD)
               Sq_Temp(ji_)=Sq_Temp(ij_)
            End Do
            ii = i*(i+1)/2
            ii_= (i-1)*nBas1+i
            Sq_Temp(ii_)=Work(ii+iOffD)
         End Do
*
         Call DGEMM_('N','T',
     &               nBas1,nBas1,nBas1,
     &               1.0d0,Sq_Temp,nBas1,
     &               Ttot_Inv,nBas1,
     &               0.0d0,Temp,nBas1)
         Call DGEMM_('N','N',
     &               nBas1,nBas1,nBas1,
     &               1.0d0,Ttot_Inv,nBas1,
     &               Temp,nBas1,
     &               0.0d0,Sq_Temp,nBas1)
*
cvv
      Call Get_iScalar('Unique atoms',nNUC)
      Call Get_cArray('Unique Atom Names',CNAME,(LENIN)*nNuc)

      Call Get_iScalar('LP_nCenter', tNuc)
c someday this code will use symmetry
*
      iPL=iPrintLevel(-1)
      If (Reduce_Prt().and.iPL.lt.3) iPL=0
*
      if(nSym.eq.1.and.tNuc.eq.nNuc.and.iPL.ge.2) then

       NBAST=nBas1
      Call Allocate_iWork(ip_center,NBAST)
      Call Get_iArray('Center Index',iWork(ip_center),NBAST)
      Call Allocate_iWork(ipNBFpA,tNUC)
      Do I = 1, tNUC
          iWork(ipNBFpA+I-1) = 0
      End Do
      Do I = 1, NBAST
          iWork(ipNBFpA+iWork(ip_center+I-1)-1) =
     &     iWork(ipNBFpA+iWork(ip_center+I-1)-1) + 1
      End Do
        Call Allocate_Work(ip_Charge,tNuc)
      Call Get_dArray('Effective nuclear charge',
     &      Work(ip_Charge),tNuc)
#ifdef VV_VALE
      Call Allocate_Work(ip_VV_QAB,tNuc)
      Call Allocate_Work(ip_VV_WAB,tNuc*tNuc)
      Call Allocate_Work(ip_VV_CAB,tNuc)
      Call Allocate_Work(ip_VV_VAB,tNuc)
      Call Allocate_iWork(ip_VV_N,tNuc+1)

       call VALE(tNuc,iWork(ipNBFpA),Work(ip_Charge),
     &           nBast,Sq_Temp,
     &     CName,LenIn,iWork(ip_VV_N),
     &     Work(ip_VV_QAB),Work(ip_VV_WAB),
     &     Work(ip_VV_CAB), Work(ip_VV_VAB))


      Call Free_Work(ip_VV_N)
      Call Free_Work(ip_VV_VAB)
      Call Free_Work(ip_VV_CAB)
      Call Free_Work(ip_VV_WAB)
      Call Free_Work(ip_VV_QAB)
#endif
      Call Free_Work(ip_Charge)
      Call Free_Work(ipNBFpA)
      Call Free_Work(ip_center)
      endif
cvv
         iMu = -1
         Do l = 0, lMax
         Do ix = l, 0, -1
         Do iy = l-ix, 0, -1
            iz = l-ix-iy
            iMu = iMu + 1
*
*....... Compute local properties
*
         Do iAtom = 1, nAtoms
            call dcopy_(3,Coor(1,iAtom),1,A,1)
            Do jAtom = 1, iAtom
               call dcopy_(3,Coor(1,jAtom),1,B,1)
*
*
*              Sum up contibutions to the domain ij
*
               Acc=Zero
               iOffO=ip_sq_mu(iMu)-1
               Do j = 1, nBas1
                  Do i = 1, nBas1
                     If ((iCenter(i).eq.iAtom .and.
     &                    iCenter(j).eq.jAtom) .or.
     &                   (iCenter(i).eq.jAtom .and.
     &                    iCenter(j).eq.iAtom)) Then
                        ij = (j-1)*nBas1+i
                        Acc = Acc + Sq_Temp(ij)*Work(ij+iOffO)
                     End If
                  End Do
               End Do
               ij=iAtom*(iAtom-1)/2+jAtom
               rMP(ij,iMu,iPert)=-Acc
*
            End Do   ! jAtom
            End Do   ! iAtom
*
         End Do   ! iy
         End Do   ! ix
         End Do   ! l
*                                                                      *
************************************************************************
*                                                                      *
         Call Free_Work(ip_D(iPert))
      End Do         ! iPert
*                                                                      *
************************************************************************
*                                                                      *
*     Modify all multipole moments
*
*     1) Modify first all to a common origin (0.0,0.0,0.0)
*     2) Modify to local moments
*
#ifdef _DEBUG_
      Call xSpot ('Middle  Local_Properties')
      Call RecPrt('rMP',' ',rMP,nij,nElem*nPert)
#endif
      Do iAtom = 1, nAtoms
         ii = iAtom*(iAtom+1)/2
*
*        Observe the loop order to make sure that element (ii) and
*        (jj) always are processed before (ij)!
*
         Do jAtom = iAtom, 1, -1
            ij = iAtom*(iAtom-1)/2+jAtom
            jj = jAtom*(jAtom+1)/2
C           Write (*,*) 'ij=',ij
*
*           Step 1
*
            Do l = 1, lMax-1
               mElem = (l+1)*(l+2)*(l+3)/6
               Do iPert = 0, nPert-1
                  Call ReExpand(rMP(1,0,iPert),nij,mElem,Origin(1,l),
     &                          Origin(1,l+1),ij,l)
               End Do
               If (ij.eq.1) Call ReExpand(rMPq,1,mElem,Origin(1,l),
     &                                    Origin(1,l+1),1,l)
            End Do
C           Write (6,*) 'End step 1a'
            mElem = (lMax+1)*(lMax+2)*(lMax+3)/6
            Do iPert = 0, nPert-1
               Call ReExpand(rMP(1,0,iPert),nij,mElem,Origin(1,lMax),
     &                       Ref,ij,lMax)
            End Do
            If (ij.eq.1) Call ReExpand(rMPq,1,mElem,Origin(1,lMax),
     &                                 C_o_C,1,lmax)
C           Write (6,*) 'End step 1b'
*
*           Step 2
*
*           Establish the initial expansion centers for each domain.
*
*           Default values
*
            If (iAtom.eq.jAtom) Then
               A(1)=Coor(1,iAtom)
               A(2)=Coor(2,iAtom)
               A(3)=Coor(3,iAtom)
            Else
               A(1) = (EC(1,ii)+EC(1,jj))*Half
               A(2) = (EC(2,ii)+EC(2,jj))*Half
               A(3) = (EC(3,ii)+EC(3,jj))*Half
            End If
            call dcopy_(3,A,1,EC(1,ij),1)
            Do iPert = 0, nPert-1
               Call ReExpand(rMP(1,0,iPert),nij,mElem,Ref,A,ij,lMax)
            End Do
C           Write (6,*) 'End step 2'
*
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
* Distributes the contributions from the bonds that doesn't fulfill the requirement
* Bond_Length <= Bond_Threshold*(Bragg_Slater(iAtom)+Bragg_Slater(jAtom) to the
* two atoms involved in the bond.
*
      Call Move_Prop(rMP,EC,lMax,nElem,nAtoms,nPert,
     &               nij,iANr,Bond_Threshold)
*
*           Modify the expansion centers
*
      Num_Warnings = 0
      Call Allocate_Work(iT_values,nij)
      Call Allocate_iWork(iT_sets,nij)
      Call Allocate_iWork(iWarnings,nij)
      Call iCopy(nij,[0],0,iWork(iWarnings),1)
      Call iCopy(nij,[0],0,iWork(iT_Sets),1)
      Call dCopy_(nij,[Zero],0,Work(iT_Values),1)
      If (.Not. Standard) Then
         Call Allocate_Work(iScratch_1,nij*(2*lMax+1))
         Call Allocate_Work(iScratch_2,nij*(2*lMax+1))
         Call Allocate_Work(iScratch_A,3*nij)
         Call Allocate_Work(iScratch_B,3*nij)
         Call Allocate_Work(ixrMP,nij*nElem)
         Call Allocate_Work(ixxrMP,nij*nElem)
         Call Allocate_Work(ixnrMP,nij*nElem)
*
         Call Move_EC(rMP,EC,Work(iScratch_1),Work(iScratch_2),
     &                Work(ixrMP),Work(ixxrMP),Work(ixnrMP),lMax,
     &                Work(iScratch_A),Work(iScratch_B),nij,
     &                nElem,Coor,nAtoms,Q_Nuc,C_o_C,nPert,
     &                Bond_Threshold,iANr,Work(iT_Values),
     &                iWork(iT_Sets),iWork(iWarnings),Num_Warnings,
     &                Opt_Method,iPlot,iPrint)
*
         Call Free_Work(ixnrMP)
         Call Free_Work(ixxrMP)
         Call Free_Work(ixrMP)
         Call Free_Work(iScratch_B)
         Call Free_Work(iScratch_A)
         Call Free_Work(iScratch_2)
         Call Free_Work(iScratch_1)
      End If
      If (iPrint .ge. 1 .OR. iPlot .ge. 1 .OR. Num_Warnings .gt. 0) Then
         Call Print_T_Values(Work(iT_Values),iWork(iT_Sets),iANr,EC,
     &                       Bond_Threshold,nAtoms,nij,Standard,
     &                       iWork(iWarnings),Num_Warnings,iPrint)
      End If
      Call Free_iWork(iWarnings)
      Call Free_iWork(iT_values)
      Call Free_iWork(iT_sets)
*
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
      Call RecPrt('EC',' ',EC,3,nij)
      Call RecPrt('rMP',' ',rMP,nij,nElem*nPert)
      Call RecPrt('rMPq',' ',rMPq,1,nElem)
      Call xSpot('Exit  Local_Properties')
#endif
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_logical(Utility)
      End
