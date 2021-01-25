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
      Subroutine Get_DEcorr(nh1,Grad,nGrad,DFTFOCK)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "mxdm.fh"
#include "addr.fh"
#include "WrkSpc.fh"
#include "infscf.fh"
      Real*8  Grad(nGrad), Ec_AB(2)
      Character*4 DFTFOCK
      Character*16  ADDC_KSDFT
      COMMON  / ADDcorr_C   / ADDC_KSDFT
      COMMON  / ADDcorr_R   / DE_KSDFT_c
      Logical Do_SpinAV
      COMMON  / SPAVE_L  / Do_SpinAV
      COMMON  / SPAVE_I  / ip_DSc
*
      Call GetMem('F-DS','Allo','Real',ipF_DFT,2*nBT)
      Call GetMem('D-DS','Allo','Real',ip_D_DS,2*nBT)
      ip_Da=ip_D_DS
      ip_Db=ip_D_DS+nBT
*
      Do iAB=1,2
       iOff=0
       jOff=0
       lOff=0
       Do iSym=1,nSym
          ipDaa=ip_Da+jOff
          If (iAB.eq.1) Then
             nXoX=nOcc(iSym,1)
             iXoX=0
          Else
             nXoX=nConstr(iSym)
             iXoX=nOcc(iSym,1)-nConstr(iSym)
          EndIf
          mAdCMOO=mAdCMO+iOff+nBas(iSym)*iXoX
          Call DGEMM_tri('N','T',nBas(iSym),nBas(iSym),nXoX,
     &                     1.0d0,Work(mAdCMOO),nBas(iSym),
     &                           Work(mAdCMOO),nBas(iSym),
     &                     0.0d0,Work(ipDaa),nBas(iSym))
          ipDbb=ip_Db+jOff
          If (iAB.eq.1) Then
             nXoX=nOcc(iSym,2)
             iXoX=0
          Else
             nXoX=nConstr(iSym)
             iXoX=nOcc(iSym,2)-nConstr(iSym)
          EndIf
          mAdCMOO=mAdCMO_ab+iOff+nBas(iSym)*iXoX
          Call DGEMM_tri('N','T',nBas(iSym),nBas(iSym),nXoX,
     &                     1.0d0,Work(mAdCMOO),nBas(iSym),
     &                           Work(mAdCMOO),nBas(iSym),
     &                     0.0d0,Work(ipDbb),nBas(iSym))
*
          If (Do_SpinAV) Then
             Do j=1,nBas(iSym)
                Do i=1,j
                   iDSc=ip_DSc-1+nBas(iSym)*(j-1)+i
                   ji=j*(j-1)/2+i
                   iDaa=ipDaa-1+ji
                   Work(iDaa)=Work(iDaa)-Work(iDSc)
                   iDbb=ipDbb-1+ji
                   Work(iDbb)=Work(iDbb)+Work(iDSc)
                End Do
             End Do
             lOff=lOff+nBas(iSym)**2
          EndIf
*
          Do j=1,nBas(iSym)
             Do i=1,j-1
                ji=j*(j-1)/2+i
                iDaa=ipDaa-1+ji
                Work(iDaa)=2.0d0*Work(iDaa)
                iDbb=ipDbb-1+ji
                Work(iDbb)=2.0d0*Work(iDbb)
             End Do
          End Do
          iOff=iOff+nBas(iSym)*nOrb(iSym)
          jOff=jOff+nBas(iSym)*(nBas(iSym)+1)/2
       End Do
*
       Call Get_Ecorr_dft(nh1,Grad,nGrad,DFTFOCK,ipF_DFT,ip_D_DS,
     &                        ADDC_KSDFT,Ec_AB(iAB))
      End Do
*----------------------------------------------------------------------*
      DE_KSDFT_c=Ec_AB(1)-Ec_AB(2)
*----------------------------------------------------------------------*
*
      Call GetMem('D-DS','Free','Real',ip_D_DS,2*nBT)
      Call GetMem('F-DS','Free','Real',ipF_DFT,2*nBT)
      Return
      End
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine Get_Ecorr_dft(nh1,Grad,nGrad,DFTFOCK,ipF_DFT,ip_D_DS,
     &                             KSDFT,Ec_AB)
      Implicit Real*8 (a-h,o-z)

#include "real.fh"
#include "WrkSpc.fh"
#include "nq_info.fh"
#include "debug.fh"
      Real*8  Grad(nGrad)
      Logical Do_MO,Do_TwoEl,Do_Grad
      Character*4 DFTFOCK
      Character*16  KSDFT
      COMMON  / OFembed_R2/ dFMD
      External VWN_III_emb,
     &         VWN_V_emb,
     &         cBLYP_emb,
     &         cPBE_emb,
     &         Checker
*
      lKSDFT=LEN(KSDFT)
      Debug=.False.
*                                                                      *
************************************************************************
*                                                                      *
*     DFT functionals, compute integrals over the potential
*
      Func            =Zero
      Dens_I          =Zero
      Grad_I          =Zero
      Tau_I           =Zero
      Do_MO           =.False.
      Do_TwoEl        =.False.
      Do_Grad=.false.
*
      nFckDim=2
      nD=2
      dFMD_=dFMD
      dFMD=1.0d0
*                                                                      *
************************************************************************
*                                                                      *
*      LSDA LDA SVWN
*
       If (KSDFT.eq.'LSDA ' .or.
     &     KSDFT.eq.'LDA '  .or.
     &     KSDFT.eq.'SVWN ') Then
c        ExFac=Get_ExFac(KSDFT)
         Functional_type=LDA_type
         Call DrvNQ(VWN_III_emb,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*      LSDA5 LDA5 SVWN5
*
       Else If (KSDFT.eq.'LSDA5' .or.
     &          KSDFT.eq.'LDA5'  .or.
     &          KSDFT.eq.'SVWN5') Then
c        ExFac=Get_ExFac(KSDFT)
         Functional_type=LDA_type
         Call DrvNQ(VWN_V_emb,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     BLYP
*
      Else If (KSDFT.eq.'BLYP') Then
c        ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Call DrvNQ(cBLYP_emb,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     PBE
*
      Else If (KSDFT.eq.'PBE') Then
c        ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Call DrvNQ(cPBE_emb,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     Checker
*
      Else If (KSDFT.eq.'CHECKER') Then
c        ExFac=Zero
         Functional_type=meta_GGA_type2
         Call DrvNQ(Checker,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
      Else
         Call WarningMessage(2,
     &               ' Get_Ecorr_dft: Unsupported functional type!')
         Write (6,*) '         Functional=',KSDFT(1:lKSDFT)
         Call Quit_OnUserError()
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Ec_AB=Func
*
#ifdef _DEBUGPRINT_
      write(6,*) ' Correlation energy: ',Ec_AB
      write(6,*)
      write(6,*) ' Correlation potentials: (itri,F_alpha,F_beta)'
      write(6,*)
      Do i=1,nh1
        Write(6,'(i4,3f22.16)') i,Work(ipF_DFT+i-1),
     &                            Work(ipF_DFT+i-1+nh1)
      End Do
#endif
*
      dFMD=dFMD_
      Return
      End
