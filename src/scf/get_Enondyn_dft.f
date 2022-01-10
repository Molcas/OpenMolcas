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
      Subroutine Get_Enondyn_dft(nh1,Grad,nGrad,DFTFOCK)
      use SCF_Arrays, only: CMO
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "mxdm.fh"
#include "WrkSpc.fh"
#include "infscf.fh"
      Real*8  Grad(nGrad)
      Character*4 DFTFOCK
#include "dcscf.fh"
*
      Erest_xc=0.0d0
      Call GetMem('F-DS','Allo','Real',ipF_DFT,2*nBT)
      Call GetMem('D-DS','Allo','Real',ip_D_DS,2*nBT)
      ip_Da=ip_D_DS
      ip_Db=ip_D_DS+nBT
*
      iOff=1
      jOff=0
      Do iSym=1,nSym
         ipDaa=ip_Da+jOff
         Call DGEMM_tri('N','T',nBas(iSym),nBas(iSym),nOcc(iSym,1),
     &                    1.0d0,CMO(iOff,1),nBas(iSym),
     &                          CMO(iOff,1),nBas(iSym),
     &                    0.0d0,Work(ipDaa),nBas(iSym))
         ipDbb=ip_Db+jOff
         Call DGEMM_tri('N','T',nBas(iSym),nBas(iSym),nOcc(iSym,2),
     &                    1.0d0,CMO(iOff,2),nBas(iSym),
     &                          CMO(iOff,2),nBas(iSym),
     &                    0.0d0,Work(ipDbb),nBas(iSym))
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
*----------------------------------------------------------------------*
      Call Get_Fmat_nondyn(Work(ip_Da),Work(ip_Db),nBT,.true.)
*----------------------------------------------------------------------*
*
*----------------------------------------------------------------------*
      Call Get_Exc_dft(nh1,Grad,nGrad,DFTFOCK,ipF_DFT,ip_D_DS,
     &                     KSDFT)
*----------------------------------------------------------------------*
*
      Call GetMem('D-DS','Free','Real',ip_D_DS,2*nBT)
      Call GetMem('F-DS','Free','Real',ipF_DFT,2*nBT)
      Return
      End
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine Get_Exc_dft(nh1,Grad,nGrad,DFTFOCK,ipF_DFT,ip_D_DS,
     &                           KSDFT)
      Implicit Real*8 (a-h,o-z)
      External LSDA, Overlap, BLYP, B3LYP, HFS, HFB,
     &         XAlpha, LSDA5, B3LYP5,TLYP,
     &         NucAtt,
     &         PBE, PBE0, M06L, M06, M062X,
     &         M06HF
#include "real.fh"
#include "WrkSpc.fh"
#include "nq_info.fh"
#include "debug.fh"
      Real*8  Grad(nGrad)
      Logical Do_MO,Do_TwoEl,Do_Grad
      Character*4 DFTFOCK
      Character*16  KSDFT
#include "dcscf.fh"
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
*
************************************************************************
*                                                                      *
       Select Case(KSDFT)
*
************************************************************************
*                                                                      *
*      LSDA LDA SVWN
*
       Case ('LSDA ','LDA ','SVWN ')
c        ExFac=Get_ExFac(KSDFT)
         Functional_type=LDA_type
         Call DrvNQ(LSDA   ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*      LSDA5 LDA5 SVWN5
*
       Case('LSDA5','LDA5','SVWN5')
c        ExFac=Get_ExFac(KSDFT)
         Functional_type=LDA_type
         Call DrvNQ(LSDA5  ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     HFB
*
       Case('HFB')
c        ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Call DrvNQ(HFB    ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*      HFS
*
       Case('HFS')
c        ExFac=Get_ExFac(KSDFT)
         Functional_type=LDA_type
         Call DrvNQ(HFS    ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*      XALPHA
*
       Case('XALPHA')
c        ExFac=Get_ExFac(KSDFT)
         Functional_type=LDA_type
         Call DrvNQ(XAlpha ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     Overlap
*
      Case('Overlap')
c        ExFac=Get_ExFac(KSDFT)
         Functional_type=LDA_type
         Call DrvNQ(Overlap,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     NucAtt
*
      Case('NucAtt')
c        ExFac=One
         Functional_type=LDA_type
         Call DrvNQ(NucAtt,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     BLYP
*
      Case('BLYP')
c        ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Call DrvNQ(BLYP   ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     TLYP
*
      Case('TLYP')
c        ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Call DrvNQ(TLYP   ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     B3LYP
*
      Case('B3LYP ')
c        ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Call DrvNQ(B3LYP  ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     B3LYP5
*
      Case('B3LYP5')
c        ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Call DrvNQ(B3LYP5 ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     PBE
*
      Case('PBE')
c        ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Call DrvNQ(PBE   ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     PBE0
*
      Case('PBE0')
c        ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Call DrvNQ(PBE0  ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     M06-L
*
      Case('M06L')
c        ExFac=Get_ExFac(KSDFT)
         Functional_type=meta_GGA_type1
         Call DrvNQ(M06L,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     M06
      Case('M06 ')
c        ExFac=Get_ExFac(KSDFT)
         Functional_type=meta_GGA_type1
         Call DrvNQ(M06,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     M06-2X
      Case('M062X')
c        ExFac=Get_ExFac(KSDFT)
         Functional_type=meta_GGA_type1
         Call DrvNQ(M062X,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     M06-HF
      Case('M06HF')
c        ExFac=Get_ExFac(KSDFT)
         Functional_type=meta_GGA_type1
         Call DrvNQ(M06HF,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
      Case Default
         Call WarningMessage(2,
     &               ' Get_Exc_dft: Undefined functional type!')
         Write (6,*) '         Functional=',KSDFT(1:lKSDFT)
         Call Quit_OnUserError()
       End Select
*                                                                      *
************************************************************************
*                                                                      *
      Erest_xc=Erest_xc-Func
*
#ifdef _DEBUGPRINT_
      write(6,*) ' XC-part of energy-restoring term : ',-Func
      write(6,*)
      write(6,*) ' XC-potentials: (itri,F_alpha,F_beta)'
      write(6,*)
      Do i=1,nh1
        Write(6,'(i4,3f22.16)') i,Work(ipF_DFT+i-1),
     &                            Work(ipF_DFT+i-1+nh1)
      End Do
#endif
*
      Return
      End
