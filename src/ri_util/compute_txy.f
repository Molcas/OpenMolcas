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
      subroutine compute_txy(DM1,nDM,Txy,nTxy,nAuxVec,nIrrep,Diag,
     &                       DMTmp,nAct)
************************************************************************
*                                                                      *
*     Compute the matrices needed for CD-CASSCF gradients              *
*                                                                      *
*     input : G2    = 2-body density matrix                           *
*             nDM    = size of the one-body DM                         *
*                                                                      *
************************************************************************
      use pso_stuff
      Implicit none
#include "real.fh"
      Integer nTxy,nAct(0:7),nCumAct(0:7),nCumAct2(0:7)
      Integer nDM,i,j,icol,iline
      Integer ista,iend,jsta,jend,ksta,kend,lsta,lend,isym,jsym,
     &        ksym,lsym,klsym,it,iu,iv,ix,itu,ivx,ituvx,ituvx2,
     &        nvx,itx,itv,iuv,iux,nkl,Txy_sta,Txy_sta2
      Integer nAuxVec,iVec,nIrrep
      Real*8 Fac,Fac2,tmp
      Real*8 DM1(nDM,nAuxVec),
     &       Txy(nTxy,nAuxVec),Diag(nDM,nAuxVec),
     &       DMtmp(nDM*(nDM+1)/2)
      Logical debug
      debug=.true.
      debug=.false.
*
      nCumAct(0)=0
      Do i=1,nIrrep-1
        nCumAct(i)=nCumAct(i-1)+nAct(i-1)
      End Do
*
      Do iVec=1,nAuxVec
************************************************************************
*                                                                      *
*     Remove Coulomb and exchange contribution from the 2DM            *
*                                                                      *
************************************************************************
         Txy_sta=1
         Txy_sta2=1
         Do klsym=0,nIrrep-1 ! compound symmetry
           nkl=0
*
           Do lSym=0,nIrrep-1
             lsta=nCumAct(lsym)+1
             lend=nCumAct(lsym)+nAct(lSym)

             ksym=iEOR(lsym,klsym)
             If (ksym.gt.lsym) Go To 100
             ksta=nCumAct(ksym)+1
             kend=nCumAct(ksym)+nAct(ksym)
             If (kSym.eq.lSym) Then
                nvx=nAct(lSym)*(nAct(lSym)+1)/2
             Else
                nvx=nAct(lSym)*nAct(ksym)
             EndIf
             nCumAct2(lSym)=nkl
*
             Do jsym=0,lSym
*
               jsta=nCumAct(jsym)+1
               jend=nCumAct(jsym)+nAct(jSym)

               isym=iEOR(jsym,klsym)
               If (isym.gt.jsym) Go To 101
               ista=nCumAct(iSym)+1
               iend=nCumAct(iSym)+nAct(iSym)
*
               iline=nkl
*
               Do ix=lsta,lend
                 If (kSym.eq.lSym) kend=ix
                 Do iv=ksta,kend
                   ivx=max(ix,iv)*(max(ix,iv)-1)/2+min(ix,iv)
                   If (jSym.eq.lSym) jend=ix
                   iline=iline+1
                   icol=nCumAct2(jSym)
                   Do iu=jsta,jend
                     iux=max(ix,iu)*(max(ix,iu)-1)/2+min(ix,iu)
                     iuv=max(iu,iv)*(max(iu,iv)-1)/2+min(iu,iv)
                     If (iSym.eq.jSym) iend=iu
                     Do it=ista,iend
                       itu=max(iu,it)*(max(iu,it)-1)/2+min(iu,it)
                       If (itu.gt.ivx) Go to 102
                       itx=max(ix,it)*(max(ix,it)-1)/2+min(ix,it)
                       itv=max(iv,it)*(max(iv,it)-1)/2+min(iv,it)
                       ituvx=max(ivx,itu)*(max(ivx,itu)-1)/2+
     &                      min(itu,ivx)
                       icol=icol+1
                       ituvx2=iline*(iline-1)/2+icol
*
                       Fac=One
                       If (ix.ne.iv) Fac=Two*Fac
                       If (it.ne.iu) Fac=Two*Fac
                       Fac2=One
                       If (it.eq.iu) Fac2=Two
*
                       DMTmp(ituvx2)=Fac*(Fac2*G2(ituvx,iVec))
                       If (.Not.lSA) Then
*For SA-CASSCF, don't remove Coulomb and exchange
                         If (iSym.eq.jSym)
     &                      DMTmp(ituvx2)=DMTmp(ituvx2)-
     &                         Fac*(DM1(itu,iVec)*DM1(ivx,iVec))
                         If (iSym.eq.kSym)
     &                      DMTmp(ituvx2)=DMTmp(ituvx2)+
     &                         Fac*(Quart*DM1(itv,iVec)*DM1(iux,iVec))
                         If (iSym.eq.lSym)
     &                      DMTmp(ituvx2)=DMTmp(ituvx2)+
     &                         Fac*(Quart*DM1(itx,iVec)*DM1(iuv,iVec))
                       EndIf
*
 102                   Continue
*
                     End Do
                   End Do
                 End Do
               End Do
*
 101           Continue
             End Do
             nkl=nkl+nvx
 100         Continue
           End Do
*
************************************************************************
*                                                                      *
*     Eigen-decompose the density                                      *
*                                                                      *
************************************************************************
*
**       Diagonalize G2
*
           Call Cho_DZero(Txy(Txy_sta2,iVec),nkl**2)
           call dcopy_(nkl,[One],0,Txy(Txy_sta2,iVec),nkl+1)
*
           Call NIdiag(DMTmp,Txy(Txy_sta2,iVec),nkl,nkl,0)
*
**       Multiply by Sqrt[eigenvalue]
*
           Do i=1,nkl
             Diag(i+Txy_sta-1,iVec)=DMTmp(i*(i+1)/2)
             tmp=Sqrt(abs(DMTmp(i*(i+1)/2)))
             Do j=1,nkl
                Txy(Txy_sta2+(i-1)*nkl+j-1,iVec)=
     &              Txy(Txy_sta2+(i-1)*nkl+j-1,iVec)*tmp
             End Do
           End Do
*
**       Since there is no screening yet
*
           nnP(klsym)=nkl
*
           Txy_sta=Txy_sta+nkl
           Txy_sta2=Txy_sta2+nkl**2
         End Do
      End Do
      end
