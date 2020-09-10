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
      subroutine Tr_prm_cnt(idbg,nBas_Cont,nBas_Prim)
      use Basis_Info
      implicit real*8(a-h,o-z)
#include "itmax.fh"
#include "info.fh"
#include "rinfo.fh"
#include "stdalloc.fh"
#include "real.fh"
      Integer icaddr(MxAO),numc(MxAO),ihelp(MxAtom,MxAng),numb(MxAO),
     &        mcaddr(MxAO), nBas_Cont(8), nBas_Prim(0:7)
      Logical New_Center,New_l,New_m, Old_Center, Old_l
      Real*8, Dimension(:), Allocatable :: Tr
*     contracted basis, atomic basis functions
*
*     symmetry info
*
*     lant(i): number of atoms in i:th symmetry bf
*     expand the coefficient matrix into symmetry basis set
*     auxiliary
*     icaddr(i): adresses in coeff for a symmetry adapted function
*
*                                                                      *
* THIS ROUTINE ONLY WORKS WITH SPHERICAL FUNCTIONS. NO CARTESIAN D:S!  *
************************************************************************
*                                                                      *
      nSym=nIrrep
      iPrint=0
      If (iprint.ge.10.or.idbg.gt.0) then
         write(idbg,*) ' in repmat', nsym
         write(idbg,*) nSym, (nBas(i),i=0,nsym-1)
         write(idbg,*) nSym, (nrBas(i),i=1,nsym)
         write(idbg,*) nSym, (nBas_Prim(i),i=0,nsym-1)
         write(idbg,*) nSym, (nBas_Cont(i),i=1,nsym)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     set up pointer
*
      k=0
      ia=0  ! center index
      ka=0  ! shell index
      Do iCnttp=1,nCnttp
         Do icnt = 1, dbsc(iCnttp)%nCntr
            ia=ia+1
            Do la=1,nAngr(ia)+1
               ka=ka+1
               ihelp(ia,la)=k
               k=k+nPrimr(ka)*nBasisr(ka)
            End Do
         End Do
      End Do
      If (iPrint.ge.10.or.idbg.gt.0) Then
         write(idbg,*) ' Help vector'
         ia=0
         Do iCnttp=1,nCnttp
            Do icnt = 1, dbsc(iCnttp)%nCntr
               ia=ia+1
               write(idbg,'(10i5)') (ihelp(ia,j),j=1,nAngr(ia)+1)
           End Do
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Loop over irreps
*
      k=0
      Do iSym = 1, nSym
         numck=1
         numcl=0
         kbias=0
*
*        Loop over basis functions in irrep
*
         Do iCont = 1, nBas_Cont(iSym)
            k=k+1
            If (iCont.gt.1) Then
               New_Center=icent(k).ne.icent(k-1)
               New_l=lnang(k).ne.lnang(k-1)
               New_m=lmag(k).ne.lmag(k-1)
               If (New_m)  kbias=kbias-numc(k-1)
               If (New_Center.or.New_l)  kbias=0
            Else
               New_Center=.True.
               New_l=.True.
            End If
            kbias = kbias + 1
            ka=0
            ia=0
            Do iCnttp = 1, nCnttp
               Do iCnt = 1, dbsc(iCnttp)%nCntr
                  ia=ia+1
                  Do la = 1, nAngr(ia)+1
                     ka=ka+1
                     Old_Center=icent(k).eq.ia
                     Old_l=lnang(k).eq.(la-1)
                     If (idbg.gt.0) write(idbg,*) ' at numck', k,ia,
     &                   icent(k),la-1,lnang(k),ia,New_Center,New_l
                     If (Old_Center.and.Old_l) Then
                        numc(k)=nBasisr(ka)
                        numb(k)=nPrimr(ka)
                        icaddr(k)=ihelp(ia,la)+(kbias-1)*nPrimr(ka)
                        If (k.gt.1.and.kbias.eq.1) numck=numcl+numck
                        mcaddr(k)=numck
                        numcl=nPrimr(ka)
                     End If
                  End Do  ! la
               End Do     ! iCnt
            End Do        ! iCnttp
         End Do           ! iCont
      End Do              ! iSym
      k=0
      If (iPrint.ge.10.or.idbg.gt.0) then
         ic=0
         ip=0
         Do iSym = 1, nSym
            Write (idbg,*)        ' symmetry',iSym
            Write (idbg,*)        ' numb'
            Write (idbg,'(20i4)') (numb(i+ic),  i=1,nBas_Cont(iSym))
            Write (idbg,*)        ' numc'
            Write (idbg,'(20i4)') (numc(i+ic),  i=1,nBas_Cont(iSym))
            Write (idbg,*)        ' Pointer to contraction vector'
            Write (idbg,'(20i4)') (icaddr(i+ic),i=1,nBas_Cont(iSym))
            Write (idbg,*)        ' mcaddr'
            Write (idbg,'(20i4)') (mcaddr(i+ic),i=1,nBas_Cont(iSym))
            ic=ic+nBas_Cont(iSym)
            ip=ip+nBas_Prim(iSym-1)
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
      nSize=0
      Do iSym = 1, nSym
         nSize=nSize+nBas_Cont(iSym)*nBas_Prim(iSym-1)
      End Do
      Call mma_allocate(Tr,nSize,label='Tr')
      Call DCopy_(nSize,[Zero],0,Tr,1)
*                                                                      *
************************************************************************
*                                                                      *
*     Generate the transformation matrix.
*
      ibasL=0
      iOff = 0
*
*---- Loop over irreps
*
      Do iSym = 1, nSym
*
*        loop over contracted
*
         Do iBas = 1, nBas_Cont(iSym)
            ibasL=ibasL+1
            ipbasL=mcaddr(ibasL) -1
*
*           loop over uncontracted
*
            Do iPrim=1,numb(ibasL)
               ipbasL = ipbasL+1
               index=iOff + (iBas-1)*nBas_Prim(iSym-1) + ipbasL
c              Write (*,*) iBas,ipbasL
               Tr(index)=rCof(icaddr(ibasL)+iPrim)
*
            End Do     ! iprim
         End Do
         iOff = iOff + nBas_Cont(iSym)*nBas_Prim(iSym-1)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Call Put_dArray('NEMO TPC',Tr,nSize)
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(Tr)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
