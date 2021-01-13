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
      SubRoutine TMOMInt(wavevector,iOpt)
************************************************************************
*                                                                      *
* Object: driver for computation of TMOM integrals                     *
*                                                                      *
************************************************************************
      Use MpmC
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      use Sizes_of_Seward, only: S
      use Basis_Info, only: nBas
      use Symmetry_Info, only: nIrrep
#endif
      Implicit Real*8 (A-H,O-Z)
      External EMFInt, EMFMem
*     ipList: list of pointers to the integrals of each component
*             of the operator
*     OperI: list which irreps a particular component of the operator
*            belongs to
*     OperC: list the character of each component of the operator
*     CoorO: list of origins of the operator, one for each component
      Integer, Dimension(:), Allocatable :: ipList, OperI, OperC
      Real*8, Dimension(:), Allocatable :: CoorO, Nuc
      Real*8 wavevector(3)
#ifdef _DEBUGPRINT_
#include "stdalloc.fh"
      Real*8, Allocatable :: Int_R(:), Int_I(:), Temp_Int(:)
      Real*8, Allocatable :: Int_R_O(:), Int_I_O(:)
      Integer IOFF(8,8)
#endif
#include "print.fh"
#include "nq_info.fh"
#include "real.fh"
#include "wldata.fh"
#include "property_label.fh"
#include "oneswi.fh"
#include "warnings.fh"
      Character*8 Label
      Dimension dum(1),idum(1)
*
#ifdef _DEBUGPRINT_
      MulTab(i,j)=iEor(i-1,j-1)+1
#endif
*
      Call Set_Basis_Mode('Valence')
      Call Setup_iSD()
*
************************************************************************
************************************************************************
*                                                                      *
*     Electromagnetic field radiation integrals.                       *
*                                                                      *
*     Note that the integral is neither symmetric nor antisymmetric!   *
*                                                                      *
************************************************************************
************************************************************************
      rHrmt=-One ! Not used
*
*     B*s Magnetic * Spin: not that this boils down to just integrals
*     over A.
*
      If (iOpt.eq.2) Then
         rHrmt=One
         nOrdOp = 0
         Label='TMOM0'
         nComp = 2
         Call Allocate_Aux()
*        Here we put in the k-vector
         Call FZero(CoorO,3*nComp)
         Call dcopy_(3,wavevector,1,CoorO,1)
*
*        The electromagnetic field operator contributes to all
*        irreducible irreps, hence OperI=255. Since the operator
*        itself is not symmetry-adapted OperC is set to a dummy value.
*
         OperI(:) = 255
         OperC(:) = 0 ! Dummy
*
         Call dcopy_(nComp,[Zero],0,Nuc,1)
         Call OneEl(EMFInt,EMFMem,Label,ipList,OperI,nComp,
     &              CoorO,nOrdOp,Nuc,rHrmt,OperC,
     &              dum,1,dum,idum,0,0,
     &              dum,1,0)
*
         Call Deallocate_Aux()
#ifdef _DEBUGPRINT_
*
      Call mma_allocate(CoorO,6,Label='CoorO')
      CoorO(:)=0.0D0
      CoorO(1:3)=wavevector
      Write (6,*) 'Wavevector=',Wavevector
*
*     This section of the code is for pure debugging and will replace
*     exact operator with truncated expansions of the operator in
*     terms of multipole integrals
*
      iOpt0=0 ! Write
      iOpt1=1 ! Read just data size and symmetry
      iOpt2=2 ! Read
      iRc=-1
      Label='TMOM0  R'
      iComp=1
*     Pick up the size and the symmetry label.
      Call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl_TMOM)
      nInts_TMOM=idum(1)
!     Write (*,*) 'nInts_TMOM=',nInts_TMOM
      Call mma_allocate(Int_R,nInts_TMOM+4,Label='Int_R')
      Call mma_allocate(Int_I,nInts_TMOM+4,Label='Int_I')
      Call mma_allocate(Int_R_O,nInts_TMOM+4,Label='Int_R_O')
      Call mma_allocate(Int_I_O,nInts_TMOM+4,Label='Int_I_O')
*
      Call RdOne(iRc,iOpt0,Label,iComp,Int_R_O,iSyLbl_TMOM)
      Label='TMOM0  I'
      iComp=1
      Call RdOne(iRc,iOpt0,Label,iComp,Int_I_O,iSyLbl_TMOM)
      Len=0
      IOFF(:,:)=-1
      Do i=1,nIrrep
         Do j=1,i
            ij=MulTab(i,j)-1
            If(iAnd(2**ij,iSyLbl_TMOM).ne.0) Then
               IOFF(i,j)=Len+1
               If(i.eq.j) Then
                  Len_=nBas(i-1)*(nBas(i-1)+1)/2
               Else
                  Len_=nBas(i-1)*nBas(j-1)
               End If
               Len = Len + Len_
            End If
         End Do
      End Do
!     Write (*,*) 'Len=',Len
!     Do I = 1, 8
!        Write (6,*) (IOFF(I,J),J=1,8)
!     End Do
*
      Int_R(:)=0.0D0
      Int_R(nInts_TMOM+1:nInts_TMOM+3)=CoorO
      Int_I(:)=0.0D0
      Int_I(nInts_TMOM+1:nInts_TMOM+3)=CoorO
*
      S%nMltpl=9
      iCase=1
      Phase=1.0D0
      Do iMltpl= 0, S%nMltpl
         Write (Label,'(A,I2)') 'Mltpl ',iMltpl
         nComp=(iMltpl+1)*(iMltpl+2)/2
         iComp=0
         Do ix = iMltpl, 0, -1
            x=CoorO(1)**ix
            Do iy = iMltpl-ix, 0, -1
               xy=x*CoorO(2)**iy
               iz = iMltpl-ix-iy
               xyz=xy*CoorO(3)**iz
*
               Fact=Phase*xyz/(Gamma(Dble(ix)+1.0D0)
     &                        *Gamma(Dble(iy)+1.0D0)
     &                        *Gamma(Dble(iz)+1.0D0))
*
*              Write (*,*) 'Fact=',Fact
*              Write (6,*) CoorO(1)**ix, ix
*              Write (6,*) CoorO(2)**iy, iy
*              Write (6,*) CoorO(3)**iz, iz
*
               iComp=iComp+1
               If (Fact.eq.0.0D0) cycle
               Call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl)
!              Write (*,*) 'iRC=',iRC
               If (iRC.ne.0) Then
                  Write (6,*) 'TMOMINT: Error reading ',Label
                  Call Abend()
               End If
               nInts=idum(1)
!              Write (6,*) 'nInts=',nInts
               Call mma_allocate(Temp_Int,nInts+4,Label='Temp_Int')
               Call RdOne(iRc,iOpt0,Label,iComp,Temp_Int,iSyLbl)
!              Write (*,*) 'Temp_Int(1)',Temp_Int(1)
*
               Len=0
               Do i=1,nIrrep
                  Do j=1,i
                     ij=MulTab(i,j)-1
                     If(iAnd(2**ij,iSyLbl).ne.0) Then
                        jOff=Len+1
                        If(i.eq.j) Then
                           Len_=nBas(i-1)*(nBas(i-1)+1)/2
                        Else
                           Len_=nBas(i-1)*nBas(j-1)
                        End If
!                       Write (*,*) 'jOff,IOFF(i,j),Len_=',
!    &                               jOff,IOFF(i,j),Len_
                        If (iCase.eq.1) Then
*                          Contribution to the real part
                           Call DaXpY_(Len_,Fact,Temp_Int(jOff),1,
     &                                 Int_R(IOFF(i,j)),1)
                        Else
*                          Contribution to the imaginary part
                           Call DaXpY_(Len_,Fact,Temp_Int(jOff),1,
     &                                 Int_I(IOFF(i,j)),1)
*                          If (i.eq.1.and.j.eq.1.and.nIrrep.eq.1) Then
*                             Write (6,*) nBas(i-1),nBas(j-1)
*                             Write (*,*) 'Fact=',Fact
*                             Write (6,*) CoorO(1)**ix, ix
*                             Write (6,*) CoorO(2)**iy, iy
*                             Write (6,*) CoorO(3)**iz, iz
*                             Write (6,*) 'Temp_Int(385)=',
*    &                                     Temp_Int(jOff+384)
*                             Write (6,*) 'Int_I(385)=',
*    &                                     Int_I(IOFF(i,j)+384)
*                             Write (6,*) 'Int_I_O(385)=',
*    &                                     Int_I_O(IOFF(i,j)+384)
*                          End If
*                          If (i.eq.2.and.j.eq.1.and.nIrrep.eq.8) Then
*                             Write (6,*) nBas(i-1),nBas(j-1)
*                             Write (*,*) 'Fact=',Fact
*                             Write (6,*) CoorO(1)**ix, ix
*                             Write (6,*) CoorO(2)**iy, iy
*                             Write (6,*) CoorO(3)**iz, iz
*                             Write (6,*) 'Temp_Int(55)=',
*    &                                     Temp_Int(jOff+54)
*                             Write (6,*) 'Int_I(55)=',
*    &                                     Int_I(IOFF(i,j)+54)
*                             Write (6,*) 'Int_I_O(55)=',
*    &                                     Int_I_O(IOFF(i,j)+54)
*                          End If
                        End If
                        Len = Len + Len_
                     End If
                  End Do
               End Do
               Call mma_deallocate(Temp_Int)
            End Do
         End Do
*
         If (iCase.eq.1) Then
            iCase=2
         Else
            iCase=1
            Phase=-Phase
         End If
      End Do
*
*     Compare exact integrals with approximated.
*
#define _COMPARE_
#ifdef _COMPARE_
      Len=0
      Do i=1,nIrrep
         Do j=1,i
            ij=MulTab(i,j)-1
            If(iAnd(2**ij,iSyLbl_TMOM).ne.0) Then
               If(i.eq.j) Then
                  Len_=nBas(i-1)*(nBas(i-1)+1)/2
               Else
                  Len_=nBas(i-1)*nBas(j-1)
               End If
               Do iLen = 1, Len_
*                    Write (6,*) 'isym,jsym,iLen=',i,j,iLen
*                    Write (*,*) 'Int_R,Int_Q=',Int_R_O(Len+iLen),
*    &                                          Int_R(Len+iLen)
*                    Write (*,*) 'Int_I,Int_J=',Int_I_O(Len+iLen),
*    &                                          Int_I(Len+iLen)
                  temp= Abs(Int_R_O(Len+iLen)-Int_R(Len+iLen))/
     &                  Max(Abs(Int_R_O(Len+iLen)),
     &                      Abs(Int_R(Len+iLen)),1.0D-8)
                  If (temp.gt.1.0D-2) Then
                     Write (6,*) 'isym,jsym,iLen=',i,j,iLen
                     Write (6,*) 'Int_R,Int_Q=',Int_R_O(Len+iLen),
     &                                          Int_R(Len+iLen)
                  End If
                  temp= Abs(Int_I_O(Len+iLen)-Int_I(Len+iLen))/
     &                  Max(Abs(Int_I_O(Len+iLen)),
     &                      Abs(Int_I(Len+iLen)),1.0D-8)
                  If (temp.gt.1.0D-2) Then
                     Write (6,*) 'isym,jsym,iLen=',i,j,iLen
                     Write (6,*) 'Int_I,Int_J=',Int_I_O(Len+iLen),
     &                                          Int_I(Len+iLen)
                  End If
               End Do
!     Write (*,*) 'i,j=',i,j
!     Write (*,*) 'Int_R,Int_Q=',Int_R_O(Len+1),Int_R(Len+1)
!     Write (*,*) 'Int_I,Int_J=',Int_I_O(Len+1),Int_I(Len+1)
!     Write (*,*) 'Int_R,Int_Q=',
!    &                     DDOT_(Len_,1.0D0,0,Int_R_O(Len+1),1),
!    &                     DDOT_(Len_,1.0D0,0,Int_R(Len+1),1)
!     Write (*,*) 'Int_I,Int_J=',
!    &                     DDOT_(Len_,1.0D0,0,Int_I_O(Len+1),1),
!    &                     DDOT_(Len_,1.0D0,0,Int_I(Len+1),1)
               Len = Len + Len_
            End If
         End Do
      End Do
#endif
*
*     Overwrite the integrals with a truncated expansion.
*
      Label='TMOM0  R'
      iComp=1
      Call WrOne(iRc,iOpt0,Label,iComp,Int_R,iSyLbl_TMOM)
      Label='TMOM0  I'
      iComp=1
      Call WrOne(iRc,iOpt0,Label,iComp,Int_I,iSyLbl_TMOM)
*
      Call mma_deallocate(Int_R_O)
      Call mma_deallocate(Int_I_O)
      Call mma_deallocate(Int_R)
      Call mma_deallocate(Int_I)
      Call mma_deallocate(CoorO)
*
#endif
      End If
*
*     A*nabla. Note that when used the numbers are multiplied with -i to
*     generate A*p.
*
      If (iOpt.le.2) Then
      nOrdOp = 1
      Label='TMOM'
      nComp = 12
      Call Allocate_Aux()
*     Here we put in the k-vector
      Call FZero(CoorO,3*nComp)
      Call dcopy_(3,wavevector,1,CoorO,1)
*
*     The electromagnetic field operator contributes to all
*     irreducible irreps, hence OperI=255. Since the operator
*     itself is not symmetry-adapted OperC is set to a dummy value.
*
      OperI(:) = 255
      OperC(:) = 0 ! Dummy
*
      Call dcopy_(nComp,[Zero],0,Nuc,1)
      Call OneEl(EMFInt,EMFMem,Label,ipList,OperI,nComp,
     &           CoorO,nOrdOp,Nuc,rHrmt,OperC,
     &           dum,1,dum,idum,0,0,
     &           dum,1,0)
*
      Call Deallocate_Aux()
      End If
*
*     The A^2 term
*
      If (iOpt.gt.2) Then
      nOrdOp = 0
      Label='TMOM2'
      nComp = 2
      Call Allocate_Aux()
*     Here we put in the k-vector
      Call FZero(CoorO,3*nComp)
      Call dcopy_(3,wavevector,1,CoorO,1)
*     Change the argument to 2xA
      Call dscal_(3,2.0D0,CoorO,1)
*
*     The electromagnetic field operator contributes to all
*     irreducible irreps, hence OperI=255. Since the operator
*     itself is not symmetry-adapted OperC is set to a dummy value.
*
      OperI(:) = 255
      OperC(:) = 0 ! Dummy
*
      Call dcopy_(nComp,[Zero],0,Nuc,1)
      Call OneEl(EMFInt,EMFMem,Label,ipList,OperI,nComp,
     &           CoorO,nOrdOp,Nuc,rHrmt,OperC,
     &           dum,1,dum,idum,0,0,
     &           dum,1,0)
*
      Call Deallocate_Aux()
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call Free_iSD()
*                                                                      *
************************************************************************
*                                                                      *
      Return
*
      Contains
      Subroutine Allocate_Aux()
      Implicit None
#include "stdalloc.fh"
*
      Call mma_Allocate(ipList,nComp,Label='ipList')
      Call mma_Allocate(OperI,nComp,Label='OperI')
      Call mma_Allocate(OperC,nComp,Label='OperC')
      Call mma_Allocate(CoorO,3*nComp,Label='CoorO')
      Call mma_Allocate(Nuc,nComp,Label='Nuc')
*
      Return
      End Subroutine Allocate_Aux
      Subroutine Deallocate_Aux()
      Implicit None
#include "stdalloc.fh"
*
      Call mma_Deallocate(OperC)
      Call mma_Deallocate(OperI)
      Call mma_Deallocate(ipList)
      Call mma_Deallocate(CoorO)
      Call mma_Deallocate(Nuc)
*
      Return
      End Subroutine Deallocate_Aux
*
      End Subroutine TMOMInt
