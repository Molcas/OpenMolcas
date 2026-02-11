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
      use definitions, only: iwp, wp, u6
      Use MpmC
      Use Integral_interfaces, only: int_kernel, int_mem
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      use definitions, only: u6
      use OneDat, only: sOpSiz
      use Sizes_of_Seward, only: S
      use Basis_Info, only: nBas
      use Symmetry_Info, only: nIrrep, MulTab=>Mul
      use stdalloc, only: mma_allocate, mma_deallocate
#endif
      use Constants, only: Zero, One, Two
      Implicit None
      Procedure(int_kernel) :: EMFInt
      Procedure(int_mem) :: EMFMem
      Real(kind=wp), intent(in):: wavevector(3)
      integer(kind=iwp), Intent(in):: iOpt

*     ipList: list of pointers to the integrals of each component
*             of the operator
*     OperI: list which irreps a particular component of the operator
*            belongs to
*     OperC: list the character of each component of the operator
*     CoorO: list of origins of the operator, one for each component
      Integer(kind=iwp), Allocatable :: ipList(:), OperI(:), OperC(:)
      Real(kind=wp), Allocatable :: CoorO(:), Nuc(:)
#ifdef _DEBUGPRINT_
      Real(kind=wp), Allocatable :: Int_R(:), Int_I(:), Temp_Int(:)
      Real(kind=wp), Allocatable :: Int_R_O(:), Int_I_O(:)
      integer(kind=iwp) IOFF(8,8)
#endif
#include "warnings.h"
      Character(LEN=8) Label
      real(kind=wp) dum(1), rHrmt
      Integer(kind=iwp) nComp,nOrdOp
*
#ifdef _DEBUGPRINT_
      Integer(kind=iwp) idum(1),i,iCase,ij,ilen,iMltpl,iOpt0,iOpt1,iRc,
     &                  iSyLbl,ix,iy,iz,j,jOff,Len,Len_,nInts,iComp,
     &                  iSyLbl_TMOM, nInts_TMOM
      real(kind=wp) Fact, Phase, Temp, x, xy, xyz
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
     &              dum,1,0)
*
         Call Deallocate_Aux()
#ifdef _DEBUGPRINT_
*
      Call mma_allocate(CoorO,6,Label='CoorO')
      CoorO(:)=Zero
      CoorO(1:3)=wavevector
      Write (u6,*) 'Wavevector=',Wavevector
*
*     This section of the code is for pure debugging and will replace
*     exact operator with truncated expansions of the operator in
*     terms of multipole integrals
*
      iOpt0=0 ! Write
      iOpt1=ibset(0,sOpSiz) ! Read just data size and symmetry
      iRc=-1
      Label='TMOM0  R'
      iComp=1
*     Pick up the size and the symmetry label.
      Call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl_TMOM)
      nInts_TMOM=idum(1)
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
*
      Int_R(:)=Zero
      Int_R(nInts_TMOM+1:nInts_TMOM+3)=CoorO
      Int_I(:)=Zero
      Int_I(nInts_TMOM+1:nInts_TMOM+3)=CoorO
*
      S%nMltpl=9
      iCase=1
      Phase=One
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
               Fact=Phase*xyz/(Gamma(Dble(ix)+One)
     &                        *Gamma(Dble(iy)+One)
     &                        *Gamma(Dble(iz)+One))
*
               iComp=iComp+1
               If (Fact.eq.Zero) cycle
               Call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl)
               If (iRC.ne.0) Then
                  Write (u6,*) 'TMOMINT: Error reading ',Label
                  Call Abend()
               End If
               nInts=idum(1)
               Call mma_allocate(Temp_Int,nInts+4,Label='Temp_Int')
               Call RdOne(iRc,iOpt0,Label,iComp,Temp_Int,iSyLbl)
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
                        If (iCase.eq.1) Then
*                          Contribution to the real part
                           Call DaXpY_(Len_,Fact,Temp_Int(jOff),1,
     &                                 Int_R(IOFF(i,j)),1)
                        Else
*                          Contribution to the imaginary part
                           Call DaXpY_(Len_,Fact,Temp_Int(jOff),1,
     &                                 Int_I(IOFF(i,j)),1)
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
                  temp= Abs(Int_R_O(Len+iLen)-Int_R(Len+iLen))/
     &                  Max(Abs(Int_R_O(Len+iLen)),
     &                      Abs(Int_R(Len+iLen)),1.0D-8)
                  If (temp.gt.1.0D-2) Then
                     Write (u6,*) 'isym,jsym,iLen=',i,j,iLen
                     Write (u6,*) 'Int_R,Int_Q=',Int_R_O(Len+iLen),
     &                                          Int_R(Len+iLen)
                  End If
                  temp= Abs(Int_I_O(Len+iLen)-Int_I(Len+iLen))/
     &                  Max(Abs(Int_I_O(Len+iLen)),
     &                      Abs(Int_I(Len+iLen)),1.0D-8)
                  If (temp.gt.1.0D-2) Then
                     Write (u6,*) 'isym,jsym,iLen=',i,j,iLen
                     Write (u6,*) 'Int_I,Int_J=',Int_I_O(Len+iLen),
     &                                          Int_I(Len+iLen)
                  End If
               End Do
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
      Call dscal_(3,Two,CoorO,1)
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
*
      Contains
      Subroutine Allocate_Aux()
      Use stdalloc, Only: mma_allocate
      Implicit None
*
      Call mma_Allocate(ipList,nComp,Label='ipList')
      Call mma_Allocate(OperI,nComp,Label='OperI')
      Call mma_Allocate(OperC,nComp,Label='OperC')
      Call mma_Allocate(CoorO,3*nComp,Label='CoorO')
      Call mma_Allocate(Nuc,nComp,Label='Nuc')
*
      End Subroutine Allocate_Aux
      Subroutine Deallocate_Aux()
      Use stdalloc, Only: mma_deallocate
      Implicit None
*
      Call mma_Deallocate(OperC)
      Call mma_Deallocate(OperI)
      Call mma_Deallocate(ipList)
      Call mma_Deallocate(CoorO)
      Call mma_Deallocate(Nuc)
*
      End Subroutine Deallocate_Aux
*
      End Subroutine TMOMInt
