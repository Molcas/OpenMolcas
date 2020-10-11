************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2000, Gunnar Karlstrom                                 *
*               2000, Roland Lindh                                     *
************************************************************************
      Subroutine ener(h1,TwoHam,D,RepNuc,nh1,First,Dff,D_Tot,Grid,
     &                nGrid_,DipMom,EField,DipEff,PolEff,
     &                Cord,MaxAto,Z_Nuc,nPolComp,nAnisopol,pField,tmpF)
************************************************************************
*                                                                      *
*     Object: to compute the needed terms from the interaction between *
*             the charge of the QM system, the dipole moments on the   *
*             grid, and the counter charge on the boundary of the      *
*             cavity. The latter can be divided into three terms,      *
*             nuclear, electronic and Langevin dipole moment terms.    *
*                                                                      *
*     Authors: G. Karlstroem                                           *
*              Dept. of Theor. Chem., Univ. of Lund, Sweden.           *
*                                                                      *
*              and                                                     *
*                                                                      *
*              R. Lindh                                                *
*              Dept. of Chem. Phys., Univ. of Lund, Sweden.            *
*                                                                      *
*              March 2000                                              *
************************************************************************
      use Symmetry_Info, only: nIrrep, iChBas
      use Basis_Info, only: nBas
      use Temporary_Parameters, only: PrPrt
      Implicit Real*8 (a-h,o-z)
      External EFInt,EFMem
#include "real.fh"
#include "rctfld.fh"
#include "stdalloc.fh"
      Real*8 h1(nh1), TwoHam(nh1), D(nh1), D_tot(nh1), EF_Grid(3),
     &       Grid(3,nGrid_), DipMom(3,nGrid_), EField(4,nGrid_),
     &       DipEff(nGrid_), PolEff(nPolComp,nGrid_), CCoor(3),
     &       Cord(3,MaxAto), Z_Nuc(MaxAto),pField(4,nGrid_),
     &     tmpF(4,nGrid_)
      Logical First, Dff, Save_tmp , NonEq
      Character*8 Label
      Integer, Allocatable:: ips(:), lOper(:), kOper(:)
      Real*8, Allocatable::  C_Coor(:,:)
      Real*8, Allocatable:: Integrals(:)
*
#include "oneel_interface.fh"
*                                                                      *
************************************************************************
*                                                                      *
*---- First compute contributions due to the interaction of the QM
*     charge distributions with the corresponding counter charges on
*     the surface of the cavity.
*
*     Here we have three terms
*     1) nuclear-nuclear, constant: added to RepNuc once
*     2) nuclear-electronic, constant: added to h1 once
*     3) electronic-electronic, linearly dependent of the density: TwoHam
*
      tmp_RepNuc=RepNuc
      NonEq=.False.
      If(lRFCav) Call RctFld(h1,TwoHam,D,RepNuc,nh1,First,Dff,NonEq)
*                                                                      *
************************************************************************
*                                                                      *
*1)   Compute interaction between the dipole moment on the grid and
*     the EF on the grid (Edip2). Eself is the self interaction term
*     of the dipole moments on the grid.
*
*     Both terms depend on D_tot and the modifies RepNuc each iteration.
*
*2)   Compute interaction term between the EF of the nuclear charges
*     and the dipole moments on the grid (EnucDip). This scalar terms
*     depends also on D_tot and modifies RepNuc at each iteration.
*
      qq0=Zero
*
      Edip2=Zero
      Eint=Zero
      Eself=Zero
      Enucdip=Zero
      Esimple=Zero
      iMax=1
      agsum=Zero

*
      Do iGrid=1,nGrid_
*------- Grid coordinates
         ghx=Grid(1,iGrid)
         ghy=Grid(2,iGrid)
         ghz=Grid(3,iGrid)
*------- Langevin dipol moment
         dx=DipMom(1,iGrid)
         dy=DipMom(2,iGrid)
         dz=DipMom(3,iGrid)
*------- Total Electric Field
         fx=EField(1,iGrid)
         fy=EField(2,iGrid)
         fz=EField(3,iGrid)
*------- Electric Field from permanent multipoles
         pfx=pField(1,iGrid)
         pfy=pField(2,iGrid)
         pfz=pField(3,iGrid)

*
*------- Twice the energy of Langevin dipole moments in EF
*
         Edip2=Edip2-fx*dx-fy*dy-fz*dz

         fDd=fx*dx+fy*dy+fz*dz   ! f dot d

*        Do NOT halve the contribution from permanent multipoles
         Eint=Eint-Half*(fDd + pfx*dx+pfy*dy+pfz*dz)

         Esimple=Esimple-Half*(tmpF(1,iGrid)*dx+tmpF(2,iGrid)*dy+
     &        tmpF(3,iGrid)*dz)

*
*------- Compute self energies of dipole moments
*
         ftot2=fx*fx+fy*fy+fz*fz
         ftot=sqrt(ftot2)
         x=ftot*DipEff(iGrid)*tk
         If (x.le.0.0000001d0) Then
            ag=Zero
            alang=Zero
         Else
            ex=exp(x)
            emx=One/ex
            alang=(ex+emx)/(ex-emx)-One/x
            ag=-log((ex-emx)/(Two*x))/tk
         End If
*

         If(iGrid.gt.nAnisopol) Then  ! isotropic
c            uind=DipEff(iGrid)*alang+ftot*PolEff(1,iGrid)
c            Eself=Eself+ftot*uind+ag-Half*ftot**2*PolEff(1,iGrid)
            Eself=Eself+fDd-Half*ftot2*PolEff(1,iGrid)
            agsum=agsum+ag
         Else ! anisotropic
            Eself=Eself+Half*(fx*dx+fy*dy+fz*dz)
         EndIf
*
*------- Compute the EF at the grid point due to the nucleus
*
         nOrdOp=1
         Call EFNuc(Grid(1,iGrid),Z_Nuc,Cord,MaxAto,
     &              EF_Grid,nOrdOp)
*
         Enucdip=Enucdip-EF_Grid(1)*dx
     &                  -EF_Grid(2)*dy
     &                  -EF_Grid(3)*dz
*
      End Do          ! iGrid
*
#ifdef _DEBUGPRINT_
      Write(6,*)'Esimple             =',Esimple
      Write(6,*)'EnucRctfld          =',RepNuc-tmp_RepNuc
      Write(6,*)'Eself               =',Eself
      Write(6,*)'Ag term             =',agsum
      Write(6,*)'Electrostatic energy=',Eint
      Write(6,*)'EnucDip (half)      =',Half*Enucdip
#endif

      RepNuc = RepNuc + Eself + agsum + Eint + Half*Enucdip
*
c     Note:Before PS made this code compatible with XFIE, the following
c     quantity was reported as Electrostatic energy: Half*Edip2
c     and: RepNuc = RepNuc + Eself + Half*Edip2 + Half*Enucdip


*                                                                      *
************************************************************************
*                                                                      *
*---- Compute interaction term between the field of the electronic
*     charge and the dipole moments. This term is modified on each
*     iteration.
*
*---- Retrieve the statatic h1
*
      Label='h1_raw  '
      Call Get_Temp(Label,h1,nh1)

*
*---- Compute the modification to h1
*
*     Evaluate electronic EF integrals at each point and add to h1
*
*
      rHrmt=One
      nOrdOp = 1
      nComp = 3
      Sig=-One
      ixyz=1
      iSymX = 2**IrrFnc(ixyz)
      ixyz=2
      iSymY = 2**IrrFnc(ixyz)
      ixyz=4
      iSymZ = 2**IrrFnc(ixyz)
      ixyz=3
      iSymXY = 2**IrrFnc(ixyz)
      ixyz=5
      iSymXZ = 2**IrrFnc(ixyz)
      ixyz=6
      iSymYZ = 2**IrrFnc(ixyz)
      ixyz=7
      iSyXYZ = 2**IrrFnc(ixyz)
*
      Call mma_allocate(ips,nComp,Label='ips')
      Call mma_allocate(lOper,nComp,Label='lOper')
      Call mma_allocate(kOper,nComp,Label='kOper')
      Call mma_allocate(C_Coor,3,nComp,Label='CCoor')
*
      Save_tmp=PrPrt
      PrPrt=.True.
*
      RepHlp=Zero
      Do iGrid = 1, nGrid_
         Write (Label,'(A,I5)') 'EF ',iGrid
         call dcopy_(3,Grid(1,iGrid),1,CCoor,1)
         iSymC = 1
         If (Ccoor(1).ne.Zero) iSymC = iOr(iSymC,iSymX)
         If (Ccoor(2).ne.Zero) iSymC = iOr(iSymC,iSymY)
         If (Ccoor(3).ne.Zero) iSymC = iOr(iSymC,iSymZ)
         If (Ccoor(1).ne.Zero .and. Ccoor(2).ne.Zero)
     &      iSymC = iOr(iSymC,iSymXY)
         If (Ccoor(1).ne.Zero .and. Ccoor(3).ne.Zero)
     &      iSymC = iOr(iSymC,iSymXZ)
         If (Ccoor(2).ne.Zero .and. Ccoor(3).ne.Zero)
     &      iSymC = iOr(iSymC,iSymYZ)
         If (Ccoor(1).ne.Zero .and. Ccoor(2).ne.Zero .and.
     &       Ccoor(3).ne.Zero) iSymC = iOr(iSymC,iSyXYZ)
*
         iComp=0
         Do ix = nOrdOp, 0, -1
            Do iy = nOrdOp-ix, 0, -1
               iComp = iComp + 1
               iz = nOrdOp-ix-iy
               ixyz=0
               If (Mod(ix,2).ne.0) ixyz=iOr(ixyz,1)
               If (Mod(iy,2).ne.0) ixyz=iOr(ixyz,2)
               If (Mod(iz,2).ne.0) ixyz=iOr(ixyz,4)
               iSym = 2**IrrFnc(ixyz)
               If (Ccoor(iComp).ne.Zero ) iSym = iOr(iSym,1)
               lOper(iComp) = MltLbl(iSymC,iSym)
               kOper(iComp) = iChBas(iComp+1)
               call dcopy_(3,Ccoor,1,C_Coor(1,iComp),1)
            End Do
         End Do
*
         Call OneEl_Integrals(EFInt,EFMem,Label,
     &                        ips,lOper,nComp,C_Coor,
     &                        nOrdOp,rHrmt,kOper,Integrals)
*
         Eeldip=Zero
         Do iComp = 1, nComp
            ip=ips(iComp)
            iSmLbl = lOper(iComp)
*                                                                      *
************************************************************************
*                                                                      *
*-----------Compute properties directly from integrals
*
            nInt=n2Tri(iSmLbl)
            If (nInt.ne.0.and.Abs(DipMom(iComp,iGrid)).gt.1.0D-20) Then
               Call CmpInt(Integrals(ip),nInt,nBas,nIrrep,iSmLbl)
               If (nInt.ne.nh1) Then
                  Call WarningMessage(2,'Ener: nInt.ne.nh1')
                  Write (6,*) 'nInt=',nInt
                  Write (6,*) 'nh1=',nh1
                  Call Abend()
               End If
*
*------------- Accumulate contribution to h1
*
               alfa=Sig*DipMom(iComp,iGrid)
               Call DaXpY_(nInt,alfa,Integrals(ip),1,h1,1)
               EelDip=EelDip-alfa*DDot_(nh1,D_tot,1,Integrals(ip),1)
*
            End If
*
*                                                                      *
************************************************************************
*                                                                      *
*
         End Do ! iComp
*
*------- Deallocate memory for integral
*
         Call mma_deallocate(Integrals)
*
         RepHlp=RepHlp+EelDip*half
!        Write (6,*) 'Eeldip=',Eeldip,RepHlp
*
      End Do    ! iGrid
      RepNuc=RepNuc+RepHlp
#ifdef _DEBUGPRINT_
      Write(6,*)'RepHlp              =',RepHlp
#endif

*
      PrPrt=Save_tmp
*
      Call mma_deallocate(C_Coor)
      Call mma_deallocate(kOper)
      Call mma_deallocate(lOper)
      Call mma_deallocate(ips)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
