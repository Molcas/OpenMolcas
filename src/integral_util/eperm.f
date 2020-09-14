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
      Subroutine eperm(D_Tot,nDens,Ravxyz,Cavxyz,nCavxyz_,dEF,
     &                 Grid,nGrid_,Cord,MaxAto,Z_Nuc,xfEF)
************************************************************************
*                                                                      *
*     Object:                                                          *
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
      use external_centers
      use Symmetry_Info, only: iChBas
      Implicit Real*8 (a-h,o-z)
      External EFInt, EFMem
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "rctfld.fh"
#include "print.fh"
#include "stdalloc.fh"
      Real*8 D_Tot(nDens), Ravxyz(nCavxyz_), Cavxyz(nCavxyz_),
     &       dEF(4,nGrid_), Grid(3,nGrid_), Origin(3), CCoor(3),
     &       Cord(3,MaxAto), Z_Nuc(MaxAto),xfEF(4,nGrid_)
      Logical Save_tmp
      Character*8 Label
      Dimension FactOp(1), l_Oper(1)
      Integer, Allocatable:: ips(:), lOper(:), kOper(:)
      Real*8, Allocatable::  C_Coor(:,:), Nuc(:)
*                                                                      *
************************************************************************
*                                                                      *
*-----Statement Functions
*
      iOff(ixyz) = ixyz*(ixyz+1)*(ixyz+2)/6
*                                                                      *
************************************************************************
*                                                                      *
      iPrint=5
      If(.not.lRFCav) Goto 99  !Skip calculation of moments if no cavity
      Call FZero(Origin,3)
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute the multipole moment of the QC system,
*     both nuclear(1) and electronic(2).
*
*     Cavxyz: Multipole moments in cartesian basis
*     Ravxyz: temporary storage
*
      Call FZero(Cavxyz,nCavxyz_)
*
*---- 1) Compute M(nuc,nl), nuclear multipole moments
*
      Do iMax = 0, lMax
         ip = 1 + iOff(iMax)
         Call RFNuc(Origin,Ravxyz(ip),iMax)
      End Do

      If (iPrint.ge.99) Call RecPrt('Nuclear Multipole Moments',
     &                              ' ',Ravxyz,1,nCavxyz_)
*
*---- 2) Compute the electronic contribution to the charge distribution.
*
*     M(el,nl) =  - Sum(p,q) Dpq <p|M(nl)|q>
*
      nOpr=1
      FactOp=One
      l_Oper=1
*-----Reset array for storage of multipole moment expansion
      Do iMltpl = 1, lMax
         Do ix = iMltpl, 0, -1
            If (Mod(ix,2).eq.0) Then
               iSymX=1
            Else
               ixyz=1
               iSymX=2**IrrFnc(ixyz)
               If (Origin(1).ne.Zero) iSymX = iOr(iSymX,1)
            End If
            Do iy = iMltpl-ix, 0, -1
               If (Mod(iy,2).eq.0) Then
                  iSymY=1
               Else
                  ixyz=2
                  iSymY=2**IrrFnc(ixyz)
                  If (Origin(2).ne.Zero) iSymY = iOr(iSymY,1)
               End If
               iz = iMltpl-ix-iy
               If (Mod(iz,2).eq.0) Then
                  iSymZ=1
               Else
                  ixyz=4
                  iSymZ=2**IrrFnc(ixyz)
                  If (Origin(3).ne.Zero) iSymZ = iOr(iSymZ,1)
               End If
*
               iTemp = MltLbl(iSymX,MltLbl(iSymY,iSymZ,
     &                            nIrrep),nIrrep)
               l_Oper=iOr(l_Oper,iTemp)
            End Do
         End Do
      End Do
      If (iPrint.ge.19) Then
         Write (6,*) '1st order total density'
         lOff = 1
         Do iIrrep = 0, nIrrep-1
            n = nBas(iIrrep)*(nBas(iIrrep)+1)/2
            Write (Label,'(A,I1)')
     &       'Diagonal Symmetry Block ',iIrrep+1
            Call Triprt(Label,' ',D_Tot(lOff),nBas(iIrrep))
            lOff = lOff + n
         End Do
      End If
*

      Call Drv1_RF(FactOp,nOpr,D_tot,nh1,Origin,l_Oper,Cavxyz,lMax)
*
      If (iPrint.ge.99) Call RecPrt('Electronic Multipole Moments',
     &                              ' ',Cavxyz,1,nCavxyz_)
*
*---- Add nuclear MM expansion to the electronic one
*
      Call DaXpY_(nCavxyz_,One,Ravxyz,1,Cavxyz,1)
*
      If (iPrint.ge.99) Call RecPrt('Electronic+Nuclear Moments',
     &                              ' ',Cavxyz,1,nCavxyz_)


      If(lXF) Then
         Call XFMoment(lMax,Cavxyz,Ravxyz,nCavxyz_,Origin)
      EndIf

      If (iPrint.ge.99) Call RecPrt('Total Multipole Moments ',
     &                              ' ',Cavxyz,1,nCavxyz_)
*                                                                      *
************************************************************************
*                                                                      *
*---- Loop over the Langevin grid, compute EF and store in dEF.
*
*     The code here is initially with the loop over the grid outermost.
*     We might have to change this later on!
*
*     Be happy, don't worry!
*
 99   Continue
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
      call mma_allocate(Nuc,nComp,Label='Nuc')
      Call mma_allocate(C_Coor,3,nComp,Label='CCoor')
*
      Save_tmp=PrPrt
      PrPrt=.True.
*
      Do iGrid = 1, nGrid_
         Write (Label,'(A,I5)') 'EF ',iGrid
         call dcopy_(3,Grid(1,iGrid),1,Ccoor,1)
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
               lOper(iComp) = MltLbl(iSymC,iSym,nIrrep)
               kOper(iComp) = iChBas(iComp+1)
               call dcopy_(3,Ccoor,1,C_Coor(1,iComp),1)
            End Do
         End Do
*

         Call EFNuc(C_Coor,Z_Nuc,Cord,MaxAto,Nuc,nOrdOp)
         Call OneEl_Property(EFInt,EFMem,Label,
     &                       ips,lOper,nComp,C_Coor,
     &                       nOrdOp,Nuc,rHrmt,kOper,
     &                       D_Tot,nDens,dEF(1,iGrid),Sig)

*        Field contribution from XF
         Call EFXF(C_Coor,XF,nXF,nOrd_XF,iXPolType,
     &        xfEF(1,iGrid),XMolnr,nXMolnr,iGrid,scal14)
*
      End Do

*     Add XF contribution to the total field
      Call DaXpY_(4*nGrid,One,xfEF,1,dEF,1)

      PrPrt=Save_tmp
*
      Call mma_deallocate(C_Coor)
      Call mma_deallocate(Nuc)
      Call mma_deallocate(kOper)
      Call mma_deallocate(lOper)
      Call mma_deallocate(ips)
*                                                                      *
************************************************************************
*                                                                      *
*---- Loop over grid points and compute square norm of EF. Store the
*     largest of the latter in fmax (used in solver).
*
      fmax=Zero
      Do iGrid = 1, nGrid_
         ftest=dEF(1,iGrid)**2
     &        +dEF(2,iGrid)**2
     &        +dEF(3,iGrid)**2
         dEF(4,iGrid)=ftest
         fmax=Max(fmax,ftest)
      End Do ! iGrid

c      Call RecPrt('eperm: dEF ',' ',dEF,4,nGrid_)
*
      Return
      End
