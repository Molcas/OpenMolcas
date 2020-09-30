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
      Subroutine edip(Ravxyz,Cavxyz,lMax_,
     &                EF,DipMom,dEF,PolEff,DipEff,Grid,nGrid_,
     &                nPolComp,nAnisopol,nXF,iXPolType,nXMolnr,XMolnr)

************************************************************************
*                                                                      *
*     Object: to solve equation system iteratively.                    *
*                                                                      *
*     Input:                                                           *
*            dEF   : the electric field of the QM system               *
*            Cavxyz: the MM expansion of the QM system                 *
*            DipEff: Effective dipole moments                          *
*            PolEff: Effective polarizabilities                        *
*            Grid  : list of grid centers                              *
*            nGird_: effective list length                             *
*                                                                      *
*     Scratch:                                                         *
*            Ravxyz: incremental charge distribution on the boundary   *
*                    of the cavity                                     *
*                                                                      *
*     Output:                                                          *
*            EF    : Total EF                                          *
*            DipMom: Langevin dipole moments on the grid               *
*                                                                      *
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
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "rctfld.fh"
*
      Real*8 Ravxyz((lMax+1)*(lMax+2)*(lMax+3)/6),
     &       Cavxyz((lMax+1)*(lMax+2)*(lMax+3)/6)
      Real*8 Grid(3,nGrid_), EF (4,nGrid_), DipMom   (3,nGrid_),
     &       dEF(4,nGrid_), PolEff(nPolComp,nGrid_), DipEff(nGrid_)
      Integer XMolnr(nXMolnr,nXF)
      Logical NonEq,lExcl
*
#ifdef _DEBUG_
      Call RecPrt('edip: dEF(permanent) ',' ',dEF,4,nGrid_)
      Call RecPrt('edip: PolEff ',' ',PolEff,nPolComp,nGrid_)
      Call RecPrt('edip: DipEff ',' ',DipEff,1,nGrid_)
      Call RecPrt('edip: Grid ',' ',Grid,3,nGrid_)
      write(6,*)
     &'nGrid_,nPolComp,nAnisopol,tk,dampIter,dipCutoff,clim,lDamping',
     & nGrid_,nPolComp,nAnisopol,tk,dampIter,dipCutoff,clim,lDamping
      do i=1,nGrid_
         write(6,*) 'EDOTr ', i,Grid(1,i)*dEF(1,i)+
     &        Grid(2,i)*dEF(2,i)+Grid(3,i)*dEF(3,i)
      EndDo
#endif

      nCavxyz_=(lMax+1)*(lMax+2)*(lMax+3)/6
      qqo = Zero

      NonEq=.False.

*
#ifdef _DEBUG_
      Write (6,*)
      Write (6,*) 'Iter fmax             testa'
#endif
      Iter=0
555   testa=fmax*afac
      Iter=Iter+1
*
*---- Loop over Langevin grid and make EF and dipol moments at the
*     grid self consistent.
*
      Do iGrid = 1, nGrid_

         fx=dEF(1,iGrid)+EF(1,iGrid)
         fy=dEF(2,iGrid)+EF(2,iGrid)
         fz=dEF(3,iGrid)+EF(3,iGrid)
         ftot=fx*fx+fy*fy+fz*fz
*------- Update EF and square norm
*
         EF(1,iGrid)=fx
         EF(2,iGrid)=fy
         EF(3,iGrid)=fz
         EF(4,iGrid)=ftot
*
*------- Reset update vector
*
         dEF(1,iGrid)=Zero
         dEF(2,iGrid)=Zero
         dEF(3,iGrid)=Zero
         dEF(4,iGrid)=Zero

      EndDo
      Do iGrid = 1, nGrid_
         fx=EF(1,iGrid)
         fy=EF(2,iGrid)
         fz=EF(3,iGrid)
         ftot=EF(4,iGrid)
*
*------- Skip if square norm below threshold
*
c         If (dEF(4,iGrid).lt.testa) Go To 666
*
         ghx=Grid(1,iGrid)
         ghy=Grid(2,iGrid)
         ghz=Grid(3,iGrid)

*
*------- Pick up dipole moment at grid point
*
         dx=DipMom(1,iGrid)
         dy=DipMom(2,iGrid)
         dz=DipMom(3,iGrid)

*

c         Dip_Eff=DipEff(iGrid)*DBLE(Min(Iter,100))/100.0D0
         Dip_Eff=DipEff(iGrid)
*
*------- Compute new dipole moment as a function of the EF, effective dipole
*        moment and effective polarizability.
*
         If (Dip_Eff.lt.1.0D-10) Then
            If(iGrid.gt.nAnisoPol) Then   ! isotropic
               DipMom(1,iGrid)=fx*PolEff(1,iGrid)
               DipMom(2,iGrid)=fy*PolEff(1,iGrid)
               DipMom(3,iGrid)=fz*PolEff(1,iGrid)
            Else  ! anisotropic
               DipMom(1,iGrid)=fx*PolEff(1,iGrid)+
     &              fy*PolEff(2,iGrid)+fz*PolEff(3,iGrid)
               DipMom(2,iGrid)=fx*PolEff(2,iGrid)+
     &              fy*PolEff(4,iGrid)+fz*PolEff(5,iGrid)
               DipMom(3,iGrid)=fx*PolEff(3,iGrid)+
     &              fy*PolEff(5,iGrid)+fz*PolEff(6,iGrid)
            EndIf
         Else   ! NB!! Only isotropic implemented for dip>0
            fftots=sqrt(ftot)
            ftots=One/fftots
            x=Dip_Eff*tk*fftots
            ex=exp(x)
            emx=One/ex
            alang=(ex+emx)/(ex-emx)-One/x
c            alang=x/Three  !Linear approximation
            i=iGrid
            radabs=sqrt(Grid(1,i)*Grid(1,i)+Grid(2,i)*Grid(2,i)
     &           +Grid(3,i)*Grid(3,i))
            uind=Dip_Eff*alang+ftot*PolEff(1,iGrid)*ftots
            DipMom(1,iGrid)=uind*fx*ftots
            DipMom(2,iGrid)=uind*fy*ftots
            DipMom(3,iGrid)=uind*fz*ftots
         End If

*Grid
*
*------- Compute the change in the dipole moment between the old (dx,dy,dz) and
*        the new (DipMom).

*------- Try damping the change in dipole moment for better convergence
*
         DipMom(1,iGrid)=(One-dampIter)*DipMom(1,iGrid)+dampIter*dx
         DipMom(2,iGrid)=(One-dampIter)*DipMom(2,iGrid)+dampIter*dy
         DipMom(3,iGrid)=(One-dampIter)*DipMom(3,iGrid)+dampIter*dz

         dx=DipMom(1,iGrid)-dx
         dy=DipMom(2,iGrid)-dy
         dz=DipMom(3,iGrid)-dz
*
*------- Given the charge (qqo=0.0) and the change of the dipole moment
*        at this point modify the multipole expansion around the origin
*        accordingly. On the first iteration we will have the MM of the
*        QM in Cavxyz too, in subsequential iterations we will only deal
*        with incremental contributions.
*

         Call qlm(ghx,ghy,ghz,qqo,dx,dy,dz,lMax,Cavxyz)
*
*------- Loop over the whole grid and update the EF due to the change of
*        the dipole moment at the grid point "iGrid".
*
         Tr1=Zero
         If(lDamping) Then
            If(iGrid.gt.nAnisopol) Then
               Tr1=PolEff(1,iGrid)
            Else
               Tr1=(PolEff(1,iGrid)+PolEff(4,iGrid)+PolEff(6,iGrid))
     &              /Three
            EndIf
         EndIf
         Do jGrid = 1, nGrid_
            If (iGrid.eq.jGrid) Go To 777
            scal=One
            If(lAmberpol.and.(iXPolType.gt.0).and.(iGrid.le.nXF)
     &           .and.(jGrid.le.nXF)) Then
               lExcl=.False.
               Do i=1,nXMolnr
                  If(XMolnr(1,jGrid).eq.XMolnr(i,iGrid)) lExcl=.True.
                  If(XMolnr(1,jGrid).eq.-XMolnr(i,iGrid)) scal=scal14
               EndDo
               If(lExcl) Then
*     exclude field from iGrid when calculating the field at jGrid
c                  Write(6,*)'EXCLUDE dip', iGrid, ' at ', jGrid
                  Goto 777
               Elseif (scal.lt.One) Then
c                  Write(6,*)'SCALE dip', iGrid, ' at ', jGrid,
c     &                 ' with ', scal

               EndIf
            EndIf
            ghx1=Grid(1,jGrid)
            ghy1=Grid(2,jGrid)
            ghz1=Grid(3,jGrid)
            rx=ghx-ghx1
            ry=ghy-ghy1
            rz=ghz-ghz1
            r2=(rx*rx+ry*ry+rz*rz)
            If(r2.lt.dipCutoff**2) Go To 777

            r2i=One/r2
            ska=dx*rx+dy*ry+dz*rz
            disti=sqrt(r2i)
            dist3=r2i*disti
            temp=Three*ska*r2i
            If(lDamping) Then
               If(jGrid.gt.nAnisopol) Then
                  Tr2=PolEff(1,jGrid)
               Else
                  Tr2=(PolEff(1,jGrid)+PolEff(4,jGrid)+PolEff(6,jGrid))
     &                 /Three
               EndIf
               s = 2.3268D0*(Tr1*Tr2)**(1.0/6.0)
               v = min(1.0D0,sqrt(r2)/s)
               d1 = 4.0*v**3 - 3.0*v**4
               d2 = v**4
c               Write(6,*)'DAMP', d1, d2, Tr1, Tr2, sqrt(r2)
               dEF(1,jGrid)=dEF(1,jGrid)-(dx*d1-temp*rx*d2)*dist3*scal
               dEF(2,jGrid)=dEF(2,jGrid)-(dy*d1-temp*ry*d2)*dist3*scal
               dEF(3,jGrid)=dEF(3,jGrid)-(dz*d1-temp*rz*d2)*dist3*scal
            Else
               dEF(1,jGrid)=dEF(1,jGrid)-(dx-temp*rx)*dist3*scal
               dEF(2,jGrid)=dEF(2,jGrid)-(dy-temp*ry)*dist3*scal
               dEF(3,jGrid)=dEF(3,jGrid)-(dz-temp*rz)*dist3*scal
            EndIf
 777        Continue
         End Do           ! jGrid
*
c666     Continue
      End Do           ! iGrid



      If(lRFCav) Then
*
*---- Compute the charge distribution on the boundary of the cavity due to the
*     MM expansion at origin.
*
         call dcopy_(nCavxyz_,Cavxyz,1,Ravxyz,1)

         Call AppFld(Ravxyz,rds,Eps,lMax,EpsInf,NonEq)


*
*---- Compute EF at the grid due to the charge distribution in MM expansion
*     for the QM system plus the dipole moments on the grid.
*
         Do iGrid = 1, nGrid_
            ghx1=Grid(1,iGrid)
            ghy1=Grid(2,iGrid)
            ghz1=Grid(3,iGrid)
            fax=Zero
            fay=Zero
            faz=Zero
*
*------- Given the charge distribution on the boundary of the cavity
*        compute EF at (ghx1,ghy1,ghz1).
*
            Call hmod(ghx1,ghy1,ghz1,v_dummy,fax,fay,faz,Ravxyz,lmax)
*
*------- Accumulate in update vector
*
            dEF(1,iGrid)=dEF(1,iGrid)+fax
            dEF(2,iGrid)=dEF(2,iGrid)+fay
            dEF(3,iGrid)=dEF(3,iGrid)+faz
         EndDo                  !iGrid

      EndIf ! if (lRFCav)

      fmax=Zero
      Do iGrid = 1, nGrid_
         ftest=dEF(1,iGrid)**2
     &        +dEF(2,iGrid)**2
     &        +dEF(3,iGrid)**2
         dEF(4,iGrid)=ftest
         fmax=Max(ftest,fmax)
      End Do          ! iGrid
*
      Call FZero(Cavxyz,nCavxyz_)
*
*---- Check convergence
*
c      Call RecPrt('DipMom ',' ',DipMom,3,nGrid_)

#ifdef _DEBUG_
      Write (6,*) Iter,fmax,testa
#endif
      If (fmax.gt.clim) Go To 555
*
*---- Now we have a MM from QM + Langevin grid which is consistent with the
*     charge distribution on the boundary of the cavity. The Langevin
*     distribution of dipole moments is also internally consistent!
*

#ifdef _DEBUG_
      Call RecPrt('edip: converged DipMom ',' ',DipMom,3,nGrid_)

*     Write out dipoles and a pointcharge representation of the dipoles
      Write(6,*)'QREP'
      do i=1,nGrid_
         dipabs=sqrt(DipMom(1,i)**2+DipMom(2,i)**2+DipMom(3,i)**2)
         del=0.01
         Write(6,*)Grid(1,i)+DipMom(1,i)/dipabs*del,
     &             Grid(2,i)+DipMom(2,i)/dipabs*del,
     &             Grid(3,i)+DipMom(3,i)/dipabs*del, dipabs/del/2.0
         Write(6,*)Grid(1,i)-DipMom(1,i)/dipabs*del,
     &             Grid(2,i)-DipMom(2,i)/dipabs*del,
     &             Grid(3,i)-DipMom(3,i)/dipabs*del, -dipabs/del/2.0
      EndDo

      do i=1,nGrid_
         ddotr=Grid(1,i)*DipMom(1,i)+
     &        Grid(2,i)*DipMom(2,i)+Grid(3,i)*DipMom(3,i)
         dipabs=sqrt(DipMom(1,i)*DipMom(1,i)+DipMom(2,i)*DipMom(2,i)
     &        +DipMom(3,i)*DipMom(3,i))
         radabs=sqrt(Grid(1,i)*Grid(1,i)+Grid(2,i)*Grid(2,i)
     &        +Grid(3,i)*Grid(3,i))
         write(6,*)'RADPOL',radabs,dipabs/(scala*scalb*scalc),
     &        ddotr/(dipabs*radabs)
      EndDo
#endif

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(lMax_)
      End
