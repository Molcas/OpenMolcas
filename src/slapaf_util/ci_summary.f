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
* From: D.R. Yarkony, J. Phys. Chem. A 105 (2001) 6277-6293
* and: J. Chem. Theory Comput. 12 (2016) 3636-3653
************************************************************************
      Subroutine CI_Summary(Lu)
      use Slapaf_Info, only: Gx, Gx0, NAC, Energy
      Implicit None
      Integer Lu, n, i
      Real*8, Dimension(:), Allocatable :: g, h, tmp
      Real*8 gg, hh, gh, sg, sh, dgh, deltagh, beta_ang, norm_g, norm_h,
     &       st, srel, shead, peaked, bif, aux
      Real*8, External :: dDot_
      Character(Len=2) LabA
      Character(Len=40) Description
#include "info_slapaf.fh"
#include "stdalloc.fh"
#include "nadc.fh"
#include "real.fh"
*
      n=3*SIZE(Gx,2)
      Call mma_Allocate(g,n)
      Call mma_Allocate(h,n)
*
*     Compute the orthogonal branching vectors
*     Note that d(E1-E0)/dx is stored, but we want to use d((E1-E0)/2)/dx
*     (and forces instead of gradients)
*
      gg=dDot_(n,Gx0(1,1,iter),1,Gx0(1,1,iter),1)*Quart
      hh=dDot_(n,NAC,1,NAC,1)
      gh=-dDot_(n,Gx0(1,1,iter),1,NAC,1) !Factor 2 included
      beta_ang=Atan2(gh,gg-hh)*Half
      Call dCopy_(n,Gx0(1,1,iter),1,g,1)
      Call dScal_(n,-Half*Cos(beta_ang),g,1)
      Call dAxpY_(n,Sin(beta_ang),NAC,1,g,1)
      Call dCopy_(n,NAC,1,h,1)
      Call dScal_(n,Cos(beta_ang),h,1)
      Call dAxpY_(n,Half*Sin(beta_ang),Gx0(1,1,iter),1,h,1)
      gg=dDot_(n,g,1,g,1)
      hh=dDot_(n,h,1,h,1)
      norm_g=Sqrt(gg)
      norm_h=Sqrt(hh)
      If (norm_g.gt.1.0D-12) Then
        Call dScal_(n,One/norm_g,g,1)
      Else
        Call dCopy_(n,[Zero],0,g,1)
      End If
      If (norm_h.gt.1.0D-12) Then
        Call dScal_(n,One/norm_h,h,1)
      Else
        Call dCopy_(n,[Zero],0,h,1)
      End If
*     Ensure that the asymmetry will be positive
*     this fixes which vector is x and which is y
      If (hh.gt.gg) Then
        Call SwapVe(g,h,n)
        aux=gg
        gg=hh
        hh=aux
        aux=norm_g
        norm_g=norm_h
        norm_h=aux
      End If
      sg=-dDot_(n,Gx(1,1,iter),1,g,1)
      sh=-dDot_(n,Gx(1,1,iter),1,h,1)
*     Ensure that the tilt heading will be in the first quadrant
*     this fixes the signs of the x and y vectors
      If (sg.lt.Zero) Then
        sg=Abs(sg)
        Call DScal_(n,-One,g,1)
      End If
      If (sh.lt.Zero) Then
        sh=Abs(sh)
        Call DScal_(n,-One,h,1)
      End If
      st=Sqrt(sg**2+sh**2)
      dgh=Sqrt((gg+hh)/Two)
      LabA=''
      deltagh=gg-hh
      If ((gg+hh).gt.1.0D-12) Then
        deltagh=deltagh/(gg+hh)
      Else
        LabA=' *'
      End If
      If (dgh.gt.1.0D-12) Then
        srel=st/dgh
      Else
        srel=Zero
      End If
      shead=Atan2(sh,sg)
*
*     peaked/sloped, bifurcating/single-path parameters
*
      peaked=srel**2/(One-deltagh**2)*(One-deltagh*Cos(Two*shead))
      bif=((One+deltagh)*Cos(shead)**2)**(One/Three)
      bif=bif+((One-deltagh)*Sin(shead)**2)**(One/Three)
      bif=(srel/(Two*deltagh))**(Two/Three)*bif
      Description=''
      If (peaked.lt.1) Then
        Description=Trim(Description)//'peaked (P<1)'
      Else If (peaked.gt.1) Then
        Description=Trim(Description)//'sloped (P>1)'
      Else
        Description=Trim(Description)//'* (P=1)'
      End If
      If (bif.lt.1) Then
        Description=Trim(Description)//' bifurcating (B<1)'
      Else If (bif.gt.1) Then
        Description=Trim(Description)//' single-path (B>1)'
      Else
        Description=Trim(Description)//' * (B=1)'
      End If
*
*     Disable Last_Energy to prevent further rotations
*
      CallLast=.False.
*
      Write(Lu,*)
      Call CollapseOutput(1,'Conical Intersection Characterization')
      Write(Lu,'(3X,A)')    '-------------------------------------'
      Write(Lu,*)
      Write(Lu,*) 'See: J. Chem. Theory Comput. 12 (2016) 3636-3653'
      Write(Lu,*)
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Call RecPrt('Gradient difference','',Gx0(1,1,iter),n,1)
      Call RecPrt('Coupling vector','',NAC,n,1)
      Call RecPrt('Average gradient','',Gx(1,1,iter),n,1)
      Write(Lu,100) 'Beta angle:',beta_ang
      Write(Lu,*)
#endif
      Write(Lu,100) 'Pitch (delta_gh):',dgh,' Eh/a0'
      Write(Lu,100) 'Asymmetry (Delta_gh):',deltagh,Trim(LabA)
#ifdef _DEBUGPRINT_
      Write(Lu,100) 'Total tilt (s):',st,' Eh/a0'
#endif
      Write(Lu,100) 'Relative tilt (sigma=s/delta_gh):',srel
      Write(Lu,100) 'Tilt heading (theta_s):',shead
      Write(Lu,*)
      Write(Lu,101) 'P:',peaked
      Write(Lu,101) 'B:',bif
      Write(Lu,101) 'Type: '//Trim(Description)
      Write(Lu,*)
      Write(Lu,*) 'Local linear representation:'
      Call mma_Allocate(tmp,n)
      Do i=1,SIZE(Gx,2)
        tmp(0*SIZE(Gx,2)+i) = g((i-1)*3+1)
        tmp(1*SIZE(Gx,2)+i) = g((i-1)*3+2)
        tmp(2*SIZE(Gx,2)+i) = g((i-1)*3+3)
      End Do
      Call RecPrt('Local x','',tmp,SIZE(Gx,2),3)
      Do i=1,SIZE(Gx,2)
        tmp(0*SIZE(Gx,2)+i) = h((i-1)*3+1)
        tmp(1*SIZE(Gx,2)+i) = h((i-1)*3+2)
        tmp(2*SIZE(Gx,2)+i) = h((i-1)*3+3)
      End Do
      Call RecPrt('Local y','',tmp,SIZE(Gx,2),3)
      Write(Lu,*)
      Write(Lu,110) Energy(iter),sg,sh
      Write(Lu,120) Two*dgh,deltagh
      Call mma_Deallocate(tmp)
      Call CollapseOutput(0,'Conical Intersection Characterization')
100   Format (5X,A,T40,ES12.5,A)
101   Format (5X,A,T11,ES12.5)
110   Format (5X,'Average energy: ',F15.8,' + ',F12.8,'*x + ',
     &        F12.8,'*y')
120   Format (5X,
     &        'Energy difference: ',F12.8,'*sqrt(r^2 + ',F12.8,'*t)',
     &        /,10X,'r^2 = x^2 + y^2',/,10X,'t = x^2 - y^2')
*
      Call mma_Deallocate(g)
      Call mma_Deallocate(h)
      Return
*
      End Subroutine CI_Summary
