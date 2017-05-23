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
      Subroutine FndSph(NAt,ICharg,ToAng,C,IAt,ITypRad,NSphInp,
     $                  Alpha,XSph,YSph,ZSph,Rad,NOrd,iPrint)
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "rctfld.fh"
      Parameter (ISAX=1000)
      Dimension C(3,*),IAt(*)
      Dimension XSph(*),YSph(*),ZSph(*),Rad(*),NOrd(*)
      Data IOut /6/, Zero /0.0d0/
*
* Assign GEPOL sphere positions and radii according to solute atoms
* nature
*
      If(ITypRad.eq.1) then
* United Atom Topological Model (UATM) radii:
        Call GetMem('Chg','Allo','Real',ip_Chg,NAt)
        call dcopy_(NAt,Zero,0,Work(ip_chg),1)
        Call UATM(IOut,ICharg,NAt,NSinit,ToAng,Rad,Alpha,C,IAt,NOrd,
     &            Work(ip_Chg),iPrint)
        Call GetMem('Chg','Free','Real',ip_Chg,NAt)
      ElseIf(ITypRad.eq.2) then
* Pauling radii on each atom:
        Do 100 I = 1, NAt
          NOrd(I) = I
          Rad(I) = Pauling(IAt(I))
          Alpha = 1.2
 100      NSinit = NAt
      ElseIf(ITypRad.eq.3) then
* Sphere radii given in the input
        Do 200 I = 1, NSphInp
          NOrd(I) = NOrdInp(I)
          Rad(I) = RadInp(I)
          Alpha = 1.2
 200      NSinit = NSphInp
      Else
        write(6,'(''Unrecognized radii type !'')')
        Call Abend()
      EndIf
      If((ITypRad.eq.2.or.ITypRad.eq.3).and.iPrint.gt.5)
     &  Call PrtCav(IOut,ITypRad,NSinit,NOrd,Alpha,Rad)

      Do 300 I = 1, NSinit
        XSph(I) = C(1,NOrd(I))
        YSph(I) = C(2,NOrd(I))
        ZSph(I) = C(3,NOrd(I))
 300    Rad(I) = Rad(I) * Alpha

      Return
      End
      Subroutine PrtCav(IOut,ITyp,NS,NOrd,Alpha,Rad)
      Implicit real*8 (A-H,O-Z)
*
* Print out sphere radii for Pauling or input cavitites
*
      Dimension NOrd(*),Rad(*)
*
      Write(iOut,*)
      Write(iOut,*)
      Write(iOut,'(6X,A)') 'Polarized Continuum Model Cavity'
      Write(iOut,'(6X,A)') '================================'
      If(ITyp.eq.2) Write(iOut,'(6X,A)') 'Pauling radii'
      If(ITyp.eq.3) Write(iOut,'(6X,A)') 'Sphere radii from input'
      Write(iOut,*)
      Write(IOut,'(6X,'' Nord  Alpha  Radius'')')
      Do 100 IS = 1, NS
        Write(IOut,'(6X,1X,I3,3X,F4.2,3X,F5.3)')
     &   NOrd(IS),Alpha,Rad(IS)
 100  Continue
      Write(IOut,'(6X,1X,78(''-''))')
      Write(IOut,*)
      Return
      End
