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
      Subroutine GenRadQuad_PAM(iNQ,nR_Eff,mr,Alpha,Process,QuadR,
     &                          nQuadR)
      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "WrkSpc.fh"
#include "real.fh"
#include "nq_info.fh"
#include "debug.fh"
      Real*8 Alpha(2), QuadR(2,nQuadR)
      Real*8 mr(2), ln_rn
      Logical Process
*
*---- Last point at infinity is eliminated
*
      If (Debug) Write (6,*) 'New Algorithm (Malmqvist)'
*
*-----Reading of the datas
      Alpha_Min=Alpha(1)
      Alpha_Max=Alpha(2)
      l_Max=2*Int(mr(1))
      If (Debug) Write (6,*) 'l_Max=',l_Max
      Relative_Max_Error=mr(2)
      If (Debug) Write (6,*) 'Relative_Max_Error=',Relative_Max_Error
*
*     Compute an approximative R_D_0
*
       Do k = 0, l_Max, l_Max-1
      R_D_0=Relative_Max_Error/(10.0D0**k)
      Dr=-Log10(R_D_0)
*
*     Starting value of h
*
      h=One/(0.47D0*Dr+0.93D0)
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute a correct h for the approximative R_D_0
*
      C1=Four*Sqrt(Two)*Pi
      C2=Pi**2/Two
 99   Continue
      h_=C2/( -Log( Ten**(-Dr)*h / C1 ) )
      If (Abs(h_-h).gt.1.0D-4) Then
         h=h_
         Go To 99
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Now find h from the correct R_D, i.e. from the highest
*     angular momentum available.
*
      Dr=-Log10(Relative_Max_Error)
 98   Continue
      h_=C2/(
     &        -Log(
     &             Ten**(-Dr)*(h/C1)*(h/Pi)**(DBLE(k)/Two)
     &             *(G((DBLE(k)+Three)/Two)/G(Three/Two))
     &            )
     &       )

      If (Debug) Write (6,*) 'h h_ ',h, h_
      If (Abs(h_-h).gt.1.0D-5) Then
         h=h_
         Go To 98
      End If
*
      If (k.eq.0) h0=h
       End Do
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute table of R_Max as a function of l
*
      Do i = l_Max, 0, -2
         D_m=-4.0D0
         If (l_Max.eq. 4) D_m=-2.3D0
         If (l_Max.eq. 2) D_m=-1.0D0
         If (l_Max.eq. 0) D_m= 1.9D0
         If (l_Max.eq.-2) D_m= 9.1D0
*
         ggg=(Two/(DBLE(i)+Three))*(D_m-Log(One/Ten**(-Dr)))
         R_Max(i)=Sqrt(Exp(ggg)/Alpha_Max)
         If (Debug) Then
            write(6,*) 'i        =',i
            write(6,*) 'l_Max    =',l_Max
            write(6,*) 'ggg      =',ggg
            write(6,*) 'R_Max(i) =',R_Max(i)
         End If
      End Do
      If (Debug) Write (6,*) 'h0,h=',h0,h
*
*     For hybrid grid use R_Max for l=0 and h for l=l_max
*     r1=R_Max(l_Max)
      r1=R_Max(0)
      ln_rn=1.7D0-log(Alpha_Min)/Two
      rn=Exp(ln_rn)
      gamma=r1/(Exp(h)-One)
      n_High=Int(Log(rn/gamma+One)/h+One)
      If (Debug) Then
        Write (6,*)
        Write (6,*) 'r1,Alpha_Min    =',r1,Alpha_Min
        Write (6,*) 'rn,Alpha_Max    =',rn,Alpha_MAx
        Write (6,*) 'h,Dr,n_High     =',h,Dr,n_High
      End If
*
*-----Store the radius and the associated weights
*
      If (Debug) Write(6,*) 'n_High',n_High
      iR = 0
      Do k = 0, n_High
         a = DBLE(k)*h
         rk = Gamma*(Exp(a)-One)
*----Note that the point at r=0 is eliminated
         If (rk.ne.Zero) Then
           iR = iR + 1
           If (Process) Then
              QuadR(1,ir)=rk
              Correction=One
*
*             Gregorious correction for points close to the nuclei
*
              If (k.eq.0) Correction= 46.D0/120.D0
              If (k.eq.1) Correction=137.D0/120.D0
              If (k.eq.2) Correction=118.D0/120.D0
              If (k.eq.3) Correction=119.D0/120.D0
*
              QuadR(2,ir) = h * ( rk + gamma ) * Correction
              QuadR(2,ir) = QuadR(1,ir)**2 * QuadR(2,ir)
           End If
         End If
*
      End Do
*-----Store the value of the maximum radius to which we should integrate
*       for the partitionning and the number of effective radii.
      nR_Eff = iR
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(iNQ)
      End
      Function G(Arg)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 G
      g=-1000.0D0
*
      Arg_=DBLE(Int(Arg))
      If (Abs(Arg-Arg_).lt.Half/Two) Then
*        Integer argument
         G=One
         rG=One
      Else
*        fractional argument
         G=Sqrt(Pi)
         rG=Half
      End If
*
 99   Continue
         If (Abs(rG-Arg).lt.Half/Two) goto 666
         G=rG*G
         rG=rG+One
      Go To 99
666   continue
      return
*
      End
      Function RD(m,h)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 RD, h
*
      RD=(Four*Sqrt(Two)*Pi/h)*Exp(-Pi**2/(Two*h))
      If (m.eq.0) Return
      RD=(G(Three/Two)/G((DBLE(m)+Three)/Two))
     &  *(Pi/h)**(DBLE(m)/Two) * RD
*
      Return
      End
