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
      Subroutine MAGN( EXCH,N,X,Y,Z,H,W,zJ,THRS,dM,sM,nT,T,sopt,
     &                 WZ,ZB,S,M, m_paranoid, DBG)
c this Subroutine is a wrapper for various MAGN subroutines

      Implicit None
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)          :: EXCH, N, nT
      Real(kind=8), intent(in)    :: X, Y, Z, H, zJ
      Real(kind=8), intent(in)    :: W(EXCH), T(nT)
      Complex(kind=8), intent(in) :: dM(3,EXCH,EXCH)
      Complex(kind=8), intent(in) :: sM(3,EXCH,EXCH)
      Logical, intent(in)          :: sopt

      Real(kind=8), intent(out)   :: ZB(nT), WZ(N)
      Real(kind=8), intent(out)   :: S(3,nT), M(3,nT)

      Real(kind=8), intent(in)    :: THRS
      Logical, intent(in)          :: m_paranoid
      Logical, intent(in)          :: DBG
c local variables:


      If( abs(zJ) .lt. tiny(0.0_wp) ) Then


         If(DBG) Write(6,*) 'Enter MAGN_NO_MF :'

         Call MAGN_NO_MF( EXCH, N, X,Y,Z, H, W, dM, sM, nT, T, sopt,
     &                    WZ, ZB, S, M, DBG )

         If(DBG) Write(6,*) 'Exit MAGN_NO_MF :'


      Else ! zJ .ne. 0.0_wp


         If(DBG) Write(6,*) 'Enter MAGN_ZJ_PAR :'

         Call MAGN_ZJ_PAR( EXCH, N, X,Y,Z, H, W, zJ, dM, sM,
     &                     nT, T, sopt, WZ, ZB, S, M, thrs,
     &                     m_paranoid, DBG )

         If(DBG) Write(6,*) 'Exit MAGN_ZJ_PAR :'


      End If



      Return
      End subroutine magn

