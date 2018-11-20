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
      Subroutine ZEEM_SA(N, H, dX,dY,dZ, W, M, sM, S, zJ,  WM,ZM,
     &                DBG, RWORK, HZEE, WORK, W_c )
c
      Implicit None
      Integer, parameter       :: wp=SELECTED_REAL_KIND(p=15,r=307)
c input variables:
      Integer,         intent(in) :: N
      Real(kind=wp),   intent(in) :: H,dX,dY,dZ,zJ
      Real(kind=wp),   intent(in) :: W(N)
      Real(kind=wp),   intent(in) :: S(3)
      Complex(kind=wp),intent(in) :: sM(3,N,N)
      Complex(kind=wp),intent(in) ::  M(3,N,N)
c output variables:
      Real(kind=wp),    intent(out) ::  WM(N)
      Complex(kind=wp), intent(out) :: ZM(N,N)
c local variables:
      Integer          :: i,j,info
      Real(kind=wp)    :: mB
      Real(kind=wp)    :: RWORK(3*N-2)
      Complex(kind=wp) :: HZEE(N*(N+1)/2)
      Complex(kind=wp) :: WORK(2*N-1), R, P, RP
      Complex(kind=wp) :: H_c, dX_c, dY_c, dZ_c, zJ_c, W_c(N), S_c(3)
      Complex(kind=wp) :: mB_c
      Logical          :: DBG
      Call qEnter('ZEEM')
      mB = 0.4668643740_wp !   in cm-1*T-1

      ! initialization

c      Call dcopy_(     N       ,  0.0_wp        , 0,    WM, 1)
c      Call dcopy_(  (3*N-2)    ,  0.0_wp        , 0, RWORK, 1)
c      Call zcopy_(   N**2      , (0.0_wp,0.0_wp), 0,    ZM, 1)
c      Call zcopy_(  (N*(N+1)/2), (0.0_wp,0.0_wp), 0,  HZEE, 1)
c      Call zcopy_(  (2*N-1)    , (0.0_wp,0.0_wp), 0,  WORK, 1)
         WM=0.0_wp
      RWORK=0.0_wp
        ZM=(0.0_wp,0.0_wp)
      HZEE=(0.0_wp,0.0_wp)
      WORK=(0.0_wp,0.0_wp)

      info=0
      If(DBG) then
        Write(6,'(A,4ES20.10)') 'dX,dY,dZ,H =',dX,dY,dZ,H
        Write(6,'(A,4ES20.10)') 'Sx,y,z=',S(1),S(2),S(3)
        Write(6,*) 'zJ = ', zJ
      End If

      ! initialize
      H_c =CMPLX( H,0.0_wp,wp)
      dX_c=CMPLX(dX,0.0_wp,wp)
      dY_c=CMPLX(dY,0.0_wp,wp)
      dZ_c=CMPLX(dZ,0.0_wp,wp)
      zJ_c=CMPLX(zJ,0.0_wp,wp)
      mB_c=CMPLX(mB,0.0_wp,wp)
      Call zcopy_(N, (0.0_wp,0.0_wp), 0,  W_c, 1)
      Call zcopy_(3, (0.0_wp,0.0_wp), 0,  S_c, 1)
      Do i=1,N
         W_c(i)=CMPLX(W(i),0.0_wp,wp)
      End Do
      Do i=1,3
         S_c(i)=CMPLX(S(i),0.0_wp,wp)
      End Do

      If(DBG) Write(6,*)' H_c = ', H_c
      If(DBG) Write(6,*)'dX_c = ',dX_c
      If(DBG) Write(6,*)'dY_c = ',dY_c
      If(DBG) Write(6,*)'dZ_c = ',dZ_c
      If(DBG) Write(6,*)'zJ_c = ',zJ_c
      If(DBG) Write(6,*)'mB_c = ',mB_c

      ! build the Zeeman Hamiltonian
      If ( zJ==0.0_wp ) Then

        Do i=1,N
          Do j=1,i
            R = (0.0_wp,0.0_wp)
            R =  dX_c * M(1,j,i)
     &         + dY_c * M(2,j,i)
     &         + dZ_c * M(3,j,i)

            HZEE( j+(i-1)*i/2 ) = HZEE( j+(i-1)*i/2 ) - mB_c * H_c * R

          End Do
        End Do

      Else ! zJ .ne. 0

        Do i=1,N
          Do j=1,i
            R = (0.0_wp,0.0_wp)
            P = (0.0_wp,0.0_wp)
            RP= (0.0_wp,0.0_wp)

            R =  dX_c *  M(1,j,i)
     &         + dY_c *  M(2,j,i)
     &         + dZ_c *  M(3,j,i)

            P =  dX_c * SM(1,j,i) * S_c(1)
     &         + dY_c * SM(2,j,i) * S_c(2)
     &         + dZ_c * SM(3,j,i) * S_c(3)

            RP = mB_c * H_c * R + zJ_c  * P

            HZEE( j+(i-1)*i/2 ) = HZEE( j+(i-1)*i/2 ) - RP
          End Do
        End Do

      End If ! zJ



      ! add diagonal energies:
      Do i=1,N
        HZEE(i+(i-1)*i/2) = HZEE(i+(i-1)*i/2) + W_c(i)
      End Do

      If(DBG) then
        Write(6,'(A)') 'HZEE:'
        Do i=1,N
          Do j=1,i
          Write(6,'(2i3,2x,100(2ES16.8,2x))') i, j, HZEE(j+(i-1)*i/2)
          End Do
        End Do
      End If
      ! diagonalization
      Call zhpev_('V','U',N, HZEE,WM,ZM, N, WORK,RWORK,INFO)

      If(DBG) then
        Do i=1,N
          Write(6,'(A,i3,A,F20.13,A,i2,A,99(2F16.10,1x))')
     &       'WM(',i,')=',WM(i)!,' ZM(j,',i,'):',(ZM(j,i),j=1,N)
        End Do
      End If
      Call qExit('ZEEM')
      Return
      End
