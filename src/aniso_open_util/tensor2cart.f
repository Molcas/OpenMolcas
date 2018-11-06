      Subroutine tensor2cart(Jt,Jc)
!     Hexch = S1.Jsph.S2

      Implicit None
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Complex(kind=wp), intent(in) :: Jt(-1:1,-1:1)
      Real(kind=wp), intent(out)   :: Jc(3,3)
      ! local variables
      Complex(kind=wp) :: i,pp,pz,pm,zp,zz,zm,mp,mz,mm

      jc(1:3,1:3)=0.0_wp
       i=(0.0_wp,1.0_wp)
      pp=(0.0_wp,0.0_wp)
      pz=(0.0_wp,0.0_wp)
      pm=(0.0_wp,0.0_wp)
      zp=(0.0_wp,0.0_wp)
      zz=(0.0_wp,0.0_wp)
      zm=(0.0_wp,0.0_wp)
      mp=(0.0_wp,0.0_wp)
      mz=(0.0_wp,0.0_wp)
      mm=(0.0_wp,0.0_wp)

      pp=jt( 1, 1)
      pz=jt( 1, 0)
      pm=jt( 1,-1)
      zp=jt( 0, 1)
      zz=jt( 0, 0)
      zm=jt( 0,-1)
      mp=jt(-1, 1)
      mz=jt(-1, 0)
      mm=jt(-1,-1)

      jc(1,1)=0.5_wp * dble(    mm -  mp -  pm +  pp )
      jc(2,2)=0.5_wp * dble( -  mm -  mp -  pm -  pp )
      jc(1,2)=0.5_wp * dble( -i*mm -i*mp +i*pm +i*pp )
      jc(2,1)=0.5_wp * dble( -i*mm +i*mp -i*pm +i*pp )

      jc(1,3)=sqrt(0.5_wp) * dble(   mz-  pz )
      jc(3,1)=sqrt(0.5_wp) * dble(   zm-  zp )

      jc(2,3)=sqrt(0.5_wp) * dble(-i*mz-i*pz )
      jc(3,2)=sqrt(0.5_wp) * dble(-i*zm-i*zp )

      jc(3,3)=dble(zz)

      Return
      End Subroutine tensor2cart






      Subroutine tensor2cart_minus(Jt,Jc)
!     Hexch = -Jsph.S1.S2

      Implicit None
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Complex(kind=wp), intent(in) :: Jt(-1:1,-1:1)
      Real(kind=wp), intent(out)   :: Jc(3,3)
      ! local variables
      Complex(kind=wp) :: i,pp,pz,pm,zp,zz,zm,mp,mz,mm

      jc(1:3,1:3)=0.0_wp
       i=(0.0_wp,1.0_wp)
      pp=(0.0_wp,0.0_wp)
      pz=(0.0_wp,0.0_wp)
      pm=(0.0_wp,0.0_wp)
      zp=(0.0_wp,0.0_wp)
      zz=(0.0_wp,0.0_wp)
      zm=(0.0_wp,0.0_wp)
      mp=(0.0_wp,0.0_wp)
      mz=(0.0_wp,0.0_wp)
      mm=(0.0_wp,0.0_wp)

      pp=jt( 1, 1)
      pz=jt( 1, 0)
      pm=jt( 1,-1)
      zp=jt( 0, 1)
      zz=jt( 0, 0)
      zm=jt( 0,-1)
      mp=jt(-1, 1)
      mz=jt(-1, 0)
      mm=jt(-1,-1)

      jc(1,1)=0.5_wp * dble( -  mm +  mp +  pm -  pp )
      jc(2,2)=0.5_wp * dble(    mm +  mp +  pm +  pp )

      jc(1,2)=0.5_wp * dble(  i*mm -i*mp +i*pm -i*pp )
      jc(2,1)=0.5_wp * dble(  i*mm +i*mp -i*pm -i*pp )

      jc(1,3)=sqrt(0.5_wp) * dble(  -zm + zp )
      jc(3,1)=sqrt(0.5_wp) * dble(  -mz + pz )

      jc(2,3)=sqrt(0.5_wp) * dble( i*zm+i*zp )
      jc(3,2)=sqrt(0.5_wp) * dble( i*mz+i*pz )

      jc(3,3)=dble(-zz)

      Return
      End Subroutine tensor2cart_minus








      Subroutine cart2tensor(Jc,Jt)
! this subroutine computed tensorial parameters for
!     Hexch= S1.Jcart.S2

      Implicit None
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Real(kind=wp), intent(in)    :: Jc(3,3)
      Complex(kind=wp), intent(out):: Jt(-1:1,-1:1)
      ! local variables
      Complex(kind=wp) :: c2, cx2, i
      Complex(kind=wp) :: xx,xy,xz,yx,yy,yz,zx,zy,zz

      jt(-1:1,-1:1)=(0.0_wp,0.0_wp)

      i=(0.0_wp,1.0_wp)
      c2=(0.5_wp,0.0_wp)
      cx2=cmplx(sqrt(0.5_wp),0.0_wp,wp)

      xx=cmplx( jc(1,1), 0.0_wp, wp )
      xy=cmplx( jc(1,2), 0.0_wp, wp )
      xz=cmplx( jc(1,3), 0.0_wp, wp )
      yx=cmplx( jc(2,1), 0.0_wp, wp )
      yy=cmplx( jc(2,2), 0.0_wp, wp )
      yz=cmplx( jc(2,3), 0.0_wp, wp )
      zx=cmplx( jc(3,1), 0.0_wp, wp )
      zy=cmplx( jc(3,2), 0.0_wp, wp )
      zz=cmplx( jc(3,3), 0.0_wp, wp )

      jt( 1, 1)=c2*(  xx -i*xy -i*yx -yy )
      jt(-1,-1)=c2*(  xx +i*xy +i*yx -yy )
      jt( 1,-1)=c2*( -xx -i*xy +i*yx -yy )
      jt(-1, 1)=c2*( -xx +i*xy -i*yx -yy )

      jt( 1, 0)=cx2*(-xz+i*yz)
      jt(-1, 0)=cx2*( xz+i*yz)

      jt( 0, 1)=cx2*(-zx+i*zy)
      jt( 0,-1)=cx2*( zx+i*zy)

      jt( 0, 0)=zz

      Return
      End Subroutine cart2tensor




      Subroutine cart2tensor_minus(Jc,Jt)
! this subroutine computed tensorial parameters for
!     Hexch= -Jcart.S1.S2

      Implicit None
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Real(kind=wp), intent(in)    :: Jc(3,3)
      Complex(kind=wp), intent(out):: Jt(-1:1,-1:1)
      ! local variables
      Complex(kind=wp) :: c2, cx2, i
      Complex(kind=wp) :: xx,xy,xz,yx,yy,yz,zx,zy,zz

      jt(-1:1,-1:1)=(0.0_wp,0.0_wp)

      i=(0.0_wp,1.0_wp)
      c2=(0.5_wp,0.0_wp)
      cx2=cmplx(sqrt(0.5_wp),0.0_wp,wp)

      xx=cmplx( jc(1,1), 0.0_wp, wp )
      xy=cmplx( jc(1,2), 0.0_wp, wp )
      xz=cmplx( jc(1,3), 0.0_wp, wp )
      yx=cmplx( jc(2,1), 0.0_wp, wp )
      yy=cmplx( jc(2,2), 0.0_wp, wp )
      yz=cmplx( jc(2,3), 0.0_wp, wp )
      zx=cmplx( jc(3,1), 0.0_wp, wp )
      zy=cmplx( jc(3,2), 0.0_wp, wp )
      zz=cmplx( jc(3,3), 0.0_wp, wp )

      jt( 1, 1)=c2*( -xx +i*xy +i*yx +yy )
      jt(-1,-1)=c2*( -xx -i*xy -i*yx +yy )
      jt( 1,-1)=c2*(  xx -i*xy +i*yx +yy )
      jt(-1, 1)=c2*(  xx +i*xy -i*yx +yy )

      jt( 1, 0)=cx2*( zx-i*zy)
      jt(-1, 0)=cx2*(-zx-i*zy)

      jt( 0, 1)=cx2*( xz-i*yz)
      jt( 0,-1)=cx2*(-xz-i*yz)

      jt( 0, 0)= -zz

      Return
      End Subroutine cart2tensor_minus


