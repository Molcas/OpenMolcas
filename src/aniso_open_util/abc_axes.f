* $ this file belongs to the Molcas repository $
      Subroutine abc_axes(cryst, coord, xyz, abc, Do_option, iReturn )

      ! this Subroutine performs a transformation of the main axes
      ! (magnetic, anisotropic etc.) from a xyz system in
      ! crystallographic "abc" system, and vice-versa;
      !    Do_option = 1 =>  transform from xyz to abc
      !    Do_option = 2 =>  transform from abc to xyz
      !    If Do_option has other value, abort
      !    coord(3)-- Cartesian coordinates of the main magnetic center

      Implicit None
      Integer, parameter          :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)         :: Do_option
      Real(kind=wp), intent(inout):: xyz(3,3), abc(3,3)
      Real(kind=wp), intent(in)   :: cryst(6), coord(3)
      Integer, intent(out)        :: iReturn
      ! local variables:
      Integer       :: i
      Real(kind=wp) :: a,b,c,al,bt,gm,cal,cbt,cgm,sal,sbt,sgm,v,pi,x,y,z
      Real(kind=wp) :: xyz2(3,3),pX(3),pY(3),pZ(3)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! initializations
      pi  = 3.1415926535897932384626433832795028841971693993751_wp
      a   = 0.0_wp
      b   = 0.0_wp
      c   = 0.0_wp
      al  = 0.0_wp
      bt  = 0.0_wp
      gm  = 0.0_wp
      cal = 0.0_wp
      cbt = 0.0_wp
      cgm = 0.0_wp
      sal = 0.0_wp
      sbt = 0.0_wp
      sgm = 0.0_wp
      X   = 0.0_wp
      Y   = 0.0_wp
      Z   = 0.0_wp
      pX  = 0.0_wp
      pY  = 0.0_wp
      pZ  = 0.0_wp
      xyz2= 0.0_wp
!
      a   = cryst(1)
      b   = cryst(2)
      c   = cryst(3)
      al  = cryst(4)*pi/180._wp
      bt  = cryst(5)*pi/180._wp
      gm  = cryst(6)*pi/180._wp
      cal = cos(al)
      cbt = cos(bt)
      cgm = cos(gm)
      sal = sin(al)
      sbt = sin(bt)
      sgm = sin(gm)

      v = sqrt(1.0_wp-cal*cal-cbt*cbt-cgm*cgm+2.0_wp*cal*cbt*cgm )

      If ( Do_option .eq. 1 ) Then
         abc=0.0_wp
         Do i=1,3
            xyz2(1,i) = 1.0_wp*xyz(1,i)+coord(1)
            xyz2(2,i) = 1.0_wp*xyz(2,i)+coord(2)
            xyz2(3,i) = 1.0_wp*xyz(3,i)+coord(3)
         End Do

         Do i=1,3
            X=xyz2(1,i)
            Y=xyz2(2,i)
            Z=xyz2(3,i)

            pX(1)= 1.0_wp/a
            pY(1)=-cgm/(a*sgm)
            pZ(1)=((cal*cgm-cbt)/(a*v*sgm))

            pX(2)= 0.0_wp
            pY(2)= 1.0_wp/(b*sgm)
            pZ(2)= (cbt*cgm-cal)/(b*v*sgm)

            pX(3)= 0.0_wp
            pY(3)= 0.0_wp
            pZ(3)= sgm/(c*v)

            abc(1,i) = pX(1)*X + pY(1)*Y + pZ(1)*Z
            abc(2,i) = pX(2)*X + pY(2)*Y + pZ(2)*Z
            abc(3,i) = pX(3)*X + pY(3)*Y + pZ(3)*Z
         End Do

      Else If (  Do_option .eq. 2 ) Then

         xyz=0.0_wp
         Do i=1,3
            X=abc(1,i)*a
            Y=abc(2,i)*b
            Z=abc(3,i)*c

            pX(1)= 1.0_wp
            pY(1)= cgm
            pZ(1)= cbt

            pX(2)= 0.0_wp
            pY(2)= sgm
            pZ(2)= (cal-cbt*cgm)/sgm

            pX(3)= 0.0_wp
            pY(3)= 0.0_wp
            pZ(3)= v/sgm

            xyz(1,i) = X+Y*cgm+Z*cbt
            xyz(2,i) =   Y*sgm+Z*((cal-cbt*cgm)/sgm)
            xyz(3,i) =         Z*v/sgm

            xyz(1,i) = pX(1)*X + pY(1)*Y + pZ(1)*Z
            xyz(2,i) = pX(2)*X + pY(2)*Y + pZ(2)*Z
            xyz(3,i) = pX(3)*X + pY(3)*Y + pZ(3)*Z
         End Do
      Else
         Write(6,'(A)') 'the Do_option is not specified. '
         Write(6,'(A)') 'the program continues without ABCC option'
         iReturn=1
         Go To 190
      End If

 190  continue
      Return
      End
