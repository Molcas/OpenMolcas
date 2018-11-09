      Real*8 function distance(N,C1,C2)
      Implicit None
      Integer, parameter        :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)       :: N
      Real(kind=wp), intent(in) :: C1(N), C2(N)
      ! local variables
      Integer       :: i
      Real(kind=wp) :: X, R
      distance=0.0_wp
      X=0.0_wp
      Do i=1,N
        R=0.0_wp
        R=C1(i)-C2(i)
        X=X+R*R
      End Do
      distance=sqrt(X)
      Return
      End
