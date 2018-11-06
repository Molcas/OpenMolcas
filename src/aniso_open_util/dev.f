      Real*8 function dev(N, Fcal, Fexp)
c this function Returns the standard deviation between experimental
c data points and the computed ones
c     N --- number of data points ( Integer, input)
c  Fcal --- calculated array of size (N), Real(kind=wp) ::, input
c  Fexp --- experimental array of size (N), Real(kind=wp) ::, input
c   dev --- standard deviation, Real(kind=wp) ::, output;
      Implicit None
      Integer, parameter        :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)       :: N
      Real(kind=wp), intent(in) :: Fcal(N), Fexp(N)
      Integer                   :: i
      Real(kind=wp)             :: diff, X
      dev=0.0_wp
      X=0.0_wp
      Do i=1,N
        diff=0.0_wp
        diff=Fcal(i)-Fexp(i)
           X=X+diff*diff/dble(N)
      End Do
      dev=sqrt(X)
      Return
      End
