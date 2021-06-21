module interpol_mod

use constants_mod

integer, parameter :: RP = selected_real_kind(15)
!real(dp), parameter :: PI = 4.0_rp * atan(1.0_rp)
!real(dp), parameter :: TOL = 4.0e-16

CONTAINS

!pure     SUBROUTINE DiffMatrix(x,N,D)
!        IMPLICIT NONE
!        INTEGER, INTENT(IN)                             :: N
!        real(dp)(KIND=RP), DIMENSION(0:N), INTENT(IN)       :: x
!        real(dp)(KIND=RP), DIMENSION(0:N,0:N), INTENT(OUT)  :: D
!!
!        ! lokale Variablen
!        INTEGER                                         :: i,j
!        real(dp)(KIND=RP), DIMENSION(0:N)                   :: w
!!
!        ! Baryzentrische Gewichte bestimmen
!        CALL BarycentricWeights(x,N,w)
!        ! DiffMatrix gemaess Aufg.1e) befuellen
!        DO i=0,N
!            DO j=0,N
!                IF (i.NE.j) THEN
!                    D(i,j) = w(j) / ( w(i) * (x(i)-x(j)) )
!                    ! auf der Diagonalen alle anderen subtrahieren
!                    D(i,i) = D(i,i) - D(i,j)
!                END IF
!            END DO
!        END DO
!!
!        RETURN
!!
!    END SUBROUTINE DiffMatrix

pure SUBROUTINE BarycentricWeights(N_in,xGP,wBary)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N_in               !< polynomial degree
real(dp),INTENT(IN)    :: xGP(0:N_in)        !< Gauss point positions for the reference interval [-1,1]
real(dp),INTENT(OUT)   :: wBary(0:N_in)      !< barycentric weights
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iGP,jGP
!==================================================================================================================================
wBary(:)=1.
DO iGP=1,N_in
  DO jGP=0,iGP-1
    wBary(jGP)=wBary(jGP)*(xGP(jGP)-xGP(iGP))
    wBary(iGP)=wBary(iGP)*(xGP(iGP)-xGP(jGP))
  END DO ! jGP
END DO ! iGP
wBary(:)=1./wBary(:)
END SUBROUTINE BarycentricWeights

pure SUBROUTINE DiffMatrix(N_in,xGP,D)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N_in              !< polynomial degree
real(dp),INTENT(IN)    :: xGP(0:N_in)       !< Gauss point positions for the reference interval [-1,1]
real(dp),INTENT(OUT)   :: D(0:N_in,0:N_in)  !< differentiation Matrix
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iGP,iLagrange
real(dp)               :: wBary(0:N_in)
!==================================================================================================================================
CALL BarycentricWeights(N_in,xGP,wBary)
D(:,:)=0.
DO iLagrange=0,N_in
  DO iGP=0,N_in
    IF(iLagrange.NE.iGP)THEN
      D(iGP,iLagrange)=wBary(iLagrange)/(wBary(iGP)*(xGP(iGP)-xGP(iLagrange)))
      D(iGP,iGP)=D(iGP,iGP)-D(iGP,iLagrange)
    END IF ! (iLagrange.NE.iGP)
  END DO ! iGP
END DO ! iLagrange
END SUBROUTINE

pure function LagrangePolynomialDerivative(N, xs, j, x) result(lp)

    integer, intent(in)     :: N
    real(dp), intent(in)    :: xs(N)
    integer, intent(in)     :: j
    real(dp), intent(in)    :: x

    real(dp)                :: lp
    real(dp)                :: acc
    integer                 :: i,m

    lp = 0.0_dp

    do i = 1,N
        if (i == j) cycle
        acc = 1.0_dp
        do m = 1,N
            if (m == i .or. m == j) cycle
            acc = acc * (x-xs(m))/(xs(j)-xs(m))
        end do
        lp = lp + acc/(xs(j)-xs(i))
    end do

end function

pure function LagrangePolynomialDerivativeMatrix(N,M, xs, ys) result(mat)

    integer, intent(in)     :: N,M
    real(dp), intent(in)    :: xs(N),ys(M)

    real(dp)                :: mat(M,N)
    integer                 :: i,j

    do i = 1,M; do j = 1,N;
        mat(i,j) = LagrangePolynomialDerivative(N,xs,j,ys(i))
    end do; end do

end function

pure function LagrangePolynomial(n, nodes, j, x) result(lp)
    
    integer, intent(in) :: n
    real(dp), intent(in) :: nodes(0:)
    integer, intent(in) :: j
    real(dp), intent(in) :: x

    real(dp) :: lp
    integer :: i

    lp = 1.0_dp
    do i = 0,n-1
        if (i == j) cycle
        lp = lp * (x - nodes(i))/(nodes(j) - nodes(i))
    end do
     
end function

pure     SUBROUTINE LagrangeInterpolatingPolynomials(y,x,w,N,l)

!(Algorithm 34)

        IMPLICIT NONE
        INTEGER, INTENT(IN)                         :: N
        real(dp), INTENT(IN)                   :: y
        real(dp), DIMENSION(0:N), INTENT(IN)   :: x,w
        real(dp), DIMENSION(0:N), INTENT(OUT)  :: l
    !
        ! lokale Variablen
        LOGICAL                                     :: stuetzstelle
        INTEGER                                     :: j
        real(dp)                               :: s,t
        real(dp), PARAMETER                    :: tol = 10.0E-16_dp
    !
        ! Falls y eine Stuetzstelle ist:
        stuetzstelle = .FALSE.
        DO j=0,N
            l(j) = 0.0_RP
            IF (ABS(y-x(j)).LT.tol) THEN
                l(j) = 1.0_RP
                stuetzstelle = .TRUE.
            END IF
        END DO
    !
!        IF (stuetzstelle) THEN
!            EXIT
!        END IF

        ! Falls y keine Stuetzstelle ist:
        IF (.NOT. stuetzstelle) THEN
            s = 0.0_RP

            DO j=0,N
                t = w(j) / (y-x(j))
                l(j) = t
                s = s + t
            END DO
            DO j=0,N
                l(j) = l(j)/s
            END DO
        END IF
!
    RETURN
!
    END SUBROUTINE LagrangeInterpolatingPolynomials
!
! ------------------------------------------------------------------------------------------------------------- !
!
pure     SUBROUTINE Evaluation(y,x,f,w,N,f_y)
        IMPLICIT NONE
        INTEGER, INTENT(IN)                         :: N
        real(dp), INTENT(IN)                   :: y
        real(dp), DIMENSION(0:N), INTENT(IN)   :: x,f,w
        real(dp), INTENT(OUT)                   :: f_y
    !
        ! lokale Variablen
        LOGICAL                                     :: stuetzstelle
        INTEGER                                     :: j
        real(dp)                               :: numerator,denominator,hilf
        real(dp), PARAMETER                    :: tol = 1.0E-16_dp
    !
        ! Falls y eine Stuetzstelle ist:
        stuetzstelle = .FALSE.
        DO j=0,N
            IF (ABS(y-x(j)).LT.tol) THEN
                f_y = f(j)
                stuetzstelle = .TRUE.
                EXIT
            END IF
        END DO
    !
        ! Falls y keine Stuetzstelle ist:
        numerator = 0
        denominator = 0
        IF (.NOT. stuetzstelle) THEN
            DO j=0,N
                hilf = w(j) / (y-x(j))
                numerator = numerator + hilf * f(j)
                denominator = denominator + hilf
            END DO
            f_y = numerator / denominator
        END IF
    !
        RETURN
    !
    END SUBROUTINE Evaluation

pure     SUBROUTINE LegendrePolynomialAndDerivativeAlbert(N,x,L_N,dL_N)
        IMPLICIT NONE
        INTEGER, INTENT(IN)         :: N
        real(dp), INTENT(IN)   :: x
        real(dp), INTENT(OUT)  :: L_N,dL_N
    !
        ! lokale Variablen
        INTEGER                     :: k
        real(dp)               :: L_N_1,L_N_2,dL_N_1,dL_N_2
    !
        ! Fuer den einfachsten Fall explizit vorgeben
        IF (N.EQ.0) THEN
            L_N = 0.0_RP
            dL_N = 0.0_RP
    !
        !  Ebenso fuer den zweit-einfachsten Fall
        ELSE IF (N.EQ.1) THEN
            L_N = x
            dL_N = 1.0_RP
    !
        ! Ansonsten auf die 3-Term-Rekursion zurueckgreifen
        ELSE
            ! Dazu wieder die beiden ersten Polynome vorgeben
            L_N_2 = 1.0_RP
            L_N_1 = x
            dL_N_2 = 0.0_RP
            dL_N_1 = 1.0_RP
    !
        DO k=2,N
            ! Rekursionsformel fuer die Legendre-Polynome verwenden
            L_N = ( (2.0_RP*real(k,dP)-1.0_RP) / real(k,dp) ) * x * L_N_1 - (real(k,dP)-1.0_RP) / real(k,dP) * L_N_2
            ! Rekursionsformel fuer die Ableitungen verwenden
            dL_N = dL_N_2 + ( 2.0_RP*real(k,dP) - 1.0_RP ) * L_N_1
            ! Werte fuer den naechsten Schritt updaten
            L_N_2 = L_N_1
            L_N_1 = L_N
            dL_N_2 = dL_N_1
            dL_N_1 = dL_N
        END DO
    !
        END IF
    !
        RETURN
    END SUBROUTINE

pure ELEMENTAL SUBROUTINE LegendrePolynomialAndDerivative(N_in,x,L,Lder)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N_in   !< (IN)  polynomial degree, (N+1) CLpoints
real(dp),INTENT(IN)    :: x      !< (IN)  coordinate value in the interval [-1,1]
real(dp),INTENT(OUT)   :: L      !< (OUT) Legedre polynomial evaluated at \f$ \xi: L_N(\xi), \partial/\partial\xi L_N(\xi) \f$
real(dp),INTENT(OUT)   :: Lder   !< (OUT) Legedre polynomial deriv. evaluated at \f$ \xi: L_N(\xi), \partial/\partial\xi L_N(\xi) \f$
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iLegendre
real(dp)    :: L_Nm1,L_Nm2 ! L_{N_in-2},L_{N_in-1}
real(dp)    :: Lder_Nm1,Lder_Nm2 ! Lder_{N_in-2},Lder_{N_in-1}
!==================================================================================================================================
IF(N_in .EQ. 0)THEN
  L=1.0_dp
  Lder=0.0_dp
ELSEIF(N_in .EQ. 1) THEN
  L=x
  Lder=1.0_dp
ELSE ! N_in > 1
  L_Nm2=1.0_dp
  L_Nm1=x
  Lder_Nm2=0.0_dp
  Lder_Nm1=1.0_dp
  DO iLegendre=2,N_in
    L=(real(2*iLegendre-1,dp)*x*L_Nm1 - real(iLegendre-1,dp)*L_Nm2)/real(iLegendre,dp)
    Lder=Lder_Nm2 + real(2*iLegendre-1,dp)*L_Nm1
    L_Nm2=L_Nm1
    L_Nm1=L
    Lder_Nm2=Lder_Nm1
    Lder_Nm1=Lder
  END DO !iLegendre=2,N_in
END IF ! N_in
!normalize
L=L*SQRT(real(N_in,dp)+0.5_dp)
Lder=Lder*SQRT(real(N_in,dp)+0.5_dp)
END SUBROUTINE LegendrePolynomialAndDerivative

pure SUBROUTINE EquidistantNodesAndWeights(N_in,xGP,wGP)

    IMPLICIT NONE

    INTEGER,INTENT(IN)              :: N_in              !< polynomial degree, (N_in+1) Gausspoints
    real(dp),INTENT(OUT)            :: xGP(0:N_in)       !< Gauss point positions for the reference interval [-1,1]
    real(dp),INTENT(OUT),OPTIONAL   :: wGP(0:N_in)       !< Gauss point weights

    integer :: i

    do i = 1,N_IN+1
        xGP(i-1) = -1.0_dp + 2.0_dp*(real(i,dp)-0.5_dp)/real(N_IN+1,dp)
    end do

    wGP = 2.0_dp/real(N_IN+1,dp)

end subroutine

pure subroutine equidistandinterfacenodes_3times(N,xs,ws)

    implicit none

    integer,intent(in)              :: N
    real(dp),intent(out)            :: xs(N)
    real(dp),intent(out),optional   :: ws(N)

    integer :: i

    do i = 2,N,3
        xs(i-1) = -1.0_dp + 2.0_dp*(real(i,dp)-0.5_dp)/real(N,dp) - 3.0_dp/real(N,dp)
        xs(i+0) = -1.0_dp + 2.0_dp*(real(i,dp)-0.5_dp)/real(N,dp)
        xs(i+1) = -1.0_dp + 2.0_dp*(real(i,dp)-0.5_dp)/real(N,dp) + 3.0_dp/real(N,dp)
    end do

    !wgp = 2.0/real(n_in+1,dp)
    ws = 0.0_dp

end subroutine

pure subroutine equidistandinterfacenodes(N,xs)

    integer,intent(in)              :: N
    real(dp),intent(out)            :: xs(N)

    integer :: i

    do i = 1,N
        xs(2*i-1) = -1.0_dp + 2.0_dp*(real(i,dp)-0.5_dp)/real(N,dp) - 1.0_dp/real(N,dp)
        xs(2*i+0) = -1.0_dp + 2.0_dp*(real(i,dp)-0.5_dp)/real(N,dp)
    end do

    xs(2*N+1) = 1.0_dp

end subroutine

pure subroutine InterfaceNodes(N,xL,xC,xR)

    implicit none

    integer,intent(in)      :: N
    real(dp),intent(out)    :: xL(N),xC(N),xR(N)

    integer :: i

    do i = 1,N
        xL(i) = -1.0_dp + 2.0_dp*(real(i,dp)-0.5_dp)/real(N,dp) - 1.0_dp/real(N,dp)
        xC(i) = -1.0_dp + 2.0_dp*(real(i,dp)-0.5_dp)/real(N,dp)
        xR(i) = -1.0_dp + 2.0_dp*(real(i,dp)-0.5_dp)/real(N,dp) + 1.0_dp/real(N,dp)
    end do

end subroutine


pure SUBROUTINE NodesAndWeights(N_in,xGP,wGP)

    IMPLICIT NONE

    INTEGER,INTENT(IN)        :: N_in              !< polynomial degree, (N_in+1) Gausspoints
    real(dp),INTENT(OUT)          :: xGP(0:N_in)       !< Gauss point positions for the reference interval [-1,1]
    real(dp),INTENT(OUT),OPTIONAL :: wGP(0:N_in)       !< Gauss point weights

    call LegendreGaussNodesAndWeights(N_in, xGP, wGP)

end subroutine


pure SUBROUTINE LegendreGaussNodesAndWeights(N_in,xGP,wGP)

IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)        :: N_in              !< polynomial degree, (N_in+1) Gausspoints
real(dp),INTENT(OUT)          :: xGP(0:N_in)       !< Gauss point positions for the reference interval [-1,1]
real(dp),INTENT(OUT),OPTIONAL :: wGP(0:N_in)       !< Gauss point weights
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER, parameter          :: nIter = 10        ! max. number of newton iterations
real(dp), parameter             :: Tol   = 1.E-15_dp    ! tolerance for Newton iteration: TODO: use variable tolerance here!
INTEGER                   :: iGP,iter
real(dp)                      :: L_Np1,Lder_Np1    ! L_{N_in+1},Lder_{N_in+1}
real(dp)                      :: dx                ! Newton step
real(dp)                      :: cheb_tmp          ! temporary variable for evaluation of chebychev node positions
!==================================================================================================================================
IF(N_in .EQ. 0) THEN
  xGP=0.0_dp
  IF(PRESENT(wGP))wGP=2.0_dp
  RETURN
ELSEIF(N_in.EQ.1)THEN
  xGP(0)=-sqrt(1.0_dp/3.0_dp)
  xGP(N_in)=-xGP(0)
  IF(PRESENT(wGP))wGP=1.0_dp
  RETURN
ELSE ! N_in>1
  cheb_tmp=2.0_dp*atan(1.0_dp)/real(N_in+1,dp) ! pi/(2N+2)
  DO iGP=0,(N_in+1)/2-1 !since points are symmetric, only left side is computed
    xGP(iGP)=-cos(cheb_tmp*real(2*iGP+1,dp)) !initial guess
    ! Newton iteration
    DO iter=0,nIter
      CALL LegendrePolynomialAndDerivative(N_in+1,xGP(iGP),L_Np1,Lder_Np1)
      dx=-L_Np1/Lder_Np1
      xGP(iGP)=xGP(iGP)+dx
      IF(abs(dx).LT.Tol*abs(xGP(iGP))) EXIT
    END DO ! iter
    IF(iter.GT.nIter) THEN
      !WRITE(*,*) 'maximum iteration steps >10 in Newton iteration for Legendre Gausspoint'
      xGP(iGP)=-cos(cheb_tmp*real(2*iGP+1,dp)) !initial guess
      ! Newton iteration
      DO iter=0,nIter
        !WRITE(*,*) iter,xGP(iGP)    !DEBUG
        CALL LegendrePolynomialAndDerivative(N_in+1,xGP(iGP),L_Np1,Lder_Np1)
        dx=-L_Np1/Lder_Np1
        xGP(iGP)=xGP(iGP)+dx
        IF(abs(dx).LT.Tol*abs(xGP(iGP))) EXIT
      END DO !iter
      !CALL abort(__STAMP__,&
      !           'ERROR: Legendre Gauss nodes could not be computed up to desired precision. Code stopped!')
    END IF ! (iter.GT.nIter)
    CALL LegendrePolynomialAndDerivative(N_in+1,xGP(iGP),L_Np1,Lder_Np1)
    xGP(N_in-iGP)=-xGP(iGP)
    IF(PRESENT(wGP))THEN
      !wGP(iGP)=2./((1.-xGP(iGP)*xGP(iGP))*Lder_Np1*Lder_Np1) !if Legendre not normalized
      wGP(iGP)=(2.0_dp*N_in+3)/((1.0_dp-xGP(iGP)*xGP(iGP))*Lder_Np1*Lder_Np1)
      wGP(N_in-iGP)=wGP(iGP)
    END IF
  END DO !iGP
END IF ! N_in
IF(mod(N_in,2) .EQ. 0) THEN
  xGP(N_in/2)=0.
  CALL LegendrePolynomialAndDerivative(N_in+1,xGP(N_in/2),L_Np1,Lder_Np1)
  !IF(PRESENT(wGP))wGP(N_in/2)=2./(Lder_Np1*Lder_Np1) !if Legendre not normalized
  IF(PRESENT(wGP))wGP(N_in/2)=(2.0_dp*N_in+3.0_dp)/(Lder_Np1*Lder_Np1)
END IF ! (mod(N_in,2) .EQ. 0)
END SUBROUTINE LegendreGaussNodesAndWeights


pure     SUBROUTINE qAndLEvaluation(N,x,q,dq,L_N)
        IMPLICIT NONE
        INTEGER, INTENT(IN)         :: N
        real(dp), INTENT(IN)   :: x
        real(dp), INTENT(OUT)  :: L_N,q,dq
    !
        ! lokale Variablen
        INTEGER                     :: k
        real(dp)               :: L_N_1,L_N_2,L_k,dL_N,dL_N_1,dL_N_2,dL_k
    !
        ! NUR FUER N>=2 !
    !
        ! Initialisieren
        L_N_2 = 1.0_RP
        L_N_1 = x
        dL_N_2 = 0.0_RP
        dL_N_1 = 1.0_RP
    !
        DO k=2,N
            ! Rekursionsformel verwenden
            L_N = (2.0_RP*real(k,dP)-1.0_RP)/real(k,dP) * x * L_N_1 - (real(k,dP)-1.0_RP)/real(k,dP) * L_N_2
            ! Rekursionsformel fuer die Ableitungen verwenden
            dL_N = dL_N_2 + (2.0_RP*real(k,dP)-1.0_RP) * L_N_1
            ! Werte fuer den naechsten Schritt updaten
            L_N_2 = L_N_1
            L_N_1 = L_N
            dL_N_2 = dL_N_1
            dL_N_1 = dL_N
        END DO
    !
        ! Einen weiteren Schritt gehen
        k = N+1
        ! Rekursionsformel verwenden
        L_k = (2.0_RP*real(k,dP)-1.0_RP)/real(k,dP) * x * L_N - (real(k,dP)-1.0_RP)/real(k,dP) * L_N_2
        ! Rekursionsformel fuer die Ableitungen verwenden
        dL_k = dL_N_2 + (2.0_RP*real(k,dP)-1.0_RP) * L_N_1
        ! Benoetigtes Polynom und Ableitung bestimmen
        q = L_k - L_N_2
        dq = dL_k - dL_N_2
    !
        RETURN
    END SUBROUTINE qAndLEvaluation
!
! ------------------------------------------------------------------------------------------------------------- !
!
pure     SUBROUTINE LegendreGaussLobattoNodesAndWeights(N,x,w)
        IMPLICIT NONE
        INTEGER, INTENT(IN)                         :: N
        real(dp), DIMENSION(0:N), INTENT(OUT)  :: x,w
    !
        ! lokale Variablen
        INTEGER, PARAMETER            :: n_it = 4
        INTEGER                       :: j,k
        real(dp), PARAMETER      :: tol = 4.0E-16_dp
        real(dp)                 :: q,dq,L_N,Delta

        L_N = 0.0_dp
    !
        ! Fuer den einfachsten Fall Knoten und Gewicht explizit angeben
        IF (N.EQ.1) THEN
            x(0) = -1.0_RP
            w(0) = 1.0_RP
            x(1) = 1.0_RP
            w(1) = w(0)
    !
        ELSE
            ! Randwerte explizit vorgeben
            x(0) = -1.0_RP
            w(0) = 2.0_RP / ( real(N,dP)*(real(N,dP)+1) )
            x(N) = 1.0_RP
            w(N) = w(0)
            ! Ansonsten: Symmetrie ausnutzen (nur fuer die erste Haelfte berechnen)
            DO j=1,((N+1)/2 - 1)
                ! Anfangsschaetzung mittels asymptotischer Relation nach Parter:
                ! Setzte negatives Vorzeichen, um die Stuetzstellen aus dem Intervall (-1,0) zu erhalten.
                x(j) = -cos( ((real(j,dP)+0.25_RP)*pi)/real(N,dP) - (3.0_RP/(8_RP*real(N,dP)*pi)) * (1.0_dp/(real(j,dP)+0.25_RP)) )
    !
                ! Diese Schaetzung wird nun mithilfe des Newton-Verfahrens praezisiert
                DO k=0,n_it
                    ! Dazu wird das (N+1)-ste Polynom sowie dessen Ableitung benoetigt
                    ! (ausgewertet an der entspr. Stuetzstelle)
                    CALL qAndLEvaluation(N,x(j),q,dq,L_N)
                    ! Newton-Korrektur berechnen und anwenden
                    Delta = -q/dq
                    x(j) = x(j) + Delta
                    ! Falls vor Abarbeiten aller festgelegten Iterationen des Newton-Verfahrens der berechnete Wert
                    ! schon nahe genug an der Nullstelle liegt, wird abgebrochen.
                    IF ( ABS(Delta).LE.(tol * ABS(x(j))) ) THEN
                        EXIT
                    END IF
    !
                END DO
    !
                CALL qAndLEvaluation(N,x(j),q,dq,L_N)
                ! Nutze Symmetrie aus (Multiplikation mit -1)
                x(N-j) = -x(j)
                ! Gewichte berechnen
                w(j) = 2.0_RP / ( real(N,dP)*(real(N,dP)+1.0_RP) * (L_N*L_N) )
                ! Auch die Gewichte sind symmetrisch
                w(N-j) = w(j)
    !
            END DO
    !
        END IF
    !
        ! Falls eine ungerade Anzahl gefordert wird, noch die "mittleren" Werte angeben
        IF (MOD(N,2).EQ.0) THEN
            CALL qAndLEvaluation(N,0.0_RP,q,dq,L_N)
            ! Berechnung der Stuetzstellen und Gewichte analog
            x(N/2) = 0.0_RP
            w(N/2) = 2.0_RP / ( real(N,dP) * (real(N,dP)+1.0_RP) * (L_N*L_N) )
        END IF
    !
        RETURN
    END SUBROUTINE LegendreGaussLobattoNodesAndWeights

pure     subroutine mkVanderMondeMatrix(n, nodes, Vdm, sVdm)

        use linalg_mod, only: invert

        integer, intent(in)  :: n
        real(dp), intent(in)  :: nodes(:)
        real(dp), intent(out) ::   Vdm(:,:)
        real(dp), intent(out) ::  sVdm(:,:)
        real(dp) :: dummy
        integer :: i,j

        do i = 1,n; do j = 1,n
            call LegendrePolynomialAnDderivative(j-1, nodes(i), Vdm(i,j), dummy)
        end do; end do

        sVdm = invert(Vdm,n)
   
    end subroutine

END MODULE
