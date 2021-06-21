module gausslobatto_mod

integer, parameter :: RP = selected_real_kind(15)
real(rp), parameter :: PI = 4.0_rp * atan(1.0_rp)
real(rp), parameter :: TOL = 4.0e-16_rp

!! CONTAINS

!! SUBROUTINE LegendrePolynomialAndDerivativeAlbert(N,x,L_N,dL_N)
!!     IMPLICIT NONE
!!     INTEGER, INTENT(IN)         :: N
!!     REAL(KIND=RP), INTENT(IN)   :: x
!!     REAL(KIND=RP), INTENT(OUT)  :: L_N,dL_N
!! !
!!     ! lokale Variablen
!!     INTEGER                     :: k
!!     REAL(KIND=RP)               :: L_N_1,L_N_2,dL_N_1,dL_N_2
!! !
!!     ! Fuer den einfachsten Fall explizit vorgeben
!!     IF (N.EQ.0) THEN
!!         L_N = 0.0_RP
!!         dL_N = 0.0_RP
!! !
!!     !  Ebenso fuer den zweit-einfachsten Fall
!!     ELSE IF (N.EQ.1) THEN
!!         L_N = x
!!         dL_N = 1.0_RP
!! !
!!     ! Ansonsten auf die 3-Term-Rekursion zurueckgreifen
!!     ELSE
!!         ! Dazu wieder die beiden ersten Polynome vorgeben
!!         L_N_2 = 1.0_RP
!!         L_N_1 = x
!!         dL_N_2 = 0.0_RP
!!         dL_N_1 = 1.0_RP
!! !
!!     DO k=2,N
!!         ! Rekursionsformel fuer die Legendre-Polynome verwenden
!!         L_N = ( (2.0_RP*REAL(k,RP)-1.0_RP) / REAL(k,RP) ) * x * L_N_1 - (REAL(k,RP)-1.0_RP) / REAL(k,RP) * L_N_2
!!         ! Rekursionsformel fuer die Ableitungen verwenden
!!         dL_N = dL_N_2 + ( 2.0_RP*REAL(k,RP) - 1.0_RP ) * L_N_1
!!         ! Werte fuer den naechsten Schritt updaten
!!         L_N_2 = L_N_1
!!         L_N_1 = L_N
!!         dL_N_2 = dL_N_1
!!         dL_N_1 = dL_N
!!     END DO
!! !
!!     END IF
!! !
!!     RETURN
!! END SUBROUTINE
!! 
!! ELEMENTAL SUBROUTINE LegendrePolynomialAndDerivative(N_in,x,L,Lder)
!! IMPLICIT NONE
!! !----------------------------------------------------------------------------------------------------------------------------------
!! ! INPUT/OUTPUT VARIABLES
!! INTEGER,INTENT(IN) :: N_in   !< (IN)  polynomial degree, (N+1) CLpoints
!! REAL,INTENT(IN)    :: x      !< (IN)  coordinate value in the interval [-1,1]
!! REAL,INTENT(OUT)   :: L      !< (OUT) Legedre polynomial evaluated at \f$ \xi: L_N(\xi), \partial/\partial\xi L_N(\xi) \f$
!! REAL,INTENT(OUT)   :: Lder   !< (OUT) Legedre polynomial deriv. evaluated at \f$ \xi: L_N(\xi), \partial/\partial\xi L_N(\xi) \f$
!! !----------------------------------------------------------------------------------------------------------------------------------
!! ! LOCAL VARIABLES
!! INTEGER :: iLegendre
!! REAL    :: L_Nm1,L_Nm2 ! L_{N_in-2},L_{N_in-1}
!! REAL    :: Lder_Nm1,Lder_Nm2 ! Lder_{N_in-2},Lder_{N_in-1}
!! !==================================================================================================================================
!! IF(N_in .EQ. 0)THEN
!!   L=1.
!!   Lder=0.
!! ELSEIF(N_in .EQ. 1) THEN
!!   L=x
!!   Lder=1.
!! ELSE ! N_in > 1
!!   L_Nm2=1.
!!   L_Nm1=x
!!   Lder_Nm2=0.
!!   Lder_Nm1=1.
!!   DO iLegendre=2,N_in
!!     L=(REAL(2*iLegendre-1)*x*L_Nm1 - REAL(iLegendre-1)*L_Nm2)/REAL(iLegendre)
!!     Lder=Lder_Nm2 + REAL(2*iLegendre-1)*L_Nm1
!!     L_Nm2=L_Nm1
!!     L_Nm1=L
!!     Lder_Nm2=Lder_Nm1
!!     Lder_Nm1=Lder
!!   END DO !iLegendre=2,N_in
!! END IF ! N_in
!! !normalize
!! L=L*SQRT(REAL(N_in)+0.5)
!! Lder=Lder*SQRT(REAL(N_in)+0.5)
!! END SUBROUTINE LegendrePolynomialAndDerivative
!! 
!! 
!! 
!! !
!! ! ------------------------------------------------------------------------------------------------------------- !
!! !
!!     SUBROUTINE LegendreGaussNodesAndWeights(N,x,w)
!!         IMPLICIT NONE
!!         INTEGER, INTENT(IN)                         :: N
!!         REAL(KIND=RP), DIMENSION(0:N), INTENT(OUT)  :: x,w
!!     !
!!         ! lokale Variablen
!!         INTEGER, PARAMETER                          :: n_it = 4
!!         INTEGER                                     :: j,k
!!         REAL(KIND=RP), PARAMETER                    :: tol = 4.0E-16
!!         REAL(KIND=RP)                               :: L_N,dL_N,Delta
!! 
!!         L_N = 0.0
!!         dL_N = 0.0
!!     !
!!         ! Fuer den einfachsten Fall Knoten und Gewicht explizit angeben
!!         IF (N.EQ.0) THEN
!!             x(0) = 0.0_RP
!!             w(0) = 2.0_RP
!!     !
!!         ! Ebenso fuer den zweit-einfachsten Fall
!!         ELSE IF (N.EQ.1) THEN
!!             x(0) = -SQRT( 1.0_RP / 3.0_RP )
!!             w(0) = 1.0_RP
!!             x(1) = -x(0)
!!             w(1) = w(0)
!!     !
!!         ELSE
!!             ! Ansonsten: Symmetrie ausnutzen (nur fuer die erste Haelfte berechnen)
!!             DO j=0,((N+1)/2 - 1)
!!                 ! Anfangsschaetzung mittels Nullstellen der Tschebyschow-Polynome:
!!                 ! Setzte negatives Vorzeichen, um die Stuetzstellen aus dem Intervall (-1,0) zu erhalten.
!!                 x(j) = -cos( (2.0_RP*REAL(j,RP)+1.0_RP) / (2.0_RP*REAL(N,RP)+2.0_RP) * pi )
!!     !
!!                 ! Diese Schaetzung wird nun mithilfe des Newton-Verfahrens praezisiert
!!                 DO k=0,n_it
!!                     ! Dazu wird das (N+1)-ste Legendre-Polynom sowie dessen Ableitung benoetigt
!!                     ! (ausgewertet an der entspr. Stuetzstelle)
!!                     CALL LegendrePolynomialAndDerivative(N+1,x(j),L_N,dL_N)
!!                     ! Newton-Korrektur berechnen und anwenden
!!                     Delta = -L_N/dL_N
!!                     x(j) = x(j) + Delta
!!                     ! Falls vor Abarbeiten aller festgelegten Iterationen des Newton-Verfahrens der berechnete Wert
!!                     ! schon nahe genug an der Nullstelle liegt, wird abgebrochen.
!!                     IF ( ABS(Delta).LE.(tol * ABS(x(j))) ) THEN
!!                         EXIT
!!                     END IF
!!     !
!!                 END DO
!!     !
!!                 CALL LegendrePolynomialAndDerivative(N+1,x(j),L_N,dL_N)
!!                 ! Nutze Symmetrie aus (Multiplikation mit -1)
!!                 x(N-j) = -x(j)
!!                 ! Gewichte berechnen
!!                 w(j) = 2.0_RP / ( (1.0_RP-x(j)*x(j)) * dL_N*dL_N )
!!                 ! Auch die Gewichte sind symmetrisch
!!                 w(N-j) = w(j)
!!     !
!!             END DO
!!     !
!!         END IF
!!     !
!!         ! Falls eine ungerade Anzahl gefordert wird, noch die "mittleren" Werte angeben
!!         IF (MOD(N,2).EQ.0) THEN
!!             CALL LegendrePolynomialAndDerivative(N+1,0.0_RP,L_N,dL_N)
!!             ! Berechnung der Stuetzstellen und Gewichte analog
!!             x(N/2) = 0.0_RP
!!             w(N/2) = 2.0_RP / (dL_N*dL_N)
!!         END IF
!!     !
!!         RETURN
!!     END SUBROUTINE LegendreGaussNodesAndWeights
!! !
!! ! ------------------------------------------------------------------------------------------------------------- !
!! !
!!     SUBROUTINE qAndLEvaluation(N,x,q,dq,L_N)
!!         IMPLICIT NONE
!!         INTEGER, INTENT(IN)         :: N
!!         REAL(KIND=RP), INTENT(IN)   :: x
!!         REAL(KIND=RP), INTENT(OUT)  :: L_N,q,dq
!!     !
!!         ! lokale Variablen
!!         INTEGER                     :: k
!!         REAL(KIND=RP)               :: L_N_1,L_N_2,L_k,dL_N,dL_N_1,dL_N_2,dL_k
!!     !
!!         ! NUR FUER N>=2 !
!!     !
!!         ! Initialisieren
!!         L_N_2 = 1.0_RP
!!         L_N_1 = x
!!         dL_N_2 = 0.0_RP
!!         dL_N_1 = 1.0_RP
!!     !
!!         DO k=2,N
!!             ! Rekursionsformel verwenden
!!             L_N = (2.0_RP*REAL(k,RP)-1.0_RP)/REAL(k,RP) * x * L_N_1 - (REAL(k,RP)-1.0_RP)/REAL(k,RP) * L_N_2
!!             ! Rekursionsformel fuer die Ableitungen verwenden
!!             dL_N = dL_N_2 + (2.0_RP*REAL(k,RP)-1.0_RP) * L_N_1
!!             ! Werte fuer den naechsten Schritt updaten
!!             L_N_2 = L_N_1
!!             L_N_1 = L_N
!!             dL_N_2 = dL_N_1
!!             dL_N_1 = dL_N
!!         END DO
!!     !
!!         ! Einen weiteren Schritt gehen
!!         k = N+1
!!         ! Rekursionsformel verwenden
!!         L_k = (2.0_RP*REAL(k,RP)-1.0_RP)/REAL(k,RP) * x * L_N - (REAL(k,RP)-1.0_RP)/REAL(k,RP) * L_N_2
!!         ! Rekursionsformel fuer die Ableitungen verwenden
!!         dL_k = dL_N_2 + (2.0_RP*REAL(k,RP)-1.0_RP) * L_N_1
!!         ! Benoetigtes Polynom und Ableitung bestimmen
!!         q = L_k - L_N_2
!!         dq = dL_k - dL_N_2
!!     !
!!         RETURN
!!     END SUBROUTINE qAndLEvaluation
!! !
!! ! ------------------------------------------------------------------------------------------------------------- !
!! !
!!     SUBROUTINE LegendreGaussLobattoNodesAndWeights(N,x,w)
!!         IMPLICIT NONE
!!         INTEGER, INTENT(IN)                         :: N
!!         REAL(KIND=RP), DIMENSION(0:N), INTENT(OUT)  :: x,w
!!     !
!!         ! lokale Variablen
!!         INTEGER, PARAMETER            :: n_it = 4
!!         INTEGER                       :: j,k
!!         REAL(KIND=RP), PARAMETER      :: tol = 4.0E-16
!!         REAL(KIND=RP)                 :: q,dq,L_N,Delta
!! 
!!         L_N = 0.0
!!     !
!!         ! Fuer den einfachsten Fall Knoten und Gewicht explizit angeben
!!         IF (N.EQ.1) THEN
!!             x(0) = -1.0_RP
!!             w(0) = 1.0_RP
!!             x(1) = 1.0_RP
!!             w(1) = w(0)
!!     !
!!         ELSE
!!             ! Randwerte explizit vorgeben
!!             x(0) = -1.0_RP
!!             w(0) = 2.0_RP / ( REAL(N,RP)*(REAL(N,RP)+1) )
!!             x(N) = 1.0_RP
!!             w(N) = w(0)
!!             ! Ansonsten: Symmetrie ausnutzen (nur fuer die erste Haelfte berechnen)
!!             DO j=1,((N+1)/2 - 1)
!!                 ! Anfangsschaetzung mittels asymptotischer Relation nach Parter:
!!                 ! Setzte negatives Vorzeichen, um die Stuetzstellen aus dem Intervall (-1,0) zu erhalten.
!!                 x(j) = -cos( ((REAL(j,RP)+0.25_RP)*pi)/REAL(N,RP) - (3.0_RP/(8_RP*REAL(N,RP)*pi)) * (1/(REAL(j,RP)+0.25_RP)) )
!!     !
!!                 ! Diese Schaetzung wird nun mithilfe des Newton-Verfahrens praezisiert
!!                 DO k=0,n_it
!!                     ! Dazu wird das (N+1)-ste Polynom sowie dessen Ableitung benoetigt
!!                     ! (ausgewertet an der entspr. Stuetzstelle)
!!                     CALL qAndLEvaluation(N,x(j),q,dq,L_N)
!!                     ! Newton-Korrektur berechnen und anwenden
!!                     Delta = -q/dq
!!                     x(j) = x(j) + Delta
!!                     ! Falls vor Abarbeiten aller festgelegten Iterationen des Newton-Verfahrens der berechnete Wert
!!                     ! schon nahe genug an der Nullstelle liegt, wird abgebrochen.
!!                     IF ( ABS(Delta).LE.(tol * ABS(x(j))) ) THEN
!!                         EXIT
!!                     END IF
!!     !
!!                 END DO
!!     !
!!                 CALL qAndLEvaluation(N,x(j),q,dq,L_N)
!!                 ! Nutze Symmetrie aus (Multiplikation mit -1)
!!                 x(N-j) = -x(j)
!!                 ! Gewichte berechnen
!!                 w(j) = 2.0_RP / ( REAL(N,RP)*(REAL(N,RP)+1.0_RP) * (L_N*L_N) )
!!                 ! Auch die Gewichte sind symmetrisch
!!                 w(N-j) = w(j)
!!     !
!!             END DO
!!     !
!!         END IF
!!     !
!!         ! Falls eine ungerade Anzahl gefordert wird, noch die "mittleren" Werte angeben
!!         IF (MOD(N,2).EQ.0) THEN
!!             CALL qAndLEvaluation(N,0.0_RP,q,dq,L_N)
!!             ! Berechnung der Stuetzstellen und Gewichte analog
!!             x(N/2) = 0.0_RP
!!             w(N/2) = 2.0_RP / ( REAL(N,RP) * (REAL(N,RP)+1.0_RP) * (L_N*L_N) )
!!         END IF
!!     !
!!         RETURN
!!     END SUBROUTINE LegendreGaussLobattoNodesAndWeights
!! 
!! subroutine mkLegendreVanderMondeMatrix(n, nodes, Vdm, sVdm)
!! 
!!     use linalg_mod, only: invert
!! 
!!     integer, intent(in)  :: n
!!     real, intent(in)  :: nodes(:)
!!     real, intent(out) ::   Vdm(:,:)
!!     real, intent(out) ::  sVdm(:,:)
!!     real :: dummy
!!     integer :: i,j
!! 
!!     do i = 1,n; do j = 1,n
!!         call LegendrePolynomialAnDderivative(j-1, nodes(i), Vdm(i,j), dummy)
!!     end do; end do
!! 
!!     sVdm = invert(Vdm,n)
!! 
!! end subroutine

end module
