!
!
 SUBROUTINE RXN(n,Kl0,vpl0,vba,lkab,lkba)
! **********************************************************************
! *                                                                    *
! * Calculates log rate constants between 2 states                     *
! *   using enhanced exponential lumping                               *
! **********************************************************************
! *                                                                    *
! * Author:     Boris Fackovec                                         *
! * Date:       03/02/2016                                             *
! * Version:    1.1                                                    *
! *                                                                    *
! **********************************************************************
!
   USE ParamsRXN, ONLY: howFine, lstartdt, thresh
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n
   DOUBLE PRECISION, INTENT(IN), DIMENSION(n) :: vpl0
   DOUBLE PRECISION, INTENT(IN), DIMENSION(n,n) :: Kl0
   LOGICAL, INTENT(IN), DIMENSION(n) :: vba
   DOUBLE PRECISION, INTENT(OUT) :: lkab, lkba
!
! -------------------------------------------------------------------
!
   INTEGER :: i, j, k
   DOUBLE PRECISION :: pl, ltemp, xl, doldd, dmax, pl1, pl2
   DOUBLE PRECISION :: lds, ldt, ltnow, ltau
   INTEGER :: nFine, na
   DOUBLE PRECISION, DIMENSION(n,n) :: Kl, L1, D1, Ls, Ds, Lp, Dp, &
        Dold, Ltem, Dtem
   LOGICAL, DIMENSION(n,n) :: Nb1, Nbs, Nbp, Nbtem
   DOUBLE PRECISION, DIMENSION(n) :: vlnow, vpla, vpl, vtemp
   DOUBLE PRECISION, DIMENSION(2) :: vUjo
   DOUBLE PRECISION, DIMENSION(3*n) :: vtr
   LOGICAL :: btemp, bl
   LOGICAL, DIMENSION(3*n) :: vbtr
   DOUBLE PRECISION, DIMENSION(9999,2) :: X ! max ratio of timescales = 2^(9999/(2^3))
   INTEGER :: sizeX, nind
!   PARAMETER (howFine=3)
!   PARAMETER (lstartdt = -6d0)
!   PARAMETER (thresh = -6d0)
!
! -------------------------------------------------------------------
!
   WRITE(*,*) "Subroutine RXN"
!
!  make sure populations sum to 1 and rate mat satisfies detailed bal
!
   CALL LogSumExp(n,vpl0,pl)
   cycNormvpl: DO i=1,n
      vpl(i) = vpl0(i) - pl
   ENDDO cycNormvpl
   Kl = Kl0
   CALL SymmetriseRateMat(n,Kl,vpl)
!
!  calculate quilibrium populations
!
   WRITE(*,*) "Calc equil populs"
   na = 0
   cycCountA: DO i=1,n
      IF (vba(i)) THEN 
         na = na + 1 
         vpla(na) = vpl(i)
      ENDIF
   ENDDO cycCountA
   CALL LogSumExp(na,vpla,pl1)
   CALL LogDiffExp(0.0d0,pl1,pl2)
!
!  calculate optimal minimum (log) time step lds
!
   nFine = 2**howFine
   cycDiagRates: DO j=1,n
      cycRows: DO i=1,n
         vlnow(i) = Kl(i,j)
      ENDDO cycRows
      CALL LogSumExp(n,vlnow,pl)
      vtemp(j) = pl
   ENDDO cycDiagRates
   lds = -MAXVAL(vtemp) + lstartdt - REAL(howFine)*LOG(2.0d0)
!
!  produce log transition matrix
!
   WRITE(*,*) "Make the log trans mat"
   cycTMCol: DO j=1,n
      cycTMRow: DO i=1,n
         Ls(i,j) = Kl(i,j) + lds
         vlnow(i) = Ls(i,j)
      ENDDO cycTMRow
      vlnow(j) = vlnow(n)
      CALL LogSumExp(n-1,vlnow,ltemp)
      CALL LogDiffExp(0.0d0,ltemp,pl)
!      vUjo(1) = pl
!      vUjo(2) = Ls(j,j)
!      CALL LogSumExp(2,vUjo,ltemp)
      Ls(j,j) = pl
   ENDDO cycTMCol
!
!  produce D and N matrices for log(1-t) formalism
!
   WRITE(*,*) "produce D and N matrices"
   cycDoRows: DO i=1,n
      cycDoCols: DO j=1,n
         IF (Ls(i,j).GT.vpl(i)) THEN
            Nbs(i,j) = .TRUE.
            CALL LogDiffExp(0.0d0,vpl(i)-Ls(i,j),ltemp)
            Ds(i,j) = ltemp + Ls(i,j) - vpl(i)
         ELSE
            Nbs(i,j) = .FALSE.
            CALL LogDiffExp(0.0d0,Ls(i,j)-vpl(i),ltemp)
            Ds(i,j) = ltemp
         ENDIF
      ENDDO cycDoCols
   ENDDO cycDoRows
!
!  log transition matrices for dt = nFine*exp(lds)
!
   WRITE(*,*) "get the matrices for nFine*lds"
   L1 = Ls
   D1 = Ds
   Nb1 = Nbs
   cycDoubling: DO i=1,howFine
      WRITE(*,*) "which fine", i, "out of", howFine
      CALL MultLogMat(n,vpl,L1,D1,Nb1,L1,D1,Nb1,L1,D1,Nb1)
   ENDDO cycDoubling
   xl = 0
   X(1,1) = -99e9
   X(1,2) = xl
   ldt = lds + REAL(howFine)*LOG(2.0d0)
   doldd = 1.0d1
   dmax = 1.0d1
   nind = 2
!
!  main cycle
!
   WRITE(*,*) "Main cycle"
   cycMain: DO WHILE( ((xl.GT.-5).OR.(doldd.GT.-3)) .AND.(dmax.GT.thresh))
      Dold = D1
      CALL LPopul(n,L1,D1,Nb1,vpl,vba,xl)
      X(nind,1) = ldt
      X(nind,2) = xl
      WRITE(*,*) ldt, xl
      nind = nind + 1
      Lp = L1
      Dp = D1
      Nbp = Nb1
      ltnow = ldt
      cycLinProp: DO i=1,nFine-1
         CALL MultLogMat(n,vpl,Ltem,Dtem,Nbtem,Lp,Dp,Nbp,Ls,Ds,Nbs)
         Lp = Ltem
         Dp = Dtem
         Nbp = Nbtem
         CALL LPopul(n,Lp,Dp,Nbp,vpl,vba,xl)
         vUjo(1) = ltnow
         vUjo(2) = lds
         CALL LogSumExp(2,vUjo,ltemp)
         ltnow = ltemp
         X(nind,1) = ltnow
         X(nind,2) = xl
         WRITE(*,*) ltnow, xl
         nind = nind + 1
      ENDDO cycLinProp
      CALL MultLogMat(n,vpl,Ltem,Dtem,Nbtem,L1,D1,Nb1,L1,D1,Nb1)
      L1 = Ltem
      D1 = Dtem
      Nb1 = Nbtem
      CALL MultLogMat(n,vpl,Ltem,Dtem,Nbtem,Ls,Ds,Nbs,Ls,Ds,Nbs)
      Ls = Ltem
      Ds = Dtem
      Nbs = Nbtem
      ldt = ldt + LOG(2.0d0)
      lds = lds + LOG(2.0d0)
      doldd = LOG(SUM(ABS(D1-Dold)))
      dmax = MAXVAL(D1)
   ENDDO cycMain
   CALL TrapzLog(nind,X,ltau)
!   WRITE(*,*) "vypoctov", nind
   lkab = pl1 - ltau
   lkba = pl2 - ltau
!
   RETURN
!
! -------------------------------------------------------------------
!
 END SUBROUTINE RXN
!
!
!
 SUBROUTINE SymmetriseRateMat(n,Kl,vpl)
! **********************************************************************
! *                                                                    *
! * Ensures the rate matrix satisfies detailed balance                 *
! *   Kl must have -99e9 for zero rates and -1e99 at its diagonal      *
! *                                                                    *
! **********************************************************************
! *                                                                    *
! * Author:     Boris Fackovec                                         *
! * Date:       02/02/2016                                             *
! * Version:    1.0                                                    *
! *                                                                    *
! **********************************************************************
!
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n
   DOUBLE PRECISION, INTENT(IN), DIMENSION(n) :: vpl
   DOUBLE PRECISION, INTENT(INOUT), DIMENSION(n,n) :: Kl
!
! -------------------------------------------------------------------
!
   INTEGER :: i, j
   DOUBLE PRECISION :: recTau, pl
   DOUBLE PRECISION, DIMENSION(2) :: vUjo
   LOGICAL, DIMENSION(3*n) :: vbtr
!
! -------------------------------------------------------------------
!
!   cycDoDiags: DO i=1,n
!      Kl(i,i) = -1e99
!   ENDDO cycDoDiags
   cycDoRows: DO i=1,n-1
      cycDoCols: DO j=i+1,n
         vUjo(1) = Kl(i,j)
         vUjo(2) = Kl(j,i)
         CALL LogSumExp(2,vUjo,recTau)
         vUjo(1) = vpl(i)
         vUjo(2) = vpl(j)
         CALL LogSumExp(2,vUjo,pl)
         Kl(i,j) = vpl(i) + recTau - pl
         Kl(j,i) = vpl(j) + recTau - pl
      ENDDO cycDoCols
   ENDDO cycDoRows
!
   RETURN
!
! -------------------------------------------------------------------
!
 END SUBROUTINE SymmetriseRateMat
!
!
!
 SUBROUTINE TrapzLog(n,X,ltau)
! **********************************************************************
! *                                                                    *
! * Integrates a vector in log formalism usind trapezoidal ruel        *
! *                                                                    *
! **********************************************************************
! *                                                                    *
! * Author:     Boris Fackovec                                         *
! * Date:       02/02/2016                                             *
! * Version:    1.0                                                    *
! *                                                                    *
! **********************************************************************
!
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n
   DOUBLE PRECISION, INTENT(IN), DIMENSION(n,2) :: X
   DOUBLE PRECISION, INTENT(OUT) :: ltau
!
! -------------------------------------------------------------------
!
   INTEGER :: i
   DOUBLE PRECISION :: ldt, pl
   DOUBLE PRECISION, DIMENSION(2) :: vUjo
!
! -------------------------------------------------------------------
!
   ltau = -9.9d99
   cycTrapez: DO i=1,n-1
      CALL LogDiffExp(X(i+1,1),X(i,1),ldt)
      vUjo(1) = X(i,2)
      vUjo(2) = X(i+1,2)
      CALL LogSumExp(2,vUjo,pl)
      pl = pl - LOG(2.0d0)
      vUjo(1) = ltau
      vUjo(2) = ldt + pl
      CALL LogSumExp(2,vUjo,ltau)
   ENDDO cycTrapez
!
   RETURN
!
! -------------------------------------------------------------------
!
 END SUBROUTINE TrapzLog
!