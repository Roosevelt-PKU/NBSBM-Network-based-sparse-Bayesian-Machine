#dyn.load("nbsbc_refine.so")

dyn.load("nbsbc.so")

library(SparseM)

nbsbc_fast <- function(X, y, N, alpha = 0, beta = 0, a0 = a, b0 = b) {
	print("testing b")
	print(a0)
	print(b0)
	
	.Call("nbsbc", X, y, nrow(X), ncol(X), N, nrow(N), alpha, beta, a0, b0,
		quote(optim(c(a,b), function(x, K1, K2) (digamma(x[1]) - digamma(x[1] + x[2]) - K1)^2 +
		(digamma(x[2]) - digamma(x[1] + x[2]) - K2)^2, function(x, K1, K2) c(2 * (digamma(x[1]) -
		digamma(x[1] + x[2]) - K1) * (trigamma(x[1]) - trigamma(x[1]+ x[2])) - 2 * (digamma(x[2]) -
		digamma(x[1] + x[2]) - K2) * trigamma(x[1] + x[2]),-2 * (digamma(x[1]) - digamma(x[1] + x[2]) - K1) *
		trigamma(x[1] + x[2]) + 2 * (digamma(x[2])- digamma(x[1] + x[2]) - K2) * (trigamma(x[2]) - trigamma(x[1] + x[2]))),
		method = "L-BFGS-B", control = list(factr = 1e-20), K1 = K1, K2 = K2)$par), quote(runif(1)), environment())
}


nbsbc_slow <- function(X, y, N, alpha = 0, beta = 0, a0 = a, b0 = b) {

	# We initialize some constants

	d  <- ncol(X)
	st <- nrow(X)
	sn <- nrow(N)

	# We initialize the posterior approximation to be uniform

	a <- aOld <- 1
	b <- bOld <- 1
	p <- pOld <- rep(0.5, d)
	m <- mOld <- rep(0, d)
	v <- vOld <- rep(Inf, d)

	# We store the initial values to check for convergence

	aBackup <- a
	bBackup <- b
	pBackup <- p
	mBackup <- m
	vBackup <- v

	# We initialize the approximate terms

	firstTermMarkov  <- list(cTilde = rep(1, d), dTilde = rep(1, d), sTilde = 0)
	secondTermMarkov <- list(c1Tilde = rep(1, sn), d1Tilde = rep(1, sn), c2Tilde = rep(1, sn), d2Tilde = rep(1, sn), sTilde = rep(0, sn))
	termSpikeAndSlab <- list(cTilde = rep(1, d), dTilde = rep(1, d), mTilde = rep(0, d), vTilde = rep(Inf, d), sTilde = 0)
	termPriorEpsilon <- list(aTilde = 0, bTilde = 0, sTilde = 0)
	termLikelihood   <- list(mTilde = matrix(0, st, d), vTilde = matrix(Inf, st, d), aTilde = rep(0, st), bTilde = rep(0, st), sTilde = rep(0, st))

	iteration <- 0
	convergence <- FALSE
	indexSecondTermMarkov <- sample(1 : sn)
	while(!convergence && iteration < 100) {

		# We process the prior for epsilon

		aOld <- a - termPriorEpsilon$aTilde
		bOld <- b - termPriorEpsilon$bTilde
		a <- a0 + aOld - 1
		b <- b0 + bOld - 1
		termPriorEpsilon$aTilde <- a - aOld
		termPriorEpsilon$bTilde <- b - bOld
		termPriorEpsilon$sTilde <- -lbeta(a0, b0)

		# We process the first terms of the Markov random field prior

		pOld <- p / firstTermMarkov$cTilde / (p / firstTermMarkov$cTilde + (1 - p) / firstTermMarkov$dTilde)
		p <- 1 / (1 + (1 / pOld - 1) * exp(-2 * alpha))

		firstTermMarkov$cTilde <- p / pOld
		firstTermMarkov$dTilde <- (1 - p) / (1 - pOld)
		firstTermMarkov$sTilde <- sum(log(exp(alpha) * pOld + exp(-alpha*0) * (1 - pOld)))

		# We process the second terms of the Markov random field prior
		
		if (beta != 0) {

			for (i in indexSecondTermMarkov) {
				j <- N[ i, 1 ]
				l <- N[ i, 3 ]
                            d_j <- N[i, 4]
                            d_l <- N[i, 5]
				w <- N[ i, 2 ]
	
				pOld_j <- p[ j ] / secondTermMarkov$c1Tilde[ i ] /
					(p[ j ] / secondTermMarkov$c1Tilde[ i ] + (1 - p[ j ]) / secondTermMarkov$d1Tilde[ i ])
				pOld_l <- p[ l ] / secondTermMarkov$c2Tilde[ i ] /
					(p[ l ] / secondTermMarkov$c2Tilde[ i ] + (1 - p[ l ]) / secondTermMarkov$d2Tilde[ i ])
	
				if (pOld_j != 0.5 && pOld_l != 0.5) {
					Z1 <- pOld_j * pOld_l * exp(beta * w *(1/sqrt(d_j)-1/sqrt(d_l))^2)
					Z2 <- (1 - pOld_j) * (1 - pOld_l) * exp(beta * w * 0)
					Z3 <- (1 - pOld_j) * pOld_l * exp(beta * w /d_l)
					Z4 <- pOld_j * (1 - pOld_l) * exp(beta * w /d_j)
					p[ j ] <- (Z1 + Z4) / (Z1 + Z2 + Z3 + Z4)
					p[ l ] <- (Z1 + Z3) / (Z1 + Z2 + Z3 + Z4)
					secondTermMarkov$sTilde[ i ] <- log(Z1 + Z2 + Z3 + Z4)
				}
	
				if (p[ j ] < 1e-15)
					p[ j ] <- 1e-15
				if (p[ j ] > 1 - 1e-15)
					p[ j ] <- 1 - 1e-15
	
				if (p[ l ] < 1e-15)
					p[ l ] <- 1e-15
				if (p[ l ] > 1 - 1e-15)
					p[ l ] <- 1 - 1e-15
	
				secondTermMarkov$c1Tilde[ i ] <- p[ j ] / pOld_j
				secondTermMarkov$d1Tilde[ i ] <- (1 - p[ j ]) / (1 - pOld_j)
				secondTermMarkov$c2Tilde[ i ] <- p[ l ] / pOld_l
				secondTermMarkov$d2Tilde[ i ] <- (1 - p[ l ]) / (1 - pOld_l)
			}
		}
	
		# We process the term of the spike and slab prior

		pOld <- p / termSpikeAndSlab$cTilde / (p / termSpikeAndSlab$cTilde + (1 - p) / termSpikeAndSlab$dTilde)
		vOld <- 1 / (1 / v - 1 / termSpikeAndSlab$vTilde)

		index <- which(is.infinite(vOld))

		mOld[ index ] <- m[ index ] <- 0
		v[ index ] <- p[ index ] <- pOld[ index ]
		termSpikeAndSlab$vTilde[ index ] <- v[ index ]
		termSpikeAndSlab$mTilde[ index ] <- 0
		termSpikeAndSlab$cTilde[ index ] <- 1
		termSpikeAndSlab$dTilde[ index ] <- 1

		index <- which(is.finite(vOld) & vOld > 0)
		if (length(index) > 0) {

			if (any(vOld < 0)) cat("*")
	
			pOld <- pOld[ index ]
			vOld <- vOld[ index ]
			mOld <- m[ index ] + vOld * 1 / termSpikeAndSlab$vTilde[ index ] * (m[ index ] - termSpikeAndSlab$mTilde[ index ])
	
			Z <- pOld * 1 / sqrt(2 * pi * (vOld + 1)) * exp(-0.5 * mOld * mOld / (vOld + 1)) + (1 - pOld) *
				1 / sqrt(2 * pi * vOld) * exp(-0.5 * mOld * mOld / vOld)
			C <- 1 / (1 + (1 - pOld) / pOld * sqrt(vOld + 1) / sqrt(vOld) * exp(-0.5 * mOld * mOld *
				(1 / vOld - 1 / (vOld + 1)))) * - mOld / (vOld + 1) +
				1 / (1 + pOld / (1 - pOld) * sqrt(vOld) / sqrt(vOld + 1) * exp(-0.5 * mOld * mOld *
				(1 / (vOld + 1 ) - 1 / vOld))) *  - mOld / vOld
			Cprima <- 1 / (1 + (1 - pOld) / pOld * sqrt(vOld + 1) / sqrt(vOld) * exp(-0.5 * mOld * mOld *
	                        (1 / vOld - 1 / (vOld + 1)))) * 0.5 * (mOld * mOld / ((vOld + 1) * (vOld + 1)) - 1 / (vOld + 1)) +
	                        1 / (1 + pOld / (1 - pOld) * sqrt(vOld) / sqrt(vOld + 1) * exp(-0.5 * mOld * mOld * (1 / (vOld + 1) - 1 / vOld))) *  
	                        0.5 * (mOld * mOld / (vOld * vOld) - 1 / vOld)
			Cprimaprima <- C * C - 2 * Cprima
	
			v[ index ] <- vOld - Cprimaprima * vOld * vOld
			m[ index ] <- mOld + C * vOld
			p[ index ] <- 1 / (1 + (1 - pOld) / pOld * sqrt(vOld + 1) / sqrt(vOld) * exp(-0.5 * mOld * mOld * (1 / vOld - 1 / (vOld + 1))))

			index2 <- which(p < 1e-15)
			p[ index2 ] <- 1e-15
			index2 <- which(p > 1- 1e-15)
			p[ index2 ] <- 1 - 1e-15
			v[ index[ which(v[ index ] < 1e-15) ] ] <- 1e-15
		
			termSpikeAndSlab$vTilde[ index ] <- (1 / Cprimaprima - vOld)
			termSpikeAndSlab$mTilde[ index ] <- mOld + C * (termSpikeAndSlab$vTilde[ index ] + vOld)
			termSpikeAndSlab$cTilde[ index ] = p[ index ] / pOld
			termSpikeAndSlab$dTilde[ index ] = (1 - p[ index ]) / (1 - pOld)
	
			termSpikeAndSlab$sTilde <- sum(log(Z) + log(sqrt(1  + vOld / termSpikeAndSlab$vTilde[ index ])) + 0.5 * C^2 / Cprimaprima)
		}

		# We process the terms for the likelihood

		for (i in 1 : st) {
			aOld <- a - termLikelihood$aTilde[ i ]
			bOld <- b - termLikelihood$bTilde[ i ]
			vOld <- 1 / (1 / v - 1 / termLikelihood$vTilde[ i, ])

			if (any(vOld < 0)) { cat("*"); next }

			mOld <- m + vOld * 1 / termLikelihood$vTilde[ i, ] * (m - termLikelihood$mTilde[ i, ])
			squareRoot <- sum(vOld * X[ i, ] * X[ i, ])

			squareRoot <- sqrt(squareRoot)
			z <- y[ i ] * sum(X[ i, ] * mOld) / squareRoot
			eHatOld <- aOld / (aOld + bOld)
			Phi <- pnorm(z)
			lambda <- 1 / squareRoot * (1 - 2 * eHatOld) * 1 / sqrt(2 * pi) * exp(-0.5 * z * z) / (eHatOld + (1 - 2 * eHatOld) * Phi)
			Z <- eHatOld + (1 - 2 * eHatOld) * Phi
	
			m <- mOld + y[ i ] * lambda * vOld * X[ i, ]
			v <- vOld - y[ i ] * lambda * sum(X[ i, ] * m) / (squareRoot * squareRoot) * vOld * X[ i, ] * vOld * X[ i, ]
			v[ which(v < 1e-15) ] <- 1e-15

			termLikelihood$vTilde[ i, ] <- 1 / (1 / v - 1 / vOld)
			termLikelihood$mTilde[ i, ] <- mOld + y[ i ] * lambda * X[ i, ] * (termLikelihood$vTilde[ i, ] + vOld)
			termLikelihood$mTilde[ i, is.infinite(termLikelihood$vTilde[ i, ]) ] <- 0

			K1 <- eHatOld * (1 - Phi) / (aOld * (eHatOld + (1 - 2 * eHatOld) * Phi)) + digamma(aOld) - digamma(aOld + bOld + 1)
			K2 <- (1 - eHatOld) * Phi / (bOld * (eHatOld + (1 - 2 * eHatOld) * Phi)) + digamma(bOld) - digamma(aOld + bOld + 1)
	
			ret <- optim(c(a, b), function(x, K1, K2) (digamma(x[1]) - digamma(x[1] + x[2]) - K1)^2 + (digamma(x[2]) -
				digamma(x[1] + x[2]) - K2)^2, function(x, K1, K2) c(2 * (digamma(x[1]) - digamma(x[1] + x[2]) -
				K1) * (trigamma(x[1]) - trigamma(x[1]+ x[2])) - 2 * (digamma(x[2]) - digamma(x[1] + x[2]) - K2) *
				trigamma(x[1] + x[2]),-2 * (digamma(x[1]) - digamma(x[1] + x[2]) - K1) * trigamma(x[1] + x[2]) + 2 *
				(digamma(x[2])- digamma(x[1] + x[2]) - K2) * (trigamma(x[2]) - trigamma(x[1] + x[2]))), method =
				"L-BFGS-B", control = list(factr = 1e-20), K1 = K1, K2 = K2)$par


			a <- ret[ 1 ]
			b <- ret[ 2 ]

			termLikelihood$aTilde[ i ] <- a - aOld
			termLikelihood$bTilde[ i ] <- b - bOld

			termLikelihood$sTilde[ i ] <- log(Z) + sum(log(sqrt(1 + vOld / termLikelihood$vTilde[ i, ]))) +
				0.5 * sum((termLikelihood$mTilde[ i, ] - mOld)^2 / (termLikelihood$vTilde[ i, ] + vOld)) + lbeta(aOld, bOld) -
				lbeta(a, b)
		}

		iteration <- iteration + 1

		maxdiff <- max(abs(m - mBackup), abs(v - vBackup), abs(p - pBackup), abs(a - aBackup), abs(b - bBackup))
		if (maxdiff < 1e-6)
			convergence <- TRUE

		cat("Iteration", iteration, "-", maxdiff, a, b, "\n")

		aBackup <- a
		bBackup <- b
		pBackup <- p
		mBackup <- m
		vBackup <- v
	}
	
	# We approximate the evidence

	logEvidence <- firstTermMarkov$sTilde + sum(secondTermMarkov$sTilde) + termSpikeAndSlab$sTilde + termPriorEpsilon$sTilde + sum(termLikelihood$sTilde)

	A <- 1 + sum(termLikelihood$aTilde) + termPriorEpsilon$aTilde
	B <- 1 + sum(termLikelihood$bTilde) + termPriorEpsilon$bTilde
	D <- sum(m^2 / v) - sum(termLikelihood$mTilde^2 / termLikelihood$vTilde) - sum(termSpikeAndSlab$mTilde^2 / termSpikeAndSlab$vTilde)

	logEvidence <- logEvidence + lbeta(A, B) + D / 2 + d / 2 * log(2 * pi) + 0.5 * sum(log(v))

	cProduct <- rep(1, d)
	dProduct <- rep(1, d)

	cProduct <- cProduct * firstTermMarkov$cTilde * termSpikeAndSlab$cTilde
	dProduct <- dProduct * firstTermMarkov$dTilde * termSpikeAndSlab$dTilde

	for (i in 1 : sn) {
		j <- N[ i, 1 ]
		l <- N[ i, 3 ]
		cProduct[ j ] <- cProduct[ j ] * secondTermMarkov$c1Tilde[ i ]
		dProduct[ j ] <- dProduct[ j ] * secondTermMarkov$d1Tilde[ i ]
		cProduct[ l ] <- cProduct[ l ] * secondTermMarkov$c2Tilde[ i ]
		dProduct[ l ] <- dProduct[ l ] * secondTermMarkov$d2Tilde[ i ]
	}

	C <- sum(log(cProduct + dProduct))

	# We estimate the partition function of the Markov Random Field prior

	partitionFunction <- approximatePartitionFunctionEP(N, alpha, beta, d)

	logEvidence <- logEvidence + C - partitionFunction

	list(m = m, v = v, p = p, a = a, b = b, logEvidence = logEvidence, error = a / (a + b))
}

##
# EP method for approximating the partition function of a Markov Random field model
#

approximatePartitionFunctionEP <- function(N, alpha, beta, d) {

	# We initialize some constants

	sn <- nrow(N)

	# We initialize the posterior approximation to be uniform

	p <- pOld <- rep(0.5, d)

	# We store the initial values to check for convergence

	pBackup <- p

	# We initialize the approximate terms

	firstTermMarkov  <- list(cTilde = rep(1, d), dTilde = rep(1, d), sTilde = 0)
	secondTermMarkov <- list(c1Tilde = rep(1, sn), d1Tilde = rep(1, sn), c2Tilde = rep(1, sn), d2Tilde = rep(1, sn), sTilde = rep(0, sn))

	iteration <- 0
	convergence <- FALSE
	indexSecondTermMarkov <- sample(1 : sn)
	while(!convergence && iteration < 100) {

		# We process the first terms of the Markov random field prior

		pOld <- p / firstTermMarkov$cTilde / (p / firstTermMarkov$cTilde + (1 - p) / firstTermMarkov$dTilde)
		p <- 1 / (1 + (1 / pOld - 1) * exp(-2 * alpha))

		firstTermMarkov$cTilde <- p / pOld
		firstTermMarkov$dTilde <- (1 - p) / (1 - pOld)
		firstTermMarkov$sTilde <- sum(log(exp(alpha) * pOld + exp(-alpha*0) * (1 - pOld)))

		# We process the second terms of the Markov random field prior
		
		if (beta != 0) {
			for (i in indexSecondTermMarkov) {
				j <- N[ i, 1 ]
				l <- N[ i, 3 ]
                            d_j <- N[i, 4]
                            d_l <- N[i, 5]
				w <- N[ i, 2 ]
	
				pOld_j <- p[ j ] / secondTermMarkov$c1Tilde[ i ] /
					(p[ j ] / secondTermMarkov$c1Tilde[ i ] + (1 - p[ j ]) / secondTermMarkov$d1Tilde[ i ])
				pOld_l <- p[ l ] / secondTermMarkov$c2Tilde[ i ] /
					(p[ l ] / secondTermMarkov$c2Tilde[ i ] + (1 - p[ l ]) / secondTermMarkov$d2Tilde[ i ])
	
				Z1 <- pOld_j * pOld_l * exp(beta * w *(1/sqrt(d_j)-1/sqrt(d_l))^2)
				Z2 <- (1 - pOld_j) * (1 - pOld_l) * exp(beta * w * 0)
				Z3 <- (1 - pOld_j) * pOld_l * exp(beta * w / d_l)
				Z4 <- pOld_j * (1 - pOld_l) * exp(beta * w / d_j)
				p[ j ] <- (Z1 + Z4) / (Z1 + Z2 + Z3 + Z4)
				p[ l ] <- (Z1 + Z3) / (Z1 + Z2 + Z3 + Z4)
				secondTermMarkov$sTilde[ i ] <- log(Z1 + Z2 + Z3 + Z4)
	
				if (p[ j ] < 1e-15)
					p[ j ] <- 1e-15
				if (p[ j ] > 1 - 1e-15)
					p[ j ] <- 1 - 1e-15
	
				if (p[ l ] < 1e-15)
					p[ l ] <- 1e-15
				if (p[ l ] > 1 - 1e-15)
					p[ l ] <- 1 - 1e-15
	
				secondTermMarkov$c1Tilde[ i ] <- p[ j ] / pOld_j
				secondTermMarkov$d1Tilde[ i ] <- (1 - p[ j ]) / (1 - pOld_j)
				secondTermMarkov$c2Tilde[ i ] <- p[ l ] / pOld_l
				secondTermMarkov$d2Tilde[ i ] <- (1 - p[ l ]) / (1 - pOld_l)
			}
		}
	
		iteration <- iteration + 1

		maxdiff <- max(abs(p - pBackup))
		if (maxdiff < 1e-6)
			convergence <- TRUE

		pBackup <- p
	}
	
	# We approximate the evidence

	logEvidence <- firstTermMarkov$sTilde + sum(secondTermMarkov$sTilde)

	cProduct <- firstTermMarkov$cTilde 
	dProduct <- firstTermMarkov$dTilde

	for (i in 1 : sn) {
		j <- N[ i, 1 ]
		l <- N[ i, 3 ]
		cProduct[ j ] <- cProduct[ j ] * secondTermMarkov$c1Tilde[ i ]
		dProduct[ j ] <- dProduct[ j ] * secondTermMarkov$d1Tilde[ i ]
		cProduct[ l ] <- cProduct[ l ] * secondTermMarkov$c2Tilde[ i ]
		dProduct[ l ] <- dProduct[ l ] * secondTermMarkov$d2Tilde[ i ]
	}

	C <- sum(log(cProduct + dProduct))

	logEvidence <- logEvidence + C

	logEvidence
}
