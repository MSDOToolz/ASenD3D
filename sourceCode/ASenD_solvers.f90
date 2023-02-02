module ASenD_solvers

use ASenD_constantVals

contains

    ! Construct the [L][D][L]T factorization of a sparsely populated symmetric matrix.  The result is over-written into
    ! the lower triangular matrix A, passed in row-skyline format.  Other inputs are the total size (number of members)
    ! of the matrix A, the row ranges of ALTri, and the matrix dimension.
    subroutine getSparseLDLFact(ALTri, ASize, ARange, n)
		implicit none
		
		integer, intent(in) :: ASize, n
		integer, intent(in) :: ARange(0:n)
		real*8, intent(out) :: ALTri(ASize)

		integer, allocatable :: rowStCol(:)	
		integer :: i1, i2, i3, i4, i5, i6, imax, stCol, endCol, maxBW, rowLen
		real*8 :: fact, diag
		real*8, allocatable :: diagvec(:), dLprod(:)

		allocate(rowStCol(n),diagvec(n))

		maxBw = 0
		do i1 = 1, n
			rowLen = ARange(i1) - ARange(i1-1)
			rowStCol(i1) = i1 - rowLen + 1
			if(rowLen .gt. maxBW) then
			maxBW = rowLen
			endif
		enddo

		allocate(dLprod(maxBW))

		do i1 = 1, n
			 !if(mod(i1,n/20) .eq. 0) then
				 !    write(*,*) i1*99/n + 1, '% complete'
			 !endif
			i3 = ARange(i1) !i3 = index in ALTri of diagonal term in row i1
			i2 = ARange(i1-1) + 1  !i2 = index in ALTri of starting term in row i1
				!Calculate the product of the i1 row of ALTri with the corresponding diagonal terms
			i5 = maxBW
			i6 = i1 - 1
			do i4 = i3 - 1, i2, -1
				dLprod(i5) = diagvec(i6)*ALTri(i4)
				i5 = i5 - 1
			    i6 = i6 - 1
			enddo
				!Calculate the diagonal term of row i1, storing diagvec/ALTri
			i5 = i1 - i3 + i2  !i5 = current column in [A], corresponding to i4
			diag = ALTri(i3)
			do i4 = i2, i3-1
				diag = diag - diagvec(i5)*ALTri(i4)*ALTri(i4)
				i5 = i5 + 1
			enddo
			diagvec(i1) = diag
			ALTri(i3) = diag
			fact = r_1/diag
			imax = min(n, i1+maxBW)
			do i2 = i1+1, imax
				stCol = max(rowStCol(i1),rowStCol(i2))
				rowLen = i1 - stCol
				if(rowLen .ge. 0) then
					i3 = maxBW - rowLen + 1 !Start index of dLprod
					i4 = ARange(i2) - (i2 - stCol) !Start index of row 2
					i5 = i4 + rowLen !Index of ALTri value being computed
					diag = ALTri(i5)
					do i6 = 0, rowLen - 1
						diag = diag - dLprod(i3+i6)*ALTri(i4+i6)
					enddo
					ALTri(i5) = diag*fact
				endif
			enddo
		enddo

		deallocate(rowStCol,diagvec,dLprod)
		
		return
    end subroutine getSparseLDLFact

    ! Solves [A]{soln} = {rhs}, given the [L][D][L]T factorization of a sparse symmetric matrix A.  Other inputs are the
    ! size (total number of entries) of the matrix factorization, the row range of the matrix factorization and the 
    ! dimension of the matrix.
    subroutine solveSparseLDLFact(soln, rhs, Lmat, LSize, LRange, Ldim, zVec)
		implicit none
		
		integer, intent(in) :: LSize, Ldim
		integer, intent(in) :: LRange(0:Ldim)
		real*8, intent(in) :: rhs(Ldim), Lmat(LSize)
		real*8, intent(out) :: soln(Ldim)
		real*8, intent(out) :: zVec(Ldim)

		real*8 :: dtProd, diag
		integer :: i1, i2, i3, imax, stCol, rowLen, maxBW
		
		maxBW = 0
		do i1 = 1, Ldim
			i2 = LRange(i1-1) + 1
			rowLen = LRange(i1) - i2
			if(rowLen .gt. maxBW) then
			    maxBW = rowLen
			endif
			stCol = i1 - rowLen
			zVec(i1) = rhs(i1)
			do i3 = 0, rowLen-1
				zVec(i1) = zVec(i1) - Lmat(i2+i3)*zVec(stCol+i3)
			enddo
		enddo

		do i1 = Ldim, 1, -1
			dtProd = r_0
			diag = Lmat(LRange(i1))
			imax = min(Ldim,i1+maxBW)
			do i2 = i1+1, imax
				rowLen = LRange(i2) - LRange(i2-1)
				stCol = i2 - rowLen + 1
				if(stCol .le. i1) then
					i3 = LRange(i2) + i1 - i2
					dtProd = dtProd + Lmat(i3)*soln(i2)
				endif
			enddo
			dtProd = dtProd*diag
			soln(i1) = (zvec(i1) - dtProd)/diag
		enddo
		
		return
    end subroutine solveSparseLDLFact
    
    subroutine gMRes(xVec,Amat,Acols,ARange,ASize,ADim,Cmat,Ccols,CRange,CSize,CDim,Pmat,PRange,PSize,bVec,numVecs,maxIt)
	    implicit none
		
		integer, intent(in) :: ASize, ADim, CSize, CDim, PSize, numVecs, maxIt
		integer, intent(in) :: ACols(ASize), ARange(0:ADim), Ccols(CSize), CRange(0:CDim), PRange(0:ADim)
		real*8, intent(in) :: Amat(ASize), Cmat(CSize), Pmat(PSize), bVec(ADim)
		real*8, intent(out) :: xVec(ADim)
		
		real*8 :: resMag, dp
		integer :: i1, i2, i3, i4, i5, i6, maxOutIt, breakSrch, numR, numC
		real*8, allocatable :: zVec(:), rVec(:), tempV1(:), tempV2(:), hMat(:,:), h2Ah(:,:)
		
		allocate(zVec(ADim),rVec(ADim),tempV1(ADim),tempV2(ADim),hMat(ADim,numVecs+1),h2Ah(numVecs+1,numVecs))
		
		tempV1(:) = r_0
		do i1 = 1, CDim
		    do i2 = CRange(i1-1)+1, CRange(i1)
			    i3 = Ccols(i2)
				if(i3 .ne. 0) then
				    tempV1(i1) = tempV1(i1) + Cmat(i2)*xVec(i3)
				endif
			enddo
		enddo
		
		tempV2(:) = r_0
		do i1 = 1, CDim
		    do i2 = CRange(i1-1)+1, CRange(i1)
			    i3 = Ccols(i2)
				if(i3 .ne. 0) then
				    tempV2(i3) = tempV2(i3) + Cmat(i2)*tempV1(i1)
				endif
			enddo
		enddo
		
		do i1 = 1, ADim
		    do i2 = ARange(i1-1)+1, ARange(i1)
			    i3 = Acols(i2)
				if(i3 .ne. 0) then
				    tempV2(i1) = tempV2(i1) + Amat(i2)*xVec(i3)
				endif
			enddo
		enddo
		
		tempV1(:) = bVec(:) - tempV2(:)
		
		call solveSparseLDLFact(rVec, tempV1, Pmat, PSize, PRange, ADim, zVec)
		!rVec(:) = tempV1(:)
		
		resMag = r_0
		do i1 = 1, ADim
		    resMag = resMag + rVec(i1)*rVec(i1)
		enddo
		resMag = sqrt(resMag)
		
		maxOutIt = maxIt/numVecs + 1
		i4 = 0
		
		write(*,*) 'GMRes initial resMag: ', resMag, 'maxOutIt: ', maxOutIt
		do while(resMag .gt. 1e-12 .and. i4 .lt. maxOutIt)
		    hMat(:,1) = (r_1/resMag)*rVec(:)
			h2Ah(:,:) = r_0
			i5 = 2
			breakSrch = 0
			do while(i5 .le. numVecs+1 .and. breakSrch .eq. 0)
			    tempV1(:) = r_0
				do i1 = 1, CDim
					do i2 = CRange(i1-1)+1, CRange(i1)
						i3 = Ccols(i2)
						if(i3 .ne. 0) then
							tempV1(i1) = tempV1(i1) + Cmat(i2)*hMat(i3,i5-1)
						endif
					enddo
				enddo
				
				tempV2(:) = r_0
				do i1 = 1, CDim
					do i2 = CRange(i1-1)+1, CRange(i1)
						i3 = Ccols(i2)
						if(i3 .ne. 0) then
							tempV2(i3) = tempV2(i3) + Cmat(i2)*tempV1(i1)
						endif
					enddo
				enddo
				
				do i1 = 1, ADim
					do i2 = ARange(i1-1)+1, ARange(i1)
						i3 = Acols(i2)
						if(i3 .ne. 0) then
							tempV2(i1) = tempV2(i1) + Amat(i2)*hMat(i3,i5-1)
						endif
					enddo
				enddo
				
				call solveSparseLDLFact(hMat(:,i5), tempV2, Pmat, PSize, PRange, ADim, zVec)
				
				do i1 = 1, i5-1
				    dp = r_0
					do i2 = 1, ADim
					    dp = dp + hMat(i2,i1)*hMat(i2,i5)
					enddo
					hMat(:,i5) = hMat(:,i5) - dp*hMat(:,i1)
					h2Ah(i1,i5-1) = dp
				enddo
				dp = r_0
				do i2 = 1, ADim
				    dp = dp + hMat(i2,i5)*hMat(i2,i5)
				enddo
				dp = sqrt(dp)
				if(dp .gt. 1e-9_8) then
				    hMat(:,i5) = (r_1/dp)*hMat(:,i5)
				    h2Ah(i5,i5-1) = dp
				    i5 = i5 + 1
				else
				    breakSrch = 1
				endif
			enddo
			
			numR = i5 - 1
			numC = i5 - 2
			call rFactorMat(h2Ah,numVecs+1,numVecs,1,numR,1,numC,2)
			tempV1(:) = r_0
			do i1 = 1, ADim
			    do i2 = 1, numR
				    tempV1(i2) = tempV1(i2) + hMat(i1,i2)*rVec(i1)
				enddo
			enddo
			i1 = numVecs + 1
			tempV2(1:numC) = r_0
			call solveRFactor(tempV2(1:numC),h2Ah,tempV1(1:numR),i1,numVecs,1,numR,1,numC,2)
			
			do i1 = 1, ADim
			    do i2 = 1, numC
				   xVec(i1) = xVec(i1) + hMat(i1,i2)*tempV2(i2)
				enddo
			enddo
			
			tempV1(:) = r_0
			do i1 = 1, CDim
				do i2 = CRange(i1-1)+1, CRange(i1)
					i3 = Ccols(i2)
					if(i3 .ne. 0) then
						tempV1(i1) = tempV1(i1) + Cmat(i2)*xVec(i3)
					endif
				enddo
			enddo
			
			tempV2(:) = r_0
			do i1 = 1, CDim
				do i2 = CRange(i1-1)+1, CRange(i1)
					i3 = Ccols(i2)
					if(i3 .ne. 0) then
						tempV2(i3) = tempV2(i3) + Cmat(i2)*tempV1(i1)
					endif
				enddo
			enddo
			
			do i1 = 1, ADim
				do i2 = ARange(i1-1)+1, ARange(i1)
					i3 = Acols(i2)
					if(i3 .ne. 0) then
						tempV2(i1) = tempV2(i1) + Amat(i2)*xVec(i3)
					endif
				enddo
			enddo
			
			tempV1(:) = bVec(:) - tempV2(:)
			
			call solveSparseLDLFact(rVec, tempV1, Pmat, PSize, PRange, ADim, zVec)
			
			resMag = r_0
			do i1 = 1, ADim
				resMag = resMag + rVec(i1)*rVec(i1)
			enddo
			resMag = sqrt(resMag)
			
			write(*,*) 'resMag: ', resMag
			
			i4 = i4 + 1
		enddo
		
		write(*,*) 'finished GMRes solve, i4 = ', i4, ' resMag = ', resMag
		
		deallocate(zVec,rVec,tempV1,tempV2,hMat,h2Ah)
		
	end subroutine gMRes

    subroutine rFactorMat(Amat,rowDim,colDim,rowSt,rowEnd,colSt,colEnd,triDiag)
        implicit none

        integer, intent(in) :: rowDim, colDim, rowSt, rowEnd, colSt, colEnd, triDiag
        real*8, intent(out) :: Amat(rowDim,colDim)


        integer :: i1, i2, i3, maxRow, maxCol, numRows, numCols, row, row1, col1, col3
        real*8 :: theta, st, ct, a1, a2, diagTerm

        numRows = rowEnd - rowSt + 1
		numCols = colEnd - colSt + 1
        do i1 = 1, numCols
		    col1 = i1 + colSt - 1
			row1 = i1 + rowSt - 1
            if(triDiag .eq. 1) then  !! Matrix is tri-diagonal
                maxRow = i1 + 1
                maxCol = i1 + 2
                if(maxRow .gt. numRows) then
                   maxRow = numRows
                endif
                if(maxCol .gt. numCols) then
                   maxCol = numCols
                endif
            elseif(triDiag .eq. 2) then  !! Matrix is semi-upper triangular
                maxRow = i1 + 1
                maxCol = numCols
                if(maxRow .gt. numRows) then
                   maxRow = numRows
                endif
                if(maxCol .gt. numCols) then
                   maxCol = numCols
                endif
            else
                maxRow = numRows
                maxCol = numCols
            endif
            do i2 = i1+1, maxRow
			    row = i2 + rowSt - 1
			    if(abs(Amat(row,col1)) .gt. 1e-10) then
				    diagTerm = Amat(row1,col1)
                    if(diagTerm .eq. r_0) then
                        theta = 1.570796326794897
                    else
                        theta = atan(Amat(row,col1)/diagTerm)
                    endif
                    st = sin(theta)
                    ct = cos(theta)
                    do i3 = i1, maxCol
					    col3 = i3 + colSt - 1
                        a1 = ct*Amat(row1,col3) + st*Amat(row,col3)
                        a2 = -st*Amat(row1,col3) + ct*Amat(row,col3)
                        Amat(row1,col3) = a1
                        Amat(row,col3) = a2
                    enddo
                    Amat(row,col1) = theta
				endif
            enddo
        enddo

        return
    end subroutine rFactorMat

    subroutine solveRFactor(soln,Amat,bVec,rowDim,colDim,rowSt,rowEnd,colSt,colEnd,triDiag)
        implicit none

        integer, intent(in) :: rowDim, colDim, rowSt, rowEnd, colSt, colEnd, triDiag
        real*8, intent(in) :: Amat(rowDim,colDim)
        real*8, intent(out) :: soln(colDim), bVec(rowDim)

        integer :: i1, i2, i3, maxRow, maxCol, numRows, numCols, row, row1, col, col1
        real*8 :: r1, r2, st, ct

        numRows = rowEnd - rowSt + 1
		numCols = colEnd - colSt + 1
        do i1 = 1, numCols
		    col = i1 + colSt - 1
			row1 = i1 + rowSt - 1
            if(triDiag .ne. 0) then
                maxRow = i1 + 1
                if(maxRow .gt. numRows) then
                    maxRow = numRows
                endif
            else
                maxRow = numRows
            endif
            do i2 = i1+1, maxRow
			    row = i2 + rowSt - 1
			    if(abs(Amat(row,col)) .gt. 1e-10) then
                    st = sin(Amat(row,col))
                    ct = cos(Amat(row,col))
                    r1 = ct*bVec(row1) + st*bVec(row)
                    r2 = -st*bVec(row1) + ct*bVec(row)
                    bVec(row1) = r1
                    bVec(row) = r2
				endif
            enddo
        enddo

        do i1 = numCols, 1, -1
		    row = i1 + rowSt - 1
			col1 = i1 + colSt - 1
            if(triDiag .eq. 1) then
                maxCol = i1 + 2
                if(maxCol .gt. numCols) then
                    maxCol = numCols
                endif
            else
                maxCol = numCols
            endif
            soln(row) = bVec(row)
            do i2 = i1+1, maxCol
			    col = i2 + colSt - 1
                soln(row) = soln(row) - Amat(row,col)*soln(col)
            enddo
            soln(row) = soln(row)/Amat(row,col1)
        enddo

        return
    end subroutine solveRFactor

    subroutine symFactor(Amat,hMat,ADim)
        implicit none
        
        integer, intent(in) :: ADim
        real*8, intent(out) :: Amat(ADim,ADim), hMat(ADim,ADim)

        integer :: i1, i2, i3, i4
        real*8 :: theta, st, ct, a1, a2

        hMat(:,:) = r_0
        do i1 = 1, ADim
            hMat(i1,i1) = r_1
        enddo

        do i1 = 1, Adim
            do i3 = i1+2, Adim
                if(Amat(i1+1,i1) .eq. r_0) then
                    theta = 1.570796326794897
                else
                    theta = atan(Amat(i3,i1)/Amat(i1+1,i1))
                endif
                st = sin(theta)
                ct = cos(theta)
                do i4 = i1, ADim
                    a1 = ct*Amat(i1+1,i4) + st*Amat(i3,i4)
                    a2 = -st*Amat(i1+1,i4) + ct*Amat(i3,i4)
                    Amat(i1+1,i4) = a1
                    Amat(i3,i4) = a2
                enddo
                do i4 = 1, ADim
                    a1 = ct*Amat(i4,i1+1) + st*Amat(i4,i3)
                    a2 = -st*Amat(i4,i1+1) + ct*Amat(i4,i3)
                    Amat(i4,i1+1) = a1
                    Amat(i4,i3) = a2
                enddo
                do i4 = 1, ADim
                    a1 = ct*hMat(i4,i1+1) + st*hMat(i4,i3)
                    a2 = -st*hMat(i4,i1+1) + ct*hMat(i4,i3)
                    hMat(i4,i1+1) = a1
                    hMat(i4,i3) = a2
                enddo
            enddo
        enddo

        return
    end subroutine symFactor

    subroutine detMat(fun,dFun,AMat,dAmat,ADim,EVals,numEvals,lam,triDiag)
        implicit none

        integer, intent(in) :: ADim, triDiag, numEVals
        real*8, intent(in) :: EVals(ADim), lam
        real*8, intent(out) :: AMat(ADim,ADim), dAmat(ADim,ADim)
        real*8, intent(out) :: fun, dFun

        integer :: i1, i2, i3, maxRow, maxCol, detExp, ddetExp, rtPExp, drtPExp
        real*8 :: theta, st, ct, a1, a2, det, rtP, maxMag, minMag
        real*8 :: dtheta, dst, dct, da1, da2, ddet, drtP

        dAmat(:,:) = r_0
        do i1 = 1, ADim
            dAmat(i1,i1) = -r_1
        enddo

        do i1 = 1, ADim - 1
            if(triDiag .eq. 1) then
                maxRow = i1 + 1
                maxCol = i1 + 2
                if(maxCol .gt. ADim) then
                   maxCol = ADim
                endif
            elseif(triDiag .eq. 2) then
                maxRow = i1 + 1
                maxCol = ADim
            else
                maxRow = ADim
                maxCol = ADim
            endif
            do i2 = i1+1, maxRow
                if(Amat(i1,i1) .eq. r_0) then
                    theta = 1.570796326794897
                    dtheta = r_0
                else
                    a1 = Amat(i2,i1)/Amat(i1,i1)
                    da1 = dAmat(i2,i1)/Amat(i1,i1) - Amat(i2,i1)*dAmat(i1,i1)/(Amat(i1,i1)*Amat(i1,i1))
                    theta = atan(a1)
                    dtheta = (r_1/(r_1 + a1*a1))*da1
                endif
                st = sin(theta)
                dst = cos(theta)*dtheta
                ct = cos(theta)
                dct = -sin(theta)*dtheta
                do i3 = i1, maxCol
                    a1 = ct*Amat(i1,i3) + st*Amat(i2,i3)
                    da1 = dct*Amat(i1,i3) + ct*dAmat(i1,i3) + dst*Amat(i2,i3) + st*dAmat(i2,i3)
                    a2 = -st*Amat(i1,i3) + ct*Amat(i2,i3)
                    da2 = -dst*Amat(i1,i3) - st*dAmat(i1,i3) + dct*Amat(i2,i3) + ct*dAmat(i2,i3)
                    Amat(i1,i3) = a1
                    dAmat(i1,i3) = da1
                    Amat(i2,i3) = a2
                    dAmat(i2,i3) = da2
                enddo
                Amat(i2,i1) = theta
            enddo
        enddo

        maxMag = 10000000000d0*r_1
        minMag = 0.0000000001d0*r_1

        det = r_1
        ddet = r_0
        detExp = 0
        do i1 = 1, ADim
            a1 = det
            da1 = ddet
            det = a1*Amat(i1,i1)
            ddet = da1*Amat(i1,i1) + a1*dAmat(i1,i1)
            do while(abs(det) .gt. abs(maxMag))
                det = det*minMag
                ddet = ddet*minMag
                detExp = detExp + 10
            enddo
            do while(abs(det) .lt. abs(minMag) .and. abs(det) .ne. r_0)
                det = det*maxMag
                ddet = ddet*maxMag
                detExp = detExp - 10
            enddo
        enddo

        rtP = r_1
        drtP = r_0
        rtPExp = 0
        do i1 = 1, numEvals
            a1 = rtP
            da1 = drtP
            rtP = a1*(eVals(i1) - lam)
            drtP = da1*(eVals(i1) - lam) - a1
            do while(abs(rtP) .gt. abs(maxMag))
                rtP = rtP*minMag
                drtP = drtP*minMag
                rtPExp = rtPExp + 10
            enddo
            do while(abs(rtP) .lt. abs(minMag) .and. abs(rtP) .ne. r_0)
                rtP = rtP*maxMag
                drtP = drtP*maxMag
                rtPExp = rtPExp - 10
            enddo
        enddo

        fun = det/rtP
        dFun = ddet/rtP - det*drtP/(rtP*rtP)

        return
    end subroutine detMat

    subroutine h2AhEvals(eVals,numEval,h2Ah,mDim,lam0,dL0,Ltol,triDiag)
        implicit none

        integer, intent(in) :: numEval, mDim, triDiag
        real*8, intent(in) :: lam0, dL0, Ltol
        real*8, intent(in) :: h2Ah(mDim,mDim)
        real*8, intent(out) :: eVals(numEval)

        integer :: i1, i2, i3, det1Exp, det2Exp, rpExp
        real*8 :: lam1, lam2, det1, det2, dL, rootProd
        real*8, allocatable :: matCopy(:,:), dMat(:,:)

        allocate(matCopy(mDim,mDim),dMat(mDim,mDim))

        lam1 = lam0
        matCopy(:,:) = h2Ah(:,:)
        do i1 = 1, mDim
            matCopy(i1,i1) = matCopy(i1,i1) - lam1
        enddo
        call detMat(det1,det2,matCopy,dMat,mDim,eVals,0,lam1,triDiag)
        
        dL = dL0

        do i2 = 1, numEval
            i3 = 1
            do while(abs(dL) .gt. abs(Ltol) .and. i3 .lt. 100*mDim)
                dL = -det1/det2
                lam1 = lam1 + dL
                matCopy(:,:) = h2Ah(:,:)
                do i1 = 1, mDim
                    matCopy(i1,i1) = matCopy(i1,i1) - lam1
                enddo
                call detMat(det1,det2,matCopy,dMat,mDim,eVals,i2-1,lam1,triDiag)
                i3 = i3 + 1
            enddo
            if(i3 .eq. 100*mDim) then
                write(*,*) 'Eval ', i2, 'Not converged'
            endif
            eVals(i2) = lam1

            lam1 = lam1 - dL0
            matCopy(:,:) = h2Ah(:,:)
            do i1 = 1, mDim
                matCopy(i1,i1) = matCopy(i1,i1) - lam1
            enddo
            call detMat(det1,det2,matCopy,dMat,mDim,eVals,i2,lam1,triDiag)
            dL = dL0
        enddo

        deallocate(matCopy,dMat)
        
        return
    end subroutine h2AhEvals

    subroutine eigenSolve(eVals,eVecs,Amat,ADim,numEVals,triDiag)
        implicit none

        integer, intent(in) :: ADim, numEVals, triDiag
        real*8, intent(out) :: Amat(ADim,ADim), eVals(numEVals), eVecs(Adim,numEVals)

        integer :: i1, i2, i3
        real*8 :: matNorm, pen, resNorm, vMag, rand
        real*8, allocatable :: h2Ah(:,:), rhs(:), pDiag(:), vNext(:), h2AhVecs(:,:), matCopy(:,:), hMat(:,:)

        allocate(h2Ah(ADim,ADim),matCopy(ADim,ADim),rhs(ADim),pDiag(ADim),vNext(ADim),h2AhVecs(ADim,numEVals))
        allocate(hMat(ADim,ADim))

        write(*,*) 'Starting E-solve'

        h2Ah = Amat
		if(triDiag .eq. 0) then
            call symFactor(h2Ah,hMat,ADim)
		else
		    hMat(:,:) = 0d0
			do i1 = 1, ADim
			    hMat(i1,i1) = 1d0
			enddo
		endif

        matNorm = r_0
        pDiag(:) = r_0
        do i1 = 1, ADim
            do i2 = 1, ADim
                matNorm = matNorm + h2Ah(i2,i1)*h2Ah(i2,i1)
                pDiag(i1) = pDiag(i1) + abs(h2Ah(i2,i1))
            enddo
            pDiag(i1) = 0.0001d0*pDiag(i1)
        enddo
        matNorm = sqrt(matNorm)

        call h2AhEvals(eVals,numEvals,h2Ah,ADim,-matNorm,1e-4*matNorm,1e-12*matNorm,triDiag)

        write(*,*) 'got E vals'

        do i1 = 1, numEVals
            matCopy(:,:) = h2Ah(:,:)
            do i2 = 1, ADim
                matCopy(i2,i2) = matCopy(i2,i2) - eVals(i1) + pDiag(i2)
            enddo
            call rFactorMat(matCopy,ADim,ADim,1,ADim,1,ADim,triDiag)
            vMag = r_0
            do i2 = 1, ADim
                h2AhVecs(i2,i1) = sin(r_1*i1*i2)
                vMag = vMag + h2AhVecs(i2,i1)*h2AhVecs(i2,i1)
            enddo
            vMag = sqrt(vMag)
            h2AhVecs(:,i1) = (1d0/vMag)*h2AhVecs(:,i1)
            do i2 = 1, ADim
                rhs(i2) = pDiag(i2)*h2AhVecs(i2,i1)
            enddo
            resNorm = r_1
            i2 = 0
            do while(abs(resNorm) .gt. 1e-10 .and. i2 .lt. 50)
                call solveRFactor(vNext,matCopy,rhs,ADim,ADim,1,ADim,1,ADim,triDiag)
                vMag = r_0
                do i3 = 1, ADim
                    vMag = vMag + vNext(i3)*vNext(i3)
                enddo
                vMag = sqrt(vMag)
                vNext = (r_1/vMag)*vNext
                resNorm = r_0
                do i3 = 1, ADim
                    resNorm = resNorm + abs(vNext(i3) - h2AhVecs(i3,i1))
                    rhs(i3) = pDiag(i3)*vNext(i3)
                enddo
                resNorm = resNorm/ADim
                h2AhVecs(:,i1) = vNext
                i2 = i2 + 1
            enddo
        enddo

        eVecs(:,:) = r_0

        do i1 = 1, ADim
            do i2 = 1, numEVals
                do i3 = 1, ADim
                    eVecs(i1,i2) = eVecs(i1,i2) + hMat(i1,i3)*h2AhVecs(i3,i2)
                enddo
            enddo
        enddo

        deallocate(h2Ah,matCopy,rhs,pDiag,vNext,h2AhVecs)
        deallocate(hMat)

        return
    end subroutine eigenSolve
	
	subroutine getEigModesLDL(eVecs,eVals,numVecs,matDim,ALT,ALTRange,ALTSize,Mmat)
	    implicit none
		
		integer, intent(in) :: numVecs, matDim, ALTSize
		integer, intent(in) :: ALTRange(0:matDim)
		real*8, intent(in) :: ALT(ALTSize), Mmat(matDim)  !! Mmat = square root of diagonal mass matrix
		real*8, intent(out) :: eVecs(matDim,numVecs), eVals(numVecs)
		
		real*8, allocatable :: basisVecs(:,:), tempV1(:), tempV2(:), tempV3(:), hAh(:,:), hAhVecs(:,:), hAhVals(:)
		real*8 :: mag, dp, tVal
		integer :: i1, i2, i3, i4, i5, i6, basisSize
		
		basisSize = 5*numVecs
		allocate(basisVecs(matDim,basisSize+1))
		allocate(tempV1(matDim))
		allocate(tempV2(matDim))
		allocate(tempV3(matDim))
		allocate(hAh(basisSize,basisSize))
		allocate(hAhVecs(basisSize,basisSize))
		allocate(hAhVals(basisSize))
		
		basisVecs(:,:) = r_0
		hAh(:,:) = r_0

		mag = r_0
		do i1 = 1, matDim
		    tVal = sin(r_1*i1)
		    basisVecs(i1,1) = tVal
			mag = mag + tVal*tVal
		enddo
		mag = sqrt(mag)
		basisVecs(:,1) = (r_1/mag)*basisVecs(:,1)
		
		do i1 = 2, basisSize + 1
		    do i2 = 1, matDim
			    tempV1(i2) = Mmat(i2)*basisVecs(i2,i1-1)
			enddo
		    call solveSparseLDLFact(tempV3, tempV1, ALT, ALTSize, ALTRange, matDim, tempV2)
			do i2 = 1, matDim
			    basisVecs(i2,i1) = Mmat(i2)*tempV3(i2)
			enddo
			do i2 = 1, i1-1
			    dp = r_0
				do i3 = 1, matDim
					dp = dp + basisVecs(i3,i1)*basisVecs(i3,i2)
				enddo
				basisVecs(:,i1) = basisVecs(:,i1) - dp*basisVecs(:,i2)
				hAh(i2,i1-1) = dp
			enddo
			mag = r_0
			do i3 = 1, matDim
				mag = mag + basisVecs(i3,i1)*basisVecs(i3,i1)
			enddo
			mag = sqrt(mag)
			basisVecs(:,i1) = (r_1/mag)*basisVecs(:,i1)
			if(i1 .le. basisSize) then
			    hAh(i1,i1-1) = mag
			endif
		enddo
		
		call eigenSolve(hAhVals,hAhVecs,hAh,basisSize,basisSize,1)
		
		i1 = 0
		i2 = basisSize
		do while(i1 .lt. numVecs .and. i2 .ge. 1)
		    tempV2(:) = r_0
			do i3 = 1, matDim
			    do i4 = 1, basisSize
				    tempV2(i3) = tempV2(i3) + basisVecs(i3,i4)*hAhVecs(i4,i2) 
				enddo
			enddo
			mag = r_0
			do i3 = 1, matDim
			    mag = mag + tempV2(i3)*tempV2(i3)
			enddo
			mag = sqrt(mag)
			tempV2(:) = (r_1/mag)*tempV2(:)
			eVecs(:,i1+1) = tempV2(:)
			do i3 = 1, matDim
			    tempV1(i3) = Mmat(i3)*tempV2(i3)
			enddo
		    call solveSparseLDLFact(tempV3, tempV1, ALT, ALTSize, ALTRange, matDim, tempV2)
			do i3 = 1, matDim
			    tempV3(i3) = Mmat(i3)*tempV3(i3)
			enddo
			dp = r_0
			mag = r_0
			do i3 = 1, matDim
			    dp = dp + tempV3(i3)*eVecs(i3,i1+1)
				mag = mag + tempV3(i3)*tempV3(i3)
			enddo
			mag = sqrt(mag)
			if(abs(dp) .gt. 0.999999d0*mag) then
			    i1 = i1 + 1
				eVals(i1) = r_1/dp
			endif
			i2 = i2 - 1
		enddo
		
		do i3 = 1, i1
		    mag = 0d0
			do i4 = 1, matDim
			    tVal = eVecs(i4,i3)/Mmat(i4)
		        eVecs(i4,i3) = tVal
				mag = mag + tVal*tVal
			enddo
			mag = sqrt(mag)
			eVecs(:,i3) = (1d0/mag)*eVecs(:,i3)
		enddo
		
		deallocate(basisVecs)
		deallocate(tempV1)
		deallocate(tempV2)
		deallocate(tempV3)
		deallocate(hAh)
		deallocate(hAhVecs)
		deallocate(hAhVals)
		
	end subroutine getEigModesLDL
	
	subroutine getSparseEigenModes(eVecs,eVals,numVecs,matDim,Amat,Acols,ARange,Asize,Mmat,Cmat,Ccols,CRange,Cdim,Csize)
	    implicit none
		
		integer, intent(in) :: numVecs, matDim, Asize, Cdim, Csize
		integer, intent(in) :: Acols(ASize), ARange(0:matdim), Ccols(Csize), CRange(0:Cdim)
		real*8, intent(in) :: Amat(Asize), Mmat(matDim), Cmat(Csize)
		real*8, intent(out) :: eVecs(matDim,numVecs), eVals(numVecs)
		
		real*8, allocatable :: searchVecs(:,:), basisMags(:), tempV1(:), tempV2(:), hAh(:,:), hAhVecs(:,:), hAhVals(:)
		
		real*8 :: mag, dp, tVal
		integer :: i1, i2, i3, i4, i5, i6, brkLp, maxOutIt, basisSize, searchSize, foundLower
		
		basisSize = 3*numVecs
		searchSize = 3*basisSize
		allocate(basisMags(basisSize))
		allocate(searchVecs(matDim,searchSize))
		allocate(tempV1(matDim))
		allocate(tempV2(matDim))
		allocate(hAh(basisSize,basisSize))
		allocate(hAhVecs(basisSize,basisSize))
		allocate(hAhVals(basisSize))
		basisMags(:) = r_0
		searchVecs(:,:) = r_0
		
		mag = r_0
		do i1 = 1, matDim
		    tVal = sin(r_1*i1)
		    searchVecs(i1,basisSize) = tVal
			mag = mag + tVal*tVal
		enddo
		mag = sqrt(mag)
		searchVecs(:,basisSize) = (r_1/mag)*searchVecs(:,basisSize)
		
		foundLower = 1
		maxOutIt = matDim/searchSize
		i5 = 0
		do while(foundLower .eq. 1 .and. i5 .le. maxOutIt)
		    foundLower = 0
		    do i1 = basisSize + 1, searchSize
			    tempV1(:) = r_0
				do i2 = 1, Cdim
				    do i3 = CRange(i2-1) + 1, CRange(i2)
					    i4 = Ccols(i3)
						tempV1(i2) = tempV1(i2) + Cmat(i3)*searchVecs(i4,i1-1)
					enddo
				enddo
				tempV2(:) = r_0
				do i2 = 1, Cdim
				    do i3 = CRange(i2-1) + 1, CRange(i2)
					    i4 = Ccols(i3)
						tempV2(i4) = tempV2(i4) + Cmat(i3)*tempV1(i2)
					enddo
				enddo
				do i2 = 1, matDim
				    do i3 = ARange(i2-1)+1, ARange(i2)
					    i4 = Acols(i3)
						tempV2(i2) = tempV2(i2) + Amat(i3)*searchVecs(i4,i1-1)
					enddo
				enddo
				do i2 = 1, matDim
				    tempV2(i2) = Mmat(i2)*tempV2(i2)
				enddo
				
				mag = r_0
				do i2 = 1, matDim
				    mag = mag + tempV2(i2)*tempV2(i2)
				enddo
				mag = sqrt(mag)
				
				i2 = basisSize
				brkLp = 0
				do while(i2 .gt. 0 .and. brkLp .eq. 0)
				    if(basisMags(i2) .eq. r_0) then
					    basisMags(i2) = mag
						searchVecs(:,i2) = searchVecs(:,i1-1)
						brkLp = 1
						foundLower = 1
					elseif(mag .lt. basisMags(i2)) then
					    do i3 = 1, i2-1
						    basisMags(i3) = basisMags(i3+1)
						    searchVecs(:,i3) = searchVecs(:,i3+1)
						enddo
						basisMags(i2) = mag
						searchVecs(:,i2) = searchVecs(:,i1-1)
						brkLp = 1
						foundLower = 1
					endif
				    i2 = i2 - 1
				enddo
				
				do i2 = 1, i1-1
				    dp = r_0
					do i3 = 1, matDim
					    dp = dp + tempV2(i3)*searchVecs(i3,i2)
					enddo
					tempV2(:) = tempV2(:) - dp*searchVecs(:,i2)
				enddo
				
				mag = r_0
				do i2 = 1, matDim
				    mag = mag + tempV2(i2)*tempV2(i2)
				enddo
				mag = sqrt(mag)
				
				searchVecs(:,i1) = (r_1/mag)*tempV2(:)
			enddo
			i5 = i5 + 1
		enddo
		
		searchVecs(:,basisSize+1:basisSize+basisSize) = r_0
		
		do i1 = 1, basisSize
			tempV1(:) = r_0
			do i2 = 1, Cdim
				do i3 = CRange(i2-1) + 1, CRange(i2)
					i4 = Ccols(i3)
					tempV1(i2) = tempV1(i2) + Cmat(i3)*searchVecs(i4,i1)
				enddo
			enddo
			tempV2(:) = r_0
			do i2 = 1, Cdim
				do i3 = CRange(i2-1) + 1, CRange(i2)
					i4 = Ccols(i3)
					tempV2(i4) = tempV2(i4) + Cmat(i3)*tempV1(i2)
				enddo
			enddo
			do i2 = 1, matDim
				do i3 = ARange(i2-1)+1, ARange(i2)
					i4 = Acols(i3)
					tempV2(i2) = tempV2(i2) + Amat(i3)*searchVecs(i4,i1)
				enddo
			enddo
			do i2 = 1, matDim
				searchVecs(i2,i1+basisSize) = Mmat(i2)*tempV2(i2)
			enddo		    
		enddo
		
		hAh(:,:) = r_0
		
		do i1 = 1, basisSize
		    do i2 = 1, basisSize
			    i4 = i2 + basisSize
			    do i3 = 1, matDim
				    hAh(i1,i2) = hAh(i1,i2) + searchVecs(i3,i1)*searchVecs(i3,i4)
				enddo
			enddo
		enddo
		
		call eigenSolve(hAhVals,hAhVecs,hAh,basisSize,basisSize,0)
		
		i1 = 0
		i5 = 0
		do while(i1 .lt. numVecs .and. i5 .lt. basisSize)
		    i5 = i5 + 1
		    tempV1(:) = r_0
			do i2 = 1, matDim
			    do i3 = 1, basisSize
				    tempV1(i2) = tempV1(i2) + searchVecs(i2,i3)*hAhVecs(i3,i5)
				enddo
			enddo
			
			mag = r_0
			do i2 = 1, matDim
			    mag = tempV1(i2)*tempV1(i2)
			enddo
			mag = sqrt(mag)
			
			eVecs(:,i1+1) = (r_1/mag)*tempV1(:)
			
			tempV1(:) = r_0
			do i2 = 1, Cdim
				do i3 = CRange(i2-1) + 1, CRange(i2)
					i4 = Ccols(i3)
					tempV1(i2) = tempV1(i2) + Cmat(i3)*eVecs(i4,i1+1)
				enddo
			enddo
			tempV2(:) = r_0
			do i2 = 1, Cdim
				do i3 = CRange(i2-1) + 1, CRange(i2)
					i4 = Ccols(i3)
					tempV2(i4) = tempV2(i4) + Cmat(i3)*tempV1(i2)
				enddo
			enddo
			do i2 = 1, matDim
				do i3 = ARange(i2-1)+1, ARange(i2)
					i4 = Acols(i3)
					tempV2(i2) = tempV2(i2) + Amat(i3)*eVecs(i4,i1+1)
				enddo
			enddo
			do i2 = 1, matDim
				tempV2(i2) = Mmat(i2)*tempV2(i2)
			enddo
			
			dp = r_0
			mag = r_0
			do i2 = 1, matDim
			    dp = dp + tempV2(i2)*eVecs(i2,i1+1)
				mag = mag + tempV2(i2)+tempV2(i2)
			enddo
			mag = sqrt(mag)
			
			if(abs(dp) .ge. 0.9999d0*mag) then
			    eVals(i1+1) = dp
				i1 = i1 + 1
			endif
		enddo
		
		deallocate(basisMags)
		deallocate(searchVecs)
		deallocate(tempV1)
		deallocate(tempV2)
		deallocate(hAh)
		deallocate(hAhVecs)
		deallocate(hAhVals)
		
	end subroutine getSparseEigenModes

end module ASenD_solvers
