module ASenD_bookKeeping
    use ASenD_globalData
	use ASenD_constantVals
	use ASenD_r_elementEqns

    contains
	
	subroutine getdToElComp()
        implicit none

        integer :: i1, i2, i3, i4, i5, i6, i7, i8, i9, di, nd, inserted

        integer, allocatable :: fullMap(:,:), dElSeti(:), dNdSeti(:)
		
		!! Build d to el and d to nd tables
		allocate(dElSeti(numDVar))
		allocate(dNdSeti(numDVar))
		allocate(dToElRange(0:numDVar))
		allocate(dToNdRange(0:numDVar))
		dElSeti(:) = 0
		dNdSeti(:) = 0
		dToElRange(:) = 0
		dToNdRange(:) = 0
		
		do i1 = 1, numDVar
		    do i2 = 1, numElSets
			    if(elSetName(i2) .eq. dElSet(i1)) then
				    i3 = elSetRange(i2) - elSetRange(i2-1)
					dToElRange(i1) = i3
					dElSeti(i1) = i2
				endif
			enddo
			do i2 = 1, numNdSets
			    if(ndSetName(i2) .eq. dNdSet(i1)) then
				    i3 = ndSetRange(i2) - ndSetRange(i2-1)
					dToNdRange(i1) = i3
					dNdSeti(i1) = i2
				endif
			enddo
		enddo
		
		do i1 = 2, numDVar
		    dToElRange(i1) = dToElRange(i1-1) + dToElRange(i1)
			dToNdRange(i1) = dToNdRange(i1-1) + dToNdRange(i1)
		enddo
		
		elToDSize = dToElRange(numDVar)
		ndToDSize = dToNdRange(numDVar)
		
		allocate(dToEl(elToDSize))
		allocate(dToElCoef(elToDSize))
		allocate(elToD(elToDSize))
		allocate(elToCoef(elToDSize))
		allocate(dToNd(ndToDSize))
		allocate(dToNdCoef(ndToDSize))
		allocate(ndToD(ndToDSize))
		allocate(ndToCoef(ndToDSize))
		
		do i1 = 1, numDVar
		    i2 = dCoefListRange(i1) - dCoefListRange(i1-1)
			if(i2 .eq. 0) then
			    if(dElSeti(i1) .ne. 0) then
				    i3 = dElSeti(i1)
					i4 = elSetRange(i3-1) + 1
					i5 = elSetRange(i3)
					i6 = dToElRange(i1-1) + 1
					i7 = dToElRange(i1)
				    dToEl(i6:i7) = elementSets(i4:i5)
					dToElCoef(i6:i7) = r_1
				else
				    i3 = dNdSeti(i1)
					i4 = ndSetRange(i3-1) + 1
					i5 = ndSetRange(i3)
					i6 = dToNdRange(i1-1) + 1
					i7 = dToNdRange(i1)
				    dToNd(i6:i7) = nodeSets(i4:i5)
					dToNdCoef(i6:i7) = r_1
				endif
			elseif(i2 .eq. 1) then
			    i8 = dCoefListRange(i1)
			    if(dElSeti(i1) .ne. 0) then
				    i3 = dElSeti(i1)
					i4 = elSetRange(i3-1) + 1
					i5 = elSetRange(i3)
					i6 = dToElRange(i1-1) + 1
					i7 = dToElRange(i1)
				    dToEl(i6:i7) = elementSets(i4:i5)
					dToElCoef(i6:i7) = dCoefList(i8)
				else
				    i3 = dNdSeti(i1)
					i4 = ndSetRange(i3-1) + 1
					i5 = ndSetRange(i3)
					i6 = dToNdRange(i1-1) + 1
					i7 = dToNdRange(i1)
				    dToNd(i6:i7) = nodeSets(i4:i5)
					dToNdCoef(i6:i7) = dCoefList(i8)
				endif
			else
			    i8 = dCoefListRange(i1-1) + 1
				i9 = dCoefListRange(i1)
			    if(dElSeti(i1) .ne. 0) then
				    i3 = dElSeti(i1)
					i4 = elSetRange(i3-1) + 1
					i5 = elSetRange(i3)
					i6 = dToElRange(i1-1) + 1
					i7 = dToElRange(i1)
				    dToEl(i6:i7) = elementSets(i4:i5)
					dToElCoef(i6:i7) = dCoefList(i8:i9)
				else
				    i3 = dNdSeti(i1)
					i4 = ndSetRange(i3-1) + 1
					i5 = ndSetRange(i3)
					i6 = dToNdRange(i1-1) + 1
					i7 = dToNdRange(i1)
				    dToNd(i6:i7) = nodeSets(i4:i5)
					dToNdCoef(i6:i7) = dCoefList(i8:i9)
				endif
			endif
		enddo
		
		!! Populate elToD from dToEl
		elToD(:) = 0
		elToDRange(:) = 0
		ndToD(:) = 0
		ndToDRange(:) = 0
        	
		do i1 = 1, elToDSize
		    i2 = dToEl(i1)
			if(i2 .gt. 0) then
			    elToDRange(i2) = elToDRange(i2) + 1
			endif
		enddo
		
		do i1 = 2, numEls
		    elToDRange(i1) = elToDRange(i1) + elToDRange(i1-1)
		enddo
		
		do i1 = 1, numDVar
		    do i2 = dToElRange(i1-1)+1, dToElRange(i1)
			    i3 = dToEl(i2)
			    if(i3 .gt. 0) then
					i4 = elToDRange(i3-1) + 1
					inserted = 0
					do while(inserted .eq. 0 .and. i4 .le. elToDRange(i3))
						if(elToD(i4) .eq. 0) then
							elToD(i4) = i1
							elToCoef(i4) = dToElCoef(i2)
							inserted = 1
						elseif(elToD(i4) .eq. i1) then
							inserted = 1
						endif
						i4 = i4 + 1
					enddo
				endif
			enddo
		enddo
		
		
		!! Populate ndToD from dToNd
		
		do i1 = 1, ndToDSize
		    i2 = dToNd(i1)
			if(i2 .gt. 0) then
			    ndToDRange(i2) = ndToDRange(i2) + 1
			endif
		enddo
		
		do i1 = 2, numNodes
		    ndToDRange(i1) = ndToDRange(i1) + ndToDRange(i1-1)
		enddo
		
		do i1 = 1, numDVar
		    do i2 = dToNdRange(i1-1)+1, dToNdRange(i1)
			    i3 = dToNd(i2)
				if(i3 .gt. 0) then
					i4 = ndToDRange(i3-1) + 1
					inserted = 0
					do while(inserted .eq. 0 .and. i4 .le. ndToDRange(i3))
						if(ndToD(i4) .eq. 0) then
							ndToD(i4) = i1
							ndToCoef(i4) = dToNdCoef(i2)
							inserted = 1
						elseif(elToD(i4) .eq. i1) then
							inserted = 1
						endif
						i4 = i4 + 1
					enddo
				endif
			enddo
		enddo
		
		!! Build comprehensive D to El tables

        allocate(dToElCRange(0:numDVar))
		allocate(elToDCRange(0:numEls))
        dToElCRange(:) = 0

        allocate(fullMap(numEls,100))

        do i1 = 1, numDVar, 100
            fullMap(:,:) = 0
            do i2 = 1, numEls
                do i3 = elToDRange(i2-1) + 1, elToDRange(i2)
                    di = elToD(i3)
                    if(di .ge. i1 .and. di .lt. i1+100) then
                        i4 = mod(di-1,100) + 1
                        fullMap(i2,i4) = 1
                    endif
                enddo
                do i3 = 1, 8
                    nd = elementList(i3,i2)
                    if(nd .ne. 0) then
                        do i4 = ndToDRange(nd-1) + 1, ndToDRange(nd)
                            di = ndToD(i4)
                            if(di .ge. i1 .and. di .lt. i1+100) then
                                i5 = mod(di-1,100) + 1
                                fullMap(i2,i5) = 1
                            endif
                        enddo
                    endif
                enddo
            enddo
            do i2 = 1, 100
                i3 = i2 + i1 - 1
                do i4 = 1, numEls
                    dToElCRange(i3) = dToElCRange(i3) + fullMap(i4,i2)
                enddo
            enddo
        enddo

        do i2 = 2, numDVar
            dToElCRange(i2) = dToElCRange(i2) + dToElCRange(i2-1)
        enddo

        i2 = dToElCRange(numDVar)
        allocate(dToElComp(i2))
		allocate(elToDComp(i2))

        dToElComp(:) = 0

        do i1 = 1, numDVar, 100
            fullMap(:,:) = 0
            do i2 = 1, numEls
                do i3 = elToDRange(i2-1) + 1, elToDRange(i2)
                    di = elToD(i3)
                    if(di .ge. i1 .and. di .lt. i1+100) then
                        i4 = mod(di-1,100) + 1
                        fullMap(i2,i4) = 1
                    endif
                enddo
                do i3 = 1, 8
                    nd = elementList(i3,i2)
                    if(nd .ne. 0) then
                        do i4 = ndToDRange(nd-1) + 1, ndToDRange(nd)
                            di = ndToD(i4)
                            if(di .ge. i1 .and. di .lt. i1+100) then
                                i5 = mod(di-1,100) + 1
                                fullMap(i2,i5) = 1
                            endif
                        enddo
                    endif
                enddo
            enddo
            do i2 = 1, 100
                i3 = i2 + i1 - 1
                i5 = dToElCRange(i3-1) + 1
                do i4 = 1, numEls
                    if(fullMap(i4,i2) .eq. 1) then
                        dToElComp(i5) = i4
                        i5 = i5 + 1
                    endif
                enddo
            enddo
        enddo
		
		elToDComp(:) = 0
		elToDCRange(:) = 0
		do i1 = 1, numDVar
		    do i2 = dToElCRange(i1-1) + 1, dToElCRange(i1)
			    i3 = dToElComp(i2)
				elToDCRange(i3) = elToDCRange(i3) + 1
			enddo
		enddo
		
		do i2 = 2, numEls
            elToDCRange(i2) = elToDCRange(i2) + elToDCRange(i2-1)
        enddo
		
		do i1 = 1, numDVar
		    do i2 = dToElCRange(i1-1) + 1, dToElCRange(i1)
			    i3 = dToElComp(i2)
				inserted = 0
				i4 = elToDCRange(i3-1) + 1
				do while(inserted .eq. 0 .and. i4 .le. elToDCRange(i3))
				    if(elToDComp(i4) .eq. 0) then
					    elToDComp(i4) = i1
						inserted = 1
					elseif(elToDComp(i4) .eq. i1) then
					    inserted = 1
					endif
					i4 = i4 + 1
				enddo
			enddo
		enddo

        deallocate(dElSeti)
		deallocate(dNdSeti)
        deallocate(fullMap)

        return
    end subroutine getdToElComp
	
	subroutine findElementSurfaces()
	    implicit none
		
		integer, allocatable :: faceEl(:), faceNum(:), faceNd1(:), faceNd2(:), faceNd3(:), faceDup(:), nextFace(:)
		integer :: faces(6,6), faceNodes(6)
		integer :: i1, i2, i3, i4, i5, i6, numFaces, lowNd, inserted, swap, nextAvail
		
		i1 = numNodes + numEls*6
		allocate(faceEl(i1))
		allocate(faceNum(i1))
		allocate(faceNd1(i1))
		allocate(faceNd2(i1))
		allocate(faceNd3(i1))
		allocate(faceDup(i1))
		allocate(nextFace(i1))
		
		faceEl(:) = 0
		faceNum(:) = 0
		faceNd1(:) = 0
		faceNd2(:) = 0
		faceNd3(:) = 0
		faceDup(:) = 0
		nextFace(:) = 0
		
		nextAvail = numNodes + 1
		do i1 = 1, numEls
		    call r_getElementFaces(faces,numFaces,i1)
		    do i2 = 1, numFaces
			    i3 = 6
				do while(i3 .gt. 0 .and. faces(i3,i2) .eq. 0)
				    i3 = i3 - 1
				enddo
				faceNodes(:) = 0
				do i4 = 1, i3
				    faceNodes(i4) = elementList(faces(i4,i2),i1)
				enddo
				do i4 = 1, i3
				    do i5 = 1, i3-1
					    if(faceNodes(i5+1) .lt. faceNodes(i5)) then
						    swap = faceNodes(i5)
							faceNodes(i5) = faceNodes(i5+1)
							faceNodes(i5+1) = swap
						endif
					enddo
				enddo
				i3 = faceNodes(1)
				inserted = 0
				do while(inserted .eq. 0)
					if(faceEl(i3) .eq. 0) then
						faceEl(i3) = i1
						faceNum(i3) = i2
						faceNd1(i3) = faceNodes(1)
						faceNd2(i3) = faceNodes(2)
						faceNd3(i3) = faceNodes(3)
						faceDup(i3) = 0
						nextFace(i3) = 0
						inserted = 1
                    elseif(faceNd2(i3) .eq. faceNodes(2) .and. faceNd3(i3) .eq. faceNodes(3)) then
						faceDup(i3) = 1
						inserted = 1
					elseif(nextFace(i3) .eq. 0) then
						nextFace(i3) = nextAvail
						faceEl(nextAvail) = i1
						faceNum(nextAvail) = i2
						faceNd1(nextAvail) = faceNodes(1)
						faceNd2(nextAvail) = faceNodes(2)
						faceNd3(nextAvail) = faceNodes(3)
						faceDup(nextAvail) = 0
						nextFace(nextAvail) = 0
						nextAvail = nextAvail + 1
						inserted = 1
					else
						i3 = nextFace(i3)
					endif						
				enddo
			enddo
		enddo
		
		i1 = numNodes + 6*numEls
		
		elementSurfaces(:,:) = 0
		do i2 = 1, i1
		    if(faceEl(i2) .ne. 0 .and. faceDup(i2) .eq. 0) then
			    elementSurfaces(faceNum(i2),faceEl(i2)) = 1
			endif
		enddo
		
		do i2 = 1, numEls
		    i3 = elementType(i2)
			if(i3 .eq. 41 .or. i3 .eq. 3) then
			    elementSurfaces(1:2,i2) = 1
			endif
		enddo
		
		deallocate(faceEl)
		deallocate(faceNum)
		deallocate(faceNd1)
		deallocate(faceNd2)
		deallocate(faceNd3)
		deallocate(faceDup)
		deallocate(nextFace)
		
	end subroutine findElementSurfaces
	
	subroutine levelSortNodes(nextAvail,stopLength,ndof,nodalConn,nCRange,nCSize,ndInserted,ndBlockNum,currentBlock,numNds)
	    implicit none
		
		integer, intent(in) :: stopLength, ndof, nCSize, numNds
		integer, intent(in) :: nodalConn(nCSize), nCRange(0:numNds)
		integer, intent(out) :: nextAvail, currentBlock
		integer, intent(out) :: ndInserted(numNds), ndBlockNum(numNds)
		
		real*8 :: vec(3), dist, minDist
		integer :: minCon, stNd, currentNd, currentInd, blockSt
		integer :: i1, i2, i3, i4
	
		minCon = 1000
		stNd = 0
		do i1 = 1, numNodes
			if(nDofIndex(1,i1) .eq. ndof) then
				i2 = nCRange(i1) - nCRange(i1-1)
				if(i2 .lt. minCon) then
					stNd = i1
					minCon = i2
				endif
			endif
		enddo
		
		originalRank(nextAvail) = stNd
		currentNd = stNd
		currentInd = nextAvail
		blockSt = nextAvail
		currentBlock = currentBlock + 1
		ndBlockNum(stNd) = currentBlock
		nextAvail = nextAvail + 1
		ndInserted(stNd) = 1
		do while(nextAvail .le. stopLength)
			i1 = nCRange(currentNd-1) + 1
			do while(i1 .le. nCRange(currentNd))
				i2 = nodalConn(i1)
				if(i2 .ne. 0) then
					if(ndInserted(i2) .eq. 0) then
						originalRank(nextAvail) = i2
						nextAvail = nextAvail + 1
						ndInserted(i2) = 1
						ndBlockNum(i2) = currentBlock
						i3 = nextAvail - blockSt
						if(i3 .eq. solverBlockDim .and. nextAvail .le. numNodes) then
							minDist = r_1*1e+100_8
							stNd = 0
							do i4 = 1, numNodes
								if(ndInserted(i4) .eq. 0 .and. nDofIndex(1,i4) .eq. ndof) then
									vec(1) = nodeList(1,i4) - nodeList(1,i2)
									vec(2) = nodeList(2,i4) - nodeList(2,i2)
									vec(3) = nodeList(3,i4) - nodeList(3,i2)
									dist = sqrt(vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3))
									if(dist .lt. minDist) then
										minDist = dist
										stNd = i4
									endif
								endif
							enddo
							originalRank(nextAvail) = stNd
							currentNd = stNd
							currentInd = nextAvail
							ndInserted(stNd) = 1
							blockSt = nextAvail
							currentBlock = currentBlock + 1
							ndBlockNum(stNd) = currentBlock
							nextAvail = nextAvail + 1
							i1 = nCRange(currentNd-1)
						endif
					endif
				endif
				i1 = i1 + 1
			enddo
			currentInd = currentInd + 1
			if(currentInd .eq. nextAvail .and. nextAvail .le. stopLength) then
				i2 = originalRank(nextAvail-1)
				minDist = r_1*1e+100_8
				stNd = 0
				do i4 = 1, numNodes
					if(ndInserted(i4) .eq. 0 .and. nDofIndex(1,i4) .eq. ndof) then
						vec(1) = nodeList(1,i4) - nodeList(1,i2)
						vec(2) = nodeList(2,i4) - nodeList(2,i2)
						vec(3) = nodeList(3,i4) - nodeList(3,i2)
						dist = sqrt(vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3))
						if(dist .lt. minDist) then
							minDist = dist
							stNd = i4
						endif
					endif
				enddo
				originalRank(nextAvail) = stNd
				currentInd = nextAvail
				ndInserted(stNd) = 1
				ndBlockNum(stNd) = currentBlock
				nextAvail = nextAvail + 1
			endif
			currentNd = originalRank(currentInd)
		enddo
	
	end subroutine levelSortNodes
	
	subroutine allocateElasticMat(nodalConn,nCRange,nCSize,ndBlockNum,numNds)
	    implicit none
		
		integer, intent(in) :: nCSize, numNds
		integer, intent(in) :: nodalConn(nCSize), nCRange(0:numNds), ndBlockNum(numNds)
		integer :: i1, i2, i3, i4, i5, i6, i7, i8
		
		allocate(elMatRange(0:elMatDim))
		allocate(elMatLTRange(0:elMatDim))
		
		elMatRange(:) = 0
		elMatLTRange(:) = 0
		
		do i1 = 1, numNodes
		    do i2 = 1, 6
			    i3 = nDofIndex(i2,i1)
			    if(i3 .ne. 0) then
				    do i4 = 1, 6
					    i5 = nDofIndex(i4,i1)
						if(i5 .ne. 0) then
						    elMatRange(i3) = elMatRange(i3) + 1
						endif
					enddo
					do i4 = nCRange(i1-1)+1, nCRange(i1)
					    i5 = nodalConn(i4)
						if(i5 .ne. 0) then
						    do i6 = 1, 6
							    i7 = nDofIndex(i6,i5)
								if(i7 .ne. 0) then
								    elMatRange(i3) = elMatRange(i3) + 1
								endif
							enddo
						endif
					enddo
				endif
			enddo
		enddo
		
		do i1 = 2, elMatDim
		    elMatRange(i1) = elMatRange(i1) + elMatRange(i1-1)
		enddo
		
		elMatSize = elMatRange(elMatDim)
		
		allocate(elasticMat(elMatSize))
		allocate(elMatCols(elMatSize))
		
		do i1 = 1, numNodes
		    do i2 = 1, 6
			    i3 = nDofIndex(i2,i1)
			    if(i3 .ne. 0) then
				    elMatLTRange(i3) = 1
					do i4 = 1, 6
					    i5 = nDofIndex(i4,i1)
						if(i5 .gt. 0 .and. i5 .lt. i3) then
						    i6 = i3 - i5 + 1
							if(i6 .gt. elMatLTRange(i3) .and. i6 .le. solverMaxBW) then
							    elMatLTRange(i3) = i6
							endif
						endif
					enddo
					do i4 = nCRange(i1-1)+1, nCRange(i1)
					    i5 = nodalConn(i4)
						if(i5 .ne. 0) then
						    do i6 = 1, 6
								i7 = nDofIndex(i6,i5)
								if(i7 .ne. 0 .and. i7 .le. i3) then
									i8 = i3 - i7 + 1
									if(i8 .gt. elMatLTRange(i3) .and. i8 .le. solverMaxBW) then
										elMatLTRange(i3) = i8
									endif
								endif
							enddo
						endif
					enddo
				endif
			enddo
		enddo
		
		do i1 = 2, elMatDim
		    elMatLTRange(i1) = elMatLTRange(i1) + elMatLTRange(i1-1)
		enddo
		
		elMatLTSize = elMatLTRange(elMatDim)
		allocate(elMatLT(elMatLTSize))
		
	end subroutine allocateElasticMat
	
	subroutine allocateThermalMat(nodalConn,nCRange,nCSize,ndBlockNum,numNds)
	    implicit none
		
		integer, intent(in) :: nCSize, numNds
		integer, intent(in) :: nodalConn(nCSize), nCRange(0:numNds), ndBlockNum(numNds)
		integer :: i1, i2, i3, i4, i5, i6, i7
		
		allocate(thermMatRange(0:numNds))
		allocate(thermMatLTRange(0:numNds))
		
		thermMatRange(:) = 0
		thermMatLTRange(:) = 0
		
		do i1 = 1, numNodes
		    i3 = currentRank(i1)
			thermMatRange(i3) = thermMatRange(i3) + 1
			do i4 = nCRange(i1-1)+1, nCRange(i1)
				i5 = nodalConn(i4)
				if(i5 .ne. 0) then
					thermMatRange(i3) = thermMatRange(i3) + 1
				endif
			enddo
		enddo
		
		do i1 = 2, numNodes
		    thermMatRange(i1) = thermMatRange(i1) + thermMatRange(i1-1)
		enddo
		
		thermMatSize = thermMatRange(numNodes)
		
		allocate(thermMat(thermMatSize))
		allocate(thermMatCols(thermMatSize))
		
		do i1 = 1, numNodes
			i3 = currentRank(i1)
			if(i3 .ne. 0) then
				thermMatLTRange(i3) = 1
				do i4 = nCRange(i1-1)+1, nCRange(i1)
					i5 = nodalConn(i4)
					if(i5 .ne. 0) then
					    i6 = currentRank(i1) - currentRank(i5)
						if(i6 .ge. 0 .and. i6 .le. solverMaxBW) then !!if(ndBlockNum(i5) .eq. ndBlockNum(i1)) then
							i7 = i6 + 1
							if(i7 .gt. thermMatLTRange(i3)) then
								thermMatLTRange(i3) = i7
							endif
						endif
					endif
				enddo
			endif
		enddo
		
		do i1 = 2, numNodes
		    thermMatLTRange(i1) = thermMatLTRange(i1) + thermMatLTRange(i1-1)
		enddo
		
		thermMatLTSize = thermMatLTRange(numNodes)
		allocate(thermMatLT(thermMatLTSize))
	
	end subroutine allocateThermalMat
	
	subroutine buildElasticMPC()
	    implicit none
		
		integer, allocatable :: eqnBlockNumEqns(:)
		integer, allocatable :: eqnBlockNumTrms(:)
		integer :: i1, i2, i3, i4, i5, i6, i7, i8, inserted
		
		allocate(eqnBlockNumEqns(0:numMPC))
		allocate(eqnBlockNumTrms(numMPC))
		eqnBlockNumEqns(:) = 0
		eqnBlockNumTrms(:) = 0
		
		do i1 = 1, mpcSize
		    if(mpcDof(i1) .le. 6) then
				i2 = mpcEqn(i1)
				eqnBlockNumTrms(i2) = eqnBlockNumTrms(i2) + 1
				read(mpcNode(i1),*,err=362) i3
				i3 = i3 + 1
				if(eqnBlockNumEqns(i2) .lt. 1) then
					eqnBlockNumEqns(i2) = 1
				endif
				goto 370
362 			do i3 = 1, numNdSets
					if(ndSetName(i3) .eq. mpcNode(i1)) then
						i4 = ndSetRange(i3) - ndSetRange(i3-1)
						if(i4 .gt. eqnBlockNumEqns(i2)) then
							eqnBlockNumEqns(i2) = i4
						endif
					endif
				enddo
370 			i3 = 0
            endif
		enddo
		
		elMPCDim = 0
		elMPCSize = 0
		
		do i1 = 1, numMPC
		    elMPCDim = elMPCDim + eqnBlockNumEqns(i1)
			elMPCSize = elMPCSize + eqnBlockNumEqns(i1)*eqnBlockNumTrms(i1)
		enddo
		
		if(elMPCDim .gt. 0) then
			allocate(elMPCMat(elMPCSize))
			allocate(elMPCMatCols(elMPCSize))
			allocate(elMPCMatRange(0:elMPCDim))
			allocate(elMPCRHS(elMPCDim))
			
			elMPCMat(:) = r_0
			elMPCMatCols(:) = 0
			elMPCMatRange(:) = 0
			elMPCRHS(:) = r_0
		endif
		
		do i1 = 2, numMPC
		    eqnBlockNumEqns(i1) = eqnBlockNumEqns(i1) + eqnBlockNumEqns(i1-1)
		enddo
		
		do i1 = 1, numMPC
		    i2 = eqnBlockNumEqns(i1-1) + 1
			i3 = eqnBlockNumEqns(i1)
		    elMPCMatRange(i2:i3) = eqnBlockNumTrms(i1)
		enddo
		
		do i1 = 2, elMPCDim
		    elMPCMatRange(i1) = elMPCMatRange(i1) + elMPCMatRange(i1-1)
		enddo
		
		do i1 = 1, mpcSize
		    if(mpcDof(i1) .le. 6) then
				i2 = mpcEqn(i1)
				i3 = eqnBlockNumEqns(i2-1) + 1
				read(mpcNode(i1),*,err=426) i4
				i4 = i4 + 1
				i5 = nDofIndex(mpcDof(i1),i4)
				do while(i3 .le. eqnBlockNumEqns(i2))
				    i6 = elMPCMatRange(i3-1) + 1
					inserted = 0
					do while(inserted .eq. 0 .and. i6 .le. elMPCMatRange(i3))
					    if(elMPCMatCols(i6) .eq. 0) then
						    elMPCMat(i6) = mpcCoef(i1)
							elMPCMatCols(i6) = i5
							inserted = 1
						elseif(elMPCMatCols(i6) .eq. i5) then
						    inserted = 1
						endif
						i6 = i6 + 1
					enddo
					elMPCRHS(i3) = mpcRHS(i2)
				    i3 = i3 + 1
				enddo
				goto 466
426 			do i4 = 1, numNdSets
					if(ndSetName(i4) .eq. mpcNode(i1)) then
						i5 = ndSetRange(i4) - ndSetRange(i4-1)
						if(i5 .eq. 1) then
							i6 = nDofIndex(mpcDof(i1),nodeSets(ndSetRange(i4)))
							do while(i3 .le. eqnBlockNumEqns(i2))
								i7 = elMPCMatRange(i3-1) + 1
								inserted = 0
								do while(inserted .eq. 0 .and. i7 .le. elMPCMatRange(i3))
									if(elMPCMatCols(i7) .eq. 0) then
										elMPCMat(i7) = mpcCoef(i1)
										elMPCMatCols(i7) = i6
										inserted = 1
									elseif(elMPCMatCols(i7) .eq. i6) then
										inserted = 1
									endif
									i7 = i7 + 1
								enddo
								elMPCRHS(i3) = mpcRHS(i2)
								i3 = i3 + 1
							enddo
                        else
                            do i6 = ndSetRange(i4-1) + 1, ndSetRange(i4)
							    i7 = nDofIndex(mpcDof(i1),nodeSets(i6))
								i8 = elMPCMatRange(i3-1) + 1
								inserted = 0
								do while(inserted .eq. 0 .and. i8 .le. elMPCMatRange(i3))
								    if(elMPCMatCols(i8) .eq. 0) then
									    elMPCMat(i8) = mpcCoef(i1)
										elMPCMatCols(i8) = i7
										inserted = 1
									elseif(elMPCMatCols(i8) .eq. i7) then
									    inserted = 1
									endif
									i8 = i8 + 1
								enddo
								elMPCRHS(i3) = mpcRHS(i2)
								i3 = i3 + 1
                            enddo							
						endif
					endif
				enddo
466 			i3 = 0
            endif
		enddo
		
		deallocate(eqnBlockNumEqns)
		deallocate(eqnBlockNumTrms)
		
	end subroutine buildElasticMPC
	
	subroutine buildThermalMPC()
	    implicit none
		
		integer, allocatable :: eqnBlockNumEqns(:)
		integer, allocatable :: eqnBlockNumTrms(:)
		integer :: i1, i2, i3, i4, i5, i6, i7, i8, inserted
		
		allocate(eqnBlockNumEqns(0:numMPC))
		allocate(eqnBlockNumTrms(numMPC))
		eqnBlockNumEqns(:) = 0
		eqnBlockNumTrms(:) = 0
		
		do i1 = 1, mpcSize
		    if(mpcDof(i1) .eq. 7) then
				i2 = mpcEqn(i1)
				eqnBlockNumTrms(i2) = eqnBlockNumTrms(i2) + 1
				read(mpcNode(i1),*,err=499) i3
				i3 = i3 + 1
				if(eqnBlockNumEqns(i2) .lt. 1) then
					eqnBlockNumEqns(i2) = 1
				endif
				goto 507
499 			do i3 = 1, numNdSets
					if(ndSetName(i3) .eq. mpcNode(i1)) then
						i4 = ndSetRange(i3) - ndSetRange(i3-1)
						if(i4 .gt. eqnBlockNumEqns(i2)) then
							eqnBlockNumEqns(i2) = i4
						endif
					endif
				enddo
507 			i3 = 0
            endif
		enddo
		
		thermMPCDim = 0
		thermMPCSize = 0
		
		do i1 = 1, numMPC
		    thermMPCDim = thermMPCDim + eqnBlockNumEqns(i1)
			thermMPCSize = thermMPCSize + eqnBlockNumEqns(i1)*eqnBlockNumTrms(i1)
		enddo
		
		if(thermMPCDim .gt. 0) then
			allocate(thermMPCMat(thermMPCSize))
			allocate(thermMPCMatCols(thermMPCSize))
			allocate(thermMPCMatRange(0:thermMPCDim))
			allocate(thermMPCRHS(thermMPCDim))
			
			thermMPCMat(:) = r_0
			thermMPCMatCols(:) = 0
			thermMPCMatRange(:) = 0
			thermMPCRHS(:) = r_0
		endif
		
		do i1 = 2, numMPC
		    eqnBlockNumEqns(i1) = eqnBlockNumEqns(i1) + eqnBlockNumEqns(i1-1)
		enddo
		
		do i1 = 1, numMPC
		    i2 = eqnBlockNumEqns(i1-1) + 1
			i3 = eqnBlockNumEqns(i1)
		    thermMPCMatRange(i2:i3) = eqnBlockNumTrms(i1)
		enddo
		
		do i1 = 2, thermMPCDim
		    thermMPCMatRange(i1) = thermMPCMatRange(i1) + thermMPCMatRange(i1-1)
		enddo
		
		do i1 = 1, mpcSize
		    if(mpcDof(i1) .eq. 7) then
				i2 = mpcEqn(i1)
				i3 = eqnBlockNumEqns(i2-1) + 1
				read(mpcNode(i1),*,err=564) i4
				i4 = i4 + 1
				i5 = currentRank(i4)
				do while(i3 .le. eqnBlockNumEqns(i2))
				    i6 = thermMPCMatRange(i3-1) + 1
					inserted = 0
					do while(inserted .eq. 0 .and. i6 .le. thermMPCMatRange(i3))
					    if(thermMPCMatCols(i6) .eq. 0) then
						    thermMPCMat(i6) = mpcCoef(i1)
							thermMPCMatCols(i6) = i5
							inserted = 1
						elseif(thermMPCMatCols(i6) .eq. i5) then
						    inserted = 1
						endif
						i6 = i6 + 1
					enddo
					thermMPCRHS(i3) = mpcRHS(i2)
				    i3 = i3 + 1
				enddo
				goto 606
564 			do i4 = 1, numNdSets
					if(ndSetName(i4) .eq. mpcNode(i1)) then
						i5 = ndSetRange(i4) - ndSetRange(i4-1)
						if(i5 .eq. 1) then
						    i6 = currentRank(nodeSets(ndSetRange(i4)))
							do while(i3 .le. eqnBlockNumEqns(i2))
								i7 = thermMPCMatRange(i3-1) + 1
								inserted = 0
								do while(inserted .eq. 0 .and. i7 .le. thermMPCMatRange(i3))
									if(thermMPCMatCols(i7) .eq. 0) then
										thermMPCMat(i7) = mpcCoef(i1)
										thermMPCMatCols(i7) = i6
										inserted = 1
									elseif(thermMPCMatCols(i7) .eq. i6) then
										inserted = 1
									endif
									i7 = i7 + 1
								enddo
								thermMPCRHS(i3) = mpcRHS(i2)
								i3 = i3 + 1
							enddo
                        else
                            do i6 = ndSetRange(i4-1) + 1, ndSetRange(i4)
							    i7 = currentRank(nodeSets(i6))
								i8 = thermMPCMatRange(i3-1) + 1
								inserted = 0
								do while(inserted .eq. 0 .and. i8 .le. thermMPCMatRange(i3))
								    if(thermMPCMatCols(i8) .eq. 0) then
									    thermMPCMat(i8) = mpcCoef(i1)
										thermMPCMatCols(i8) = i7
										inserted = 1
									elseif(thermMPCMatCols(i8) .eq. i7) then
									    inserted = 1
									endif
									i8 = i8 + 1
								enddo
								thermMPCRHS(i3) = mpcRHS(i2)
								i3 = i3 + 1
                            enddo			
						endif
					endif
				enddo
606 			i3 = 0
            endif
		enddo		
		
		deallocate(eqnBlockNumEqns)
		deallocate(eqnBlockNumTrms)		
		
	end subroutine buildThermalMPC
	
	subroutine convertToLTri(lTri,lTRange,lTSize,mainMat,mainCols,mainRange,mainSize,mainDim,cMat,cCols,cRange,cSize,cDim)
	    implicit none
		
		integer, intent(in) :: lTSize, mainSize, mainDim, cSize, cDim
		real*8, intent(in) :: mainMat(mainSize), cMat(cSize)
		integer, intent(in) :: lTRange(0:mainDim), mainCols(mainSize), mainRange(0:mainDim), cCols(cSize), cRange(0:cDim)
		real*8, intent(out) :: lTri(lTSize)
		
		integer :: i1, i2, i3, i4, i5, i6, minCol
		
		lTri(:) = r_0
		do i1 = 1, mainDim
		    minCol = i1 - lTRange(i1) + lTRange(i1-1) + 1
			do i2 = mainRange(i1-1)+1, mainRange(i1)
			    i3 = mainCols(i2)
			    if(i3 .ge. minCol .and. i3 .le. i1) then
				    i4 = lTRange(i1) - i1 + i3
					lTri(i4) = lTri(i4) + mainMat(i2)
				endif
			enddo
		enddo
		
		do i1 = 1, cDim
		    do i2 = cRange(i1-1)+1, cRange(i1)
			    i4 = cCols(i2)
				if(i4 .ne. 0) then
					minCol = i4 - lTRange(i4) + lTRange(i4-1) + 1
					do i3 = cRange(i1-1)+1, cRange(i1)
						i5 = cCols(i3)
						if(i5 .ge. minCol .and. i5 .le. i4) then
						    i6 = lTRange(i4) - i4 + i5
							lTri(i6) = lTri(i6) + cMat(i2)*cMat(i3)
						endif
					enddo
				endif
			enddo
		enddo
		
	end subroutine convertToLTri
	
	subroutine convertToFullPop(aFull,aSparse,aCols,aRange,aSize,aDim,cMat,cCols,cRange,cSize,cDim)
	    implicit none
		
		integer, intent(in) :: aSize, aDim, cSize, cDim
		integer, intent(in) :: aCols(aSize), aRange(0:aDim), cCols(cSize), cRange(0:cDim)
		real*8, intent(in) :: aSparse(aSize), cMat(cSize)
		real*8, intent(out) :: aFull(aDim,aDim)
		
		integer :: i1, i2, i3, i4, i5
		
		aFull(:,:) = r_0
		
		do i1 = 1, aDim
		    do i2 = aRange(i1-1)+1, aRange(i1)
			    i3 = aCols(i2)
				if(i3 .ne. 0) then
				    aFull(i1,i3) = aFull(i1,i3) + aSparse(i2)
				endif
			enddo
		enddo
		
		do i1 = 1, cDim
		    do i2 = cRange(i1-1)+1, cRange(i1)
			    i4 = cCols(i2)
				if(i4 .ne. 0) then
					do i3 = cRange(i1-1)+1, cRange(i1)
						i5 = cCols(i3)
						if(i5 .ne. 0) then
						    aFull(i4,i5) = aFull(i4,i5) +cMat(i2)*cMat(i3)
						endif
					enddo
				endif
			enddo
		enddo
		
	end subroutine convertToFullPop
	
	subroutine analysisPrep()
	    implicit none
		
		integer, allocatable :: nodalConn(:)
		integer, allocatable :: nCRange(:)
		integer, allocatable :: ndInserted(:)
		integer, allocatable :: ndBlockNum(:)
		integer :: nCSize
		integer :: nextAvail, currentBlock
		integer :: eT, ndof, numSolidNds, numSectionNds, minCon, stNd
		integer :: numNds, dofPerNd, numIntDof, numIntPts, dofTable(2,33)
		real*8 :: intPts(3,8), ipWt(8)
		integer :: i1, i2, i3, i4, i5, i6, inserted
		
	    allocate(currentRank(numNodes))
		allocate(originalRank(numNodes))
		allocate(nDofIndex(6,numNodes))
		
		allocate(nCRange(0:numNodes))
		allocate(ndInserted(numNodes))
		allocate(ndBlockNum(numNodes))
		
		currentRank(:) = 0
		originalRank(:) = 0
		nDofIndex(:,:) = 0
		nCRange(:) = 0
		ndInserted(:) = 0
		ndBlockNum(:) = 0
		
		!! Reorder Nodes
		do i1 = 1, numEls
		    eT = elementType(i1)
		    if(eT .eq. 41 .or. eT .eq. 3 .or. eT .eq. 2) then
			    ndof = 6
			else
			    ndof = 3
            endif			
		    do i2 = 1, 8
			    i3 = elementList(i2,i1)
			    if(i3 .ne. 0) then
				    do i4 = 1, 8
					    i5 = elementList(i4,i1)
					    if(i5 .ne. 0 .and. i5 .ne. i3) then
						    nCRange(i3) = nCRange(i3) + 1
						endif
					enddo
					nDofIndex(1,i3) = ndof
				endif
			enddo
		enddo
		
		numSolidNds = 0		
		do i1 = 1, numNodes
		    if(nDofIndex(1,i1) .eq. 3) then
			    numSolidNds = numSolidNds + 1
			endif
		enddo
		numSectionNds = numNodes - numSolidNds
		
		do i1 = 2, numNodes
		    nCRange(i1) = nCRange(i1) + nCRange(i1-1)
		enddo
		
		nCSize = nCRange(numNodes)
		allocate(nodalConn(nCSize))
		
		nodalConn(:) = 0

        do i1 = 1, numEls		
		    do i2 = 1, 8
			    i3 = elementList(i2,i1)
			    if(i3 .ne. 0) then
				    do i4 = 1, 8
					    i5 = elementList(i4,i1)
					    if(i5 .ne. 0 .and. i5 .ne. i3) then
						    i6 = nCRange(i3-1) + 1
							inserted = 0
							do while(i6 .le. nCRange(i3) .and. inserted .eq. 0)
							    if(nodalConn(i6) .eq. 0) then
								    nodalConn(i6) = i5
									inserted = 1
								elseif(nodalConn(i6) .eq. i5) then
								    inserted = 1
								endif
								i6 = i6 + 1
							enddo
						endif
					enddo
				endif
			enddo
		enddo
		
		nextAvail = 1
		currentBlock = 0
		i1 = nCRange(numNodes)
		if(numSolidNds .gt. 0) then
            call levelSortNodes(nextAvail,numSolidNds,3,nodalConn,nCRange,i1,ndInserted,ndBlockNum,currentBlock,numNodes)
		endif
		
		if(nextAvail .le. numNodes) then
            call levelSortNodes(nextAvail,numNodes,6,nodalConn,nCRange,i1,ndInserted,ndBlockNum,currentBlock,numNodes)
		endif
		
		do i1 = 1, numNodes
		    currentRank(originalRank(i1)) = i1
		enddo
		
		!! Populate the global degree of freedom index
		do i1 = 1, numNodes
			! if(solverMeth .eq. 'iterative') then
				! if(nDofIndex(1,i1) .eq. 3) then
					! do i2 = 1, 3
						! nDofIndex(i2,i1) = currentRank(i1) + (i2-1)*numSolidNds
					! enddo
				! elseif(nDofIndex(1,i1) .eq. 6) then
					! do i2 = 1, 6
						! nDofIndex(i2,i1) = currentRank(i1) + (i2-1)*numSectionNds + 2*numSolidNds
					! enddo
				! endif
			! else
				if(nDofIndex(1,i1) .eq. 3) then
					do i2 = 1, 3
						nDofIndex(i2,i1) = 3*(currentRank(i1) - 1) + i2
					enddo
				elseif(nDofIndex(1,i1) .eq. 6) then
					do i2 = 1, 6
						nDofIndex(i2,i1) = 6*(currentRank(i1) - 1) + i2 - 3*numSolidNds
					enddo
				endif
			!endif
		enddo
		
		elMatDim = 3*numSolidNds + 6*numSectionNds
		
		if(.not. allocated(nodeTemp)) then
		    allocate(nodeTemp(numNodes))
			allocate(nodeTdot(numNodes))
			allocate(nodeDisp(elMatDim))
			allocate(nodeVel(elMatDim))
			allocate(nodeAcc(elMatDim))
			allocate(prevTemp(numNodes))
			allocate(prevTdot(numNodes))
			allocate(prevDisp(elMatDim))
			allocate(prevVel(elMatDim))
			allocate(prevAcc(elMatDim))
			allocate(elasticLoad(elMatDim))
			allocate(thermalLoad(numNodes))
			allocate(intVecRange(0:numEls))
			allocate(intMatRange(0:numEls))
			allocate(delDisp(elMatDim))
			allocate(swapVec(elMatDim))
			allocate(dLdu(elMatDim))
			allocate(dLdv(elMatDim))
			allocate(dLda(elMatDim))
			allocate(dLdt(numNodes))
			allocate(dLdtdot(numNodes))
			allocate(dsumTermdT(numNodes))
			allocate(dsumTermdU(elMatDim))
			allocate(tempAdj(numNodes))
			allocate(tdotAdj(numNodes))
			allocate(dispAdj(elMatDim))
			allocate(velAdj(elMatDim))
			allocate(accAdj(elMatDim))
			allocate(dRudD(elMatDim))
			allocate(dRtdD(numNodes))
			
			nodeTemp(:) = r_0
			nodeTdot(:) = r_0
			nodeDisp(:) = r_0
			nodeVel(:) = r_0
			nodeAcc(:) = r_0
			prevTemp(:) = r_0
			prevTdot(:) = r_0
			prevDisp(:) = r_0
			prevVel(:) = r_0
			prevAcc(:) = r_0
			elasticLoad(:) = r_0
			thermalLoad(:) = r_0
			intVecRange(:) = r_0
			intMatRange(:) = r_0
			delDisp(:) = r_0
			swapVec(:) = r_0
			dLdu(:) = r_0
			dLdv(:) = r_0
			dLda(:) = r_0
			dLdt(:) = r_0
			dLdtdot(:) = r_0
			dsumTermdT(:) = r_0
			dsumTermdU(:) = r_0
			tempAdj(:) = r_0
			tdotAdj(:) = r_0
			dispAdj(:) = r_0
			velAdj(:) = r_0
			accAdj(:) = r_0
			dRudD(:) = r_0
			dRtdD(:) = r_0
		endif
		
		intVecRange(0) = 0
		intMatRange(0) = 0
		do i1 = 1, numEls
			call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i1)
			i2 = numIntDof*(numNds*dofPerNd + numIntDof)
			intVecRange(i1) = numIntDof
			intMatRange(i1) = i2
		enddo
		
		do i1 = 2, numEls
		    intVecRange(i1) = intVecRange(i1) + intVecRange(i1-1)
			intMatRange(i1) = intMatRange(i1) + intMatRange(i1-1)
		enddo
		
		intVecSize = intVecRange(numEls)
		intMatSize = intMatRange(numEls)
		
        if(.not. allocated(internalDisp) .and. intVecSize .gt. 0) then
			allocate(internalDisp(intVecSize))
			allocate(prevIntDisp(intVecSize))
			allocate(intElasticLoad(intVecSize))
	        allocate(intswapVec(intVecSize))
			allocate(intdLdu(intVecSize))
			allocate(intdsumTermdU(intVecSize))
			allocate(intDispAdj(intVecSize))
			allocate(intdRudD(intVecSize))
			allocate(intElasticMat(intMatSize))
			
			internalDisp(:) = r_0
			prevIntDisp(:) = r_0
			intElasticLoad(:) = r_0
	        intswapVec(:) = r_0
			intdLdu(:) = r_0
			intdsumTermdU(:) = r_0
			intDispAdj(:) = r_0
			intdRudD(:) = r_0
			intElasticMat(:) = r_0
		endif
		
		call allocateElasticMat(nodalConn,nCRange,nCSize,ndBlockNum,numNodes)
		call allocateThermalMat(nodalConn,nCRange,nCSize,ndBlockNum,numNodes)
		call buildElasticMPC()
		call buildThermalMPC()
		
		call findElementSurfaces()
		
		deallocate(ndInserted)
		deallocate(ndBlockNum)
		deallocate(nodalConn)
		deallocate(nCRange)
		
	end subroutine analysisPrep
	
	subroutine deallocateModelData()
	    implicit none
		
		if(allocated(nodeList)) then
		    deallocate(nodeList)
			deallocate(elementList)
			deallocate(elementType)
			deallocate(elementSurfaces)
			deallocate(elementSection)
			deallocate(elementSets)
			deallocate(elSetName)
			deallocate(elSetRange)
			deallocate(nodeSets)
			deallocate(ndSetName)
			deallocate(ndSetRange)
			deallocate(elDepends)
			deallocate(elToDRange)
			deallocate(ndToDRange)
			deallocate(sectionType)
			deallocate(sectionMatId)
			deallocate(sectionMatName)
			deallocate(sectionOrient)
			deallocate(secLayupRange)
			deallocate(sectionZOffset)
			deallocate(beamProperties)
			deallocate(beamStiffness)
			deallocate(beamExpLoadCoef)
			deallocate(beamMass)
			deallocate(beamThermCond)
			deallocate(beamSpecHeat)
			deallocate(initialDisp)
			deallocate(initialVel)
			deallocate(initialAcc)
			deallocate(initialTemp)
			deallocate(initialTdot)
		endif
		
		if(allocated(layupMatName)) then
			deallocate(layupMatName)
			deallocate(layupMatID)
			deallocate(layupThickness)
			deallocate(layupAngle)
		endif
		
		if(allocated(materialName)) then
		    deallocate(materialName)
			deallocate(materialDensity)
			deallocate(materialElastic)
			deallocate(materialStiffMat)
			deallocate(materialThermCond)
			deallocate(materialThermExp)
			deallocate(materialSpecHeat)
			deallocate(materialMaxStress)
			deallocate(materialMaxStrain)
			deallocate(materialMaxStrnEngy)
		endif
		
		call deallocateLoadData()
		call deallocateConstraintData()
		
	end subroutine deallocateModelData
	
	subroutine deallocateLoadData()
	    implicit none
		
		if(allocated(loadNodes)) then
		    deallocate(loadNodes)
			deallocate(inputLoads)
			deallocate(loadsRange)
			deallocate(loadsActTime)
		    deallocate(loadType)
		endif
	
	end subroutine deallocateLoadData
	
	subroutine deallocateConstraintData()
	    implicit none
		
		if(allocated(mpcEqn)) then
		    deallocate(mpcEqn)
			deallocate(mpcNode)
			deallocate(mpcDof)
			deallocate(mpcCoef)
			deallocate(mpcRHS)
		endif
	end subroutine deallocateConstraintData
	
	subroutine deallocateDesignVarData()
	    implicit none
		
		if(allocated(r_dVec)) then
			deallocate(r_dVec)
			deallocate(c_dVec)
			deallocate(dLdD)
			deallocate(dCategory)
			deallocate(dComponent)
			deallocate(dLayer)
			deallocate(dActTime)
			deallocate(dElSet)
			deallocate(dNdSet)
			deallocate(dCoefList)
			deallocate(dCoefListRange)
			deallocate(dsumTermdD)
			deallocate(dtotVoldD)
			deallocate(dToEl)
			deallocate(dToElCoef)
			deallocate(dToElRange)
			deallocate(elToD)
			deallocate(elToCoef)
			deallocate(ndToD)
			deallocate(ndToCoef)
			deallocate(dToNd)
			deallocate(dToNdCoef)
			deallocate(dToNdRange)
		endif
		
	end subroutine deallocateDesignVarData
	
	subroutine deallocateObjectiveData()
	    implicit none
		
		if(allocated(objVal)) then
		    deallocate(objVal)
		    deallocate(objCategory)
			deallocate(objOperator)
			deallocate(objActTime)
			deallocate(objComponent)
			deallocate(objLayer)
			deallocate(objCoef)
			deallocate(objExp)
			deallocate(objElSet)
			deallocate(objElSetId)
			deallocate(objNdSet)
			deallocate(objNdSetId)
			deallocate(objTgtTag)
			deallocate(objTgtVal)
			deallocate(objTgtRange)
		endif
		
	end subroutine deallocateObjectiveData
	
	subroutine deallocateModalData()
	    implicit none
		
		if(allocated(eigenModes)) then
	        deallocate(eigenModes)
	        deallocate(eigenVals)
			deallocate(eigenFactors)
	        deallocate(diagMassMat)
	    endif
	
	end subroutine deallocateModalData
	
	subroutine deallocateAnalysisData()
	    implicit none

        if(allocated(currentRank)) then
			deallocate(currentRank)
			deallocate(originalRank)
			deallocate(nDofIndex)
            deallocate(nodeTemp)
			deallocate(nodeTdot)
			deallocate(nodeDisp)
			deallocate(nodeVel)
			deallocate(nodeAcc)
			deallocate(prevTemp)
			deallocate(prevTdot)
			deallocate(prevDisp)
			deallocate(prevVel)
			deallocate(prevAcc)
			deallocate(elasticLoad)
			deallocate(thermalLoad)
			deallocate(intVecRange)
			deallocate(intMatRange)
			deallocate(delDisp)
			deallocate(swapVec)
			deallocate(dLdu)
			deallocate(dLdv)
			deallocate(dLda)
			deallocate(dLdt)
			deallocate(dLdtdot)
			deallocate(dsumTermdT)
			deallocate(dsumTermdU)
			deallocate(tempAdj)
			deallocate(tdotAdj)
			deallocate(dispAdj)
			deallocate(velAdj)
			deallocate(accAdj)
			deallocate(dRudD)
			deallocate(dRtdD)
		endif
		
		if(allocated(internalDisp)) then
			deallocate(internalDisp)
			deallocate(prevIntDisp)
			deallocate(intElasticLoad)
	        deallocate(intswapVec)
			deallocate(intdLdu)
			deallocate(intdsumTermdU)
			deallocate(intDispAdj)
			deallocate(intdRudD)
			deallocate(intElasticMat)
		endif
		
		if(allocated(thermMat)) then
			deallocate(thermMatRange)
			deallocate(thermMatLTRange)
			deallocate(thermMat)
			deallocate(thermMatCols)
			deallocate(thermMatLT)
		endif
		
		if(allocated(thermMPCMat)) then
			deallocate(thermMPCMat)
			deallocate(thermMPCMatCols)
			deallocate(thermMPCMatRange)
			deallocate(thermMPCRHS)
		endif
		
		if(allocated(elasticMat)) then
			deallocate(elMatRange)
			deallocate(elMatLTRange)
			deallocate(elasticMat)
			deallocate(elMatCols)
			deallocate(elMatLT)
        endif
		
		if(allocated(elMPCMat)) then
			deallocate(elMPCMat)
			deallocate(elMPCMatCols)
			deallocate(elMPCMatRange)
			deallocate(elMPCRHS)
		endif
		
	
	end subroutine deallocateAnalysisData
	
	subroutine deallocateAll()
	    implicit none
		
		call deallocateModelData()
		call deallocateLoadData()
		call deallocateConstraintData()
		call deallocateDesignVarData()
		call deallocateObjectiveData()
		call deallocateAnalysisData()
		call deallocateModalData()
		
	end subroutine deallocateAll

end module ASenD_bookKeeping