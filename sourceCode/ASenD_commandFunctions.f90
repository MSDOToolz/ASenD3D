module ASenD_commandFunctions
    use ASenD_globalData
	use ASenD_constantVals
	use ASenD_input
	use ASenD_output
	use ASenD_r_elementEqns
	use ASenD_c_elementEqns
	use ASenD_r_designPropertyFunctions
	use ASenD_c_designPropertyFunctions
	use ASenD_solvers
	use ASenD_bookKeeping
	use ASenD_objective
	
	contains
	
	subroutine readAnyInput(fileName,commandTag)
	    implicit none
		
		character(len=128), intent(in) :: fileName
		character(len=32), intent(in) :: commandTag
		
		if(commandTag(1:15) .eq. '*readModelInput') then
		    call readModelInput(fileName)
			call readConstraints(fileName)
			call readLoads(fileName)
			call readInitialState(fileName)
		elseif(commandTag(1:10) .eq. '*readLoads') then
		    call readLoads(fileName)
		elseif(commandTag(1:16) .eq. '*readConstraints') then
		    call readConstraints(fileName)
		elseif(commandTag(1:17) .eq. '*readInitialState') then
		    call readInitialState(fileName)
		elseif(commandTag(1:19) .eq. '*readDesignVarInput') then
		    call readDesignVarInput(fileName)
			call getdToElComp()
		elseif(commandTag(1:20) .eq. '*readDesignVarValues') then
		    call readDesignVarValues(fileName)
		elseif(commandTag(1:19) .eq. '*readObjectiveInput') then
		    call readObjectiveInput(fileName)
		elseif(commandTag(1:19) .eq. '*readNodeResults') then
		    call readNodeResults(fileName)
		endif
		
	end subroutine readAnyInput
	
	subroutine processReadCommand(jobUnit,fileLine,iosVal)
	    implicit none
		
		integer, intent(in) :: jobUnit
		character(len=256), intent(out) :: fileLine
		integer, intent(out) :: iosVal
	
	    character(len=32) :: commandTag
		character(len=128) :: inputFileName
	    integer :: i1, i2
	
	    i1 = index(fileLine,'*')
		commandTag = fileLine(i1:i1+31)
	    read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
		i1 = index(fileLine,'*')
		do while(i1 .eq. 0 .and. iosVal .eq. 0)
			i2 = index(fileLine,':')
			if(i2 .gt. 0) then
				if(fileLine(i2-8:i2) .eq. 'fileName:') then
					read(fileLine(i2+1:i2+128),'(A)') inputFileName
					inputFileName = adjustl(inputFileName)
					inputFileName = trim(inputFileName)
					write(lfUnit,*) 'calling ', commandTag, 'file: ', inputFileName
					call readAnyInput(inputFileName,commandTag)
				endif
			endif
			read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
			i1 = index(fileLine,'*')
		enddo
	
	end subroutine processReadCommand
	
	subroutine getThermalSolnLoad(buildMat)
	    implicit none
		
		integer, intent(in) :: buildMat
		
		integer :: numNds, dofPerNd, numIntDof, numIntPts
		integer :: dofTable(2,33)
		real*8 :: intPts(3,8), ipWt(8)
		real*8 :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), den, tExp(6), tCond(3,3), sHeat
		real*8 :: ABD(9,9), sMass(6,6), stExp(6), stCond(3,3), ssHeat
		real*8 :: bStiff(6,6), bMass(6,6), btExp(6), btCond(3,3), bsHeat
		real*8 :: temp(10), Tdot(10), pTemp(10), pTdot(10)
		real*8 :: disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)
		
		real*8 :: Rt(10), dRdT(10,10)
		integer :: i1, i2, i3, i4, i5, i6, Kgi, Kgj, inserted, eType
		
		if(buildMat .eq. 1) then
		    thermMat(:) = r_0
		endif
		do i1 = 1, numEls
		    call r_getElementData(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
				locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, & 
				ABD, sMass, stExp, stCond, ssHeat, &
				bStiff, bMass, btExp, btCond, bsHeat, &
				temp,Tdot,disp,vel,acc,pTemp,pTdot,pDisp,pVel,pAcc,i1)
			eType = elementType(i1)
			call r_getElRt(Rt,dRdT,buildMat,numNds,numIntPts,intPts,ipWt, &
				locNds,orient,tCond,stCond,btCond,den,sHeat,ssHeat,bsHeat, &
				temp,pTemp,pTDot,eType)
			do i2 = 1, numNds
			    i3 = currentRank(elementList(i2,i1))
				thermalLoad(i3) = thermalLoad(i3) - Rt(i2)
			enddo
			if(buildMat .eq. 1) then
				do i2 = 1, numNds
					Kgi = currentRank(elementList(i2,i1))
					do i3 = 1, numNds
						Kgj = currentRank(elementList(i3,i1))
						i4 = thermMatRange(Kgi-1)+1
						inserted = 0
						do while(inserted .eq. 0 .and. i4 .lt. thermMatRange(Kgi))
							i6 = thermMatCols(i4)
							if(i6 .eq. 0) then
								thermMat(i4) = dRdT(i2,i3)
								thermMatCols(i4) = Kgj
								inserted = 1
							elseif(i6 .eq. Kgj) then
								thermMat(i4) = thermMat(i4) + dRdT(i2,i3)
								inserted = 1
							endif
							i4 = i4 + 1
						enddo
					enddo
				enddo			
			endif
		enddo
		
	end subroutine getThermalSolnLoad
	
	subroutine getThermalAppliedLoad(time)
	    implicit none
		
		real*8, intent(in) :: time
		
		integer :: dofTable(2,33), ndDof
		real*8 :: ndFrc(33), bdyFrc(6)
		real*8 :: trac(6), nDir(3), tTol	
		
		real*8 :: ndLd(7)
		character(len=16) :: lt
		integer :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10
		
		do i1 = 1, numLds
		    if(time .ge. loadsActTime(1,i1) .and. time .le. loadsActTime(2,i1)) then
			    lt = loadType(i1)
			    if(lt .eq. 'nodal') then
					do i2 = loadsRange(i1-1) + 1, loadsRange(i1)
						read(loadNodes(i2),*,err=55) i3
						i5 = currentRank(i3)
						thermalLoad(i5) = thermalLoad(i5) + inputLoads(7,i2)
						goto 63
55                      do i4 = 1, numNdSets
							if(ndSetName(i4) .eq. loadNodes(i2)) then
								do i5 = ndSetRange(i4-1)+1, ndSetRange(i4)
									i8 = currentRank(nodeSets(i5))
									thermalLoad(i8) = thermalLoad(i8) + inputLoads(7,i2)
								enddo
							endif
						enddo
63                      i4 = 1
				    enddo
				elseif(lt .eq. 'bodyHeatGen') then
					do i2 = loadsRange(i1-1) + 1, loadsRange(i1)
						read(loadNodes(i2),*,err=117) i3
						bdyFrc(1) = inputLoads(1,i2)
						bdyFrc(2:6) = r_0
						call r_getElBodyForce(ndFrc,dofTable,ndDof,bdyFrc,0,i3)
						do i4 = 1, ndDof
						    i5 = dofTable(1,i4)
							if(i5 .eq. 1) then
								i6 = currentRank(elementList(dofTable(2,i4),i3))
								thermalLoad(i6) = thermalLoad(i6) + ndFrc(i4)
							endif
						enddo
						goto 134
117                     do i4 = 1, numElSets
							if(elSetName(i4) .eq. loadNodes(i2)) then
								do i5 = elSetRange(i4-1)+1, elSetRange(i4)
									i6 = elementSets(i5)
								    bdyFrc(1) = inputLoads(1,i2)
									bdyFrc(2:6) = r_0
									call r_getElBodyForce(ndFrc,dofTable,ndDof,bdyFrc,0,i6)
									do i7 = 1, ndDof
										i8 = dofTable(1,i7)
										if(i8 .eq. 1) then
											i9 = currentRank(elementList(dofTable(2,i7),i6))
											thermalLoad(i9) = thermalLoad(i9) + ndFrc(i7)
										endif
									enddo
								enddo
							endif
						enddo
134                     i4 = 1
					enddo
                elseif(lt .eq. 'surfaceFlux') then
					do i2 = loadsRange(i1-1) + 1, loadsRange(i1)
						read(loadNodes(i2),*,err=147) i3
						trac(1) = inputLoads(5,i2)
						trac(2:6) = r_0
						nDir(:) = inputLoads(1:3,i2)
						tTol = inputLoads(4,i2)
						call r_getElSurfaceTraction(ndFrc,dofTable,ndDof,trac,r_0,i3,nDir,tTol)
						do i4 = 1, ndDof
						    i5 = dofTable(1,i4)
							if(i5 .eq. 1) then
								i6 = currentRank(elementList(dofTable(2,i4),i3))
								thermalLoad(i6) = thermalLoad(i6) - ndFrc(i4)
							endif
						enddo
						goto 170
147                     do i4 = 1, numElSets
							if(elSetName(i4) .eq. loadNodes(i2)) then
								do i5 = elSetRange(i4-1)+1, elSetRange(i4)
									i6 = elementSets(i5)
									trac(1) = inputLoads(5,i2)
						            trac(2:6) = r_0
									nDir(:) = inputLoads(1:3,i2)
									tTol = inputLoads(4,i2)
									call r_getElSurfaceTraction(ndFrc,dofTable,ndDof,trac,r_0,i6,nDir,tTol)
									do i7 = 1, ndDof
										i8 = dofTable(1,i7)
										if(i8 .eq. 1) then
											i9 = currentRank(elementList(dofTable(2,i7),i6))
											thermalLoad(i9) = thermalLoad(i9) - ndFrc(i7)
										endif
									enddo
								enddo
							endif
						enddo
170                     i4 = 1
					enddo				
				endif
			endif
		enddo
		
		do i2 = 1, numNodes
		    call r_getNodeDLoad(time,ndLd,i2)
			i3 = currentRank(i2)
			thermalLoad(i3) = thermalLoad(i3) + ndLd(7)
		enddo
	
	end subroutine getThermalAppliedLoad
	
	subroutine getThermalConstraintLoad()
	    implicit none
		
	    integer :: i1, i2, i3
		real*8, allocatable :: tempV(:)
		
		allocate(tempV(thermMPCDim))
		
		tempV(:) = thermMPCRHS(:)
		do i1 = 1, thermMPCDim
		    do i2 = thermMPCMatRange(i1-1)+1,thermMPCMatRange(i1)
			    i3 = thermMPCMatCols(i2)
				if(i3 .ne. 0) then
				    tempV(i1) = tempV(i1) - thermMPCMat(i2)*nodeTemp(i3)
				endif
			enddo
		enddo
		
		do i1 = 1, thermMPCDim
		    do i2 = thermMPCMatRange(i1-1)+1,thermMPCMatRange(i1)
			    i3 = thermMPCMatCols(i2)
				if(i3 .ne. 0) then
				    thermalLoad(i3) = thermalLoad(i3) + thermMPCMat(i2)*tempV(i1)
				endif
			enddo
		enddo
		
		deallocate(tempV)
		
	end subroutine getThermalConstraintLoad
	
	subroutine scaleThermalMPC()
	    implicit none
		
		real*8 :: maxMat, rowMag
		integer :: i1, i2, i3, i4
		
		maxMat = r_0
		do i1 = 1, thermMatSize
		    rowMag = abs(thermMat(i1))
		    if(rowMag .gt. maxMat) then
			    maxMat = rowMag
			endif
		enddo
		
		do i1 = 1, thermMPCDim
		    i3 = thermMPCMatRange(i1-1)+1
			i4 = thermMPCMatRange(i1)
		    rowMag = r_0
			do i2 = i3, i4
			    rowMag = rowMag + thermMPCMat(i2)*thermMPCMat(i2)
			enddo
			rowMag = sqrt(rowMag)
			thermMPCMat(i3:i4) = (r_1/rowMag)*thermMPCMat(i3:i4)
			thermMPCRHS(i1) = (r_1/rowMag)*thermMPCRHS(i1)
		enddo
		
		thermMPCMat = 100d0*maxMat*thermMPCMat
		thermMPCRHS = 100d0*maxMat*thermMPCRHS
	
	end subroutine scaleThermalMPC
	
	subroutine getEldRudUComplex(Ru,dRdU,eType,numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
					    locNds,globNds,orient,cMat,tExp,ABD,stExp,bStiff,btExp,temp,disp,vel,acc,pDisp,pVel,pAcc, &
					    statInOri,den,sMass,bMass)
	    implicit none
		
		real*8, intent(out) :: dRdU(33,33)
		real*8, intent(in) :: Ru(33), intPts(3,8), ipWt(8)
		real*8, intent(in) :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), tExp(6), ABD(9,9), stExp(6)
		real*8, intent(in) :: bStiff(6,6), btExp(6), temp(10), disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)
		real*8, intent(in) :: statInOri(3,240), den, sMass(6,6), bMass(6,6)
		integer, intent(in) :: eType, numNds, dofPerNd, numIntDof, numIntPts, dofTable(2,33)
		
		complex*16 :: c_dRdU(33,33)
		complex*16 :: c_Ru(33), c_intPts(3,8), c_ipWt(8)
		complex*16 :: c_locNds(3,10), c_globNds(3,10), c_orient(3,3), c_cMat(6,6), c_tExp(6), c_ABD(9,9), c_stExp(6)
		complex*16 :: c_bStiff(6,6), c_btExp(6), c_temp(10), c_disp(6,11), c_vel(6,11), c_acc(6,11)
		complex*16 :: c_pDisp(6,11), c_pVel(6,11), c_pAcc(6,11)
		complex*16 :: c_statInOri(3,240), c_den, c_sMass(6,6), c_bMass(6,6)
		
		integer :: i1, i2, i3, i4, i5
		
		c_dRdU = c_1*dRdU
		c_Ru = c_1*Ru
		c_intPts = c_1*intPts
		c_ipWt = c_1*ipWt
		c_locNds = c_1*locNds
		c_globNds = c_1*globNds
		c_orient = c_1*orient
		c_cMat = c_1*cMat
		c_tExp = c_1*tExp
		c_ABD = c_1*ABD
		c_stExp = c_1*stExp
		c_bStiff = c_1*bStiff
		c_btExp = c_1*btExp
		c_temp = c_1*temp
		c_disp = c_1*disp
		c_vel = c_1*vel
		c_acc = c_1*acc
		c_pDisp = c_1*pDisp
		c_pVel = c_1*pVel 
		c_pAcc = c_1*pAcc
		c_statInOri = c_1*statInOri
		c_den = c_1*den
		c_sMass = c_1*sMass
		c_bMass = c_1*bMass
		
		i1 = numNds*dofPerNd + numIntDof
		do i2 = 1, i1
		    i3 = dofTable(1,i2)
			i4 = dofTable(2,i2)
			c_disp(i3,i4) = c_disp(i3,i4) + compStep
			call c_getElRu(c_Ru,c_dRdU,0,eType,numNds,dofPerNd,numIntDof,numIntPts,dofTable,c_intPts,c_ipWt, &
				    c_locNds,c_globNds,c_orient,c_cMat,c_tExp,c_ABD,c_stExp,c_bStiff,c_btExp,c_temp,c_disp, &
					c_vel,c_acc,c_pDisp,c_pVel,c_pAcc, &
				    c_statInOri,c_den,c_sMass,c_bMass)
			do i5 = 1, i1
			    dRdU(i5,i2) = compStepInv*imag(c_Ru(i5))
			enddo
			c_disp(i3,i4) = c_disp(i3,i4) - compStep
		enddo
		
		
	end subroutine getEldRudUComplex
	
	subroutine getElasticSolnLoad(buildMat)
	    implicit none
		
		integer, intent(in) :: buildMat
		
		real*8 :: Ru(33), dRdU(33,33), intPts(3,8), ipWt(8)
		real*8 :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), tExp(6), tCond(3,3), sHeat, ABD(9,9), stExp(6)
		real*8 :: stCond(3,3), ssHeat, btCond(3,3), bsHeat
		real*8 :: bStiff(6,6), btExp(6), temp(10), disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)
		real*8 :: Tdot(10), pTemp(10), pTdot(10)
		real*8 :: statInOri(3,240), den, sMass(6,6), bMass(6,6)
		real*8 :: nnInv
		real*8 :: resVec(33), bVec(33), zVec(33)
		integer :: eType, numNds, dofPerNd, numIntDof, numIntPts, dofTable(2,33)
		integer :: i1, i2, i3, i4, i5, i6, kGi, kGj, inserted
		real*8 :: dRdUCopy(33,33)
		
		if(buildMat .eq. 1) then
		    elasticMat(:) = r_0
		endif
		do i1 = 1, numEls
		    call r_getElementData(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
				locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, & 
				ABD, sMass, stExp, stCond, ssHeat, &
				bStiff, bMass, btExp, btCond, bsHeat, &
				temp,Tdot,disp,vel,acc,pTemp,pTdot,pDisp,pVel,pAcc,i1)
			!! -------------------
			! if(i1 .eq. 200) then
			    ! write(lfUnit,*) 'element ', i1, ' data'
				! write(lfUnit,*) 'numNds: ', numNds, 'dofPerNd: ', dofPerNd, 'numIntDof: ', numIntDof, 'numIntPts: ', numIntPts
				! write(lfUnit,*) 'intPts: '
				! write(lfUnit,*) intPts
				! write(lfUnit,*) 'ipWt'
				! write(lfUnit,*) ipWt
			    ! write(lfUnit,*) 'locNds: '
				! write(lfUnit,*) locNds
				! write(lfUnit,*) 'globNds: '
				! write(lfUnit,*) globNds
				! write(lfUnit,*) 'orient: '
				! write(lfUnit,*) orient
				! write(lfUnit,*) 'cMat'
				! write(lfUnit,*) cMat
				! write(lfUnit,*) 'ABD: '
				! write(lfUnit,*) ABD
				! write(lfUnit,*) 'stExp: '
				! write(lfUnit,*) stExp
				! write(lfUnit,*) 'temp: '
				! write(lfUnit,*) temp
				! write(lfUnit,*) 'disp: '
				! write(lfUnit,*) disp
			! endif
			!! ----------------------------
			eType = elementType(i1)
			if(eType .eq. 41 .or. eType .eq. 3 .or. eType .eq. 2) then
			    call r_getInstOrient(statInOri,orient,disp,numNds,1)
			endif
		    call r_getElRu(Ru,dRdU,buildMat,eType,numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
				locNds,globNds,orient,cMat,tExp,ABD,stExp,bStiff,btExp,temp,disp,vel,acc,pDisp,pVel,pAcc, &
				statInOri,den,sMass,bMass)
			! if(buildMat .eq. 1) then
			    ! call getEldRudUComplex(Ru,dRdU,eType,numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
					    ! locNds,globNds,orient,cMat,tExp,ABD,stExp,bStiff,btExp,temp,disp,vel,acc,pDisp,pVel,pAcc, &
					    ! statInOri,den,sMass,bMass)
				! if(i1 .eq. 5) then
					! do i2 = 1, 28
						! do i3 = 1, 28
							! if(abs(dRdU(i2,i3) - dRdUCopy(i2,i3)) .gt. 1e-3) then
								! write(lfUnit,*) 'discrepancy: ', i2, i3, dRdU(i2,i3), dRdUCopy(i2,i3)
							! endif
						! enddo
					! enddo
				! endif
			! endif
			! write(lfUnit,*) 'Element ', i1, 'Ru:'
			! write(lfUnit,*) Ru(1:24)
			i3 = numNds*dofPerNd
			do i2 = 1, i3
			    i4 = dofTable(1,i2)
				i5 = elementList(dofTable(2,i2),i1)
				i6 = nDofIndex(i4,i5)
				if(i6 .gt. 0) then
				    elasticLoad(i6) = elasticLoad(i6) - Ru(i2)
				endif
			enddo
			if(numIntDof .ne. 0) then
			    i4 = intVecRange(i1-1) + 1
				i5 = intVecRange(i1)
				intElasticLoad(i4:i5) = intElasticLoad(i4:i5) - Ru(i3+1:i3+numIntDof)
			endif
			if(buildMat .eq. 1) then
			    ! -------------------
			    ! if(i1 .eq. 1) then
					! write(lfUnit,*) 'first element elastic matrix:'
					! i3 = numNds*dofPerNd + numIntDof
					! do i2 = 1, i3
						! write(lfUnit,*) dRdU(1:i3,i2)
					! enddo
					! close(lfUnit)
					! open(unit=lfUnit, file='jobLogFile.txt', action='write', access='append')
				! endif
				! ---------------------
				if(numIntDof .gt. 0) then
					i2 = numNds*dofPerNd + 1
					i3 = numNds*dofPerNd + numIntDof
					call r_condenseMat(dRdU,33,resVec,bVec,zVec,i2,i3)
					i4 = intMatRange(i1-1) + 1
					do i5 = 1, i3
						do i6 = i2, i3
							intElasticMat(i4) = dRdU(i6,i5)
							i4 = i4 + 1
						enddo
					enddo
				endif
				!! ---------------
				! if(i1 .eq. 200) then
					! dRdU(1,1) = 10000d0*dRdU(1,1)
					! dRdU(2,2) = 10000d0*dRdU(2,2)
					! dRdU(3,3) = 10000d0*dRdU(3,3)
					! dRdU(4,4) = 10000d0*dRdU(4,4)
					! dRdU(9,9) = 10000d0*dRdU(9,9)
					! dRdU(17,17) = 10000d0*dRdU(17,17)
					! dRdU(18,18) = 10000d0*dRdU(18,18)
					! bVec(:) = r_0
					! bVec(5) = 0.25d0
					! bVec(6) = 0.25d0
					! bVec(7) = 0.25d0
					! bVec(8) = 0.25d0
					! resVec(:) = r_0
					! write(lfUnit,*) 'single el mat:'
					! do i2 = 1, 24
					    ! write(lfUnit,*) dRdU(1:24,i2)
					! enddo
					! call rFactorMat(dRdU(1:24,1:24),24,24,1,24,1,24,0)
					! write(lfUnit,*) 'factored single el:'
					! do i2 = 1, 24
					    ! write(lfUnit,*) dRdU(1:24,i2)
					! enddo
					! call solveRFactor(resVec(1:24),dRdU(1:24,1:24),bVec(1:24),24,24,1,24,1,24,0)
					! write(lfUnit,*) 'single el result:'
					! do i2 = 1, numNds
					    ! write(lfUnit,*) resVec(i2), resVec(i2+8), resVec(i2+16)
					! enddo
				! endif
				!! -------------------
				i2 = numNds*dofPerNd
				do i3 = 1, i2
					i5 = dofTable(1,i3)
					i6 = elementList(dofTable(2,i3),i1)
					kGi = nDofIndex(i5,i6)
					do i4 = 1, i2
						i5 = dofTable(1,i4)
						i6 = elementList(dofTable(2,i4),i1)
						kGj = nDofIndex(i5,i6)
						i5 = elMatRange(kGi-1) + 1
						inserted = 0
						do while(inserted .eq. 0 .and. i5 .le. elMatRange(kGi))
							i6 = elMatCols(i5)
							if(i6 .eq. 0) then
								elasticMat(i5) = dRdU(i3,i4)
								elMatCols(i5) = kGj
								inserted = 1
							elseif(i6 .eq. kGj) then
								elasticMat(i5) = elasticMat(i5) + dRdU(i3,i4)
								inserted = 1
							endif
							i5 = i5 + 1
						enddo
					enddo
				enddo
			endif
		enddo
		
	end subroutine getElasticSolnLoad
	
	subroutine getElasticAppliedLoad(time)
	    implicit none
		
		real*8, intent(in) :: time
		
		integer :: isGrav, dofTable(2,33), ndDof
		real*8 :: ndLd(7), nnInv, ndFrc(33), bdyFrc(6)
		real*8 :: trac(6), pressure, nDir(3), tTol
		real*8 :: center(3), axis(3), angVel
		character(len=16) :: lt
		integer :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10
		
		write(lfUnit,*) 'building elastic applied load'
		
		do i1 = 1, numLds
		    if(time .ge. loadsActTime(1,i1) .and. time .le. loadsActTime(2,i1)) then
			    lt = loadType(i1)
			    if(lt .eq. 'nodal') then
					do i2 = loadsRange(i1-1) + 1, loadsRange(i1)
						read(loadNodes(i2),*,err=67) i3
						do i4 = 1, 6
							i5 = nDofIndex(i4,i3)
							if(i5 .ne. 0) then
								elasticLoad(i5) = elasticLoad(i5) + inputLoads(i4,i2)
							endif
						enddo
						goto 75
67                      do i4 = 1, numNdSets
							if(ndSetName(i4) .eq. loadNodes(i2)) then
								do i5 = ndSetRange(i4-1)+1, ndSetRange(i4)
									i6 = nodeSets(i5)
									do i7 = 1, 6
										i8 = nDofIndex(i7,i6)
										if(i8 .ne. 0) then
											elasticLoad(i8) = elasticLoad(i8) + inputLoads(i7,i2)
										endif
									enddo
								enddo
							endif
						enddo
75                      i4 = 1
					enddo
				elseif(lt .eq. 'bodyForce' .or. lt .eq. 'gravitational' .or. lt .eq. 'centrifugal') then
				    if(lt .eq. 'bodyForce') then
					    isGrav = 0
				    else
					    isGrav = 1
					endif
					do i2 = loadsRange(i1-1) + 1, loadsRange(i1)
						read(loadNodes(i2),*,err=308) i3
						if(lt .eq. 'centrifugal') then
						    center(:) = inputLoads(1:3,i2)
							axis(:) = inputLoads(4:6,i2)
							angVel = inputLoads(7,i2)
						    call r_getElCentAcc(bdyFrc(1:3),center,axis,angVel,i3)
							bdyFrc(4:6) = r_0
						else
						    bdyFrc(:) = inputLoads(1:6,i2)
						endif
						call r_getElBodyForce(ndFrc,dofTable,ndDof,bdyFrc,isGrav,i3)
						do i4 = 1, ndDof
						    i5 = dofTable(1,i4)
							i6 = elementList(dofTable(2,i4),i3)
							i7 = nDofIndex(i5,i6)
							elasticLoad(i7) = elasticLoad(i7) + ndFrc(i4)
						enddo
						goto 328
308                     do i4 = 1, numElSets
							if(elSetName(i4) .eq. loadNodes(i2)) then
								do i5 = elSetRange(i4-1)+1, elSetRange(i4)
									i6 = elementSets(i5)
									if(lt .eq. 'centrifugal') then
									    center(:) = inputLoads(1:3,i2)
							            axis(:) = inputLoads(4:6,i2)
							            angVel = inputLoads(7,i2)
										call r_getElCentAcc(bdyFrc(1:3),center,axis,angVel,i6)
										bdyFrc(4:6) = r_0
									else
										bdyFrc(:) = inputLoads(1:6,i2)
									endif
									call r_getElBodyForce(ndFrc,dofTable,ndDof,bdyFrc,isGrav,i6)
									do i7 = 1, ndDof
										i8 = dofTable(1,i7)
										i9 = elementList(dofTable(2,i7),i6)
										i10 = nDofIndex(i8,i9)
										elasticLoad(i10) = elasticLoad(i10) + ndFrc(i7)
									enddo
								enddo
							endif
						enddo
328                     i4 = 1
					enddo
                elseif(lt .eq. 'surfacePressure' .or. lt .eq. 'surfaceTraction') then
					do i2 = loadsRange(i1-1) + 1, loadsRange(i1)
					    if(lt .eq. 'surfacePressure') then
						    pressure = inputLoads(5,i2)
							trac(:) = r_0
						else
						    pressure = r_0
							trac(1:3) = inputLoads(5:7,i2)
							trac(4:6) = r_0
						endif
						read(loadNodes(i2),*,err=367) i3
						nDir(:) = inputLoads(1:3,i2)
						tTol = inputLoads(4,i2)
						call r_getElSurfaceTraction(ndFrc,dofTable,ndDof,trac,pressure,i3,nDir,tTol)
						do i4 = 1, ndDof
						    i5 = dofTable(1,i4)
							i6 = elementList(dofTable(2,i4),i3)
							i7 = nDofIndex(i5,i6)
							elasticLoad(i7) = elasticLoad(i7) + ndFrc(i4)
						enddo
						goto 387
367                     do i4 = 1, numElSets
							if(elSetName(i4) .eq. loadNodes(i2)) then
								do i5 = elSetRange(i4-1)+1, elSetRange(i4)
									i6 = elementSets(i5)
									nDir(:) = inputLoads(1:3,i2)
									tTol = inputLoads(4,i2)
									call r_getElSurfaceTraction(ndFrc,dofTable,ndDof,trac,pressure,i6,nDir,tTol)
									do i7 = 1, ndDof
										i8 = dofTable(1,i7)
										i9 = elementList(dofTable(2,i7),i6)
										i10 = nDofIndex(i8,i9)
										elasticLoad(i10) = elasticLoad(i10) + ndFrc(i7)
									enddo
								enddo
							endif
						enddo
387                     i4 = 1
					enddo			
				endif
			endif
		enddo
		
		do i2 = 1, numNodes
		    call r_getNodeDLoad(time,ndLd,i2)
			do i3 = 1, 6
			    i4 = nDofIndex(i3,i2)
				if(i4 .gt. 0) then
			        elasticLoad(i4) = elasticLoad(i4) + ndLd(i3)
				endif
			enddo
		enddo
		
	end subroutine getElasticAppliedLoad
	
	subroutine getElasticConstraintLoad()
	    implicit none
		
		integer :: i1, i2, i3
		real*8, allocatable :: tempV(:)
		
		allocate(tempV(elMPCDim))
		
		tempV(:) = elMPCRHS(:)
		do i1 = 1, elMPCDim
		    do i2 = elMPCMatRange(i1-1)+1,elMPCMatRange(i1)
			    i3 = elMPCMatCols(i2)
				if(i3 .ne. 0) then
				    tempV(i1) = tempV(i1) - elMPCMat(i2)*nodeDisp(i3)
				endif
			enddo
		enddo
		
		do i1 = 1, elMPCDim
		    do i2 = elMPCMatRange(i1-1)+1,elMPCMatRange(i1)
			    i3 = elMPCMatCols(i2)
				if(i3 .ne. 0) then
				    elasticLoad(i3) = elasticLoad(i3) + elMPCMat(i2)*tempV(i1)
				endif
			enddo
		enddo
		
		deallocate(tempV)
		
	end subroutine getElasticConstraintLoad
	
	subroutine scaleElasticMPC()
	    implicit none
		
		real*8 :: maxMat, rowMag
		integer :: i1, i2, i3, i4
		
		maxMat = r_0
		do i1 = 1, elMatSize
		    rowMag = abs(elasticMat(i1))
		    if(rowMag .gt. maxMat) then
			    maxMat = rowMag
			endif
		enddo
		
		do i1 = 1, elMPCDim
		    i3 = elMPCMatRange(i1-1)+1
			i4 = elMPCMatRange(i1)
		    rowMag = r_0
			do i2 = i3, i4
			    rowMag = rowMag + elMPCMat(i2)*elMPCMat(i2)
			enddo
			rowMag = sqrt(rowMag)
			elMPCMat(i3:i4) = (r_1/rowMag)*elMPCMat(i3:i4)
			elMPCRHS(i1) = (r_1/rowMag)*elMPCRHS(i1)
		enddo
		
		elMPCMat = 100d0*maxMat*elMPCMat
		elMPCRHS = 100d0*maxMat*elMPCRHS
		
	end subroutine scaleElasticMPC
	
	subroutine getDiagMassMat()
	    implicit none
		
		integer :: numNds, dofPerNd, numIntDof, numIntPts
		integer :: dofTable(2,33)
		real*8 :: intPts(3,8), ipWt(8)
		real*8 :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), den, tExp(6), tCond(3,3), sHeat
		real*8 :: ABD(9,9), sMass(6,6), sTELd(6), stCond(3,3), ssHeat
		real*8 :: bStiff(6,6), bMass(6,6), bTELd(6), btCond(3,3), bsHeat
		real*8 :: temp(10), Tdot(10), pTemp(10), pTdot(10)
		real*8 :: disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)		
		
		integer :: eType
		real*8 :: Rum(33), dRdA(33,33)
		real*8 :: statInOri(3,240)
		real*8 :: nnInv
		integer :: i1, i2, i3, i4, i5, i6

        if(.not. allocated(diagMassMat)) then
		    allocate(diagMassMat(elMatDim))
		endif

        diagMassMat(:) = r_0		
		do i1 = 1, numEls
		    call r_getElementData(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
				locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, & 
				ABD, sMass, sTELd, stCond, ssHeat, &
				bStiff, bMass, bTELd, btCond, bsHeat, &
				temp,Tdot,disp,vel,acc,pTemp,pTdot,pDisp,pVel,pAcc,i1)
			eType = elementType(i1)
			if(eType .eq. 41 .or. eType .eq. 3 .or. eType .eq. 2) then
			    call r_getInstOrient(statInOri,orient,disp,numNds,1)
			endif
			acc(:,:) = r_0
			i3 = numNds*dofPerNd
			do i2 = 1, i3
			    i4 = dofTable(1,i2)
				i5 = dofTable(2,i2)
				acc(i4,i5) = r_1
			enddo
		    call r_getElRum(Rum,dRdA,0,eType,numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
	            locNds, globNds, orient, statInOri, den, sMass, bMass, disp, acc)
			do i2 = 1, i3
			    i4 = dofTable(1,i2)
				i5 = elementList(dofTable(2,i2),i1)
				i6 = nDofIndex(i4,i5)
				diagMassMat(i6) = diagMassMat(i6) + Rum(i2)
			enddo
		enddo
		
	end subroutine getDiagMassMat
	
	subroutine updateExternalRHS(extVec,extDim,intVec,intDim)
	    implicit none
		
		integer, intent(in) :: extDim, intDim
		real*8, intent(in) :: intVec(intDim)
		real*8, intent(out) :: extVec(extDim) 
		
		real*8 :: elIntMat(33,33), temp1(33), temp2(33), zVec(33)
		real*8 :: intPts(3,8), ipWt(8)
		integer :: numNds, dofPerNd, numIntDof, numIntPts, dofTable(2,33)
		integer :: i1, i2, i3, i4, i5, i6, nodalDof, totalDof, eType
	    
		do i1 = 1, numEls
		    call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i1)
			if(numIntDof .ne. 0) then
				nodalDof = numNds*dofPerNd
				totalDof = nodalDof + numIntDof
				i4 = intMatRange(i1-1) + 1
				do i2 = 1, totalDof
				    do i3 = nodalDof+1, totalDof
					    elIntMat(i3,i2) = intElasticMat(i4)
						i4 = i4 + 1
					enddo
				enddo
				i4 = intVecRange(i1-1) + 1
				i5 = intVecRange(i1)
				temp2(nodalDof+1:totalDof) = intVec(i4:i5)
				call r_SolveLUxb(temp1, elIntMat, temp2, zVec, 33, nodalDof+1, totalDof)
				do i2 = 1, nodalDof
				    i4 = elementList(dofTable(2,i2),i1)
					i5 = dofTable(1,i2)
					i6 = nDofIndex(i5,i4)
				    do i3 = nodalDof+1, totalDof
					    extVec(i6) = extVec(i6) - elIntMat(i3,i2)*temp1(i3)
					enddo
				enddo
			endif
		enddo
	
	end subroutine updateExternalRHS
	
	subroutine updateInternalSoln(extSoln,extDim,intSoln,intRHS,intDim)
	    implicit none
		
		integer, intent(in) :: extDim, intDim
		real*8, intent(in) :: extSoln(extDim), intRHS(intDim)
		real*8, intent(out) :: intSoln(intDim)
		
		real*8 :: elIntMat(33,33), temp1(33), temp2(33), zVec(33)
		real*8 :: intPts(3,8), ipWt(8)
		integer :: numNds, dofPerNd, numIntDof, numIntPts, dofTable(2,33)
		integer :: i1, i2, i3, i4, i5, i6, nodalDof, totalDof, eType
	    
		do i1 = 1, numEls
		    call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i1)
			if(numIntDof .ne. 0) then
				nodalDof = numNds*dofPerNd
				totalDof = nodalDof + numIntDof
				i4 = intMatRange(i1-1) + 1
				do i2 = 1, totalDof
				    do i3 = nodalDof+1, totalDof
					    elIntMat(i3,i2) = intElasticMat(i4)
						i4 = i4 + 1
					enddo
				enddo
				i4 = intVecRange(i1-1) + 1
				i5 = intVecRange(i1)
				temp2(nodalDof+1:totalDof) = intRHS(i4:i5)
				do i2 = 1, nodalDof
				    i4 = elementList(dofTable(2,i2),i1)
					i5 = dofTable(1,i2)
					i6 = nDofIndex(i5,i4)
				    do i3 = nodalDof+1, totalDof
					    temp2(i3) = temp2(i3) - elIntMat(i3,i2)*extSoln(i6)
					enddo
				enddo
				call r_SolveLUxb(temp1, elIntMat, temp2, zVec, 33, nodalDof+1, totalDof)
				i4 = intVecRange(i1-1) + 1
				i5 = intVecRange(i1)
				intSoln(i4:i5) = intSoln(i4:i5) + temp1(nodalDof+1:totalDof)
			endif
		enddo
		
	end subroutine updateInternalSoln
	
	subroutine loadInitialState()
	    implicit none
		
		integer :: i1, i2, i3
		
		prevTemp(:) = r_0
		prevTdot(:) = r_0
		prevDisp(:) = r_0
		prevVel(:) = r_0
		prevAcc(:) = r_0
		if(intVecSize .gt. 0) then
		    prevIntDisp(:) = r_0
		endif
		
		do i1 = 1, numNodes
		    i3 = currentRank(i1)
			prevTemp(i3) = initialTemp(i1)
			prevTdot(i3) = initialTdot(i1)
		    do i2 = 1, 6
			    i3 = nDofIndex(i2,i1)
				if(i3 .ne. 0) then
				    prevDisp(i3) = initialDisp(i2,i1)
					prevVel(i3) = initialVel(i2,i1)
					prevAcc(i3) = initialAcc(i2,i1)
				endif
			enddo
		enddo
		
	end subroutine loadInitialState
	
	subroutine getModelDimension(modDim)
	    implicit none
		
		real*8, intent(out) :: modDim
		
		real*8 :: xMin, xMax, yMin, yMax, zMin, zMax
		real*8 :: xLen, yLen, zLen, md2
		real*8 :: ndCrd(3)
		integer :: i1, i2, i3, i4
		
		xMin = 1e+100_8
		yMin = xMin
		zMin = xMin
		xMax = -xMin
		yMax = xMax
		zMax = xMax
		
		do i1 = 1, numNodes
		    ndCrd(:) = nodeList(:,i1)
			do i2 = ndToDRange(i1-1)+1, ndToDRange(i1)
				i3 = ndToD(i2)
				if(dCategory(i3) .eq. 'nodeCoord') then
					i4 = dComponent(i3)
					ndCrd(i4) = ndCrd(i4) + ndToCoef(i2)*r_dVec(i3)
				endif
			enddo
			if(ndCrd(1) .gt. xMax) then
			    xMax = ndCrd(1)
			endif
			if(ndCrd(1) .lt. xMin) then
			    xMin = ndCrd(1)
			endif
			if(ndCrd(2) .gt. yMax) then
			    yMax = ndCrd(2)
			endif
			if(ndCrd(2) .lt. yMin) then
			    yMin = ndCrd(2)
			endif
			if(ndCrd(3) .gt. zMax) then
			    zMax = ndCrd(3)
			endif
			if(ndCrd(3) .lt. zMin) then
			    zMin = ndCrd(3)
			endif
		enddo
		
		xLen = xMax - xMin
		yLen = yMax - yMin
		zLen = zMax - zMin
		
		md2 = xLen*xLen + yLen*yLen + zLen*zLen
		modDim = sqrt(md2)
		
	end subroutine getModelDimension
	
	subroutine setSolnToMode(solnField,modeNum,maxAmp)
	    implicit none
		
		character(len=16), intent(in) :: solnField
		integer, intent(in) :: modeNum
		real*8, intent(in) :: maxAmp
		
		real*8 :: modDim, eVMax, multFact, absVi
		integer :: i1, i2
		
		eVMax = r_0
		do i1 = 1, elMatDim
		    absVi = abs(eigenModes(i1,modeNum))
		    if(absVi .gt. eVMax) then
			    eVMax = absVi
			endif
		enddo
		
		multFact = maxAmp/eVMax
		
		if(solnField .eq. 'velocity') then
		    nodeVel(:) = multFact*eigenModes(:,modeNum)
		elseif(solnField .eq. 'acceleration') then
		    nodeAcc(:) = multFact*eigenModes(:,modeNum)
		else
		    nodeDisp(:) = multFact*eigenModes(:,modeNum)
		endif
		
	end subroutine setSolnToMode
	
	subroutine getBucklingLoadFactors(bFact,numFact)
	    implicit none
		
		integer, intent(in) :: numFact
		real*8, intent(out) :: bFact(numFact)
		
		real*8, allocatable :: intSwap(:)
		real*8 :: modDim, numer, denom, K0VMag
		real*8 :: alpha(2)
		integer :: i1, i2, i3, i4, dynCopy
		
		if(intVecSize .gt. 0) then
		    allocate(intSwap(intVecSize))
			intSwap(:) = internalDisp(:)
			internalDisp(:) = r_0
		endif
		
		swapVec(:) = nodeDisp(:)
		nodeDisp(:) = r_0
		
		dynCopy = dynamic
		dynamic = 0
		
		call getElasticSolnLoad(1)
		
		do i1 = 1, numFact
		    numer = r_0
		    do i2 = 1, elMatDim
			    do i3 = elMatRange(i2-1)+1, elMatRange(i2)
				    i4 = elMatCols(i3)
					numer = numer + eigenModes(i2,i1)*elasticMat(i3)*eigenModes(i4,i1)
				enddo
			enddo
			denom = numer - eigenVals(i1)
			if(abs(denom) .gt. r_0) then
			    bFact(i1) = numer/denom
			else
			    bFact(i1) = 1e+100_8
			endif
		enddo
		
		dynamic = dynCopy

        nodeDisp(:) = swapVec(:)
        if(intVecSize .gt. 0) then
		    internalDisp(:) = intSwap(:)
            deallocate(intSwap)
		endif
		
	end subroutine getBucklingLoadFactors
	
	subroutine stepSolve(time,appLdFact)
	    implicit none
		
		real*8, intent(in) :: time, appLdFact
		
		real*8 :: delUNorm, mag
		integer :: i1, i2, i3, i4, maxNLIt, numVecs
		
		if(solveThermal .eq. 1) then
		    numVecs = solverMaxBW
			thermalLoad(:) = r_0
			call getThermalSolnLoad(0)
			call getThermalAppliedLoad(time)
			call getThermalConstraintLoad()
			if(solverMeth .eq. 'direct') then
			    call solveSparseLDLFact(delDisp(1:numNodes),thermalLoad,thermMatLT,thermMatLTSize, &
				    thermMatLTRange,numNodes,swapVec(1:numNodes))
			else
				call gMRes(delDisp(1:numNodes),thermMat,thermMatCols,thermMatRange,thermMatSize,thermMatDim, &
					thermMPCMat,thermMPCMatCols,thermMPCMatRange,thermMPCSize,thermMPCDim,thermMatLT,thermMatLTRange, &
					thermMatLTSize,thermalLoad,numVecs,thermMatDim)
			endif
			nodeTemp(:) = nodeTemp(:) + delDisp(1:numNodes)		
		endif
	
	    if(solveElastic .eq. 1) then
			if(nLGeom .eq. 0) then
				maxNLIt = 1
			else
			    maxNLIt = 20
			endif
			numVecs = solverMaxBW
			i1 = 0
			delUNorm = r_1
			do while(i1 .lt. maxNLIt .and. delUNorm .gt. 1e-12)
			    elasticLoad(:) = r_0
				if(intVecSize .gt. 0) then
				    intElasticLoad(:) = r_0
					! internalDisp(:) = r_0
				endif
				call getElasticAppliedLoad(time)
				elasticLoad(:) = appLdFact*elasticLoad(:)
				if(nLGeom .ne. 0) then
			        call getElasticSolnLoad(1)
					call convertToLTri(elMatLT,elMatLTRange,elMatLTSize,elasticMat,elMatCols,elMatRange, &
						 elMatSize,elMatDim,elMPCMat,elMPCMatCols,elMPCMatRange,elMPCSize,elMPCDim)
					call getSparseLDLFact(elMatLT,elMatLTSize,elMatLTRange,elMatDim)
				else
				    call getElasticSolnLoad(0)
				endif
				call getElasticConstraintLoad()
				if(intVecSize .gt. 0) then
				    call updateExternalRHS(elasticLoad,elMatDim,intElasticLoad,intVecSize)
				endif
				write(lfUnit,*) 'calling solver, elastic solution'
				!----
				! write(lfUnit,*) 'applied load:'
				! do i2 = 1, numNodes
				    ! do i3 = 1, 3
					    ! i4 = nDofIndex(i3,i2)
						! write(lfUnit,*) i2, i3, elasticLoad(i4)
					! enddo
				! enddo
				!----
				if(solverMeth .eq. 'direct') then
				    call solveSparseLDLFact(delDisp, elasticLoad, elMatLT, elMatLTSize, elMatLTRange, elMatDim, swapVec)
				else
					call gMRes(delDisp,elasticMat,elMatCols,elMatRange,elMatSize,elMatDim,elMPCMat,elMPCMatCols, &
						elMPCMatRange,elMPCSize,elMPCDim,elMatLT,elMatLTRange,elMatLTSize,elasticLoad,numVecs,4*elMatDim)
				endif
				!call solveRFactor(delDisp,aFull,elasticLoad,elMatDim,elMatDim,0)
				if(intVecSize .gt. 0) then
				    call updateInternalSoln(delDisp,elMatDim,internalDisp,intElasticLoad,intVecSize)
				endif
				nodeDisp(:) = nodeDisp(:) + delDisp(:)
				delUNorm = r_0
				do i2 = 1, elMatDim
				    mag = abs(delDisp(i2))
				    if(mag .gt. delUNorm) then
					    delUNorm = mag
					endif
				enddo
				i1 = i1 + 1
				write(lfUnit,*) 'solver iteration: ', i1, ', step magnitude: ', delUNorm
			enddo
		endif
	
	end subroutine stepSolve
	
	subroutine solve()
	    implicit none
		
		character(len=128) :: outfileName
		real*8 :: time, c1, c2, appLdFact
		integer :: i1, i2
		
		! ---
		integer :: i3, i4, i5
		real*8 :: x, z, P_6EI, M_EI
		! ---
		
		if(.not. allocated(currentRank)) then
		    call analysisPrep()
		endif
		
		write(lfUnit,*) 'finished analysisPrep'
		
		call loadInitialState()
		
		if(solveThermal .eq. 1) then
		    call getThermalSolnLoad(1)
			call scaleThermalMPC()
			call convertToLTri(thermMatLT,thermMatLTRange,thermMatLTSize,thermMat,thermMatCols,thermMatRange, &
			    thermMatSize,thermMatDim,thermMPCMat,thermMPCMatCols,thermMPCMatRange,thermMPCSize,thermMPCDim)
			call getSparseLDLFact(thermMatLT,thermMatLTSize,thermMatLTRange,thermMatDim)
		endif
		
		if(solveElastic .eq. 1) then
		    call getElasticSolnLoad(1)
			! --------------
			! P_6EI = 0.002d0
			! M_EI = 0.012d0
			! nodeDisp(:) = r_0
			! do i1 = 1, numNodes
			    ! x = nodeList(1,i1)
				! z = nodeList(3,i1)
				! i3 = nDofIndex(1,i1)
				! nodeDisp(i3) = -1.000035d0*P_6EI*(60d0*x - 3d0*x*x)*(z-0.05d0)
				! !nodeDisp(i3) = -M_EI*x*(z-0.05d0)
				! i3 = nDofIndex(3,i1)
				! nodeDisp(i3) = P_6EI*(30d0*x*x - x*x*x)
				! !nodeDisp(i3) = 0.5d0*M_EI*x*x
			! enddo
			
			! elasticLoad(:) = r_0
			! do i1 = 1, elMatDim
			    ! do i2 = elMatRange(i1-1)+1, elMatRange(i1)
				    ! i3 = elMatCols(i2)
					! if(i3 .ne. 0) then
					    ! elasticLoad(i1) = elasticLoad(i1) + elasticMat(i2)*nodeDisp(i3)
					! endif
				! enddo
			! enddo
			
			! write(lfUnit,*) 'reactionForce:'
			! do i1 = 1, numNodes
			    ! write(lfUnit,*) i1, elasticLoad(nDofIndex(1,i1)), elasticLoad(nDofIndex(2,i1)), elasticLoad(nDofIndex(3,i1))
			! enddo
			
			! nodeDisp(:) = r_0
			! elasticLoad(:) = r_0
			! --------------
			call scaleElasticMPC()
			!---------------
			! write(lfUnit,*) 'MPC Matrix:'
			! do i1 = 1, elMPCDim
			    ! write(lfUnit,*) 'equation', i1
				! do i2 = elMPCMatRange(i1-1) + 1, elMPCMatRange(i1)
				    ! i3 = elMPCMatCols(i2)
					! if(i3 .ne. 0) then
					    ! do i4 = 1, numNodes
						    ! do i5 = 1, 3
							    ! if(nDofIndex(i5,i4) .eq. i3) then
								    ! write(lfUnit,*) '    node', i4, 'dof', i5, 'value', elMPCMat(i2)
								! endif
							! enddo
						! enddo
					! endif
				! enddo
			! enddo
			!--------------------
			if(nLGeom .eq. 0) then
			    !! --------------------------
				! outfileName = 'sparseStiffness.txt'
				! call writeSparseMatrix(outfileName,elasticMat,elMatCols,elMatRange,elMatSize,elMatDim)
				! outfileName = 'elasticMPC.txt'
				! call writeSparseMatrix(outfileName,elMPCMat,elMPCMatCols,elMPCMatRange,elMPCSize,elMPCDim)
				!! ---------------------------------
				write(*,*) 'converting to LTri'
				call convertToLTri(elMatLT,elMatLTRange,elMatLTSize,elasticMat,elMatCols,elMatRange, &
					 elMatSize,elMatDim,elMPCMat,elMPCMatCols,elMPCMatRange,elMPCSize,elMPCDim)
				
				!! ------------------------
			    ! allocate(aFull(elMatDim,elMatDim))
			    ! call convertToFullPop(aFull,elasticMat,elMatCols,elMatRange,elMatSize,elMatDim,elMPCMat, &
				    ! elMPCMatCols,elMPCMatRange,elMPCSize,elMPCDim)
				! write(lfUnit,*) 'matrix discrepancies:'
				! do i1 = 1, elMatDim
				    ! i3 = i1 - (elMatLTRange(i1) - elMatLTRange(i1-1)) + 1
				    ! do i2 = elMatLTRange(i1-1) + 1, elMatLTRange(i1)
					    ! if(abs(elMatLT(i2) - aFull(i1,i3)) .gt. 0.00001d0) then
						    ! write(lfUnit,*) 'row ', i1, 'col', i3, 'LT:', elMatLT(i2), 'aFull:', aFull(i1,i3)
						! endif
						! i3 = i3 + 1
					! enddo
					! i3 = i1 - (elMatLTRange(i1) - elMatLTRange(i1-1)) + 1
					! do i2 = 1, i1
					    ! if(abs(aFull(i1,i2)) .gt. 0.0000001d0 .and. i2 .lt. i3) then
						    ! write(lfUnit,*) 'missing term, row', i1, 'col', i2, 'val', aFull(i1,i2)
						! endif
					! enddo
				! enddo
				! elasticLoad(:) = r_0
				! call getElasticAppliedLoad(loadTime)
				! call rFactorMat(aFull,elMatDim,elMatDim,1,elMatDim,1,elMatDim,0)
				! call solveRFactor(nodeDisp,aFull,elasticLoad,elMatDim,elMatDim,1,elMatDim,1,elMatDim,0)
				! write(lfUnit,*) 'full pop solution:'
				! do i1 = 1, numNodes
				    ! write(lfUnit,*) i1, nodeDisp(nDofIndex(1,i1)), nodeDisp(nDofIndex(2,i1)), nodeDisp(nDofIndex(3,i1))
				! enddo
				
				!! ------------------------
				!outfileName = 'elasticLTStiffess.txt'
				!call writeLowerTriMatrix(outfileName,elMatLT,elMatLTRange,elMatLTSize,elMatDim)
				!! -------------------------
				write(*,*) 'getting matrix factorization'
				call getSparseLDLFact(elMatLT,elMatLTSize,elMatLTRange,elMatDim)
				
				!! ------------------------
				!outfileName = 'elasticLTFactored.txt'
				!call writeLowerTriMatrix(outfileName,elMatLT,elMatLTRange,elMatLTSize,elMatDim)
				!! -------------------------
			endif
		endif
		
		write(*,*) 'beginning step solve'
		
		if(dynamic .eq. 0) then
		    if(nLGeom .eq. 0) then
		        call stepSolve(loadTime,1d0)
			else
			    do i1 = 1, ldRampSteps
				    appLdFact = -r_1*i1*(i1-2*ldRampSteps)/(ldRampSteps*ldRampSteps)
					call stepSolve(loadTime,appLdFact)
					if(intVecSize .gt. 0) then
						prevIntDisp(:) = internalDisp(:)
					endif
				enddo
			endif
		else
		    c1 = r_p5/nMBeta
			c2 = r_2/(delT*delT)
		    i1 = 0
			time = r_0
			if(writeSolnHist .eq. 1) then
				call writeBinarySolution(i1)
			endif
		    do while(time .lt. simPeriod)
			    time = time + delT
				i1 = i1 + 1
				call stepSolve(time,1d0)
				do i2 = 1, numNodes
				    nodeTdot(i2) = (r_1/nMGamma)*((r_1/delT)*(nodeTemp(i2) - prevTemp(i2)) - (r_1-nMGamma)*prevTdot(i2))
				enddo
				do i2 = 1, elMatDim
					nodeAcc(i2) = c1*(c2*(nodeDisp(i2)-prevDisp(i2)-delT*prevVel(i2)) - (r_1 - r_2*nMBeta)*prevAcc(i2))
					nodeVel(i2) = prevVel(i2) + delT*((r_1-nMGamma)*prevAcc(i2) + nMGamma*nodeAcc(i2))
				enddo
				prevTemp(:) = nodeTemp(:)
				prevTdot(:) = nodeTdot(:)
				prevDisp(:) = nodeDisp(:)
				if(intVecSize .gt. 0) then
				    prevIntDisp(:) = internalDisp(:)
				endif
				prevVel(:) = nodeVel(:)
				prevAcc(:) = nodeAcc(:)
				if(writeSolnHist .eq. 1) then
					call writeBinarySolution(i1)
				endif
			enddo
			numTSteps = i1
		endif
		
		! deallocate(aFull)
		
	end subroutine solve
	
	subroutine getCompleteObjective()
	    implicit none
		
		real*8 :: time
		integer :: i1, i2, i3, errFlag
		
		objVal = r_0
		
		if(dynamic .eq. 0) then
		    call getObjective(loadTime)
		else
		    i1 = numTSteps
		    call readBinarySolution(numTSteps,errFlag)
			if(errFlag .eq. 1) then
			    goto 896
			endif
			do i1 = numTSteps-1, 0, -1
			    nodeTemp(:) = prevTemp(:)
				nodeTdot(:) = prevTdot(:)
			    nodeDisp(:) = prevDisp(:)
				if(intVecSize .gt. 0) then
				    internalDisp(:) = prevIntDisp(:)
				endif
				nodeVel(:) = prevVel(:)
				nodeAcc(:) = prevAcc(:)
				call readBinarySolution(i1,errFlag)
				if(errFlag .eq. 1) then
				   goto 896
				endif
				time = delT*(i1+1)
				call getObjective(time)
			enddo
			return
896			write(lfUnit,*) 'Error: could not read the solution history file for time step ', i1
			write(lfUnit,*) 'Be sure to run a dynamic analysis with *solve with the saveSolnHistory option set to yes'
			write(lfUnit,*) 'if objective and sensitivities are desired for dynamic problems.'
		endif
		
	end subroutine getCompleteObjective
	
	subroutine augmentdLdu()
	    implicit none
		
		real*8 :: multFact
		integer :: i1
		
		multFact = 1e-20_8
		
		nodeTemp(:) = multFact*nodeTemp(:)
		nodeTdot(:) = multFact*nodeTdot(:)
		nodeDisp(:) = multFact*nodeDisp(:)
		nodeVel(:) = multFact*nodeVel(:)
		nodeAcc(:) = multFact*nodeAcc(:)
		prevTemp(:) = multFact*prevTemp(:)
		prevTdot(:) = multFact*prevTdot(:)
		prevDisp(:) = multFact*prevDisp(:)
		prevVel(:) = multFact*prevVel(:)
		prevAcc(:) = multFact*prevAcc(:)
		if(intVecSize .gt. 0) then
		    internalDisp(:) = multFact*internalDisp(:)
		    prevIntDisp(:) = multFact*prevIntDisp(:)
		endif
		
		if(solveThermal .eq. 1) then
		    !! Block 1
		    swapVec(1:numNodes) = prevTemp(:)
			prevTemp(:) = tempAdj(:)
			thermalLoad(:) = r_0
			call getThermalSolnLoad(0)
			dLdt(:) = dLdt(:) + thermalLoad(:) + tdotAdj(:)
			prevTemp(:) = swapVec(1:numNodes)
			
			!! Block 2
			swapVec(1:numNodes) = prevTdot(:)
			prevTdot(:) = tempAdj(:)
			thermalLoad(:) = r_0
			call getThermalSolnLoad(0)
			dLdtdot(:) = dLdtdot(:) + thermalLoad(:) + delT*(r_1-nMGamma)*tdotAdj(:)
			prevTdot(:) = swapVec(1:numNodes)
		endif
		
		if(solveElastic .eq. 1) then
		    !! Block 3
			swapVec(:) = prevDisp(:)
			prevDisp(:) = dispAdj(:)
			if(intVecSize .gt. 0) then
			    intswapVec(:) = prevIntDisp(:)
			    prevIntDisp(:) = intDispAdj(:)
			endif
			elasticLoad(:) = r_0
			intElasticLoad(:) = r_0
			call getElasticSolnLoad(0)
			dLdu(:) = dLdu(:) + elasticLoad(:) + accAdj(:)
			prevDisp(:) = swapVec(:)
			if(intVecSize .gt. 0) then
			    intdLdu(:) = intdLdu(:) + intElasticLoad(:)
			    prevIntDisp(:) = intswapVec(:)
			endif
			
			!! Block 4
			swapVec(:) = prevAcc(:)
			prevAcc(:) = dispAdj(:)
			elasticLoad(:) = r_0
			intElasticLoad(:) = r_0
			call getElasticSolnLoad(0)
			
			do i1 = 1, elMatDim
			    dLda(i1) = dLda(i1) + elasticLoad(i1) + r_p5*delT*delT*(r_1-r_2*nMBeta)*accAdj(i1) + delT*(r_1-nMGamma)*velAdj(i1)
			enddo
			prevAcc(:) = swapVec(:)
			
			!! Block 5
			swapVec(:) = prevVel(:)
			prevVel(:) = dispAdj(:)
			elasticLoad(:) = r_0
			intElasticLoad(:) = r_0
			call getElasticSolnLoad(0)
			dLdv(:) = dLdv(:) + elasticLoad(:) + delT*accAdj(:) + velAdj(:)
			prevVel(:) = swapVec(:)
		endif
		
		multFact = r_1/multFact
		
		nodeTemp(:) = multFact*nodeTemp(:)
		nodeTdot(:) = multFact*nodeTdot(:)
		nodeDisp(:) = multFact*nodeDisp(:)
		nodeVel(:) = multFact*nodeVel(:)
		nodeAcc(:) = multFact*nodeAcc(:)
		prevTemp(:) = multFact*prevTemp(:)
		prevTdot(:) = multFact*prevTdot(:)
		prevDisp(:) = multFact*prevDisp(:)
		prevVel(:) = multFact*prevVel(:)
		prevAcc(:) = multFact*prevAcc(:)
		internalDisp(:) = multFact*internalDisp(:)
		prevIntDisp(:) = multFact*prevIntDisp(:)
		
	end subroutine augmentdLdu
	
	subroutine solveDynamicAdjoint()
	    implicit none
		
		integer :: elNum
		integer :: numNds, dofPerNd, numIntDof, numIntPts
		integer :: dofTable(2,33)
		real*8 :: intPts(3,8), ipWt(8)
		real*8 :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), den, tExp(6), tCond(3,3), sHeat
		real*8 :: ABD(9,9), sMass(6,6), stExp(6), stCond(3,3), ssHeat
		real*8 :: bStiff(6,6), bMass(6,6), btExp(6), btCond(3,3), bsHeat
		real*8 :: temp(10), Tdot(10), pTemp(10), pTdot(10)
		real*8 :: disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)
		
		real*8 :: Ruk(33), dRdU(33,33), dRdT(33,10), stEn, dSEdT(10)
		
		integer :: i1, i2, i3, i4, i5, i6, i7, eType, numVecs
		
		if(solveElastic .eq. 1) then
		    if(dynamic .eq. 1) then
				velAdj(:) = dLdv(:)
				accAdj(:) = (-r_1/(delT*delT*nMBeta))*(dLda(:) + delT*nMGamma*velAdj(:))			
				dLdu(:) = dLdu(:) - accAdj(:)
			endif
			if(nLGeom .ne. 0) then
				call getElasticSolnLoad(1)
				call scaleElasticMPC()
				call convertToLTri(elMatLT,elMatLTRange,elMatLTSize,elasticMat,elMatCols,elMatRange, &
					 elMatSize,elMatDim,elMPCMat,elMPCMatCols,elMPCMatRange,elMPCSize,elMPCDim)
				call getSparseLDLFact(elMatLT,elMatLTSize,elMatLTRange,elMatDim)
			endif
			if(intVecSize .gt. 0) then
			    call updateExternalRHS(dLdu,elMatDim,intdLdu,intVecSize)
			endif
			numVecs = solverMaxBW
			if(solverMeth .eq. 'direct') then
			    call solveSparseLDLFact(dispAdj,dLdu,elMatLT,elMatLTSize,elMatLTRange,elMatDim,swapVec)
			else
				call gMRes(dispAdj,elasticMat,elMatCols,elMatRange,elMatSize,elMatDim,elMPCMat,elMPCMatCols, &
						elMPCMatRange,elMPCSize,elMPCDim,elMatLT,elMatLTRange,elMatLTSize,dLdu,numVecs,elMatDim)
			endif
			if(intVecSize .gt. 0) then
			    intDispAdj(:) = r_0
			    call updateInternalSoln(dispAdj,elMatDim,intDispAdj,intdLdu,intVecSize)
			endif
        endif

        if(solveThermal .eq. 1) then
		    if(dynamic .eq. 1) then
		        tdotAdj(:) = (-r_1/(delT*nMGamma))*dLdtdot(:)
			endif
			
			if(solveElastic .eq. 1) then
			    do i1 = 1, numEls
				    call r_getElementData(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
						locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, & 
						ABD, sMass, stExp, stCond, ssHeat, &
						bStiff, bMass, btExp, btCond, bsHeat, &
						temp,Tdot,disp,vel,acc,pTemp,pTdot,pDisp,pVel,pAcc,i1)
					eType = elementType(i1)
					call r_getElRuk(Ruk,dRdU,dRdT,stEn,dSEdT,1,eType,numNds,dofPerNd,numIntDof,numIntPts,dofTable, &
	                    intPts,ipWt,locNds,globNds,orient,cMat,tExp,ABD,stExp,bStiff,btExp,temp,disp)
					i2 = numNds*dofPerNd
					do i3 = 1, i2
					    i4 = dofTable(1,i3)
						i5 = elementList(dofTable(2,i3),i1)
						i6 = nDofIndex(i4,i5)
						do i4 = 1, numNds
						    i5 = currentRank(elementList(i4,i1))
						    dLdt(i5) = dLdt(i5) - dRdT(i3,i4)*dispAdj(i6)
						enddo
					enddo
					if(numIntDof .gt. 0) then
					    i2 = numNds*dofPerNd + 1
						i3 = numNds*dofPerNd + numIntDof
						i4 = intVecRange(i1-1) + 1
						do i5 = i2, i3
						    do i6 = 1, numNds
							    i7 = currentRank(elementList(i6,i1))
							    dLdt(i7) = dLdt(i7) - dRdT(i5,i6)*intDispAdj(i4)
							enddo
							i4 = i4 + 1
						enddo
					endif
				enddo
			endif
			
			if(dynamic .eq. 1) then
			    dLdt(:) = dLdt(:) - tdotAdj(:)
			endif
			
			call getThermalSolnLoad(1)
			call scaleThermalMPC()
			call convertToLTri(thermMatLT,thermMatLTRange,thermMatLTSize,thermMat,thermMatCols,thermMatRange, &
			    thermMatSize,thermMatDim,thermMPCMat,thermMPCMatCols,thermMPCMatRange,thermMPCSize,thermMPCDim)
			call getSparseLDLFact(thermMatLT,thermMatLTSize,thermMatLTRange,thermMatDim)
			numVecs = solverMaxBW
			if(solverMeth .eq. 'direct') then
			    call solveSparseLDLFact(tempAdj,dLdt,thermMatLT,thermMatLTSize,thermMatLTRange,numNodes,swapVec)
			else
				call gMRes(tempAdj(:),thermMat,thermMatCols,thermMatRange,thermMatSize,thermMatDim, &
					thermMPCMat,thermMPCMatCols,thermMPCMatRange,thermMPCSize,thermMPCDim,thermMatLT,thermMatLTRange, &
					thermMatLTSize,dLdt,numVecs,thermMatDim)
			endif
        endif		
		
	end subroutine solveDynamicAdjoint
	
	subroutine getdRthermaldD(time,dVarNum)
	    implicit none
		
		integer, intent(in) :: dVarNum
		real*8, intent(in) :: time
		
		integer :: eType
		integer :: numNds, dofPerNd, numIntDof, numIntPts
		integer :: dofTable(2,33)
		complex*16 :: intPts(3,8), ipWt(8)
		complex*16 :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), den, tExp(6), tCond(3,3), sHeat
		complex*16 :: ABD(9,9), sMass(6,6), sTELd(6), stCond(3,3), ssHeat
		complex*16 :: bStiff(6,6), bMass(6,6), bTELd(6), btCond(3,3), bsHeat
		complex*16 :: temp(10), Tdot(10), pTemp(10), pTdot(10)
		complex*16 :: disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)

		integer :: ndDof
		complex*16 :: ndFrc(33), bdyFrc(6)
		complex*16 :: trac(6), nDir(3), tTol
        character(len=16) :: lt		
		
		complex*16 :: Rt(10), dRdt(10,10), ndLd(7)
		integer :: i1, i2, i3, i4, i5, i6, i7, i8, i9
		
		dRtdD(:) = r_0
		
	    elDepends(:) = 0
		do i1 = dToElCRange(dVarNum-1)+1, dToElCRange(dVarNum)
		    i2 = dToElComp(i1)
			if(i2 .ne. 0) then
			    elDepends(i2) = 1
			endif
		enddo
		
		c_dVec(:) = c_1*r_dVec(:)
		c_dVec(dVarNum) = c_dVec(dVarNum) + compStep
		
		do i1 = dToElCRange(dVarNum-1) + 1, dToElCRange(dVarNum)
		    i2 = dToElComp(i1)
			call c_getElementData(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
				locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, & 
				ABD, sMass, sTELd, stCond, ssHeat, &
				bStiff, bMass, bTELd, btCond, bsHeat, &
	            temp,Tdot,disp,vel,acc,pTemp,pTdot,pDisp,pVel,pAcc,i2)
			eType = elementType(i2)
			call c_getElRt(Rt,dRdt,0,numNds,numIntPts,intPts,ipWt, &
				locNds,orient,tCond,stCond,btCond,den,sHeat,ssHeat,bsHeat, &
				temp,pTemp,pTDot,eType)
			do i3 = 1, numNds
			    i4 = currentRank(elementList(i3,i2))
				dRtdD(i4) = dRtdD(i4) + compStepInv*imag(Rt(i3))
			enddo
		enddo
		
		do i1 = 1, numLds
		    if(time .ge. loadsActTime(1,i1) .and. time .le. loadsActTime(2,i1)) then
			    lt = loadType(i1)
				if(lt .eq. 'bodyHeatGen') then
					do i2 = loadsRange(i1-1) + 1, loadsRange(i1)
						read(loadNodes(i2),*,err=1082) i3
						if(elDepends(i3) .eq. 1) then
							bdyFrc(1) = c_1*inputLoads(1,i2)
							bdyFrc(2:6) = c_0
							call c_getElBodyForce(ndFrc,dofTable,ndDof,bdyFrc,0,i3)
							do i4 = 1, ndDof
								i5 = dofTable(1,i4)
								if(i5 .eq. 1) then
									i6 = currentRank(elementList(dofTable(2,i4),i3))
									dRtdD(i6) = dRtdD(i6) - compStepInv*imag(ndFrc(i4))
								endif
							enddo
						endif
						goto 1099
1082                    do i4 = 1, numElSets
							if(elSetName(i4) .eq. loadNodes(i2)) then
								do i5 = elSetRange(i4-1)+1, elSetRange(i4)
									i6 = elementSets(i5)
									if(elDepends(i6) .eq. 1) then
										bdyFrc(1) = c_1*inputLoads(1,i2)
										bdyFrc(2:6) = c_0
										call c_getElBodyForce(ndFrc,dofTable,ndDof,bdyFrc,0,i6)
										do i7 = 1, ndDof
											i8 = dofTable(1,i7)
											if(i8 .eq. 1) then
												i9 = currentRank(elementList(dofTable(2,i7),i6))
												dRtdD(i9) = dRtdD(i9) - compStepInv*imag(ndFrc(i7))
											endif
										enddo
									endif
								enddo
							endif
						enddo
1099                    i4 = 1
					enddo
                elseif(lt .eq. 'surfaceFlux') then
					do i2 = loadsRange(i1-1) + 1, loadsRange(i1)
						read(loadNodes(i2),*,err=1121) i3
						if(elDepends(i3) .eq. 1) then
							trac(1) = c_1*inputLoads(5,i2)
							trac(2:6) = c_0
							nDir(:) = c_1*inputLoads(1:3,i2)
							tTol = c_1*inputLoads(4,i2)
							call c_getElSurfaceTraction(ndFrc,dofTable,ndDof,trac,c_0,i3,nDir,tTol)
							do i4 = 1, ndDof
								i5 = dofTable(1,i4)
								if(i5 .eq. 1) then
									i6 = currentRank(elementList(dofTable(2,i4),i3))
									dRtdD(i6) = dRtdD(i6) + compStepInv*imag(ndFrc(i4))
								endif
							enddo
						endif
						goto 1144
1121                    do i4 = 1, numElSets
							if(elSetName(i4) .eq. loadNodes(i2)) then
								do i5 = elSetRange(i4-1)+1, elSetRange(i4)
									i6 = elementSets(i5)
									if(elDepends(i6) .eq. 1) then
										trac(1) = c_1*inputLoads(5,i2)
										trac(2:6) = c_0
										nDir(:) = c_1*inputLoads(1:3,i2)
										tTol = c_1*inputLoads(4,i2)
										call c_getElSurfaceTraction(ndFrc,dofTable,ndDof,trac,c_0,i6,nDir,tTol)
										do i7 = 1, ndDof
											i8 = dofTable(1,i7)
											if(i8 .eq. 1) then
												i9 = currentRank(elementList(dofTable(2,i7),i6))
												dRtdD(i9) = dRtdD(i9) + compStepInv*imag(ndFrc(i7))
											endif
										enddo
									endif
								enddo
							endif
						enddo
1144                    i4 = 1
					enddo				
				endif
			endif
		enddo	
		
		do i1 = dToNdRange(dVarNum-1) + 1, dToNdRange(dVarNum)
		    i2 = dToNd(i1)
		    call c_getNodeDLoad(c_1*time,ndLd,i2)
			i3 = currentRank(i2)
			dRtdD(i3) = dRtdD(i3) - compStepInv*imag(ndLd(7))
		enddo
		
	end subroutine getdRthermaldD
	
	subroutine getdRelasticdD(time,dVarNum)
	    implicit none
		
		integer, intent(in) :: dVarNum
		real*8, intent(in) :: time
		
		integer :: eType
		integer :: numNds, dofPerNd, numIntDof, numIntPts
		integer :: dofTable(2,33)
		complex*16 :: intPts(3,8), ipWt(8)
		complex*16 :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), den, tExp(6), tCond(3,3), sHeat
		complex*16 :: ABD(9,9), sMass(6,6), stExp(6), stCond(3,3), ssHeat
		complex*16 :: bStiff(6,6), bMass(6,6), btExp(6), btCond(3,3), bsHeat
		complex*16 :: temp(10), Tdot(10), pTemp(10), pTdot(10)
		complex*16 :: disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)

  		integer :: isGrav, ndDof
		complex*16 :: ndFrc(33), bdyFrc(6)
		complex*16 :: trac(6), pressure, nDir(3), tTol
		complex*16 :: center(3), axis(3), angVel
		character(len=16) :: lt
		
		complex*16 :: Ru(33), dRdU(33,33), ndLd(7)
		complex*16 :: nnInv, statInOri(3,240)
		integer :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10
		
	    dRudD(:) = r_0
		if(intVecSize .gt. 0) then
	        intdRudD(:) = r_0
		endif
		elDepends(:) = 0
		do i1 = dToElCRange(dVarNum-1)+1, dToElCRange(dVarNum)
		    i2 = dToElComp(i1)
			if(i2 .ne. 0) then
			    elDepends(i2) = 1
			endif
		enddo
		c_dVec(:) = c_1*r_dVec(:)
		c_dVec(dVarNum) = c_dVec(dVarNum) + compStep
		
		do i1 = dToElCRange(dVarNum-1) + 1, dToElCRange(dVarNum)
		    i2 = dToElComp(i1)
			call c_getElementData(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
				locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, & 
				ABD, sMass, stExp, stCond, ssHeat, &
				bStiff, bMass, btExp, btCond, bsHeat, &
	            temp,Tdot,disp,vel,acc,pTemp,pTdot,pDisp,pVel,pAcc,i2)
			eType = elementType(i2)
			if(eType .eq. 41 .or. eType .eq. 3 .or. eType .eq. 2) then
			    call c_getInstOrient(statInOri,orient,disp,numNds,1)
			endif
			!! -------------------
			! if(i2 .eq. 1) then
			    ! write(lfUnit,*) 'dVarNum: ', dVarNum
			    ! write(lfUnit,*) 'element ', i2, ' data'
				! write(lfUnit,*) 'numNds: ', numNds, 'dofPerNd: ', dofPerNd, 'numIntDof: ', numIntDof, 'numIntPts: ', numIntPts
				! write(lfUnit,*) 'intPts: '
				! write(lfUnit,*) intPts
				! write(lfUnit,*) 'ipWt'
				! write(lfUnit,*) ipWt
			    ! write(lfUnit,*) 'locNds: '
				! write(lfUnit,*) locNds
				! write(lfUnit,*) 'globNds: '
				! write(lfUnit,*) globNds
				! write(lfUnit,*) 'orient: '
				! write(lfUnit,*) orient
				! write(lfUnit,*) 'cMat: '
				! write(lfUnit,*) cMat
				! write(lfUnit,*) 'ABD: '
				! write(lfUnit,*) ABD
				! write(lfUnit,*) 'stExp: '
				! write(lfUnit,*) stExp
				! write(lfUnit,*) 'temp: '
				! write(lfUnit,*) temp
				! write(lfUnit,*) 'disp: '
				! write(lfUnit,*) disp
			! endif
			!! ----------------------------
			call c_getElRu(Ru,dRdU,0,eType,numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
				locNds,globNds,orient,cMat,tExp,ABD,stExp,bStiff,btExp,temp,disp,vel,acc,pDisp,pVel,pAcc, &
				statInOri,den,sMass,bMass)
			i3 = numNds*dofPerNd
			do i4 = 1, i3
			    i5 = dofTable(1,i4)
				i6 = elementList(dofTable(2,i4),i2)
			    i7 = nDofIndex(i5,i6)
				dRudD(i7) = dRudD(i7) + compStepInv*imag(Ru(i4))
			enddo
			if(numIntDof .gt. 0) then
			    i3 = numNds*dofPerNd + 1
				i4 = numNds*dofPerNd + numIntDof
				i5 = intVecRange(i2-1) + 1
				do i6 = i3, i4
				    intdRudD(i5) = intdRudD(i5) + compStepInv*imag(Ru(i6))
					i5 = i5 + 1
				enddo
			endif
		enddo
		
		do i1 = 1, numLds
		    if(time .ge. loadsActTime(1,i1) .and. time .le. loadsActTime(2,i1)) then
			    lt = loadType(i1)
				if(lt .eq. 'bodyForce' .or. lt .eq. 'gravitational' .or. lt .eq. 'centrifugal') then
				    if(lt .eq. 'bodyForce') then
					    isGrav = 0
				    else
					    isGrav = 1
					endif
					do i2 = loadsRange(i1-1) + 1, loadsRange(i1)
						read(loadNodes(i2),*,err=1087) i3
						if(elDepends(i3) .eq. 1) then
							if(lt .eq. 'centrifugal') then
								center(:) = c_1*inputLoads(1:3,i2)
								axis(:) = c_1*inputLoads(4:6,i2)
								angVel = c_1*inputLoads(7,i2)
								call c_getElCentAcc(bdyFrc(1:3),center,axis,angVel,i3)
								bdyFrc(4:6) = c_0
							else
								bdyFrc(:) = c_1*inputLoads(1:6,i2)
							endif
							call c_getElBodyForce(ndFrc,dofTable,ndDof,bdyFrc,isGrav,i3)
							do i4 = 1, ndDof
								i5 = dofTable(1,i4)
								i6 = elementList(dofTable(2,i4),i3)
								i7 = nDofIndex(i5,i6)
								dRudD(i7) = dRudD(i7) - compStepInv*imag(ndFrc(i4))
							enddo
						endif
						goto 1110
1087                    do i4 = 1, numElSets
							if(elSetName(i4) .eq. loadNodes(i2)) then
								do i5 = elSetRange(i4-1)+1, elSetRange(i4)
									i6 = elementSets(i5)
									if(elDepends(i6) .eq. 1) then
										if(lt .eq. 'centrifugal') then
											center(:) = c_1*inputLoads(1:3,i2)
											axis(:) = c_1*inputLoads(4:6,i2)
											angVel = c_1*inputLoads(7,i2)
											call c_getElCentAcc(bdyFrc(1:3),center,axis,angVel,i6)
											bdyFrc(4:6) = c_0
										else
											bdyFrc(:) = c_1*inputLoads(1:6,i2)
										endif
										call c_getElBodyForce(ndFrc,dofTable,ndDof,bdyFrc,isGrav,i6)
										do i7 = 1, ndDof
											i8 = dofTable(1,i7)
											i9 = elementList(dofTable(2,i7),i6)
											i10 = nDofIndex(i8,i9)
											dRudD(i10) = dRudD(i10) - compStepInv*imag(ndFrc(i7))
										enddo
									endif
								enddo
							endif
						enddo
1110                    i4 = 1
					enddo
                elseif(lt .eq. 'surfacePressure' .or. lt .eq. 'surfaceTraction') then
					do i2 = loadsRange(i1-1) + 1, loadsRange(i1)
					    if(lt .eq. 'surfacePressure') then
						    pressure = c_1*inputLoads(5,i2)
							trac(:) = c_0
						else
						    pressure = c_0
							trac(1:3) = c_1*inputLoads(5:7,i2)
							trac(4:6) = c_0
						endif
						read(loadNodes(i2),*,err=1137) i3
						if(elDepends(i3) .eq. 1) then
						    nDir(:) = c_1*inputLoads(1:3,i2)
							tTol = c_1*inputLoads(4,i2)
							call c_getElSurfaceTraction(ndFrc,dofTable,ndDof,trac,pressure,i3,nDir,tTol)
							do i4 = 1, ndDof
								i5 = dofTable(1,i4)
								i6 = elementList(dofTable(2,i4),i3)
								i7 = nDofIndex(i5,i6)
								dRudD(i7) = dRudD(i7) - compStepInv*imag(ndFrc(i4))
							enddo
						endif
						goto 1157
1137                    do i4 = 1, numElSets
							if(elSetName(i4) .eq. loadNodes(i2)) then
								do i5 = elSetRange(i4-1)+1, elSetRange(i4)
									i6 = elementSets(i5)
									if(elDepends(i6) .eq. 1) then
									    nDir(:) = c_1*inputLoads(1:3,i2)
										tTol = c_1*inputLoads(4,i2)
										call c_getElSurfaceTraction(ndFrc,dofTable,ndDof,trac,pressure,i6,nDir,tTol)
										do i7 = 1, ndDof
											i8 = dofTable(1,i7)
											i9 = elementList(dofTable(2,i7),i6)
											i10 = nDofIndex(i8,i9)
											dRudD(i10) = dRudD(i10) - compStepInv*imag(ndFrc(i7))
										enddo
									endif
								enddo
							endif
						enddo
1157                    i4 = 1
					enddo			
				endif
			endif
		enddo
		
		do i1 = dToNdRange(dVarNum-1) + 1, dToNdRange(dVarNum)
		    i2 = dToNd(i1)
		    call c_getNodeDLoad(c_1*time,ndLd,i2)
			do i3 = 1, 6
			    i4 = nDofIndex(i3,i2)
				if(i4 .gt. 0) then
			        dRudD(i4) = dRudD(i4) - compStepInv*imag(ndLd(i3))
				endif
			enddo
		enddo
		
		c_dVec(dVarNum) = c_dVec(dVarNum) - compStep
		
	end subroutine getdRelasticdD
	
	subroutine getTotaldLdD()
	    implicit none
		
		real*8 :: time
		integer :: i1, i2, i3, i4, i5, i6, errFlag
		
		if(dynamic .eq. 0) then
		    dLdt(:) = r_0
			dLdtdot(:) = r_0
			dLdu(:) = r_0
			if(intVecSize .gt. 0) then
			    intdLdu(:) = r_0
			endif
			dLdv(:) = r_0
			dLda(:) = r_0
			call getdObjdU(loadTime)
			call solveDynamicAdjoint()
			dLdD(:) = r_0
			call getdObjdD(loadTime)
			do i1 = 1, numDVar
			    if(solveThermal .eq. 1) then
				    call getdRthermaldD(loadTime,i1)
					do i2 = 1, numNodes
					    dLdD(i1) = dLdD(i1) - tempAdj(i2)*dRtdD(i2)
					enddo
				endif
				if(solveElastic .eq. 1) then
				    call getdRelasticdD(loadTime,i1)
					do i2 = 1, elMatDim
					    dLdD(i1) = dLdD(i1) - dispAdj(i2)*dRudD(i2)
					enddo
					if(intVecSize .gt. 0) then
					    do i2 = 1, intVecSize
						    dLdD(i1) = dLdD(i1) - intDispAdj(i2)*intdRudD(i2)
						enddo
					endif
				endif
			enddo
		else
		    dLdD(:) = r_0
		    call readBinarySolution(numTSteps,errFlag)
			time = delT*numTSteps
			do i1 = numTSteps, 1, -1
			    dLdt(:) = r_0
				dLdtdot(:) = r_0
				dLdu(:) = r_0
				if(intVecSize .gt. 0) then
					intdLdu(:) = r_0
				endif
				dLdv(:) = r_0
				dLda(:) = r_0
				if(i1 .lt. numTSteps) then
				    call augmentdLdu()
				endif
				nodeTemp(:) = prevTemp(:)
				nodeTdot(:) = prevTdot(:)
				nodeDisp(:) = prevDisp(:)
				nodeVel(:) = prevVel(:)
				nodeAcc(:) = prevAcc(:)
				if(intVecSize .gt. 0) then
					internalDisp(:) = prevIntDisp(:)
				endif
				call readBinarySolution(i1-1,errFlag)
				call getdObjdU(time)
				call solveDynamicAdjoint()
				call getdObjdD(time)
				do i2 = 1, numDVar
					if(solveThermal .eq. 1) then
						call getdRthermaldD(time,i2)
						do i3 = 1, numNodes
							dLdD(i2) = dLdD(i2) - tempAdj(i3)*dRtdD(i3)
						enddo
					endif
					if(solveElastic .eq. 1) then
						call getdRelasticdD(time,i2)
						do i3 = 1, elMatDim
							dLdD(i2) = dLdD(i2) - dispAdj(i3)*dRudD(i3)
						enddo
						if(intVecSize .gt. 0) then
							do i3 = 1, intVecSize
								dLdD(i2) = dLdD(i2) - intDispAdj(i3)*intdRudD(i3)
							enddo
						endif
					endif
				enddo
				time = time - delT
			enddo
		endif
		
	end subroutine getTotaldLdD
	
	subroutine processWriteNodeEl(jobUnit,fileLine,iosVal)
	    implicit none
		
		integer, intent(in) :: jobUnit
		character(len=256), intent(out) :: fileLine
		integer, intent(out) :: iosVal
		
		integer, allocatable :: allNdSet(:), allElSet(:), writeTSteps(:)
		
		character(len=32) :: commandTag
		character(len=128) :: outputFileName, fullFileName
		
		character(len=16) :: writeFields(10)
	    character(len=16) :: extension, stepStr
		real*8 :: time
		
		integer :: nSInd, eSInd, numFields, readInt(3)
		integer :: i1, i2, i3, i4, i5, i6, i7, errFlag
		
		i1 = index(fileLine,'*')
		commandTag = fileLine(i1:i1+31)
		read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
		i1 = index(fileLine,'*')
		numFields = 0
		nSInd = 0
		eSInd = 0
		do while(i1 .eq. 0 .and. iosVal .eq. 0)
			i2 = index(fileLine,':')
			if(i2 .gt. 0) then
				if(fileLine(i2-8:i2) .eq. 'fileName:') then
					read(fileLine(i2+1:i2+128),'(A)') outputFileName
					outputFileName = adjustl(outputFileName)
					outputFileName = trim(outputFileName)
					read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
					i1 = index(fileLine,'*')
				elseif(fileLine(i2-7:i2) .eq. 'nodeSet:') then
					do i3 = 1, numNdSets
						i4 = index(fileLine,ndSetName(i3))
						if(i4 .gt. 0) then
							nSInd = i3
						endif
					enddo
					if(nSInd .eq. 0) then
						write(lfUnit,*) 'Warning: no node set named ', fileLine(i2+1:i2+64), 'was found'
						write(lfUnit,*) 'Writing results for all nodes'
					endif
					read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
					i1 = index(fileLine,'*')
				elseif(fileLine(i2-10:i2) .eq. 'elementSet:') then
				    do i3 = 1, numElSets
						i4 = index(fileLine,elSetName(i3))
						if(i4 .gt. 0) then
							eSInd = i3
						endif
					enddo
					if(eSInd .eq. 0) then
						write(lfUnit,*) 'Warning: no element set named ', fileLine(i2+1:i2+64), 'was found'
						write(lfUnit,*) 'Writing results for all elements'
					endif
					read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
					i1 = index(fileLine,'*')
				elseif(fileLine(i2-6:i2) .eq. 'fields:') then
					read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
					i2 = index(fileLine,':')
					do while(i2 .eq. 0 .and. iosVal .eq. 0)
						i3 = index(fileLine,'-')
						if(i3 .gt. 0) then
							numFields = numFields + 1
							read(fileLine(i3+1:i3+64),*) writeFields(numFields)
							if(writeFields(numFields) .eq. 'reactionForce') then
								elasticLoad(:) = r_0
								if(intVecSize .gt. 0) then
									intElasticLoad(:) = r_0
								endif
								call getElasticSolnLoad(0)
							endif
							if(writeFields(numFields) .eq. 'reactionHeatGen') then
								thermalLoad(:) = r_0
								call getThermalSolnLoad(0)
							endif
						endif
						read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
						i2 = index(fileLine,':')
					enddo
					i1 = index(fileLine,'*')
				elseif(fileLine(i2-9:i2) .eq. 'timeSteps:') then
					if(numTSteps .le. 0) then
						write(lfUnit,*) 'Error: A dynamic analysis must be run prior to writing time step dependent results.'
						write(lfUnit,*) 'Aborting ', commandTag
						goto 251
					endif
					if(.not. allocated(writeTSteps)) then
						allocate(writeTSteps(numTSteps))
					endif
					writeTSteps(:) = 0
					i3 = index(fileLine,'all')
					if(i3 .gt. 0) then
						writeTSteps(:) = 1
						read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
						i1 = index(fileLine,'*')
					else
						read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
						i2 = index(fileLine,':')
						do while(i2 .eq. 0 .and. iosVal .eq. 0)
							i3 = index(fileLine, '-')
							if(i3 .gt. 0) then
								i4 = index(fileLine,'[')
								if(i4 .gt. 0) then
									i5 = index(fileLine,']')
									readInt(3) = 1
									read(fileLine(i4+1:i5-1),*,end=254) readInt(1:3)
254                                 do i6 = readInt(1), readInt(2), readInt(3)
										if(i6 .le. numTSteps .and. i6 .gt. 0) then
											writeTSteps(i6) = 1
										else
											write(lfUnit,*) 'Error: a time step specified for writeNodeResults is out of'
											write(lfUnit,*) 'the simulated range.  Aborting ', commandTag
											goto 251
										endif
									enddo
								else
									read(fileLine(i3+1:i3+64),*) i5
									if(i5 .le. numTSteps .and. i5 .gt. 0) then
										writeTSteps(i5) = 1
									else
										write(lfUnit,*) 'Error: a time step specified for writeNodeResults is out of'
										write(lfUnit,*) 'the simulated range.  Aborting ', commandTag
										goto 251
									endif
								endif
							endif
							read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
							i1 = index(fileLine,'*')
							i2 = index(fileLine,':')
						enddo
					endif
				else
					read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
					i1 = index(fileLine,'*')
				endif
			endif
		enddo
		if(.not. allocated(writeTSteps)) then
		    if(commandTag(1:17) .eq. '*writeNodeResults') then
				if(nSInd .eq. 0) then
					if(.not. allocated(allNdSet)) then
						allocate(allNdSet(numNodes))
					endif
					do i3 = 1, numNodes
						allNdSet(i3) = i3
					enddo
					time = numTSteps*delT
					call writeNodeResults(outputFileName,writeFields,numFields,allNdSet,numNodes,time,numTSteps)
				else
					i3 = ndSetRange(nSInd-1) + 1
					i4 = ndSetRange(nSInd)
					i5 = i4 - i3 + 1
					time = numTSteps*delT
					call writeNodeResults(outputFileName,writeFields,numFields,nodeSets(i3:i4),i5,time,numTSteps)
				endif
			else
				if(eSInd .eq. 0) then
					if(.not. allocated(allElSet)) then
						allocate(allElSet(numEls))
					endif
					do i3 = 1, numEls
						allElSet(i3) = i3
					enddo
					time = numTSteps*delT
					call writeElementResults(outputFileName,writeFields,numFields,allElSet,numEls,time,numTSteps)
				else
					i3 = elSetRange(eSInd-1) + 1
					i4 = elSetRange(eSInd)
					i5 = i4 - i3 + 1
					time = numTSteps*delT
					call writeElementResults(outputFileName,writeFields,numFields,elementSets(i3:i4),i5,time,numTSteps)
				endif
			endif
		else
			do i6 = 1, numTSteps
				if(writeTSteps(i6) .eq. 1) then
					call readBinarySolution(i6,errFlag)
					if(errFlag .eq. 0) then
						nodeTemp = prevTemp(:)
						nodeTdot = prevTdot(:)
						nodeDisp = prevDisp(:)
						if(intVecSize .ge. 1) then
							internalDisp(:) = prevIntDisp(:)
						endif
						nodeVel(:) = prevVel(:)
						nodeAcc(:) = prevAcc(:)
						write(stepStr,'(I0)') i6
						i7 = index(outputFileName,'.')
						if(i7 .gt. 0) then
							read(outputFileName(i7:i7+15),*) extension
							read(outputFileName(1:i7-1),*) fullFileName
							outputFileName = fullFileName
						else
							extension = ''
						endif
						fullFileName = outputFileName // '_step' // stepStr // extension
						if(commandTag(1:17) .eq. '*writeNodeResults') then
							if(nSInd .eq. 0) then
								if(.not. allocated(allNdSet)) then
									allocate(allNdSet(numNodes))
								endif
								do i3 = 1, numNodes
									allNdSet(i3) = i3
								enddo
								time = numTSteps*delT
								call writeNodeResults(fullFileName,writeFields,numFields,allNdSet,numNodes,time,numTSteps)
							else
								i3 = ndSetRange(nSInd-1) + 1
								i4 = ndSetRange(nSInd)
								i5 = i4 - i3 + 1
								time = numTSteps*delT
								call writeNodeResults(fullFileName,writeFields,numFields,nodeSets(i3:i4),i5,time,numTSteps)
							endif
						else
							if(eSInd .eq. 0) then
								if(.not. allocated(allElSet)) then
									allocate(allElSet(numEls))
								endif
								do i3 = 1, numEls
									allElSet(i3) = i3
								enddo
								time = numTSteps*delT
								call writeElementResults(fullFileName,writeFields,numFields,allElSet,numEls,time,numTSteps)
							else
								i3 = elSetRange(eSInd-1) + 1
								i4 = elSetRange(eSInd)
								i5 = i4 - i3 + 1
								time = numTSteps*delT
								call writeElementResults(fullFileName,writeFields,numFields,elementSets(i3:i4),i5,time,numTSteps)
							endif
						endif
					else
						write(lfUnit,*) 'Error: time step ', i6, ' was never recorded in the solution history.'
						write(lfUnit,*) 'Be sure to set the saveSolnHist option to yes under the *solve command'
						write(lfUnit,*) 'if the objective gradient or time-specific intermediate results are desired'
						write(lfUnit,*) 'for a dynamic analysis.'
					endif
				endif
			enddo
			deallocate(writeTSteps)
		endif
		if(allocated(allNdSet)) then
			deallocate(allNdSet)
		endif
		if(allocated(allElSet)) then
		    deallocate(allElSet)
		endif
251		i1 = index(fileLine,'*')
		
	end subroutine processWriteNodeEl
	
	subroutine processWriteNdElProps(jobUnit,fileLine,iosVal)
	    implicit none
		
		integer, intent(in) :: jobUnit
		character(len=256), intent(out) :: fileLine
		integer, intent(out) :: iosVal
		
		integer, allocatable :: allNdSet(:), allElSet(:)
		character(len=16) :: propList(32)
		character(len=128) :: fileName
		character(len=64) :: setName
		integer :: i1, i2, i3, i4, i5, nSInd, eSInd, numProps, wtNd
		
		i1 = index(fileLine,'*')
		if(fileLine(i1:i1+20) .eq. '*writeNodeCoordinates') then
		    wtNd = 1
		else
		    wtNd = 0
		endif
		
		numProps = 0
		nSInd = 0
		eSInd = 0
		
		read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
		i1 = index(fileLine,'*')
		do while(i1 .eq. 0 .and. iosVal .eq. 0)
		    i2 = index(fileLine,':')
			if(i2 .gt. 0) then
			    if(fileLine(i2-8:i2) .eq. 'fileName:') then
				    read(fileLine(i2+1:i2+128),'(A)') fileName
					fileName = adjustl(fileName)
					fileName = trim(fileName)
					read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
		            i1 = index(fileLine,'*')
				elseif(fileLine(i2-7:i2) .eq. 'nodeSet:') then
				    read(fileLine(i2+1:i2+64),*) setName
					do i3 = 1, numNdSets
					    if(setName .eq. ndSetName(i3)) then
						    nSInd = i3
						endif
					enddo
					if(nSInd .eq. 0) then
					    write(lfUnit,*) 'Warning: there is no node set with the name ', setName
						write(lfUnit,*) 'Writing coordinates of all nodes.'
					endif
					read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
		            i1 = index(fileLine,'*')
				elseif(fileLine(i2-10:i2) .eq. 'elementSet:') then
				    read(fileLine(i2+1:i2+64),*) setName
					do i3 = 1, numElSets
					    if(setName .eq. elSetName(i3)) then
						    eSInd = i3
						endif
					enddo
					if(eSInd .eq. 0) then
					    write(lfUnit,*) 'Warning: there is no element set with the name ', setName
						write(lfUnit,*) 'Writing properties of all elements.'
					endif
					read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
		            i1 = index(fileLine,'*')
				elseif(fileLine(i2-10:i2) .eq. 'properties:') then
				    read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
		            i1 = index(fileLine,'*')
					i2 = index(fileLine,':')
					do while(i2 .eq. 0 .and. iosVal .eq. 0)
					    i3 = index(fileLine,'-')
						if(i3 .gt. 0) then
						    numProps = numProps + 1
							read(fileLine(i3+1:i3+16),*) propList(numProps)
						endif
						read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
						i2 = index(fileLine,':')
					enddo
					i1 = index(fileLine,'*')
				else
				    read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
		            i1 = index(fileLine,'*')
				endif
			else
			    read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
		        i1 = index(fileLine,'*')
			endif
		enddo
		if(wtNd .eq. 1) then
		    if(nSInd .eq. 0) then
			    allocate(allNdSet(numNodes))
				do i3 = 1, numNodes
				    allNdSet(i3) = i3
				enddo
				call writeNodeCoord(fileName,allNdSet,numNodes)
				deallocate(allNdSet)
			else
			    i3 = ndSetRange(nSInd-1) + 1
				i4 = ndSetRange(nSInd)
				i5 = i4 - i3 + 1
				call writeNodeCoord(fileName,nodeSets(i3:i4),i5)
			endif
		else
		    if(eSInd .eq. 0) then
			    allocate(allElSet(numEls))
				do i3 = 1, numEls
				    allElSet(i3) = i3
				enddo
				call writeElProperties(fileName,allElSet,numEls,propList,numProps)
				deallocate(allElSet)
			else
			    i3 = elSetRange(eSInd-1) + 1
				i4 = elSetRange(eSInd)
				i5 = i4 - i3 + 1
				call writeElProperties(fileName,elementSets(i3:i4),i5,propList,numProps)
			endif
		endif
		
	end subroutine processWriteNdElProps
	
end module ASenD_commandFunctions