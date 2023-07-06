module ASenD_objective
    use ASenD_globalData
	use ASenD_constantVals
	use ASenD_r_elementEqns
	use ASenD_c_elementEqns

    contains
	
	subroutine getObjective(time)
	    implicit none
		
		real*8, intent(in) :: time
		
		real*8 :: c1, c2, constTgt, sumTerm, sPt(3)
		real*8 :: stress(6), strain(6), mass, vol, totVol
		real*8 :: frcMom(9), shDef(9), ptT, bfrcMom(6), bDef(6)
		real*8 :: stEn, dSEdU(33), dSEdT(10), flux(3), dfluxdT(3,10), Tx(3), dTxdT(3,10)
		real*8 :: disp(6), ddispdU(6,33)
		integer :: i1, i2, i3, i4, i5, i6, ndSetNum, elSetNum, tgtListLen, layer, comp
		
		if(objElSetId(1) .eq. 0 .and. objNdSetId(1) .eq. 0) then		
			do i1 = 1, numObjTerms
				do i2 = 1, numNdSets
					if(ndSetName(i2) .eq. objNdSet(i1)) then
						objNdSetId(i1) = i2
					endif
				enddo
				do i2 = 1, numElSets
					if(elSetName(i2) .eq. objElSet(i1)) then
						objElSetId(i1) = i2
					endif
				enddo
			enddo
		endif
		
		if(.not. allocated(objTgtExpanded)) then
		    i1 = numNodes + numEls
		    allocate(objTgtExpanded(i1))
		endif
		
		do i1 = 1, numObjTerms
		    if(time .ge. objActTime(1,i1) .and. time .le. objActTime(2,i1)) then
			    if(objNdSetId(i1) .ne. 0) then
				    i3 = objNdSetId(i1)
				    i2 = ndSetRange(i3) - ndSetRange(i3-1)
				elseif(objElSetId(i1) .ne. 0) then
				    i3 = objElSetId(i1)
					i2 = elSetRange(i3) - elSetRange(i3-1)
				endif
				
				if(objTgtTag(i1) .eq. 'none') then
					objTgtExpanded(1:i2) = r_0
				elseif(objTgtTag(i1) .eq. 'single') then
					i3 = objTgtRange(i1)
					objTgtExpanded(1:i2) = objTgtVal(i3)
				elseif(objTgtTag(i1) .eq. 'vector') then
					i3 = objTgtRange(i1-1) + 1
					i4 = objTgtRange(i1)
					objTgtExpanded(1:i2) = objTgtVal(i3:i4)
				elseif(objTgtTag(i1) .eq. 'maxTensile') then
					if(objCategory(i1) .eq. 'stress') then
						elSetNum = objElSetId(i1)
						i5 = 0
						do i3 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
							i5 = i5 + 1
							i4 = sectionMatId(elementSection(elementSets(i3)))
							objTgtExpanded(i5) = materialMaxStress(objComponent(i1),i4)
						enddo
					else
						elSetNum = objElSetId(i1)
						i5 = 0
						do i3 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
							i5 = i5 + 1
							i4 = sectionMatId(elementSection(elementSets(i3)))
							objTgtExpanded(i5) = materialMaxStrain(objComponent(i1),i4)
						enddo
					endif
				elseif(objTgtTag(i1) .eq. 'maxCompressive') then
					if(objCategory(i1) .eq. 'stress') then
						elSetNum = objElSetId(i1)
						i5 = 0
						do i3 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
							i5 = i5 + 1
							i4 = sectionMatId(elementSection(elementSets(i3)))
							objTgtExpanded(i5) = -abs(materialMaxStress(objComponent(i1)+3,i4))
						enddo
					else
						elSetNum = objElSetId(i1)
						i5 = 0
						do i3 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
							i5 = i5 + 1
							i4 = sectionMatId(elementSection(elementSets(i3)))
							objTgtExpanded(i5) = -abs(materialMaxStrain(objComponent(i1)+3,i4))
						enddo
					endif
				elseif(objTgtTag(i1) .eq. 'maxShear') then
					if(objCategory(i1) .eq. 'stress') then
						elSetNum = objElSetId(i1)
						i5 = 0
						do i3 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
							i5 = i5 + 1
							i4 = sectionMatId(elementSection(elementSets(i3)))
							objTgtExpanded(i5) = materialMaxStress(objComponent(i1)+3,i4)
						enddo
					else
						elSetNum = objElSetId(i1)
						i5 = 0
						do i3 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
							i5 = i5 + 1
							i4 = sectionMatId(elementSection(elementSets(i3)))
							objTgtExpanded(i5) = materialMaxStrain(objComponent(i1)+3,i4)
						enddo
					endif
				elseif(objTgtTag(i1) .eq. 'maxStrainEnergy') then
					elSetNum = objElSetId(i1)
					i5 = 0
					do i3 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						i5 = i5 + 1
						i4 = sectionMatId(elementSection(elementSets(i3)))
						objTgtExpanded(i5) = materialMaxStrnEngy(i4)
					enddo
				endif
			
			    sumTerm = r_0
			    if(objCategory(i1) .eq. 'temperature') then
					ndSetNum = objNdSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
					i2 = 1
					do i3 = ndSetRange(ndSetNum-1) + 1, ndSetRange(ndSetNum)
						i4 = currentRank(nodeSets(i3))
						sumTerm = sumTerm + (nodeTemp(i4) - objTgtExpanded(i2))**c2
						i2 = i2 + 1
					enddo
					objVal(i1) = objVal(i1) + c1*sumTerm
				elseif(objCategory(i1) .eq. 'tDot') then
					ndSetNum = objNdSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
					i2 = 1
					do i3 = ndSetRange(ndSetNum-1) + 1, ndSetRange(ndSetNum)
						i4 = currentRank(nodeSets(i3))
						sumTerm = sumTerm + (nodeTdot(i4) - objTgtExpanded(i2))**c2
						i2 = i2 + 1
					enddo
					objVal(i1) = objVal(i1) + c1*sumTerm
				elseif(objCategory(i1) .eq. 'displacement') then
					ndSetNum = objNdSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
					i2 = 1
					do i3 = ndSetRange(ndSetNum-1) + 1, ndSetRange(ndSetNum)
						i4 = nodeSets(i3)
						i5 = nDofIndex(objComponent(i1),i4)
						sumTerm = sumTerm + (nodeDisp(i5) - objTgtExpanded(i2))**c2
						i2 = i2 + 1
					enddo
					objVal(i1) = objVal(i1) + c1*sumTerm
				elseif(objCategory(i1) .eq. 'velocity') then
					ndSetNum = objNdSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
					i2 = 1
					do i3 = ndSetRange(ndSetNum-1) + 1, ndSetRange(ndSetNum)
						i4 = nodeSets(i3)
						i5 = nDofIndex(objComponent(i1),i4)
						sumTerm = sumTerm + (nodeVel(i5) - objTgtExpanded(i2))**c2
						i2 = i2 + 1
					enddo
					objVal(i1) = objVal(i1) + c1*sumTerm
				elseif(objCategory(i1) .eq. 'acceleration') then
					ndSetNum = objNdSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
					i2 = 1
					do i3 = ndSetRange(ndSetNum-1) + 1, ndSetRange(ndSetNum)
						i4 = nodeSets(i3)
						i5 = nDofIndex(objComponent(i1),i4)
						sumTerm = sumTerm + (nodeAcc(i5) - objTgtExpanded(i2))**c2
						i2 = i2 + 1
					enddo
					objVal(i1) = objVal(i1) + c1*sumTerm
				elseif(objCategory(i1) .eq. 'stress') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
					layer = objLayer(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElStress(stress,strain,i3,layer,sPt)
							sumTerm = sumTerm + (stress(comp) - objTgtExpanded(i4))**c2
							i4 = i4 + 1
						enddo
						objVal(i1) = objVal(i1) + c1*sumTerm
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElStress(stress,strain,i3,layer,sPt)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + stress(comp)*vol
							totVol = totVol + vol
						enddo
						objVal(i1) = objVal(i1) + c1*((sumTerm/totVol) - objTgtExpanded(1))**c2
					endif
				elseif(objCategory(i1) .eq. 'strain') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
					layer = objLayer(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElStress(stress,strain,i3,layer,sPt)
							sumTerm = sumTerm + (strain(comp) - objTgtExpanded(i4))**c2
							i4 = i4 + 1
						enddo
						objVal(i1) = objVal(i1) + c1*sumTerm
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElStress(stress,strain,i3,layer,sPt)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + strain(comp)*vol
							totVol = totVol + vol
						enddo
						objVal(i1) = objVal(i1) + c1*((sumTerm/totVol) - objTgtExpanded(1))**c2
					endif
				elseif(objCategory(i1) .eq. 'strainEnergy') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElStrainEnergy(stEn,dSEdU,dSEdT,i3)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + ((stEn/vol) - objTgtExpanded(i4))**c2
							i4 = i4 + 1
						enddo
						objVal(i1) = objVal(i1) + c1*sumTerm
					elseif(objOperator(i1) .eq. 'volumeIntegral') then
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElStrainEnergy(stEn,dSEdU,dSEdT,i3)
							sumTerm = sumTerm + stEn
						enddo
						objVal(i1) = objVal(i1) + c1*(sumTerm - objTgtExpanded(1))**c2
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElStrainEnergy(stEn,dSEdU,dSEdT,i3)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + stEn
							totVol = totVol + vol
						enddo
						objVal(i1) = objVal(i1) + c1*((sumTerm/totVol) - objTgtExpanded(1))**c2
					endif
				elseif(objCategory(i1) .eq. 'shellFrcMom') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getShellFrcMom(frcMom,shDef,ptT,i3,sPt)
							sumTerm = sumTerm + (frcMom(comp) - objTgtExpanded(i4))**c2
							i4 = i4 + 1
						enddo
						objVal(i1) = objVal(i1) + c1*sumTerm
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getShellFrcMom(frcMom,shDef,ptT,i3,sPt)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + frcMom(comp)*vol
							totVol = totVol + vol
						enddo
						objVal(i1) = objVal(i1) + c1*((sumTerm/totVol) - objTgtExpanded(1))**c2
					endif
				elseif(objCategory(i1) .eq. 'shellDef') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getShellFrcMom(frcMom,shDef,ptT,i3,sPt)
							sumTerm = sumTerm + (shDef(comp) - objTgtExpanded(i4))**c2
							i4 = i4 + 1
						enddo
						objVal(i1) = objVal(i1) + c1*sumTerm
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getShellFrcMom(frcMom,shDef,ptT,i3,sPt)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + shDef(comp)*vol
							totVol = totVol + vol
						enddo
						objVal(i1) = objVal(i1) + c1*((sumTerm/totVol) - objTgtExpanded(1))**c2
					endif
				elseif(objCategory(i1) .eq. 'beamFrcMom') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getBeamFrcMom(bfrcMom,bDef,i3,sPt)
							sumTerm = sumTerm + (bfrcMom(comp) - objTgtExpanded(i4))**c2
							i4 = i4 + 1
						enddo
						objVal(i1) = objVal(i1) + c1*sumTerm
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getBeamFrcMom(bfrcMom,bDef,i3,sPt)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + bfrcMom(comp)*vol
							totVol = totVol + vol
						enddo
						objVal(i1) = objVal(i1) + c1*((sumTerm/totVol) - objTgtExpanded(1))**c2
					endif
				elseif(objCategory(i1) .eq. 'beamDef') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getBeamFrcMom(bfrcMom,bDef,i3,sPt)
							sumTerm = sumTerm + (bDef(comp) - objTgtExpanded(i4))**c2
							i4 = i4 + 1
						enddo
						objVal(i1) = objVal(i1) + c1*sumTerm
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getBeamFrcMom(bfrcMom,bDef,i3,sPt)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + bDef(comp)*vol
							totVol = totVol + vol
						enddo
						objVal(i1) = objVal(i1) + c1*((sumTerm/totVol) - objTgtExpanded(1))**c2
					endif
				elseif(objCategory(i1) .eq. 'flux') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElFlux(flux,dfluxdT,Tx,dTxdT,i3,sPt)
							sumTerm = sumTerm + (flux(comp) - objTgtExpanded(i4))**c2
							i4 = i4 + 1
						enddo
						objVal(i1) = objVal(i1) + c1*sumTerm
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElFlux(flux,dfluxdT,Tx,dTxdT,i3,sPt)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + flux(comp)*vol
							totVol = totVol + vol
						enddo
						objVal(i1) = objVal(i1) + c1*((sumTerm/totVol) - objTgtExpanded(1))**c2
					endif
				elseif(objCategory(i1) .eq. 'tempGradient') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElFlux(flux,dfluxdT,Tx,dTxdT,i3,sPt)
							sumTerm = sumTerm + (Tx(comp) - objTgtExpanded(i4))**c2
							i4 = i4 + 1
						enddo
						objVal(i1) = objVal(i1) + c1*sumTerm
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElFlux(flux,dfluxdT,Tx,dTxdT,i3,sPt)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + Tx(comp)*vol
							totVol = totVol + vol
						enddo
						objVal(i1) = objVal(i1) + c1*((sumTerm/totVol) - objTgtExpanded(1))**c2
					endif
				elseif(objCategory(i1) .eq. 'mass') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
					if(objOperator(i1) .eq. 'volumeIntegral') then
					    do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + mass
						enddo
						objVal(i1) = objVal(i1) + c1*(sumTerm - objTgtExpanded(1))**c2
					endif
				elseif(objCategory(i1) .eq. 'volume') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
					if(objOperator(i1) .eq. 'volumeIntegral') then
					    do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + vol
						enddo
						objVal(i1) = objVal(i1) + c1*(sumTerm - objTgtExpanded(1))**c2
					endif
				elseif(objCategory(i1) .eq. 'massDisp') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElMassVol(mass,vol,i3)
							call r_getElementDisp(disp,ddispdU,i3,sPt)
							sumTerm = sumTerm + (mass*disp(comp) - objTgtExpanded(i4))**c2
							i4 = i4 + 1
						enddo
						objVal(i1) = objVal(i1) + c1*sumTerm
					endif
				endif
			endif
		enddo
		
	end subroutine getObjective
	
	subroutine getdObjdU(time)
	    implicit none
		
		real*8, intent(in) :: time
		
		real*8 :: c1, c2, constTgt, sumTerm, tVal, sPt(3)
		real*8 :: stress(6), strain(6), mass, vol, totVol
		real*8 :: dsdU(6,33), dsdT(6,10), dedU(6,33)
		real*8 :: frcMom(9), shDef(9), ptT, bfrcMom(6), bDef(6)
		real*8 :: dfMdU(9,33), dfMdT(9,33), dsDdU(9,33), dPtTdT(10)
		real*8 :: dbfrcMomdU(6,33), dbfrcMomdT(6,10), dbDefdU(6,33)
		real*8 :: stEn, dSEdU(33), dSEdT(10), flux(3), dfluxdT(3,10), Tx(3), dTxdT(3,10)
		real*8 :: disp(6), ddispdU(6,33)
		real*8 :: intPts(3,8), ipWt(8)
		integer :: numNds, dofPerNd, numIntDof, numIntPts, dofTable(2,33)
		integer :: i1, i2, i3, i4, i5, i6, i7, i8, i9, ndSetNum, elSetNum, tgtListLen, layer, comp
		
		if(objElSetId(1) .eq. 0 .and. objNdSetId(1) .eq. 0) then		
			do i1 = 1, numObjTerms
				do i2 = 1, numNdSets
					if(ndSetName(i2) .eq. objNdSet(i1)) then
						objNdSetId(i1) = i2
					endif
				enddo
				do i2 = 1, numElSets
					if(elSetName(i2) .eq. objElSet(i1)) then
						objElSetId(i1) = i2
					endif
				enddo
			enddo
		endif
		
		if(.not. allocated(objTgtExpanded)) then
		    i1 = numNodes + numEls
		    allocate(objTgtExpanded(i1))
		endif
		
		do i1 = 1, numObjTerms
		    if(time .ge. objActTime(1,i1) .and. time .le. objActTime(2,i1)) then
			    if(objNdSetId(i1) .ne. 0) then
				    i3 = objNdSetId(i1)
				    i2 = ndSetRange(i3) - ndSetRange(i3-1)
				elseif(objElSetId(i1) .ne. 0) then
				    i3 = objElSetId(i1)
					i2 = elSetRange(i3) - elSetRange(i3-1)
				endif
				
				if(objTgtTag(i1) .eq. 'none') then
					objTgtExpanded(1:i2) = r_0
				elseif(objTgtTag(i1) .eq. 'single') then
					i3 = objTgtRange(i1)
					objTgtExpanded(1:i2) = objTgtVal(i3)
				elseif(objTgtTag(i1) .eq. 'vector') then
					i3 = objTgtRange(i1-1) + 1
					i4 = objTgtRange(i1)
					objTgtExpanded(1:i2) = objTgtVal(i3:i4)
				elseif(objTgtTag(i1) .eq. 'maxTensile') then
					if(objCategory(i1) .eq. 'stress') then
						elSetNum = objElSetId(i1)
						i5 = 0
						do i3 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
							i5 = i5 + 1
							i4 = sectionMatId(elementSection(elementSets(i3)))
							objTgtExpanded(i5) = materialMaxStress(objComponent(i1),i4)
						enddo
					else
						elSetNum = objElSetId(i1)
						i5 = 0
						do i3 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
							i5 = i5 + 1
							i4 = sectionMatId(elementSection(elementSets(i3)))
							objTgtExpanded(i5) = materialMaxStrain(objComponent(i1),i4)
						enddo
					endif
				elseif(objTgtTag(i1) .eq. 'maxCompressive') then
					if(objCategory(i1) .eq. 'stress') then
						elSetNum = objElSetId(i1)
						i5 = 0
						do i3 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
							i5 = i5 + 1
							i4 = sectionMatId(elementSection(elementSets(i3)))
							objTgtExpanded(i5) = -abs(materialMaxStress(objComponent(i1)+3,i4))
						enddo
					else
						elSetNum = objElSetId(i1)
						i5 = 0
						do i3 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
							i5 = i5 + 1
							i4 = sectionMatId(elementSection(elementSets(i3)))
							objTgtExpanded(i5) = -abs(materialMaxStrain(objComponent(i1)+3,i4))
						enddo
					endif
				elseif(objTgtTag(i1) .eq. 'maxShear') then
					if(objCategory(i1) .eq. 'stress') then
						elSetNum = objElSetId(i1)
						i5 = 0
						do i3 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
							i5 = i5 + 1
							i4 = sectionMatId(elementSection(elementSets(i3)))
							objTgtExpanded(i5) = materialMaxStress(objComponent(i1)+3,i4)
						enddo
					else
						elSetNum = objElSetId(i1)
						i5 = 0
						do i3 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
							i5 = i5 + 1
							i4 = sectionMatId(elementSection(elementSets(i3)))
							objTgtExpanded(i5) = materialMaxStrain(objComponent(i1)+3,i4)
						enddo
					endif
				elseif(objTgtTag(i1) .eq. 'maxStrainEnergy') then
					elSetNum = objElSetId(i1)
					i5 = 0
					do i3 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						i5 = i5 + 1
						i4 = sectionMatId(elementSection(elementSets(i3)))
						objTgtExpanded(i5) = materialMaxStrnEngy(i4)
					enddo
				endif
			
			    sumTerm = r_0
				dsumTermdT(:) = r_0
				dsumTermdU(:) = r_0
				if(intVecSize .gt. 0) then
				    intdsumTermdU(:) = r_0
				endif
			    if(objCategory(i1) .eq. 'temperature') then
					ndSetNum = objNdSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
					i2 = 1
					do i3 = ndSetRange(ndSetNum-1) + 1, ndSetRange(ndSetNum)
						i4 = currentRank(nodeSets(i3))
						sumTerm = sumTerm + (nodeTemp(i4) - objTgtExpanded(i2))**c2
						dsumTermdT(i4) = dsumTermdT(i4) + c2*(nodeTemp(i4) - objTgtExpanded(i2))**(c2-r_1)
						i2 = i2 + 1
					enddo
					dLdt(:) = dLdt(:) + c1*dsumTermdT(:)
				elseif(objCategory(i1) .eq. 'tDot') then
					ndSetNum = objNdSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
					i2 = 1
					do i3 = ndSetRange(ndSetNum-1) + 1, ndSetRange(ndSetNum)
						i4 = currentRank(nodeSets(i3))
						sumTerm = sumTerm + (nodeTdot(i4) - objTgtExpanded(i2))**c2
						dsumTermdT(i4) = dsumTermdT(i4) + c2*(nodeTdot(i4) - objTgtExpanded(i2))**(c2-r_1)
						i2 = i2 + 1
					enddo
					dLdtdot(:) = dLdtdot(:) + c1*dsumTermdT(:)
				elseif(objCategory(i1) .eq. 'displacement') then
					ndSetNum = objNdSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
					i2 = 1
					do i3 = ndSetRange(ndSetNum-1) + 1, ndSetRange(ndSetNum)
						i4 = nodeSets(i3)
						i5 = nDofIndex(objComponent(i1),i4)
						sumTerm = sumTerm + (nodeDisp(i5) - objTgtExpanded(i2))**c2
						dsumTermdU(i5) = dsumTermdU(i5) + c2*(nodeDisp(i5) - objTgtExpanded(i2))**(c2-r_1)
						i2 = i2 + 1
					enddo
					dLdu(:) = dLdu(:) + c1*dsumTermdU(:)
				elseif(objCategory(i1) .eq. 'velocity') then
					ndSetNum = objNdSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
					i2 = 1
					do i3 = ndSetRange(ndSetNum-1) + 1, ndSetRange(ndSetNum)
						i4 = nodeSets(i3)
						i5 = nDofIndex(objComponent(i1),i4)
						sumTerm = sumTerm + (nodeVel(i5) - objTgtExpanded(i2))**c2
						dsumTermdU(i5) = dsumTermdU(i5) + c2*(nodeVel(i5) - objTgtExpanded(i2))**(c2-r_1)
						i2 = i2 + 1
					enddo
					dLdv(:) = dLdv(:) + c1*dsumTermdU(:)
				elseif(objCategory(i1) .eq. 'acceleration') then
					ndSetNum = objNdSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
					i2 = 1
					do i3 = ndSetRange(ndSetNum-1) + 1, ndSetRange(ndSetNum)
						i4 = nodeSets(i3)
						i5 = nDofIndex(objComponent(i1),i4)
						sumTerm = sumTerm + (nodeAcc(i5) - objTgtExpanded(i2))**c2
						dsumTermdU(i5) = dsumTermdU(i5) + c2*(nodeAcc(i5) - objTgtExpanded(i2))**(c2-r_1)
						i2 = i2 + 1
					enddo
					dLda(:) = dLda(:) + c1*dsumTermdU(:)
				elseif(objCategory(i1) .eq. 'stress') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
					layer = objLayer(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElStress(stress,strain,i3,layer,sPt)
							sumTerm = sumTerm + (stress(comp) - objTgtExpanded(i4))**c2
							call r_getdElStressdU(dsdU,dsdT,dedU,i3,layer,sPt)
							call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
							i5 = numNds*dofPerNd + numIntDof
							tVal = c2*(stress(comp) - objTgtExpanded(i4))**(c2-r_1)
							do i6 = 1, i5
							    i7 = dofTable(1,i6)
								i8 = dofTable(2,i6)
								if(i8 .le. numNds) then
								    i9 = nDofIndex(i7,elementList(i8,i3))
									dsumTermdU(i9) = dsumTermdU(i9) + tVal*dsdU(comp,i6)
								else
								    i9 = intVecRange(i3-1) + i6 - numNds*dofPerNd
									intdsumTermdU(i9) = intdsumTermdU(i9) + tVal*dsdU(comp,i6)
								endif
							enddo
							do i6 = 1, numNds
							    i7 = currentRank(elementList(i6,i3))
								dsumTermdT(i7) = dsumTermdT(i7) + tVal*dsdT(comp,i6)
							enddo
							i4 = i4 + 1
						enddo
						dLdu(:) = dLdu(:) + c1*dsumTermdU(:)
						if(intVecSize .gt. 0) then
						    intdLdu(:) = intdLdu(:) + c1*intdsumTermdU(:)
						endif
						dLdt(:) = dLdt(:) + c1*dsumTermdT(:)
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElStress(stress,strain,i3,layer,sPt)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + stress(comp)*vol
							call r_getdElStressdU(dsdU,dsdT,dedU,i3,layer,sPt)
							call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
							i5 = numNds*dofPerNd + numIntDof
							do i6 = 1, i5
							    i7 = dofTable(1,i6)
								i8 = dofTable(2,i6)
								if(i8 .le. numNds) then
								    i9 = nDofIndex(i7,elementList(i8,i3))
									dsumTermdU(i9) = dsumTermdU(i9) + dsdU(comp,i6)*vol
								else
								    i9 = intVecRange(i3-1) + i6 - numNds*dofPerNd
									intdsumTermdU(i9) = intdsumTermdU(i9) + dsdU(comp,i6)*vol
								endif
							enddo
							do i6 = 1, numNds
							    i7 = currentRank(elementList(i6,i3))
								dsumTermdT(i7) = dsumTermdT(i7) + dsdT(comp,i6)*vol
							enddo
							totVol = totVol + vol
						enddo
						tVal = (c2*c1/totVol)*((sumTerm/totVol) - objTgtExpanded(1))**(c2-r_1)
						dLdu(:) = dLdu(:) + tVal*dsumTermdU(:)
						if(intVecSize .gt. 0) then
						    intdLdu(:) = intdLdu(:) + tVal*intdsumTermdU(:)
						endif
						dLdt(:) = dLdt(:) + tVal*dsumTermdT(:)
					endif
				elseif(objCategory(i1) .eq. 'strain') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
					layer = objLayer(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElStress(stress,strain,i3,layer,sPt)
							sumTerm = sumTerm + (strain(comp) - objTgtExpanded(i4))**c2
							call r_getdElStressdU(dsdU,dsdT,dedU,i3,layer,sPt)
							call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
							i5 = numNds*dofPerNd + numIntDof
							tVal = c2*(strain(comp) - objTgtExpanded(i4))**(c2-r_1)
							do i6 = 1, i5
							    i7 = dofTable(1,i6)
								i8 = dofTable(2,i6)
								if(i8 .le. numNds) then
								    i9 = nDofIndex(i7,elementList(i8,i3))
									dsumTermdU(i9) = dsumTermdU(i9) + tVal*dedU(comp,i6)
								else
								    i9 = intVecRange(i3-1) + i6 - numNds*dofPerNd
									intdsumTermdU(i9) = intdsumTermdU(i9) + tVal*dedU(comp,i6)
								endif
							enddo
							i4 = i4 + 1
						enddo
						dLdu(:) = dLdu(:) + c1*dsumTermdU(:)
						if(intVecSize .gt. 0) then
						    intdLdu(:) = intdLdu(:) + c1*intdsumTermdU(:)
						endif
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElStress(stress,strain,i3,layer,sPt)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + strain(comp)*vol
							call r_getdElStressdU(dsdU,dsdT,dedU,i3,layer,sPt)
							call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
							i5 = numNds*dofPerNd + numIntDof
							do i6 = 1, i5
							    i7 = dofTable(1,i6)
								i8 = dofTable(2,i6)
								if(i8 .le. numNds) then
								    i9 = nDofIndex(i7,elementList(i8,i3))
									dsumTermdU(i9) = dsumTermdU(i9) + dedU(comp,i6)*vol
								else
								    i9 = intVecRange(i3-1) + i6 - numNds*dofPerNd
									intdsumTermdU(i9) = intdsumTermdU(i9) + dedU(comp,i6)*vol
								endif
							enddo
							totVol = totVol + vol
						enddo
						tVal = (c2*c1/totVol)*((sumTerm/totVol) - objTgtExpanded(1))**(c2-r_1)
						dLdu(:) = dLdu(:) + tVal*dsumTermdU(:)
						if(intVecSize .gt. 0) then
						    intdLdu(:) = intdLdu(:) + tVal*intdsumTermdU(:)
						endif
					endif
				elseif(objCategory(i1) .eq. 'strainEnergy') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElStrainEnergy(stEn,dSEdU,dSEdT,i3)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + ((stEn/vol) - objTgtExpanded(i4))**c2
							call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
							i5 = numNds*dofPerNd + numIntDof
							tVal = (c2/vol)*((stEn/vol) - objTgtExpanded(i4))**(c2-r_1)
							do i6 = 1, i5
							    i7 = dofTable(1,i6)
								i8 = dofTable(2,i6)
								if(i8 .le. numNds) then
								    i9 = nDofIndex(i7,elementList(i8,i3))
									dsumTermdU(i9) = dsumTermdU(i9) + tVal*dSEdU(i6)
								else
								    i9 = intVecRange(i3-1) + i6 - numNds*dofPerNd
									intdsumTermdU(i9) = intdsumTermdU(i9) + tVal*dSEdU(i6)
								endif
							enddo
							do i6 = 1, numNds
							    i7 = currentRank(elementList(i6,i3))
								dsumTermdT(i7) = dsumTermdT(i7) + tVal*dSEdT(i6)
							enddo
							i4 = i4 + 1
						enddo
						dLdu(:) = dLdu(:) + c1*dsumTermdU(:)
						if(intVecSize .gt. 0) then
						    intdLdu(:) = intdLdu(:) + c1*intdsumTermdU(:)
						endif
						dLdt(:) = dLdt(:) + c1*dsumTermdT(:)
					elseif(objOperator(i1) .eq. 'volumeIntegral') then
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElStrainEnergy(stEn,dSEdU,dSEdT,i3)
							sumTerm = sumTerm + stEn
							call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
							i5 = numNds*dofPerNd + numIntDof
							do i6 = 1, i5
							    i7 = dofTable(1,i6)
								i8 = dofTable(2,i6)
								if(i8 .le. numNds) then
								    i9 = nDofIndex(i7,elementList(i8,i3))
									dsumTermdU(i9) = dsumTermdU(i9) + dSEdU(i6)
								else
								    i9 = intVecRange(i3-1) + i6 - numNds*dofPerNd
									intdsumTermdU(i9) = intdsumTermdU(i9) + dSEdU(i6)
								endif
							enddo
							do i6 = 1, numNds
							    i7 = currentRank(elementList(i6,i3))
								dsumTermdT(i7) = dsumTermdT(i7) + dSEdT(i6)
							enddo
						enddo
						tVal = c2*c1*(sumTerm - objTgtExpanded(1))**(c2-r_1)
						dLdu(:) = dLdu(:) + tVal*dsumTermdU(:)
						if(intVecSize .gt. 0) then
						    intdLdu(:) = intdLdu(:) + tVal*intdsumTermdU(:)
						endif
						dLdt(:) = dLdt(:) + tVal*dsumTermdT(:)
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElStrainEnergy(stEn,dSEdU,dSEdT,i3)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + stEn
							call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
							i5 = numNds*dofPerNd + numIntDof
							do i6 = 1, i5
							    i7 = dofTable(1,i6)
								i8 = dofTable(2,i6)
								if(i8 .le. numNds) then
								    i9 = nDofIndex(i7,elementList(i8,i3))
									dsumTermdU(i9) = dsumTermdU(i9) + dSEdU(i6)
								else
								    i9 = intVecRange(i3-1) + i6 - numNds*dofPerNd
									intdsumTermdU(i9) = intdsumTermdU(i9) + dSEdU(i6)
								endif
							enddo
							do i6 = 1, numNds
							    i7 = currentRank(elementList(i6,i3))
								dsumTermdT(i7) = dsumTermdT(i7) + dSEdT(i6)
							enddo
							totVol = totVol + vol
						enddo
						tVal = (c2*c1/totVol)*((sumTerm/totVol) - objTgtExpanded(1))**(c2-r_1)
						dLdu(:) = dLdu(:) + tVal*dsumTermdU(:)
						if(intVecSize .gt. 0) then
						    intdLdu(:) = intdLdu(:) + tVal*intdsumTermdU(:)
						endif
						dLdt(:) = dLdt(:) + tVal*dsumTermdT(:)
					endif
				elseif(objCategory(i1) .eq. 'shellFrcMom') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
					layer = objLayer(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getShellFrcMom(frcMom,shDef,ptT,i3,sPt)
							sumTerm = sumTerm + (frcMom(comp) - objTgtExpanded(i4))**c2
							call r_getdShellFrcMomdU(dfMdU,dfMdT,dsDdU,dPtTdT,i3,sPt)
							call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
							i5 = numNds*dofPerNd + numIntDof
							tVal = c2*(frcMom(comp) - objTgtExpanded(i4))**(c2-r_1)
							do i6 = 1, i5
							    i7 = dofTable(1,i6)
								i8 = dofTable(2,i6)
								if(i8 .le. numNds) then
								    i9 = nDofIndex(i7,elementList(i8,i3))
									dsumTermdU(i9) = dsumTermdU(i9) + tVal*dfMdU(comp,i6)
								else
								    i9 = intVecRange(i3-1) + i6 - numNds*dofPerNd
									intdsumTermdU(i9) = intdsumTermdU(i9) + tVal*dfMdU(comp,i6)
								endif
							enddo
							do i6 = 1, numNds
							    i7 = currentRank(elementList(i6,i3))
								dsumTermdT(i7) = dsumTermdT(i7) + tVal*dfMdT(comp,i6)
							enddo
							i4 = i4 + 1
						enddo
						dLdu(:) = dLdu(:) + c1*dsumTermdU(:)
						if(intVecSize .gt. 0) then
						    intdLdu(:) = intdLdu(:) + c1*intdsumTermdU(:)
						endif
						dLdt(:) = dLdt(:) + c1*dsumTermdT(:)
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getShellFrcMom(frcMom,shDef,ptT,i3,sPt)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + frcMom(comp)*vol
							call r_getdShellFrcMomdU(dfMdU,dfMdT,dsDdU,dPtTdT,i3,sPt)
							call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
							i5 = numNds*dofPerNd + numIntDof
							do i6 = 1, i5
							    i7 = dofTable(1,i6)
								i8 = dofTable(2,i6)
								if(i8 .le. numNds) then
								    i9 = nDofIndex(i7,elementList(i8,i3))
									dsumTermdU(i9) = dsumTermdU(i9) + dfMdU(comp,i6)*vol
								else
								    i9 = intVecRange(i3-1) + i6 - numNds*dofPerNd
									intdsumTermdU(i9) = intdsumTermdU(i9) + dfMdU(comp,i6)*vol
								endif
							enddo
							do i6 = 1, numNds
							    i7 = currentRank(elementList(i6,i3))
								dsumTermdT(i7) = dsumTermdT(i7) + dfMdT(comp,i6)*vol
							enddo
							totVol = totVol + vol
						enddo
						tVal = (c2*c1/totVol)*((sumTerm/totVol) - objTgtExpanded(1))**(c2-r_1)
						dLdu(:) = dLdu(:) + tVal*dsumTermdU(:)
						if(intVecSize .gt. 0) then
						    intdLdu(:) = intdLdu(:) + tVal*intdsumTermdU(:)
						endif
						dLdt(:) = dLdt(:) + tVal*dsumTermdT(:)
					endif
				elseif(objCategory(i1) .eq. 'shellDef') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getShellFrcMom(frcMom,shDef,ptT,i3,sPt)
							sumTerm = sumTerm + (shDef(comp) - objTgtExpanded(i4))**c2
							call r_getdShellFrcMomdU(dfMdU,dfMdT,dsDdU,dPtTdT,i3,sPt)
							call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
							i5 = numNds*dofPerNd + numIntDof
							tVal = c2*(shDef(comp) - objTgtExpanded(i4))**(c2-r_1)
							do i6 = 1, i5
							    i7 = dofTable(1,i6)
								i8 = dofTable(2,i6)
								if(i8 .le. numNds) then
								    i9 = nDofIndex(i7,elementList(i8,i3))
									dsumTermdU(i9) = dsumTermdU(i9) + tVal*dsDdU(comp,i6)
								else
								    i9 = intVecRange(i3-1) + i6 - numNds*dofPerNd
									intdsumTermdU(i9) = intdsumTermdU(i9) + tVal*dsDdU(comp,i6)
								endif
							enddo
							i4 = i4 + 1
						enddo
						dLdu(:) = dLdu(:) + c1*dsumTermdU(:)
						if(intVecSize .gt. 0) then
						    intdLdu(:) = intdLdu(:) + c1*intdsumTermdU(:)
						endif
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getShellFrcMom(frcMom,shDef,ptT,i3,sPt)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + shDef(comp)*vol
							call r_getdShellFrcMomdU(dfMdU,dfMdT,dsDdU,dPtTdT,i3,sPt)
							call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
							i5 = numNds*dofPerNd + numIntDof
							do i6 = 1, i5
							    i7 = dofTable(1,i6)
								i8 = dofTable(2,i6)
								if(i8 .le. numNds) then
								    i9 = nDofIndex(i7,elementList(i8,i3))
									dsumTermdU(i9) = dsumTermdU(i9) + dsDdU(comp,i6)*vol
								else
								    i9 = intVecRange(i3-1) + i6 - numNds*dofPerNd
									intdsumTermdU(i9) = intdsumTermdU(i9) + dsDdU(comp,i6)*vol
								endif
							enddo
							totVol = totVol + vol
						enddo
						tVal = (c2*c1/totVol)*((sumTerm/totVol) - objTgtExpanded(1))**(c2-r_1)
						dLdu(:) = dLdu(:) + tVal*dsumTermdU(:)
						if(intVecSize .gt. 0) then
						    intdLdu(:) = intdLdu(:) + tVal*intdsumTermdU(:)
						endif
					endif
				elseif(objCategory(i1) .eq. 'beamFrcMom') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getBeamFrcMom(bfrcMom,bDef,i3,sPt)
							sumTerm = sumTerm + (bfrcMom(comp) - objTgtExpanded(i4))**c2
							call r_getdBeamFrcMomdU(dbfrcMomdU,dbfrcMomdT,dbDefdU,i3,sPt)
							call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
							i5 = numNds*dofPerNd + numIntDof
							tVal = c2*(bfrcMom(comp) - objTgtExpanded(i4))**(c2-r_1)
							do i6 = 1, i5
							    i7 = dofTable(1,i6)
								i8 = dofTable(2,i6)
								if(i8 .le. numNds) then
								    i9 = nDofIndex(i7,elementList(i8,i3))
									dsumTermdU(i9) = dsumTermdU(i9) + tVal*dbfrcMomdU(comp,i6)
								else
								    i9 = intVecRange(i3-1) + i6 - numNds*dofPerNd
									intdsumTermdU(i9) = intdsumTermdU(i9) + tVal*dbfrcMomdU(comp,i6)
								endif
							enddo
							do i6 = 1, numNds
							    i7 = currentRank(elementList(i6,i3))
								dsumTermdT(i7) = dsumTermdT(i7) + tVal*dbfrcMomdT(comp,i6)
							enddo
							i4 = i4 + 1
						enddo
						dLdu(:) = dLdu(:) + c1*dsumTermdU(:)
						if(intVecSize .gt. 0) then
						    intdLdu(:) = intdLdu(:) + c1*intdsumTermdU(:)
						endif
						dLdt(:) = dLdt(:) + c1*dsumTermdT(:)
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getBeamFrcMom(bfrcMom,bDef,i3,sPt)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + bfrcMom(comp)*vol
							call r_getdBeamFrcMomdU(dbfrcMomdU,dbfrcMomdT,dbDefdU,i3,sPt)
							call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
							i5 = numNds*dofPerNd + numIntDof
							do i6 = 1, i5
							    i7 = dofTable(1,i6)
								i8 = dofTable(2,i6)
								if(i8 .le. numNds) then
								    i9 = nDofIndex(i7,elementList(i8,i3))
									dsumTermdU(i9) = dsumTermdU(i9) + dbfrcMomdU(comp,i6)*vol
								else
								    i9 = intVecRange(i3-1) + i6 - numNds*dofPerNd
									intdsumTermdU(i9) = intdsumTermdU(i9) + dbfrcMomdU(comp,i6)*vol
								endif
							enddo
							do i6 = 1, numNds
							    i7 = currentRank(elementList(i6,i3))
								dsumTermdT(i7) = dsumTermdT(i7) + dbfrcMomdT(comp,i6)*vol
							enddo
							totVol = totVol + vol
						enddo
						tVal = (c2*c1/totVol)*((sumTerm/totVol) - objTgtExpanded(1))**(c2-r_1)
						dLdu(:) = dLdu(:) + tVal*dsumTermdU(:)
						if(intVecSize .gt. 0) then
						    intdLdu(:) = intdLdu(:) + tVal*intdsumTermdU(:)
						endif
						dLdt(:) = dLdt(:) + tVal*dsumTermdT(:)
					endif
				elseif(objCategory(i1) .eq. 'beamDef') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getBeamFrcMom(bfrcMom,bDef,i3,sPt)
							sumTerm = sumTerm + (bDef(comp) - objTgtExpanded(i4))**c2
							call r_getdBeamFrcMomdU(dbfrcMomdU,dbfrcMomdT,dbDefdU,i3,sPt)
							call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
							i5 = numNds*dofPerNd + numIntDof
							tVal = c2*(bDef(comp) - objTgtExpanded(i4))**(c2-r_1)
							do i6 = 1, i5
							    i7 = dofTable(1,i6)
								i8 = dofTable(2,i6)
								if(i8 .le. numNds) then
								    i9 = nDofIndex(i7,elementList(i8,i3))
									dsumTermdU(i9) = dsumTermdU(i9) + tVal*dbDefdU(comp,i6)
								else
								    i9 = intVecRange(i3-1) + i6 - numNds*dofPerNd
									intdsumTermdU(i9) = intdsumTermdU(i9) + tVal*dbDefdU(comp,i6)
								endif
							enddo
							i4 = i4 + 1
						enddo
						dLdu(:) = dLdu(:) + c1*dsumTermdU(:)
						if(intVecSize .gt. 0) then
						    intdLdu(:) = intdLdu(:) + c1*intdsumTermdU(:)
						endif
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getBeamFrcMom(bfrcMom,bDef,i3,sPt)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + bDef(comp)*vol
							call r_getdBeamFrcMomdU(dbfrcMomdU,dbfrcMomdT,dbDefdU,i3,sPt)
							call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
							i5 = numNds*dofPerNd + numIntDof
							do i6 = 1, i5
							    i7 = dofTable(1,i6)
								i8 = dofTable(2,i6)
								if(i8 .le. numNds) then
								    i9 = nDofIndex(i7,elementList(i8,i3))
									dsumTermdU(i9) = dsumTermdU(i9) + dbDefdU(comp,i6)*vol
								else
								    i9 = intVecRange(i3-1) + i6 - numNds*dofPerNd
									intdsumTermdU(i9) = intdsumTermdU(i9) + dbDefdU(comp,i6)*vol
								endif
							enddo
							totVol = totVol + vol
						enddo
						tVal = (c2*c1/totVol)*((sumTerm/totVol) - objTgtExpanded(1))**(c2-r_1)
						dLdu(:) = dLdu(:) + tVal*dsumTermdU(:)
						if(intVecSize .gt. 0) then
						    intdLdu(:) = intdLdu(:) + tVal*intdsumTermdU(:)
						endif
					endif
				elseif(objCategory(i1) .eq. 'flux') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElFlux(flux,dfluxdT,Tx,dTxdT,i3,sPt)
							sumTerm = sumTerm + (flux(comp) - objTgtExpanded(i4))**c2
							call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
							tVal = c2*(flux(comp) - objTgtExpanded(i4))**(c2-r_1)
							do i6 = 1, numNds
							    i7 = currentRank(elementList(i6,i3))
								dsumTermdT(i7) = dsumTermdT(i7) + tVal*dfluxdT(comp,i6)
							enddo
							i4 = i4 + 1
						enddo
						dLdt(:) = dLdt(:) + c1*dsumTermdT(:)
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElFlux(flux,dfluxdT,Tx,dTxdT,i3,sPt)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + flux(comp)*vol
							call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
							do i6 = 1, numNds
							    i7 = currentRank(elementList(i6,i3))
								dsumTermdT(i7) = dsumTermdT(i7) + dfluxdT(comp,i6)*vol
							enddo
							totVol = totVol + vol
						enddo
						tVal = (c2*c1/totVol)*((sumTerm/totVol) - objTgtExpanded(1))**(c2-r_1)
						dLdt(:) = dLdt(:) + c1*tVal*dsumTermdT(:)
					endif
				elseif(objCategory(i1) .eq. 'tempGradient') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElFlux(flux,dfluxdT,Tx,dTxdT,i3,sPt)
							sumTerm = sumTerm + (Tx(comp) - objTgtExpanded(i4))**c2
							call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
							tVal = c2*(Tx(comp) - objTgtExpanded(i4))**(c2-r_1)
							do i6 = 1, numNds
							    i7 = currentRank(elementList(i6,i3))
								dsumTermdT(i7) = dsumTermdT(i7) + tVal*dTxdT(comp,i6)
							enddo
							i4 = i4 + 1
						enddo
						objVal(i1) = objVal(i1) + c1*sumTerm
						dLdt(:) = dLdt(:) + c1*dsumTermdT(:)
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElFlux(flux,dfluxdT,Tx,dTxdT,i3,sPt)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + Tx(comp)*vol
							call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
							do i6 = 1, numNds
							    i7 = currentRank(elementList(i6,i3))
								dsumTermdT(i7) = dsumTermdT(i7) + dTxdT(comp,i6)*vol
							enddo
							totVol = totVol + vol
						enddo
						tVal = (c2*c1/totVol)*((sumTerm/totVol) - objTgtExpanded(1))**(c2-r_1)
						dLdt(:) = dLdt(:) + tVal*dsumTermdT(:)
					endif
				elseif(objCategory(i1) .eq. 'massDisp') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElMassVol(mass,vol,i3)
							call r_getElementDisp(disp,ddispdU,i3,sPt)
							sumTerm = sumTerm + (mass*disp(comp) - objTgtExpanded(i4))**c2
							call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
							i5 = numNds*dofPerNd + numIntDof
							tVal = c2*mass*(mass*disp(comp) - objTgtExpanded(i4))**(c2-r_1)
							do i6 = 1, i5
							    i7 = dofTable(1,i6)
								i8 = dofTable(2,i6)
								if(i8 .le. numNds) then
								    i9 = nDofIndex(i7,elementList(i8,i3))
									dsumTermdU(i9) = dsumTermdU(i9) + tVal*ddispdU(comp,i6)
								else
								    i9 = intVecRange(i3-1) + i6 - numNds*dofPerNd
									intdsumTermdU(i9) = intdsumTermdU(i9) + tVal*ddispdU(comp,i6)
								endif
							enddo
							i4 = i4 + 1
						enddo
						objVal(i1) = objVal(i1) + c1*sumTerm
						dLdu(:) = dLdu(:) + c1*dsumTermdU(:)
						if(intVecSize .gt. 0) then
						    intdLdu(:) = intdLdu(:) + c1*intdsumTermdU(:)
						endif
					endif
				endif
			endif
		enddo
		
	end subroutine getdObjdU
	
	subroutine getdObjdD(time)
	    implicit none
		
		real*8, intent(in) :: time
		
		real*8 :: c1, c2, constTgt, sumTerm, tVal, sPt(3)
		real*8 :: stress(6), strain(6), mass, vol, totVol
		real*8 :: frcMom(9), shDef(9), ptT, bfrcMom(6), bDef(6)
		real*8 :: stEn, dSEdU(33), dSEdT(10), flux(3), dfluxdT(3,10), Tx(3), dTxdT(3,10)
		real*8 :: disp(6), ddispdU(6,33)
		
		complex*16 :: c_stress(6), c_strain(6), c_mass, c_vol, c_totVol
		complex*16 :: c_frcMom(9), c_shDef(9), c_ptT, c_bfrcMom(6), c_bDef(6)
		complex*16 :: c_stEn, c_dSEdU(33), c_dSEdT(10), c_flux(3), c_dfluxdT(3,10), c_Tx(3), c_dTxdT(3,10)
		complex*16 :: c_disp(6), c_ddispdU(6,33)
		
		integer :: i1, i2, i3, i4, i5, i6, ndSetNum, elSetNum, tgtListLen, layer, comp
		
		c_dVec(:) = c_1*r_dVec(:)
		
		if(objElSetId(1) .eq. 0 .and. objNdSetId(1) .eq. 0) then		
			do i1 = 1, numObjTerms
				do i2 = 1, numNdSets
					if(ndSetName(i2) .eq. objNdSet(i1)) then
						objNdSetId(i1) = i2
					endif
				enddo
				do i2 = 1, numElSets
					if(elSetName(i2) .eq. objElSet(i1)) then
						objElSetId(i1) = i2
					endif
				enddo
			enddo
		endif
		
		if(.not. allocated(objTgtExpanded)) then
		    i1 = numNodes + numEls
		    allocate(objTgtExpanded(i1))
		endif
		
		do i1 = 1, numObjTerms
		    if(time .ge. objActTime(1,i1) .and. time .le. objActTime(2,i1)) then
			    if(objNdSetId(i1) .ne. 0) then
				    i3 = objNdSetId(i1)
				    i2 = ndSetRange(i3) - ndSetRange(i3-1)
				elseif(objElSetId(i1) .ne. 0) then
				    i3 = objElSetId(i1)
					i2 = elSetRange(i3) - elSetRange(i3-1)
				endif
				
				if(objTgtTag(i1) .eq. 'none') then
					objTgtExpanded(1:i2) = r_0
				elseif(objTgtTag(i1) .eq. 'single') then
					i3 = objTgtRange(i1)
					objTgtExpanded(1:i2) = objTgtVal(i3)
				elseif(objTgtTag(i1) .eq. 'vector') then
					i3 = objTgtRange(i1-1) + 1
					i4 = objTgtRange(i1)
					objTgtExpanded(1:i2) = objTgtVal(i3:i4)
				elseif(objTgtTag(i1) .eq. 'maxTensile') then
					if(objCategory(i1) .eq. 'stress') then
						elSetNum = objElSetId(i1)
						i5 = 0
						do i3 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
							i5 = i5 + 1
							i4 = sectionMatId(elementSection(elementSets(i3)))
							objTgtExpanded(i5) = materialMaxStress(objComponent(i1),i4)
						enddo
					else
						elSetNum = objElSetId(i1)
						i5 = 0
						do i3 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
							i5 = i5 + 1
							i4 = sectionMatId(elementSection(elementSets(i3)))
							objTgtExpanded(i5) = materialMaxStrain(objComponent(i1),i4)
						enddo
					endif
				elseif(objTgtTag(i1) .eq. 'maxCompressive') then
					if(objCategory(i1) .eq. 'stress') then
						elSetNum = objElSetId(i1)
						i5 = 0
						do i3 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
							i5 = i5 + 1
							i4 = sectionMatId(elementSection(elementSets(i3)))
							objTgtExpanded(i5) = -abs(materialMaxStress(objComponent(i1)+3,i4))
						enddo
					else
						elSetNum = objElSetId(i1)
						i5 = 0
						do i3 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
							i5 = i5 + 1
							i4 = sectionMatId(elementSection(elementSets(i3)))
							objTgtExpanded(i5) = -abs(materialMaxStrain(objComponent(i1)+3,i4))
						enddo
					endif
				elseif(objTgtTag(i1) .eq. 'maxShear') then
					if(objCategory(i1) .eq. 'stress') then
						elSetNum = objElSetId(i1)
						i5 = 0
						do i3 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
							i5 = i5 + 1
							i4 = sectionMatId(elementSection(elementSets(i3)))
							objTgtExpanded(i5) = materialMaxStress(objComponent(i1)+3,i4)
						enddo
					else
						elSetNum = objElSetId(i1)
						i5 = 0
						do i3 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
							i5 = i5 + 1
							i4 = sectionMatId(elementSection(elementSets(i3)))
							objTgtExpanded(i5) = materialMaxStrain(objComponent(i1)+3,i4)
						enddo
					endif
				elseif(objTgtTag(i1) .eq. 'maxStrainEnergy') then
					elSetNum = objElSetId(i1)
					i5 = 0
					do i3 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						i5 = i5 + 1
						i4 = sectionMatId(elementSection(elementSets(i3)))
						objTgtExpanded(i5) = materialMaxStrnEngy(i4)
					enddo
				endif
			
			    sumTerm = r_0
				dsumTermdD(:) = r_0
				if(objCategory(i1) .eq. 'stress') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
					layer = objLayer(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElStress(stress,strain,i3,layer,sPt)
							sumTerm = sumTerm + (stress(comp) - objTgtExpanded(i4))**c2
							do i5 = elToDCRange(i3-1)+1, elToDCRange(i3)
							    i6 = elToDComp(i5)
								c_dVec(i6) = c_dVec(i6) + compStep
								call c_getElStress(c_stress,c_strain,i3,layer,c_1*sPt)
								tVal = c2*(stress(comp) - objTgtExpanded(i4))**(c2-r_1)
								dsumTermdD(i6) = dsumTermdD(i6) + tVal*imag(c_stress(comp))*compStepInv
								c_dVec(i6) = c_dVec(i6) - compStep
							enddo
							i4 = i4 + 1
						enddo
						dLdD(:) = dLdD(:) + c1*dsumTermdD(:)
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						dtotVoldD = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElStress(stress,strain,i3,layer,sPt)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + stress(comp)*vol
							totVol = totVol + vol
							do i5 = elToDCRange(i3-1)+1, elToDCRange(i3)
							    i6 = elToDComp(i5)
								c_dVec(i6) = c_dVec(i6) + compStep
								call c_getElStress(c_stress,c_strain,i3,layer,c_1*sPt)
								call c_getElMassVol(c_mass,c_vol,i3)
								tVal = compStepInv*(imag(c_stress(comp))*vol + stress(comp)*imag(c_vol))
								dsumTermdD(i6) = dsumTermdD(i6) + tVal
								dtotVoldD(i6) = dtotVoldD(i6) + imag(c_vol)*compStepInv
								c_dVec(i6) = c_dVec(i6) - compStep
							enddo
						enddo
						tVal = c2*c1*((sumTerm/totVol) - objTgtExpanded(1))**(c2-r_1)
						dLdD(:) = dLdD(:) + (tVal/totVol)*(dsumTermdD(:) - (sumTerm/totVol)*dtotVoldD(:))
					endif
				elseif(objCategory(i1) .eq. 'strain') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
					layer = objLayer(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElStress(stress,strain,i3,layer,sPt)
							sumTerm = sumTerm + (strain(comp) - objTgtExpanded(i4))**c2
							do i5 = elToDCRange(i3-1)+1, elToDCRange(i3)
							    i6 = elToDComp(i5)
								c_dVec(i6) = c_dVec(i6) + compStep
								call c_getElStress(c_stress,c_strain,i3,layer,c_1*sPt)
								tVal = c2*(strain(comp) - objTgtExpanded(i4))**(c2-r_1)
								dsumTermdD(i6) = dsumTermdD(i6) + tVal*imag(c_strain(comp))*compStepInv
								c_dVec(i6) = c_dVec(i6) - compStep
							enddo
							i4 = i4 + 1
						enddo
						dLdD(:) = dLdD(:) + c1*dsumTermdD(:)
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						dtotVoldD = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElStress(stress,strain,i3,layer,sPt)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + strain(comp)*vol
							totVol = totVol + vol
							do i5 = elToDCRange(i3-1)+1, elToDCRange(i3)
							    i6 = elToDComp(i5)
								c_dVec(i6) = c_dVec(i6) + compStep
								call c_getElStress(c_stress,c_strain,i3,layer,c_1*sPt)
								call c_getElMassVol(c_mass,c_vol,i3)
								tVal = compStepInv*(imag(c_strain(comp))*vol + strain(comp)*imag(c_vol))
								dsumTermdD(i6) = dsumTermdD(i6) + tVal
								dtotVoldD(i6) = dtotVoldD(i6) + imag(c_vol)*compStepInv
								c_dVec(i6) = c_dVec(i6) - compStep
							enddo
						enddo
						tVal = c2*c1*((sumTerm/totVol) - objTgtExpanded(1))**(c2-r_1)
						dLdD(:) = dLdD(:) + (tVal/totVol)*(dsumTermdD(:) - (sumTerm/totVol)*dtotVoldD(:))
					endif
				elseif(objCategory(i1) .eq. 'strainEnergy') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElStrainEnergy(stEn,dSEdU,dSEdT,i3)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + ((stEn/vol) - objTgtExpanded(i4))**c2
							do i5 = elToDCRange(i3-1)+1, elToDCRange(i3)
							    i6 = elToDComp(i5)
								c_dVec(i6) = c_dVec(i6) + compStep
								call c_getElStrainEnergy(c_stEn,c_dSEdU,c_dSEdT,i3)
								call c_getElMassVol(c_mass,c_vol,i3)
								tVal = c2*((stEn/vol) - objTgtExpanded(i4))**(c2-r_1)
								dsumTermdD(i6) = dsumTermdD(i6) + tVal*compStepInv*(imag(c_stEn) - stEn*imag(c_vol)/vol)/vol
								c_dVec(i6) = c_dVec(i6) - compStep
							enddo
							i4 = i4 + 1
						enddo
						dLdD(:) = dLdD(:) + c1*dsumTermdD(:)
					elseif(objOperator(i1) .eq. 'volumeIntegral') then
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElStrainEnergy(stEn,dSEdU,dSEdT,i3)
							sumTerm = sumTerm + stEn
							do i5 = elToDCRange(i3-1)+1, elToDCRange(i3)
							    i6 = elToDComp(i5)
								c_dVec(i6) = c_dVec(i6) + compStep
								call c_getElStrainEnergy(c_stEn,c_dSEdU,c_dSEdT,i3)
								dsumTermdD(i6) = dsumTermdD(i6) + compStepInv*imag(c_stEn)
								c_dVec(i6) = c_dVec(i6) - compStep
							enddo
						enddo
						tVal = c2*c1*(sumTerm - objTgtExpanded(1))**(c2-r_1)
						dLdD(:) = dLdD(:) + tVal*dsumTermdD(:)
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElStrainEnergy(stEn,dSEdU,dSEdT,i3)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + stEn
							totVol = totVol + vol
							do i5 = elToDCRange(i3-1)+1, elToDCRange(i3)
							    i6 = elToDComp(i5)
								c_dVec(i6) = c_dVec(i6) + compStep
								call c_getElStrainEnergy(c_stEn,c_dSEdU,c_dSEdT,i3)
								call c_getElMassVol(c_mass,c_vol,i3)
								dsumTermdD(i6) = dsumTermdD(i6) + compStepInv*imag(c_stEn)
								dtotVoldD(i6) = dtotVoldD(i6) + compStepInv*imag(c_vol)
								c_dVec(i6) = c_dVec(i6) - compStep
							enddo
						enddo
						tVal = c2*c1*((sumTerm/totVol) - objTgtExpanded(1))**(c2-r_1)
						dLdD(:) = dLdD(:) + (tVal/totVol)*(dsumTermdD(:) - (sumTerm/totVol)*dTotVoldD(:))
					endif
				elseif(objCategory(i1) .eq. 'shellFrcMom') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getShellFrcMom(frcMom,shDef,ptT,i3,sPt)
							sumTerm = sumTerm + (frcMom(comp) - objTgtExpanded(i4))**c2
							do i5 = elToDCRange(i3-1)+1, elToDCRange(i3)
							    i6 = elToDComp(i5)
								c_dVec(i6) = c_dVec(i6) + compStep
								call c_getShellFrcMom(c_frcMom,c_shDef,c_ptT,i3,c_1*sPt)
								tVal = c2*(frcMom(comp) - objTgtExpanded(i4))**(c2-r_1)
								dsumTermdD(i6) = dsumTermdD(i6) + tVal*imag(c_frcMom(comp))*compStepInv
								c_dVec(i6) = c_dVec(i6) - compStep
							enddo
							i4 = i4 + 1
						enddo
						dLdD(:) = dLdD(:) + c1*dsumTermdD(:)
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getShellFrcMom(frcMom,shDef,ptT,i3,sPt)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + frcMom(comp)*vol
							totVol = totVol + vol
							do i5 = elToDCRange(i3-1)+1, elToDCRange(i3)
							    i6 = elToDComp(i5)
								c_dVec(i6) = c_dVec(i6) + compStep
								call c_getShellFrcMom(c_frcMom,c_shDef,c_ptT,i3,c_1*sPt)
								call c_getElMassVol(c_mass,c_vol,i3)
								tVal = compStepInv*(imag(c_frcMom(comp))*vol + frcMom(comp)*imag(c_vol))
								dsumTermdD(i6) = dsumTermdD(i6) + tVal
								dtotVoldD(i6) = dtotVoldD(i6) + compStepInv*imag(c_vol)
								c_dVec(i6) = c_dVec(i6) - compStep
							enddo
						enddo
						tVal = c2*c1*((sumTerm/totVol) - objTgtExpanded(1))**(c2-r_1)
						dLdD(:) = dLdD(:) + (tVal/totVol)*(dsumTermdD(:) - (sumTerm/totVol)*dtotVoldD(:))
					endif
				elseif(objCategory(i1) .eq. 'shellDef') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getShellFrcMom(frcMom,shDef,ptT,i3,sPt)
							sumTerm = sumTerm + (shDef(comp) - objTgtExpanded(i4))**c2
							do i5 = elToDCRange(i3-1)+1, elToDCRange(i3)
							    i6 = elToDComp(i5)
								c_dVec(i6) = c_dVec(i6) + compStep
								call c_getShellFrcMom(c_frcMom,c_shDef,c_ptT,i3,c_1*sPt)
								tVal = c2*(shDef(comp) - objTgtExpanded(i4))**(c2-r_1)
								dsumTermdD(i6) = dsumTermdD(i6) + tVal*imag(c_shDef(comp))*compStepInv
								c_dVec(i6) = c_dVec(i6) - compStep
							enddo
							i4 = i4 + 1
						enddo
						dLdD(:) = dLdD(:) + c1*dsumTermdD(:)
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getShellFrcMom(frcMom,shDef,ptT,i3,sPt)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + shDef(comp)*vol
							totVol = totVol + vol
							do i5 = elToDCRange(i3-1)+1, elToDCRange(i3)
							    i6 = elToDComp(i5)
								c_dVec(i6) = c_dVec(i6) + compStep
								call c_getShellFrcMom(c_frcMom,c_shDef,c_ptT,i3,c_1*sPt)
								call c_getElMassVol(c_mass,c_vol,i3)
								tVal = compStepInv*(imag(c_shDef(comp))*vol + shDef(comp)*imag(c_vol))
								dsumTermdD(i6) = dsumTermdD(i6) + tVal
								dtotVoldD(i6) = dtotVoldD(i6) + compStepInv*imag(c_vol)
								c_dVec(i6) = c_dVec(i6) - compStep
							enddo
						enddo
						tVal = c2*c1*((sumTerm/totVol) - objTgtExpanded(1))**(c2-r_1)
						dLdD(:) = dLdD(:) + (tVal/totVol)*(dsumTermdD(:) - (sumTerm/totVol)*dtotVoldD(:))
					endif
				elseif(objCategory(i1) .eq. 'beamFrcMom') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getBeamFrcMom(bfrcMom,bDef,i3,sPt)
							sumTerm = sumTerm + (bfrcMom(comp) - objTgtExpanded(i4))**c2
							do i5 = elToDCRange(i3-1)+1, elToDCRange(i3)
							    i6 = elToDComp(i5)
								c_dVec(i6) = c_dVec(i6) + compStep
								call c_getBeamFrcMom(c_bfrcMom,c_bDef,i3,c_1*sPt)
								tVal = c2*(bfrcMom(comp) - objTgtExpanded(i4))**(c2-r_1)
								dsumTermdD(i6) = dsumTermdD(i6) + tVal*imag(c_bfrcMom(comp))*compStepInv
								c_dVec(i6) = c_dVec(i6) - compStep
							enddo
							i4 = i4 + 1
						enddo
						dLdD(:) = dLdD(:) + c1*dsumTermdD(:)
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getBeamFrcMom(bfrcMom,bDef,i3,sPt)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + bfrcMom(comp)*vol
							totVol = totVol + vol
							do i5 = elToDCRange(i3-1)+1, elToDCRange(i3)
							    i6 = elToDComp(i5)
								c_dVec(i6) = c_dVec(i6) + compStep
								call c_getBeamFrcMom(c_bfrcMom,c_bDef,i3,c_1*sPt)
								call c_getElMassVol(c_mass,c_vol,i3)
								tVal = compStepInv*(imag(c_bfrcMom(comp))*vol + bfrcMom(comp)*imag(c_vol))
								dsumTermdD(i6) = dsumTermdD(i6) + tVal
								dtotVoldD(i6) = dtotVoldD(i6) + compStepInv*imag(c_vol)
								c_dVec(i6) = c_dVec(i6) - compStep
							enddo
						enddo
						tVal = c2*c1*((sumTerm/totVol) - objTgtExpanded(1))**(c2-r_1)
						dLdD(:) = dLdD(:) + (tVal/totVol)*(dsumTermdD(:) - (sumTerm/totVol)*dtotVoldD(:))
					endif
				elseif(objCategory(i1) .eq. 'beamDef') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getBeamFrcMom(bfrcMom,bDef,i3,sPt)
							sumTerm = sumTerm + (bDef(comp) - objTgtExpanded(i4))**c2
							do i5 = elToDCRange(i3-1)+1, elToDCRange(i3)
							    i6 = elToDComp(i5)
								c_dVec(i6) = c_dVec(i6) + compStep
								call c_getBeamFrcMom(c_bfrcMom,c_bDef,i3,c_1*sPt)
								tVal = c2*(bDef(comp) - objTgtExpanded(i4))**(c2-r_1)
								dsumTermdD(i6) = dsumTermdD(i6) + tVal*imag(c_bDef(comp))*compStepInv
								c_dVec(i6) = c_dVec(i6) - compStep
							enddo
							i4 = i4 + 1
						enddo
						dLdD(:) = dLdD(:) + c1*dsumTermdD(:)
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getBeamFrcMom(bfrcMom,bDef,i3,sPt)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + bDef(comp)*vol
							totVol = totVol + vol
							do i5 = elToDCRange(i3-1)+1, elToDCRange(i3)
							    i6 = elToDComp(i5)
								c_dVec(i6) = c_dVec(i6) + compStep
								call c_getBeamFrcMom(c_bfrcMom,c_bDef,i3,c_1*sPt)
								call c_getElMassVol(c_mass,c_vol,i3)
								tVal = compStepInv*(imag(c_bDef(comp))*vol + bDef(comp)*imag(c_vol))
								dsumTermdD(i6) = dsumTermdD(i6) + tVal
								dtotVoldD(i6) = dtotVoldD(i6) + compStepInv*imag(c_vol)
								c_dVec(i6) = c_dVec(i6) - compStep
							enddo
						enddo
						tVal = c2*c1*((sumTerm/totVol) - objTgtExpanded(1))**(c2-r_1)
						dLdD(:) = dLdD(:) + (tVal/totVol)*(dsumTermdD(:) - (sumTerm/totVol)*dtotVoldD(:))
					endif
				elseif(objCategory(i1) .eq. 'flux') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElFlux(flux,dfluxdT,Tx,dTxdT,i3,sPt)
							sumTerm = sumTerm + (flux(comp) - objTgtExpanded(i4))**c2
							do i5 = elToDCRange(i3-1)+1, elToDCRange(i3)
							    i6 = elToDComp(i5)
								c_dVec(i6) = c_dVec(i6) + compStep
								call c_getElFlux(c_flux,c_dfluxdT,c_Tx,c_dTxdT,i3,c_1*sPt)
								tVal = c2*(flux(comp) - objTgtExpanded(i4))**(c2-r_1)
								dsumTermdD(i6) = dsumTermdD(i6) + tVal*compStepInv*imag(c_flux(comp))
								c_dVec(i6) = c_dVec(i6) - compStep
							enddo
							i4 = i4 + 1
						enddo
						dLdD(:) = dLdD(:) + c1*dsumTermdD(:)
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElFlux(flux,dfluxdT,Tx,dTxdT,i3,sPt)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + flux(comp)*vol
							totVol = totVol + vol
							do i5 = elToDCRange(i3-1)+1, elToDCRange(i3)
							    i6 = elToDComp(i5)
								c_dVec(i6) = c_dVec(i6) + compStep
								call c_getElFlux(c_flux,c_dfluxdT,c_Tx,c_dTxdT,i3,c_1*sPt)
								call c_getElMassVol(c_mass,c_vol,i3)
								tVal = compStepInv*(imag(c_flux(comp))*vol + flux(comp)*imag(c_vol))
								dsumTermdD(i6) = dsumTermdD(i6) + tVal
								dtotVoldD(i6) = dtotVoldD(i6) + compStepInv*imag(c_vol)
								c_dVec(i6) = c_dVec(i6) - compStep
							enddo
						enddo
						tVal = c2*c1*((sumTerm/totVol) - objTgtExpanded(1))**(c2-r_1)
						dLdD(:) = dLdD(:) + (tVal/totVol)*(dsumTermdD(:) - (sumTerm/totVol)*dtotVoldD(:))
					endif
				elseif(objCategory(i1) .eq. 'tempGradient') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElFlux(flux,dfluxdT,Tx,dTxdT,i3,sPt)
							sumTerm = sumTerm + (Tx(comp) - objTgtExpanded(i4))**c2
							do i5 = elToDCRange(i3-1)+1, elToDCRange(i3)
							    i6 = elToDComp(i5)
								c_dVec(i6) = c_dVec(i6) + compStep
								call c_getElFlux(c_flux,c_dfluxdT,c_Tx,c_dTxdT,i3,c_1*sPt)
								tVal = c2*(Tx(comp) - objTgtExpanded(i4))**(c2-r_1)
								dsumTermdD(i6) = dsumTermdD(i6) + tVal*compStepInv*imag(c_Tx(comp))
								c_dVec(i6) = c_dVec(i6) - compStep
							enddo
							i4 = i4 + 1
						enddo
						dLdD(:) = dLdD(:) + c1*dsumTermdD(:)
					elseif(objOperator(i1) .eq. 'volumeAverage') then
						totVol = r_0
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElFlux(flux,dfluxdT,Tx,dTxdT,i3,sPt)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + Tx(comp)*vol
							totVol = totVol + vol
							do i5 = elToDCRange(i3-1)+1, elToDCRange(i3)
							    i6 = elToDComp(i5)
								c_dVec(i6) = c_dVec(i6) + compStep
								call c_getElFlux(c_flux,c_dfluxdT,c_Tx,c_dTxdT,i3,c_1*sPt)
								call c_getElMassVol(c_mass,c_vol,i3)
								tVal = compStepInv*(imag(c_Tx(comp))*vol + Tx(comp)*imag(c_vol))
								dsumTermdD(i6) = dsumTermdD(i6) + tVal
								dtotVoldD(i6) = dtotVoldD(i6) + compStepInv*imag(c_vol)
								c_dVec(i6) = c_dVec(i6) - compStep
							enddo
						enddo
						tVal = c2*c1*((sumTerm/totVol) - objTgtExpanded(1))**(c2-r_1)
						dLdD(:) = dLdD(:) + (tVal/totVol)*(dsumTermdD(:) - (sumTerm/totVol)*dtotVoldD(:))
					endif
				elseif(objCategory(i1) .eq. 'mass') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
					if(objOperator(i1) .eq. 'volumeIntegral') then
					    do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + mass
							do i5 = elToDCRange(i3-1)+1, elToDCRange(i3)
							    i6 = elToDComp(i5)
								c_dVec(i6) = c_dVec(i6) + compStep
								call c_getElMassVol(c_mass,c_vol,i3)
								dsumTermdD(i6) = dsumTermdD(i6) + compStepInv*imag(c_mass)
								c_dVec(i6) = c_dVec(i6) - compStep
							enddo
						enddo
						tVal = c2*c1*(sumTerm - objTgtExpanded(1))**(c2-r_1)
						dLdD(:) = dLdD(:) + tVal*dsumTermdD(:)
					endif
				elseif(objCategory(i1) .eq. 'volume') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
					if(objOperator(i1) .eq. 'volumeIntegral') then
					    do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElMassVol(mass,vol,i3)
							sumTerm = sumTerm + vol
							do i5 = elToDCRange(i3-1)+1, elToDCRange(i3)
							    i6 = elToDComp(i5)
								c_dVec(i6) = c_dVec(i6) + compStep
								call c_getElMassVol(c_mass,c_vol,i3)
								dsumTermdD(i6) = dsumTermdD(i6) + compStepInv*imag(c_vol)
								c_dVec(i6) = c_dVec(i6) - compStep
							enddo
						enddo
						objVal(i1) = objVal(i1) + c1*(sumTerm - objTgtExpanded(1))**c2
						tVal = c2*c1*(sumTerm - objTgtExpanded(1))**(c2-r_1)
						dLdD(:) = dLdD(:) + tVal*dsumTermdD(:)
					endif
				elseif(objCategory(i1) .eq. 'massDisp') then
				    elSetNum = objElSetId(i1)
					c1 = objCoef(i1)
					c2 = objExp(i1)
				    comp = objComponent(i1)
					sPt(:) = r_0
					if(objOperator(i1) .eq. 'powerNorm') then
					    i4 = 1
						do i2 = elSetRange(elSetNum-1)+1, elSetRange(elSetNum)
						    i3 = elementSets(i2)
							call r_getElMassVol(mass,vol,i3)
							call r_getElementDisp(disp,ddispdU,i3,sPt)
							sumTerm = sumTerm + (mass*disp(comp) - objTgtExpanded(i4))**c2
							do i5 = elToDCRange(i3-1)+1, elToDCRange(i3)
							    i6 = elToDComp(i5)
								c_dVec(i6) = c_dVec(i6) + compStep
								call c_getElMassVol(c_mass,c_vol,i3)
								tVal = c2*(mass*disp(comp) - objTgtExpanded(i4))**(c2-r_1)
								dsumTermdD(i6) = dsumTermdD(i6) + compStepInv*imag(c_mass)*disp(comp)
								c_dVec(i6) = c_dVec(i6) - compStep
							enddo
							i4 = i4 + 1
						enddo
						dLdD(:) = dLdD(:) + c1*dsumTermdD(:)
					endif
				endif
			endif
		enddo
		
	end subroutine getdObjdD
	
end module ASenD_objective