module ASenD_output
    use ASenD_globalData
	use ASenD_constantVals
	use ASenD_input
	use ASenD_r_elementEqns
	use ASenD_r_designPropertyFunctions
	
	contains
	
!! General output
	
	subroutine writeNodeResults(fileName,fields,numFields,ndSet,nSLen,time,tStep)
	    implicit none
		
		integer, intent(in) :: numFields, nSlen, tStep
		real*8, intent(in) :: time
		character(len=128), intent(in) :: fileName
		character(len=16), intent(in) :: fields(10)
		integer, intent(in) :: ndSet(nSLen)
		
		character(len=16) :: intStr
		character(len=128) :: fltStr
		integer :: i1, i2, i3, i4, i5, i6
		
		open(unit=15, file=fileName, status='replace', action='write')
		
		write(15,*) 'time: ', time
		write(15,*) 'timeStep: ', tStep
		
		do i1 = 1, numFields
		    if(fields(i1) .eq. 'temperature') then
			    write(15,*) 'temperature:'
			    do i2 = 1, nSLen
				    i3 = ndSet(i2)
					i4 = currentRank(i3)
					write(intStr,'(I0)') i3 - 1
					write(fltStr,'(G18.12)') nodeTemp(i4)
					write(15,*) '    - [', trim(intStr), ', ', trim(fltStr), ']'
				enddo
			elseif(fields(i1) .eq. 'tdot') then
			    write(15,*) 'tdot:'
			    do i2 = 1, nSLen
				    i3 = ndSet(i2)
					i4 = currentRank(i3)
					write(intStr,'(I0)') i3 - 1
					write(fltStr,'(G18.12)') nodeTdot(i4)
					write(15,*) '    - [', trim(intStr), ', ', trim(fltStr), ']'
				enddo
			elseif(fields(i1) .eq. 'reactionHeatGen') then
			    write(15,*) 'reactionHeatGen:'
			    do i2 = 1, nSLen
				    i3 = ndSet(i2)
					i4 = currentRank(i3)
					write(intStr,'(I0)') i3 - 1
					write(fltStr,'(G18.12)') -thermalLoad(i4)
					write(15,*) '    - [', trim(intStr), ', ', trim(fltStr), ']'
				enddo
			elseif(fields(i1) .eq. 'displacement') then
			    write(15,*) 'displacement:'
			    do i2 = 1, nSLen
				    i3 = ndSet(i2)
					write(intStr,'(I0)') i3 - 1
					i6 = 1
					do i4 = 1, 6
					    i5 = nDofIndex(i4,i3)
						if(i5 .gt. 0) then
						    write(fltStr(i6:i6+17),'(G18.12)') nodeDisp(i5)
							i6 = i6 + 18
							fltStr(i6:i6) = ','
							i6 = i6 + 1
						else
						    write(fltStr(i6:i6+17),'(G18.12)') r_0
							i6 = i6 + 18
							fltStr(i6:i6) = ','
							i6 = i6 + 1
						endif
					enddo
					write(15,*) '    - [', trim(intStr), ', ', fltStr(1:i6-2), ']'
				enddo
			elseif(fields(i1) .eq. 'velocity') then
			    write(15,*) 'velocity:'
			    do i2 = 1, nSLen
				    i3 = ndSet(i2)
					write(intStr,'(I0)') i3 - 1
					i6 = 1
					do i4 = 1, 6
					    i5 = nDofIndex(i4,i3)
						if(i5 .gt. 0) then
						    write(fltStr(i6:i6+17),'(G18.12)') nodeVel(i5)
							i6 = i6 + 18
							fltStr(i6:i6) = ','
							i6 = i6 + 1
						else
						    write(fltStr(i6:i6+17),'(G18.12)') r_0
							i6 = i6 + 18
							fltStr(i6:i6) = ','
							i6 = i6 + 1
						endif
					enddo
					write(15,*) '    - [', trim(intStr), ', ', fltStr(1:i6-2), ']'
				enddo
			elseif(fields(i1) .eq. 'acceleration') then
			    write(15,*) 'acceleration:'
			    do i2 = 1, nSLen
				    i3 = ndSet(i2)
					write(intStr,'(I0)') i3 - 1
					i6 = 1
					do i4 = 1, 6
					    i5 = nDofIndex(i4,i3)
						if(i5 .gt. 0) then
						    write(fltStr(i6:i6+17),'(G18.12)') nodeAcc(i5)
							i6 = i6 + 18
							fltStr(i6:i6) = ','
							i6 = i6 + 1
						else
						    write(fltStr(i6:i6+17),'(G18.12)') r_0
							i6 = i6 + 18
							fltStr(i6:i6) = ','
							i6 = i6 + 1
						endif
					enddo
					write(15,*) '    - [', trim(intStr), ', ', fltStr(1:i6-2), ']'
				enddo
			elseif(fields(i1) .eq. 'reactionForce') then
			    write(15,*) 'reactionForce:'
			    do i2 = 1, nSLen
				    i3 = ndSet(i2)
					write(intStr,'(I0)') i3 - 1
					i6 = 1
					do i4 = 1, 6
					    i5 = nDofIndex(i4,i3)
						if(i5 .gt. 0) then
						    write(fltStr(i6:i6+17),'(G18.12)') -elasticLoad(i5)
							i6 = i6 + 18
							fltStr(i6:i6) = ','
							i6 = i6 + 1
						else
						    write(fltStr(i6:i6+17),'(G18.12)') r_0
							i6 = i6 + 18
							fltStr(i6:i6) = ','
							i6 = i6 + 1
						endif
					enddo
					write(15,*) '    - [', trim(intStr), ', ', fltStr(1:i6-2), ']'
				enddo
			else
			    write(lfUnit,*) 'Warning: Unrecognized field name ', fields(i1), ' listed in writeNodeResults'
				write(lfUnit,*) ' '
			endif
		enddo
		
		close(15)
		
	end subroutine writeNodeResults
	
	subroutine writeElementResults(fileName,fields,numFields,elSet,eSLen,time,tStep)
	    implicit none
		
		integer, intent(in) :: numFields, eSLen, tStep
		real*8, intent(in) :: time
		character(len=128), intent(in) :: fileName
		character(len=16), intent(in) :: fields(10)
		integer, intent(in) :: elSet(eSLen)
		
		integer :: numNds, dofPerNd, numIntDof, numIntPts
		integer :: dofTable(2,33)
		real*8 :: intPts(3,8), ipWt(8)
		
		real*8 :: flux(3), dfluxdT(3,10), Tx(3), dTxdT(3,10)
		
		real*8 :: frcMom(9), shDef(9), ptT
		
		character(len=16) :: elNumStr, secNumStr, ipStr, layNumStr
        character(len=32) :: intStr
		character(len=256) :: fltStr
		real*8 :: stress(6), strain(6), sEDen, strength, fIndex, maxFInd
		integer :: i1, i2, i3, i4, i5, i6, i7, eType, secNum, matID, layerInd, numLayers
		
		open(unit=138, file=fileName, status='replace', action='write')
		
		do i1 = 1, numFields
		    if(fields(i1) .eq. 'stress') then
			    write(138,*) 'stress:'
				write(138,*) '## element, section, integration point, layer, S11, S22, S33, S12, S13, S23, max index '
				do i2 = 1, eSLen
				    i3 = elSet(i2)
					eType = elementType(i3)
					secNum = elementSection(i3)
					numLayers = secLayupRange(secNum) - secLayupRange(secNum-1)
					call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
					if(eType .eq. 41 .or. eType .eq. 3) then
					    do i4 = 1, numIntPts
						    do i5 = 1, numLayers
							    call r_getElStress(stress,strain,i3,i5,intPts(:,i4))
								write(elNumStr,'(I0)') i3 - 1
								write(secNumStr,'(I0)') secNum - 1
								write(ipStr,'(I0)') i4 - 1
								write(layNumStr,'(I0)') i5 - 1
								intStr = trim(elNumStr) // ', ' // trim(secNumStr) // ', ' // trim(ipStr) // ', ' // trim(layNumStr)
                                i6 = 1
                                do i7 = 1, 6
								    write(fltStr(i6:i6+17),'(G18.12)') stress(i7)
									i6 = i6 + 18
									fltStr(i6:i6) = ','
									i6 = i6 + 1
                                enddo
								layerInd = secLayupRange(secNum-1) + i5
								matID = layupMatId(layerInd)
								maxFInd = r_0
								do i7 = 1, 3
								    strength = materialMaxStress(i7,matID)
								    if(strength .ne. r_0) then
									    fIndex = stress(i7)/strength
										if(fIndex .gt. maxFInd) then
										    maxFInd = fIndex
										endif
									endif
									strength = materialMaxStress(i7+3,matID)
								    if(strength .ne. r_0) then
									    fIndex = stress(i7)/strength
										if(fIndex .gt. maxFInd) then
										    maxFInd = fIndex
										endif
									endif
									strength = materialMaxStress(i7+6,matID)
								    if(strength .ne. r_0) then
									    fIndex = abs(stress(i7+3))/strength
										if(fIndex .gt. maxFInd) then
										    maxFInd = fIndex
										endif
									endif
								enddo
								if(maxFInd .eq. r_0) then
								    fltStr(i6:i6+3) = ' NA,'
									i6 = i6 + 4
								else
								    write(fltStr(i6:i6+17),'(G18.12)') maxFInd
									i6 = i6 + 18
									fltStr(i6:i6) = ','
									i6 = i6 + 1
								endif
								write(138,*) '    - [', trim(intStr), ', ', fltStr(1:i6-2), ']'
							enddo
						enddo
					elseif(eType .eq. 2) then
					else
					    do i4 = 1, numIntPts
							call r_getElStress(stress,strain,i3,0,intPts(:,i4))
							write(elNumStr,'(I0)') i3 - 1
							write(secNumStr,'(I0)') secNum - 1
							write(ipStr,'(I0)') i4 - 1
							layNumStr = ' NA'
							intStr = trim(elNumStr) // ', ' // secNumStr// ', ' // trim(ipStr) // ', ' // trim(layNumStr)
							i6 = 1
							do i7 = 1, 6
								write(fltStr(i6:i6+17),'(G18.12)') stress(i7)
								i6 = i6 + 18
								fltStr(i6:i6) = ','
								i6 = i6 + 1
							enddo
							matID = sectionMatId(secNum)
							maxFInd = r_0
							do i7 = 1, 3
								strength = materialMaxStress(i7,matID)
								if(strength .ne. r_0) then
									fIndex = stress(i7)/strength
									if(fIndex .gt. maxFInd) then
										maxFInd = fIndex
									endif
								endif
								strength = materialMaxStress(i7+3,matID)
								if(strength .ne. r_0) then
									fIndex = stress(i7)/strength
									if(fIndex .gt. maxFInd) then
										maxFInd = fIndex
									endif
								endif
								strength = materialMaxStress(i7+6,matID)
								if(strength .ne. r_0) then
									fIndex = abs(stress(i7+3))/strength
									if(fIndex .gt. maxFInd) then
										maxFInd = fIndex
									endif
								endif
							enddo
							if(maxFInd .eq. r_0) then
								fltStr(i6:i6+3) = ' NA,'
								i6 = i6 + 4
							else
								write(fltStr(i6:i6+17),'(G18.12)') maxFInd
								i6 = i6 + 18
								fltStr(i6:i6) = ','
								i6 = i6 + 1
							endif
							write(138,*) '    - [', trim(intStr), ', ', fltStr(1:i6-2), ']'
						enddo
					endif
				enddo
			elseif(fields(i1) .eq. 'strain') then
			    write(138,*) 'strain:'
				write(138,*) '## element, section, integration point, layer, E11, E22, E33, E12, E13, E23, max index '
				do i2 = 1, eSLen
				    i3 = elSet(i2)
					eType = elementType(i3)
					secNum = elementSection(i3)
					numLayers = secLayupRange(secNum) - secLayupRange(secNum-1)
					call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
					if(eType .eq. 41 .or. eType .eq. 3) then
					    do i4 = 1, numIntPts
						    do i5 = 1, numLayers
							    call r_getElStress(stress,strain,i3,i5,intPts(:,i4))
								write(elNumStr,'(I0)') i3 - 1
								write(secNumStr,'(I0)') secNum - 1
								write(ipStr,'(I0)') i4 - 1
								write(layNumStr,'(I0)') i5 - 1
								intStr = trim(elNumStr) // ', ' // trim(secNumStr) // ', ' // trim(ipStr) // ', ' // trim(layNumStr)
                                i6 = 1
                                do i7 = 1, 6
								    write(fltStr(i6:i6+17),'(G18.12)') strain(i7)
									i6 = i6 + 18
									fltStr(i6:i6) = ','
									i6 = i6 + 1
                                enddo
								layerInd = secLayupRange(secNum-1) + i5
								matID = layupMatId(layerInd)
								maxFInd = r_0
								do i7 = 1, 3
								    strength = materialMaxStrain(i7,matID)
								    if(strength .ne. r_0) then
									    fIndex = strain(i7)/strength
										if(fIndex .gt. maxFInd) then
										    maxFInd = fIndex
										endif
									endif
									strength = materialMaxStrain(i7+3,matID)
								    if(strength .ne. r_0) then
									    fIndex = strain(i7)/strength
										if(fIndex .gt. maxFInd) then
										    maxFInd = fIndex
										endif
									endif
									strength = materialMaxStrain(i7+6,matID)
								    if(strength .ne. r_0) then
									    fIndex = abs(strain(i7+3))/strength
										if(fIndex .gt. maxFInd) then
										    maxFInd = fIndex
										endif
									endif
								enddo
								if(maxFInd .eq. r_0) then
								    fltStr(i6:i6+3) = ' NA,'
									i6 = i6 + 4
								else
								    write(fltStr(i6:i6+17),'(G18.12)') maxFInd
									i6 = i6 + 18
									fltStr(i6:i6) = ','
									i6 = i6 + 1
								endif
								write(138,*) '    - [', trim(intStr), ', ', fltStr(1:i6-2), ']'
							enddo
						enddo
					elseif(eType .eq. 2) then
					else
					    do i4 = 1, numIntPts
							call r_getElStress(stress,strain,i3,0,intPts(:,i4))
							write(elNumStr,'(I0)') i3 - 1
							write(secNumStr,'(I0)') secNum - 1
							write(ipStr,'(I0)') i4 - 1
							layNumStr = ' NA'
							intStr = trim(elNumStr) // ', ' // secNumStr// ', ' // trim(ipStr) // ', ' // trim(layNumStr)
							i6 = 1
							do i7 = 1, 6
								write(fltStr(i6:i6+17),'(G18.12)') strain(i7)
								i6 = i6 + 18
								fltStr(i6:i6) = ','
								i6 = i6 + 1
							enddo
							matID = sectionMatId(secNum)
							maxFInd = r_0
							do i7 = 1, 3
								strength = materialMaxStrain(i7,matID)
								if(strength .ne. r_0) then
									fIndex = strain(i7)/strength
									if(fIndex .gt. maxFInd) then
										maxFInd = fIndex
									endif
								endif
								strength = materialMaxStrain(i7+3,matID)
								if(strength .ne. r_0) then
									fIndex = strain(i7)/strength
									if(fIndex .gt. maxFInd) then
										maxFInd = fIndex
									endif
								endif
								strength = materialMaxStrain(i7+6,matID)
								if(strength .ne. r_0) then
									fIndex = abs(strain(i7+3))/strength
									if(fIndex .gt. maxFInd) then
										maxFInd = fIndex
									endif
								endif
							enddo
							if(maxFInd .eq. r_0) then
								fltStr(i6:i6+3) = ' NA,'
								i6 = i6 + 4
							else
								write(fltStr(i6:i6+17),'(G18.12)') maxFInd
								i6 = i6 + 18
								fltStr(i6:i6) = ','
								i6 = i6 + 1
							endif
							write(138,*) '    - [', trim(intStr), ', ', fltStr(1:i6-2), ']'
						enddo
					endif
				enddo
            elseif(fields(i1) .eq. 'strainEnergyDen') then
			    write(138,*) 'strainEnergyDen:'
				write(138,*) '## element, section, integration point, layer, strain energy density, max Index '
				do i2 = 1, eSLen
				    i3 = elSet(i2)
					eType = elementType(i3)
					secNum = elementSection(i3)
					numLayers = secLayupRange(secNum) - secLayupRange(secNum-1)
					call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
					if(eType .eq. 41 .or. eType .eq. 3) then
					    do i4 = 1, numIntPts
						    do i5 = 1, numLayers
							    call r_getElStress(stress,strain,i3,i5,intPts(:,i4))
								write(elNumStr,'(I0)') i3 - 1
								write(secNumStr,'(I0)') secNum - 1
								write(ipStr,'(I0)') i4 - 1
								write(layNumStr,'(I0)') i5 - 1
								intStr = trim(elNumStr)//', '//trim(secNumStr)//', '//trim(ipStr)//', '//trim(layNumStr)
                                sEDen = r_0
                                do i7 = 1, 6
								    sEDen = sEDen + r_p5*stress(i7)*strain(i7)
                                enddo
								write(fltStr(1:18),'(G18.12)') sEDen
								fltStr(19:19) = ','
								layerInd = secLayupRange(secNum-1) + i5
								matID = layupMatId(layerInd)
								strength = materialMaxStrnEngy(matID)
								if(strength .ne. r_0) then
								    maxFInd = sEDen/strength
									write(fltStr(20:37),'(G18.12)') maxFInd
								else
								    fltStr(20:37) = ' NA               '
								endif
								write(138,*) '    - [', trim(intStr), ', ', fltStr(1:37), ']'
							enddo
						enddo
					elseif(eType .eq. 2) then
					else
					    do i4 = 1, numIntPts
							call r_getElStress(stress,strain,i3,0,intPts(:,i4))
							write(elNumStr,'(I0)') i3 - 1
							write(secNumStr,'(I0)') secNum - 1
							write(ipStr,'(I0)') i4 - 1
							layNumStr = ' NA'
							intStr = trim(elNumStr) // ', ' // trim(secNumStr) // ', ' // trim(ipStr) // ', ' // trim(layNumStr)
							sEDen = r_0
							do i7 = 1, 6
								sEDen = sEDen + r_p5*stress(i7)*strain(i7)
							enddo
							write(fltStr(1:18),'(G18.12)') sEDen
							fltStr(19:19) = ','
							matID = sectionMatId(secNum)
							strength = materialMaxStrnEngy(matID)
							if(strength .ne. r_0) then
								maxFInd = sEDen/strength
								write(fltStr(20:37),'(G18.12)') maxFInd
							else
								fltStr(20:37) = ' NA               '
							endif							
						    write(138,*) '    - [', trim(intStr), ', ', fltStr(1:37), ']'
						enddo
					endif
				enddo
			elseif(fields(i1) .eq. 'heatFlux') then
			    write(138,*) 'heatFlux:'
				write(138,*) '## element, section, integration point, q1, q2, q3 '
				do i2 = 1, eSLen
				    i3 = elSet(i2)
					eType = elementType(i3)
					secNum = elementSection(i3)
					call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
					do i4 = 1, numIntPts
						call r_getElFlux(flux,dfluxdT,Tx,dTxdT,i3,intPts(:,i4))
						write(elNumStr,'(I0)') i3 - 1
						write(secNumStr,'(I0)') secNum - 1
						write(ipStr,'(I0)') i4 - 1
						intStr = trim(elNumStr) // ', ' // trim(secNumStr) // ', ' // trim(ipStr)
						i6 = 1
						do i7 = 1, 3
							write(fltStr(i6:i6+17),'(G18.12)') flux(i7)
							i6 = i6 + 18
							fltStr(i6:i6) = ','
							i6 = i6 + 1
						enddo
						write(138,*) '    - [', trim(intStr), ', ', fltStr(1:i6-2), ']'
					enddo
				enddo
			elseif(fields(i1) .eq. 'tempGradient') then
			    write(138,*) 'tempGradient:'
				write(138,*) '## element, section, integration point, Tg1, Tg2, Tg3 '
				do i2 = 1, eSLen
				    i3 = elSet(i2)
					eType = elementType(i3)
					secNum = elementSection(i3)
					call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
					do i4 = 1, numIntPts
						call r_getElFlux(flux,dfluxdT,Tx,dTxdT,i3,intPts(:,i4))
						write(elNumStr,'(I0)') i3 - 1
						write(secNumStr,'(I0)') secNum - 1
						write(ipStr,'(I0)') i4 - 1
						intStr = trim(elNumStr) // ', ' // trim(secNumStr) // ', ' // trim(ipStr)
						i6 = 1
						do i7 = 1, 3
							write(fltStr(i6:i6+17),'(G18.12)') Tx(i7)
							i6 = i6 + 18
							fltStr(i6:i6) = ','
							i6 = i6 + 1
						enddo
						write(138,*) '    - [', trim(intStr), ', ', fltStr(1:i6-2), ']'
					enddo
				enddo
			elseif(fields(i1) .eq. 'shellDeformation') then
			    write(138,*) 'shellDeformation:'
				write(138,*) '## element, section, integration point, E11, E22, E12, K11, K22, K12 '
				do i2 = 1, eSLen
				    i3 = elSet(i2)
					eType = elementType(i3)
					if(eType .eq. 41 .or. eType .eq. 3) then
						secNum = elementSection(i3)
						call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
						do i4 = 1, numIntPts
							call r_getShellFrcMom(frcMom,shDef,ptT,i3,intPts(:,i4))
							write(elNumStr,'(I0)') i3 - 1
							write(secNumStr,'(I0)') secNum - 1
							write(ipStr,'(I0)') i4 - 1
							intStr = trim(elNumStr) // ', ' // trim(secNumStr) // ', ' // trim(ipStr)
							i6 = 1
							do i7 = 1, 6
								write(fltStr(i6:i6+17),'(G18.12)') shDef(i7)
								i6 = i6 + 18
								fltStr(i6:i6) = ','
								i6 = i6 + 1
							enddo
							write(138,*) '    - [', trim(intStr), ', ', fltStr(1:i6-2), ']'
						enddo
					endif
				enddo
			elseif(fields(i1) .eq. 'shellForceMoment') then
			    write(138,*) 'shellForceMoment:'
				write(138,*) '## element, section, integration point, F11, F22, F12, M11, M22, M12 '
				do i2 = 1, eSLen
				    i3 = elSet(i2)
					eType = elementType(i3)
					if(eType .eq. 41 .or. eType .eq. 3) then
						secNum = elementSection(i3)
						call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
						do i4 = 1, numIntPts
							call r_getShellFrcMom(frcMom,shDef,ptT,i3,intPts(:,i4))
							write(elNumStr,'(I0)') i3 - 1
							write(secNumStr,'(I0)') secNum - 1
							write(ipStr,'(I0)') i4 - 1
							intStr = trim(elNumStr) // ', ' // trim(secNumStr) // ', ' // trim(ipStr)
							i6 = 1
							do i7 = 1, 6
								write(fltStr(i6:i6+17),'(G18.12)') frcMom(i7)
								i6 = i6 + 18
								fltStr(i6:i6) = ','
								i6 = i6 + 1
							enddo
							write(138,*) '    - [', trim(intStr), ', ', fltStr(1:i6-2), ']'
						enddo
					endif
				enddo
			elseif(fields(i1) .eq. 'beamDeformation') then
			    write(138,*) 'beamDeformation:'
				write(138,*) '## element, section, integration point, E11, E21, E31, K11, K21, K31 '
				do i2 = 1, eSLen
				    i3 = elSet(i2)
					eType = elementType(i3)
					if(eType .eq. 2) then
						secNum = elementSection(i3)
						call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
						do i4 = 1, numIntPts
							call r_getBeamFrcMom(frcMom,shDef,i3,intPts(:,i4))
							write(elNumStr,'(I0)') i3 - 1
							write(secNumStr,'(I0)') secNum - 1
							write(ipStr,'(I0)') i4 - 1
							intStr = trim(elNumStr) // ', ' // trim(secNumStr) // ', ' // trim(ipStr)
							i6 = 1
							do i7 = 1, 6
								write(fltStr(i6:i6+17),'(G18.12)') shDef(i7)
								i6 = i6 + 18
								fltStr(i6:i6) = ','
								i6 = i6 + 1
							enddo
							write(138,*) '    - [', trim(intStr), ', ', fltStr(1:i6-2), ']'
						enddo
					endif
				enddo
			elseif(fields(i1) .eq. 'beamForceMoment') then
			    write(138,*) 'beamForceMoment:'
				write(138,*) '## element, section, integration point, F11, F21, F31, M11, M21, M31 '
				do i2 = 1, eSLen
				    i3 = elSet(i2)
					eType = elementType(i3)
					if(eType .eq. 2) then
						secNum = elementSection(i3)
						call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
						do i4 = 1, numIntPts
							call r_getBeamFrcMom(frcMom,shDef,i3,intPts(:,i4))
							write(elNumStr,'(I0)') i3 - 1
							write(secNumStr,'(I0)') secNum - 1
							write(ipStr,'(I0)') i4 - 1
							intStr = trim(elNumStr) // ', ' // trim(secNumStr) // ', ' // trim(ipStr)
							i6 = 1
							do i7 = 1, 6
								write(fltStr(i6:i6+17),'(G18.12)') frcMom(i7)
								i6 = i6 + 18
								fltStr(i6:i6) = ','
								i6 = i6 + 1
							enddo
							write(138,*) '    - [', trim(intStr), ', ', fltStr(1:i6-2), ']'
						enddo
					endif
				enddo
			else
			    write(lfUnit,*) 'Warning: Unrecognized field name ', fields(i1), ' listed in writeNodeResults'
				write(lfUnit,*) ' '
			endif
		enddo
		
		close(138)
		
	end subroutine writeElementResults
	
	subroutine writeModalResults(fileName,writeModes)
	    implicit none

        integer, intent(in) :: writeModes
		character(len=128), intent(in) :: fileName

        character(len=32) :: intStr		
		character(len=128) :: fltStr
		character(len=256) :: writeLine
		integer :: i1, i2, i3, i4, i5, i6
		
		open(unit=900, file=fileName, status='replace', action='write')
		
		if(modalType .eq. 0) then
		    write(900,*) '# Modal buckling analysis results'
		else
		    write(900,*) '# Modal frequency analysis results'
		endif
		
		write(900,*) 'eigenValues:'
		if(modalType .eq. 0) then
			write(900,*) '# mode number, eigenvalue, load factor'
		else
		    write(900,*) '# mode number, eigenvalue, natural frequency'
		endif
		do i1 = 1, numEigModes
		    write(intStr,'(I0)') i1 - 1
			intStr = adjustl(intStr)
			write(fltStr(1:18),'(G18.12)') eigenVals(i1)
			fltStr(19:20) = ', '
			write(fltStr(21:38),'(G18.12)') eigenFactors(i1)
			writeLine = '    -[' // trim(intStr) // ', ' // fltStr(1:38) // ']'
			write(900,*) trim(writeLine)
		enddo
		
		if(writeModes .eq. 1) then
		    write(900,*) 'eigenModes:'
			do i1 = 1, numEigModes
			    write(900,*) '    - mode:', i1 - 1
				write(900,*) '      displacement:'
				do i2 = 1, numNodes
				    write(intStr,'(I0)') i2 - 1
					i5 = 1
					do i3 = 1, 6
					    i4 = nDofIndex(i3,i2)
						if(i4 .gt. 0) then
						    write(fltStr(i5:i5+17),'(G18.12)') eigenModes(i4,i1)
							i5 = i5+18
							fltStr(i5:i5) = ','
							i5 = i5 + 1
						else
						    write(fltStr(i5:i5+17),'(G18.12)') r_0
							i5 = i5+18
							fltStr(i5:i5) = ','
							i5 = i5 + 1
						endif
					enddo
					writeLine = '          - [' // trim(intStr) // ', ' // fltStr(1:i5-2) // ']'
					write(900,*) trim(writeLine)
				enddo
			enddo
		endif
		
		close(900)
		
	end subroutine writeModalResults
	
	subroutine writeNodeCoord(fileName,nSet,nSLen)
	    implicit none
		
		integer, intent(in) :: nSLen
		integer, intent(in) :: nSet(nSLen)
		character(len=128), intent(in) :: fileName
		
		character(len=32) :: intStr
		character(len=64) :: fltStr
		character(len=128) :: writeLine
		real*8 :: ndCrd(3)
		integer :: i1, i2, i3, i4, i5
		
		open(unit=965, file=fileName, status='replace', action='write')
		
		write(965,*) 'nodes:'
		
		do i5 = 1, nSLen
		    i1 = nSet(i5)
		    ndCrd(:) = nodeList(:,i1)
			do i2 = ndToDRange(i1-1)+1, ndToDRange(i1)
				i3 = ndToD(i2)
				if(dCategory(i3) .eq. 'nodeCoord') then
					i4 = dComponent(i3)
					ndCrd(i4) = ndCrd(i4) + ndToCoef(i2)*r_dVec(i3)
				endif
			enddo
			write(intStr,'(I0)') i1 - 1
			i3 = 1
			do i2 = 1, 3
			    write(fltStr(i3:i3+17),'(G18.12)') ndCrd(i2)
				i3 = i3 + 18
				fltStr(i3:i3) = ','
				i3 = i3 + 1
			enddo
			writeLine = '    - [' // trim(intStr) // ', ' // fltStr(1:i3-2) // ']'
			write(965,*) trim(writeLine)
		enddo
		
		close(965)
		
	end subroutine writeNodeCoord
	
	subroutine writeElProperties(fileName,elSet,eSLen,propList,numProps)
	    implicit none
		
		integer, intent(in) :: eSLen, numProps
		integer, intent(in) :: elSet(eSLen)
		character(len=128), intent(in) :: fileName
		character(len=32), intent(in) :: propList(numProps)
		
		character(len=512) :: writeLine, fltStr
		character(len=16) :: intStr
		real*8 :: elProps(21), tCondMat(3,3), orient(3,3), nodes(3,10), CMat(6,6), ABDMat(9,9)
		integer :: i1, i2, i3, i4, i5, i6, i7, i8, eType, eSec, matID, layerNum
		
		open(unit=1010, file=fileName, status='replace', action='write')
		
		write(1010,*) 'elementProperties:'
		
		do i1 = 1, eSLen
		    i2 = elSet(i1)
			write(1010,*) '    - element: ', (i2-1)
			eType = elementType(i2)
			eSec = elementSection(i2)
			do i3 = 1, numProps
			    if(propList(i3) .eq. 'density') then
				    if(eType .eq. 41 .or. eType .eq. 3) then
					    write(1010,*) '      density:'
						write(1010,*) '      #   layer, value'
						layerNum = 0
						do i4 = secLayupRange(eSec-1)+1, secLayupRange(eSec)
						    layerNum = layerNum + 1
							matID = layupMatId(i4)
							elProps(1) = materialDensity(matID)
							do i5 = elToDRange(i2-1) + 1, elToDRange(i2)
								i6 = elToD(i5)
								if(dCategory(i6) .eq. 'density' .and. dLayer(i6) .eq. layerNum) then
									elProps(1) = elProps(1) + elToCoef(i5)*r_dVec(i6)
								endif
							enddo
							write(1010,*) '          - [', (layerNum-1), ',', elProps(1), ']'
						enddo
					else
					    matID = sectionMatId(eSec)
						if(matID .gt. 0) then
							elProps(1) = materialDensity(matID)
							do i4 = elToDRange(i2-1) + 1, elToDRange(i2)
								i5 = elToD(i4)
								if(dCategory(i5) .eq. 'density') then
									elProps(1) = elProps(1) + elToCoef(i4)*r_dVec(i5)
								endif
							enddo
							write(1010,*) '      density: ', elProps(1)
						else
						    write(lfUnit,*) 'Warning: element ', i2, ' in section ', eSec, ' has no material definition.'
							write(lfUnit,*) 'Cannont write density to element property output file'
						endif
					endif
				elseif(propList(i3) .eq. 'specHeat') then
				    if(eType .eq. 41 .or. eType .eq. 3) then
					    write(1010,*) '      specHeat:'
						write(1010,*) '      #   layer, value'
						layerNum = 0
						do i4 = secLayupRange(eSec-1)+1, secLayupRange(eSec)
						    layerNum = layerNum + 1
							matID = layupMatId(i4)
							elProps(1) = materialSpecHeat(matID)
							do i5 = elToDRange(i2-1) + 1, elToDRange(i2)
								i6 = elToD(i5)
								if(dCategory(i6) .eq. 'specHeat' .and. dLayer(i6) .eq. layerNum) then
									elProps(1) = elProps(1) + elToCoef(i5)*r_dVec(i6)
								endif
							enddo
							write(1010,*) '          - [', (layerNum-1), ',', elProps(1), ']'
						enddo
					elseif(eType .eq. 2 .and. abs(beamSpecHeat(eSec)) .gt. r_0) then
						call r_getBeamSpecHeat(elProps(1),i2)
						write(1010,*) '      sectnSpecHeat: ', elProps(1)
					else
					    matID = sectionMatId(eSec)
						if(matID .gt. 0) then
							elProps(1) = materialSpecHeat(matID)
							do i4 = elToDRange(i2-1) + 1, elToDRange(i2)
								i5 = elToD(i4)
								if(dCategory(i5) .eq. 'specHeat') then
									elProps(1) = elProps(1) + elToCoef(i4)*r_dVec(i5)
								endif
							enddo
							write(1010,*) '      specHeat: ', elProps(1)
						else
							write(lfUnit,*) 'Warning: element ', i2, ' in section ', eSec, ' has no material definition.'
							write(lfUnit,*) 'Cannot write specific heat capacity to element property output file.'
						endif
					endif
				elseif(propList(i3) .eq. 'thickness') then
				    if(eType .eq. 41 .or. eType .eq. 3) then
					    write(1010,*) '      thickness:'
						write(1010,*) '      #   layer, value'
						layerNum = 0
						do i4 = secLayupRange(eSec-1)+1, secLayupRange(eSec)
						    layerNum = layerNum + 1
							elProps(1) = layupThickness(i4)
							do i5 = elToDRange(i2-1) + 1, elToDRange(i2)
								i6 = elToD(i5)
								if(dCategory(i6) .eq. 'thickness' .and. dLayer(i6) .eq. layerNum) then
									elProps(1) = elProps(1) + elToCoef(i5)*r_dVec(i6)
								endif
							enddo
							write(1010,*) '          - [', (layerNum-1), ',', elProps(1), ']'
						enddo
					endif
				elseif(propList(i3) .eq. 'angle') then
				    if(eType .eq. 41 .or. eType .eq. 3) then
					    write(1010,*) '      angle:'
						write(1010,*) '      #   layer, value'
						layerNum = 0
						do i4 = secLayupRange(eSec-1)+1, secLayupRange(eSec)
						    layerNum = layerNum + 1
							elProps(1) = layupAngle(i4)
							do i5 = elToDRange(i2-1) + 1, elToDRange(i2)
								i6 = elToD(i5)
								if(dCategory(i6) .eq. 'thickness' .and. dLayer(i6) .eq. layerNum) then
									elProps(1) = elProps(1) + elToCoef(i5)*r_dVec(i6)
								endif
							enddo
							write(1010,*) '          - [', (layerNum-1), ',', elProps(1), ']'
						enddo
					endif
				elseif(propList(i3) .eq. 'zOffset') then
				    if(eType .eq. 41 .or. eType .eq. 3) then
						elProps(1) = sectionZOffset(eSec)
						do i5 = elToDRange(i2-1) + 1, elToDRange(i2)
							i6 = elToD(i5)
							if(dCategory(i6) .eq. 'zOffset') then
								elProps(1) = elProps(1) + elToCoef(i5)*r_dVec(i6)
							endif
						enddo
						write(1010,*) '      zOffset: ', elProps(1)
					endif
				elseif(propList(i3) .eq. 'area') then
				    if(eType .eq. 2) then
						elProps(1) = beamProperties(1,eSec)
						do i5 = elToDRange(i2-1) + 1, elToDRange(i2)
							i6 = elToD(i5)
							if(dCategory(i6) .eq. 'area') then
								elProps(1) = elProps(1) + elToCoef(i5)*r_dVec(i6)
							endif
						enddo
						write(1010,*) '      area: ', elProps(1)
					endif
				elseif(propList(i3) .eq. 'polarMoment') then
				    if(eType .eq. 2) then
						elProps(1) = beamProperties(7,eSec)
						do i5 = elToDRange(i2-1) + 1, elToDRange(i2)
							i6 = elToD(i5)
							if(dCategory(i6) .eq. 'polarMoment') then
								elProps(1) = elProps(1) + elToCoef(i5)*r_dVec(i6)
							endif
						enddo
						write(1010,*) '      polarMoment: ', elProps(1)
					endif
				elseif(propList(i3) .eq. 'modulus') then
				    if(eType .eq. 41 .or. eType .eq. 3) then
					    write(1010,*) '      modulus:'
						write(1010,*) '      #   layer, value'
						layerNum = 0
						do i4 = secLayupRange(eSec-1)+1, secLayupRange(eSec)
						    layerNum = layerNum + 1
							matID = layupMatId(i4)
							elProps(1:3) = materialElastic(1:3,matID)
							do i5 = elToDRange(i2-1) + 1, elToDRange(i2)
								i6 = elToD(i5)
								if(dCategory(i6) .eq. 'modulus' .and. dLayer(i6) .eq. layerNum) then
								    i7 = dComponent(i6)
									elProps(i7) = elProps(i7) + elToCoef(i5)*r_dVec(i6)
								endif
							enddo
							write(intStr,'(I0)') layerNum - 1
							i8 = 1
							do i5 = 1, 3
							    write(fltStr(i8:i8+17),'(G18.12)') elProps(i5)
								i8 = i8 + 18
								fltStr(i8:i8) = ','
								i8 = i8 + 1
							enddo
							writeLine = '          - [' // trim(intStr) // ', ' // fltStr(1:i8-2) // ']'
							write(1010,*) trim(writeLine)
						enddo
					else
					    matID = sectionMatId(eSec)
						if(matID .gt. 0) then
							elProps(1:3) = materialElastic(1:3,matID)
							do i4 = elToDRange(i2-1) + 1, elToDRange(i2)
								i5 = elToD(i4)
								if(dCategory(i5) .eq. 'modulus') then
								    i7 = dComponent(i5)
									elProps(i7) = elProps(i7) + elToCoef(i4)*r_dVec(i5)
								endif
							enddo
							i8 = 1
							do i5 = 1, 3
							    write(fltStr(i8:i8+17),'(G18.12)') elProps(i5)
								i8 = i8 + 18
								fltStr(i8:i8) = ','
								i8 = i8 + 1
							enddo
							writeLine = '      modulus: [' // fltStr(1:i8-2) // ']'
							write(1010,*) trim(writeLine)
						else
						    write(lfUnit,*) 'Warning: element ', i2, ' in section ', eSec, ' has no material definition.'
							write(lfUnit,*) 'Cannont write modulus to element property output file'
						endif
					endif
				elseif(propList(i3) .eq. 'shearModulus') then
				    if(eType .eq. 41 .or. eType .eq. 3) then
					    write(1010,*) '      shearModulus:'
						write(1010,*) '      #   layer, value'
						layerNum = 0
						do i4 = secLayupRange(eSec-1)+1, secLayupRange(eSec)
						    layerNum = layerNum + 1
							matID = layupMatId(i4)
							elProps(1:3) = materialElastic(7:9,matID)
							do i5 = elToDRange(i2-1) + 1, elToDRange(i2)
								i6 = elToD(i5)
								if(dCategory(i6) .eq. 'shearModulus' .and. dLayer(i6) .eq. layerNum) then
								    i7 = dComponent(i6)
									elProps(i7) = elProps(i7) + elToCoef(i5)*r_dVec(i6)
								endif
							enddo
							write(intStr,'(I0)') layerNum - 1
							i8 = 1
							do i5 = 1, 3
							    write(fltStr(i8:i8+17),'(G18.12)') elProps(i5)
								i8 = i8 + 18
								fltStr(i8:i8) = ','
								i8 = i8 + 1
							enddo
							writeLine = '          - [' // trim(intStr) // ', ' // fltStr(1:i8-2) // ']'
							write(1010,*) trim(writeLine)
						enddo
					else
					    matID = sectionMatId(eSec)
						if(matID .gt. 0) then
							elProps(1:3) = materialElastic(7:9,matID)
							do i4 = elToDRange(i2-1) + 1, elToDRange(i2)
								i5 = elToD(i4)
								if(dCategory(i5) .eq. 'shearModulus') then
								    i7 = dComponent(i5)
									elProps(i7) = elProps(i7) + elToCoef(i4)*r_dVec(i5)
								endif
							enddo
							i8 = 1
							do i5 = 1, 3
							    write(fltStr(i8:i8+17),'(G18.12)') elProps(i5)
								i8 = i8 + 18
								fltStr(i8:i8) = ','
								i8 = i8 + 1
							enddo
							writeLine = '      shearModulus: [' // fltStr(1:i8-2) // ']'
							write(1010,*) trim(writeLine)
						else
						    write(lfUnit,*) 'Warning: element ', i2, ' in section ', eSec, ' has no material definition.'
							write(lfUnit,*) 'Cannont write shearModulus to element property output file'
						endif
					endif
				elseif(propList(i3) .eq. 'poissonRatio') then
				    if(eType .eq. 41 .or. eType .eq. 3) then
					    write(1010,*) '      poissonRatio:'
						write(1010,*) '      #   layer, value'
						layerNum = 0
						do i4 = secLayupRange(eSec-1)+1, secLayupRange(eSec)
						    layerNum = layerNum + 1
							matID = layupMatId(i4)
							elProps(1:3) = materialElastic(4:6,matID)
							do i5 = elToDRange(i2-1) + 1, elToDRange(i2)
								i6 = elToD(i5)
								if(dCategory(i6) .eq. 'poissonRatio' .and. dLayer(i6) .eq. layerNum) then
								    i7 = dComponent(i6)
									elProps(i7) = elProps(i7) + elToCoef(i5)*r_dVec(i6)
								endif
							enddo
							write(intStr,'(I0)') layerNum - 1
							i8 = 1
							do i5 = 1, 3
							    write(fltStr(i8:i8+17),'(G18.12)') elProps(i5)
								i8 = i8 + 18
								fltStr(i8:i8) = ','
								i8 = i8 + 1
							enddo
							writeLine = '          - [' // trim(intStr) // ', ' // fltStr(1:i8-2) // ']'
							write(1010,*) trim(writeLine)
						enddo
					else
					    matID = sectionMatId(eSec)
						if(matID .gt. 0) then
							elProps(1:3) = materialElastic(4:6,matID)
							do i4 = elToDRange(i2-1) + 1, elToDRange(i2)
								i5 = elToD(i4)
								if(dCategory(i5) .eq. 'poissonRatio') then
								    i7 = dComponent(i5)
									elProps(i7) = elProps(i7) + elToCoef(i4)*r_dVec(i5)
								endif
							enddo
							i8 = 1
							do i5 = 1, 3
							    write(fltStr(i8:i8+17),'(G18.12)') elProps(i5)
								i8 = i8 + 18
								fltStr(i8:i8) = ','
								i8 = i8 + 1
							enddo
							writeLine = '      poissonRatio: [' // fltStr(1:i8-2) // ']'
							write(1010,*) trim(writeLine)
						else
						    write(lfUnit,*) 'Warning: element ', i2, ' in section ', eSec, ' has no material definition.'
							write(lfUnit,*) 'Cannont write poissonRatio to element property output file.'
						endif
					endif
				elseif(propList(i3) .eq. 'thermalCond') then
				    if(eType .eq. 41 .or. eType .eq. 3) then
					    write(1010,*) '      thermalCond:'
						write(1010,*) '      #   layer, value'
						layerNum = 0
						do i4 = secLayupRange(eSec-1)+1, secLayupRange(eSec)
						    layerNum = layerNum + 1
							matID = layupMatId(i4)
							elProps(1:6) = materialThermCond(1:6,matID)
							do i5 = elToDRange(i2-1) + 1, elToDRange(i2)
								i6 = elToD(i5)
								if(dCategory(i6) .eq. 'thermalCond' .and. dLayer(i6) .eq. layerNum) then
								    i7 = dComponent(i6)
									elProps(i7) = elProps(i7) + elToCoef(i5)*r_dVec(i6)
								endif
							enddo
							write(intStr,'(I0)') layerNum - 1
							i8 = 1
							do i5 = 1, 6
							    write(fltStr(i8:i8+17),'(G18.12)') elProps(i5)
								i8 = i8 + 18
								fltStr(i8:i8) = ','
								i8 = i8 + 1
							enddo
							writeLine = '          - [' // trim(intStr) // ', ' // fltStr(1:i8-2) // ']'
							write(1010,*) trim(writeLine)
						enddo
					elseif(eType .eq. 2 .and. abs(beamThermCond(1,eSec)) .gt. r_0) then
					    call r_getBeamTCond(tCondMat,i2)
						elProps(1) = tCondMat(1,1)
						elProps(2) = tCondMat(2,2)
						elProps(3) = tCondMat(3,3)
						elProps(4) = tCondMat(1,2)
						elProps(5) = tCondMat(1,3)
						elProps(6) = tCondMat(2,3)
						i8 = 1
						do i5 = 1, 6
							write(fltStr(i8:i8+17),'(G18.12)') elProps(i5)
							i8 = i8 + 18
							fltStr(i8:i8) = ','
							i8 = i8 + 1
						enddo
						writeLine = '      sectnthermalCond: [' // fltStr(1:i8-2) // ']'
						write(1010,*) trim(writeLine)
					else
					    matID = sectionMatId(eSec)
						if(matID .gt. 0) then
							elProps(1:6) = materialThermCond(1:6,matID)
							do i4 = elToDRange(i2-1) + 1, elToDRange(i2)
								i5 = elToD(i4)
								if(dCategory(i5) .eq. 'thermalCond') then
								    i7 = dComponent(i5)
									elProps(i7) = elProps(i7) + elToCoef(i4)*r_dVec(i5)
								endif
							enddo
							i8 = 1
							do i5 = 1, 6
							    write(fltStr(i8:i8+17),'(G18.12)') elProps(i5)
								i8 = i8 + 18
								fltStr(i8:i8) = ','
								i8 = i8 + 1
							enddo
							writeLine = '      thermalCond: [' // fltStr(1:i8-2) // ']'
							write(1010,*) trim(writeLine)
						else
						    write(lfUnit,*) 'Warning: element ', i2, ' in section ', eSec, ' has no material definition.'
							write(lfUnit,*) 'Cannont write thermalCond to element property output file.'
						endif
					endif
				elseif(propList(i3) .eq. 'thermalExp') then
				    if(eType .eq. 41 .or. eType .eq. 3) then
					    write(1010,*) '      thermalExp:'
						write(1010,*) '      #   layer, value'
						layerNum = 0
						do i4 = secLayupRange(eSec-1)+1, secLayupRange(eSec)
						    layerNum = layerNum + 1
							matID = layupMatId(i4)
							elProps(1:6) = materialThermExp(1:6,matID)
							do i5 = elToDRange(i2-1) + 1, elToDRange(i2)
								i6 = elToD(i5)
								if(dCategory(i6) .eq. 'thermalExp' .and. dLayer(i6) .eq. layerNum) then
								    i7 = dComponent(i6)
									elProps(i7) = elProps(i7) + elToCoef(i5)*r_dVec(i6)
								endif
							enddo
							write(intStr,'(I0)') layerNum - 1
							i8 = 1
							do i5 = 1, 6
							    write(fltStr(i8:i8+17),'(G18.12)') elProps(i5)
								i8 = i8 + 18
								fltStr(i8:i8) = ','
								i8 = i8 + 1
							enddo
							writeLine = '          - [' // trim(intStr) // ', ' // fltStr(1:i8-2) // ']'
							write(1010,*) trim(writeLine)
						enddo
					elseif(eType .eq. 2 .and. abs(beamExpLoadCoef(1,eSec)) .gt. r_0) then
					    call r_getBeamExpLoad(elProps(1:6),i2)
						i8 = 1
						do i5 = 1, 6
							write(fltStr(i8:i8+17),'(G18.12)') elProps(i5)
							i8 = i8 + 18
							fltStr(i8:i8) = ','
							i8 = i8 + 1
						enddo
						writeLine = '      sectnThermalExp: [' // fltStr(1:i8-2) // ']'
						write(1010,*) trim(writeLine)
					else
					    matID = sectionMatId(eSec)
						if(matID .gt. 0) then
							elProps(1:6) = materialThermExp(1:6,matID)
							do i4 = elToDRange(i2-1) + 1, elToDRange(i2)
								i5 = elToD(i4)
								if(dCategory(i5) .eq. 'thermalExp') then
								    i7 = dComponent(i5)
									elProps(i7) = elProps(i7) + elToCoef(i4)*r_dVec(i5)
								endif
							enddo
							i8 = 1
							do i5 = 1, 6
							    write(fltStr(i8:i8+17),'(G18.12)') elProps(i5)
								i8 = i8 + 18
								fltStr(i8:i8) = ','
								i8 = i8 + 1
							enddo
							writeLine = '      thermalExp: [' // fltStr(1:i8-2) // ']'
							write(1010,*) trim(writeLine)
						else
						    write(lfUnit,*) 'Warning: element ', i2, ' in section ', eSec, ' has no material definition.'
							write(lfUnit,*) 'Cannont write thermalExp to element property output file.'
						endif
					endif
				elseif(propList(i3) .eq. 'areaMoment') then
				    if(eType .eq. 2) then
					    call r_getBeamSecProps(elProps(1:7),i2)
						i8 = 1
						do i5 = 2, 6
							write(fltStr(i8:i8+17),'(G18.12)') elProps(i5)
							i8 = i8 + 18
							fltStr(i8:i8) = ','
							i8 = i8 + 1
						enddo
						writeLine = '      areaMoment: [' // fltStr(1:i8-2) // ']'
						write(1010,*) trim(writeLine)
					endif
				elseif(propList(i3) .eq. 'orientation') then
				    call r_getOrientation(orient,i2)
					if(eType .eq. 41 .or. eType .eq. 3 .or. eType .eq. 2) then
					    call r_getNodeCoord(nodes,i2)
					    call correctAlpha(orient,nodes,eType)
					endif
					i8 = 1
					do i5 = 1, 3
					    do i6 = 1, 3
						    write(fltStr(i8:i8+17),'(G18.12)') orient(i5,i6)
							i8 = i8 + 18
							fltStr(i8:i8) = ','
							i8 = i8 + 1
						enddo
					enddo
					writeLine = '      orientation: [' // fltStr(1:i8-2) // ']'
				    write(1010,*) trim(writeLine)
				elseif(propList(i3) .eq. 'stiffnessMat') then
				    if(eType .eq. 2) then
					    call r_getBeamStiffness(CMat,i2)
					elseif(eType .eq. 41 .or. eType .eq. 3) then
					    call r_getShellStiffness(ABDMat,i2)
						CMat(:,:) = ABDMat(1:6,1:6)
					else
					    call r_getMaterialStiffness(CMat,i2)
					endif
					i8 = 1
					do i5 = 1, 6
					    do i6 = i5, 6
						    write(fltStr(i8:i8+17),'(G18.12)') CMat(i5,i6)
							i8 = i8 + 18
							fltStr(i8:i8) = ','
							i8 = i8 + 1
						enddo
					enddo
					writeLine = '      stiffnessMat: [' // fltStr(1:i8-2) // ']'
				    write(1010,*) trim(writeLine)
				elseif(propList(i3) .eq. 'massMat') then
				    if(eType .eq. 2) then
					    call r_getBeamMass(CMat,i2)
					elseif(eType .eq. 41 .or. eType .eq. 3) then
					    call r_getShellMass(CMat,i2)
					endif
					i8 = 1
					do i5 = 1, 6
					    do i6 = i5, 6
						    write(fltStr(i8:i8+17),'(G18.12)') CMat(i5,i6)
							i8 = i8 + 18
							fltStr(i8:i8) = ','
							i8 = i8 + 1
						enddo
					enddo
					writeLine = '      massMat: [' // fltStr(1:i8-2) // ']'
				    write(1010,*) trim(writeLine)
				endif
			enddo
		enddo
		
		close(1010)
		
	end subroutine writeElProperties
	
	subroutine writeDVarValues(fileName,fields,numFields)
	    implicit none
		
		character(len=128), intent(in) :: fileName
		integer, intent(in) :: numFields
		character(len=16), intent(in) :: fields(numFields)
		
		character(len=128) :: writeLine
		integer :: i1, i2, i3
		
		open(unit=610, file=fileName, status='replace', action='write')
		
		write(610,*) 'designVariableValues:'
		
		do i1 = 1, numDVar
		    writeLine(1:13) = '    - value: '
			write(writeLine(14:31),'(G18.12)') r_dVec(i1)
			write(610,*) writeLine(1:31)
			do i2 = 1, numFields
			    if(fields(i2) .eq. 'category') then
					write(610,*) '      category: ', dCategory(i1)
				elseif(fields(i2) .eq. 'component') then
					write(610,*) '      component: ', dComponent(i1)
				elseif(fields(i2) .eq. 'layer') then
					write(610,*) '      layer: ', dLayer(i1)
				elseif(fields(i2) .eq. 'activeTime') then
				    writeLine(1:19) = '      activeTime: ['
					write(writeLine(20:37),'(G18.12)') dActTime(1,i1)
					writeLine(38:39) = ', '
					write(writeLine(40:57),'(G18.12)') dActTime(2,i1)
					writeLine(58:58) = ']'
					write(610,*) writeLine(1:58)
				endif
			enddo
		enddo
		
		close(610)
		
	end subroutine writeDVarValues
	
	subroutine writeObjective(fileName,fields,numFields,writeGrad)
	    implicit none
		
		character(len=128), intent(in) :: fileName
		integer, intent(in) :: numFields
		character(len=16), intent(in) :: fields(numFields)
		integer, intent(in) :: writeGrad
		
		character(len=128) :: writeLine
		character(len=32) :: intStr, fltStr
		
		real*8 :: totalObj
		integer :: i1, i2, i3
		
		open(unit=653, file=fileName, status='replace', action='write')
		
		write(653,*) 'objective:'
		
		totalObj = r_0
		do i1 = 1, numObjTerms
		    totalObj = totalObj + objVal(i1)
		enddo
		
        writeLine = '    totalValue: '
        write(writeLine(19:36),'(G18.12)') totalObj
		write(653,*) trim(writeLine)

        write(653,*) '    termValues:'
        do i1 = 1, numObjTerms
		    writeLine = '        - value: '
            write(writeLine(18:35),'(G18.12)') objVal(i1)
		    write(653,*) writeLine(1:35)
			do i2 = 1, numFields
			    if(fields(i2) .eq. 'category') then
				    write(653,*) '          category: ', objCategory(i1)
				elseif(fields(i2) .eq. 'operator') then
				    write(653,*) '          operator: ', objOperator(i1)
				elseif(fields(i2) .eq. 'component') then
				    write(653,*) '          component: ', objComponent(i1)
				elseif(fields(i2) .eq. 'layer') then
				    write(653,*) '          layer: ', objLayer(i1)
				elseif(fields(i2) .eq. 'elementSet') then
				    write(653,*) '          elementSet: ', objElSet(i1)
				elseif(fields(i2) .eq. 'nodeSet') then
				    write(653,*) '          nodeSet: ', objNdSet(i1)
				elseif(fields(i2) .eq. 'coefficient') then
				    writeLine = '          coefficient: '
					write(writeLine(23:40),'(G18.12)') objCoef(i1)
				    write(653,*) writeLine(1:40)
				elseif(fields(i2) .eq. 'exponent') then
				    writeLine = '          exponent: '
					write(writeLine(21:38),'(G18.12)') objExp(i1)
				    write(653,*) writeLine(1:38)
				elseif(fields(i2) .eq. 'activeTime') then
				    writeLine(1:19) = '      activeTime: ['
					write(writeLine(20:37),'(G18.12)') dActTime(1,i1)
					writeLine(38:39) = ', '
					write(writeLine(40:57),'(G18.12)') dActTime(2,i1)
					writeLine(58:58) = ']'
					write(653,*) writeLine(1:58)
				endif
			enddo
		enddo
		
		if(writeGrad .eq. 1) then
			write(653,*) 'objectiveGradient:'
			
			do i1 = 1, numDVar
				write(intStr,'(I0)') i1 - 1
				write(fltStr,'(G18.12)') dLdD(i1)
				writeLine = '    - [' // trim(intStr) // ', ' // trim(fltStr) // ']'
				write(653,*) trim(writeLine)
			enddo
		endif
		
		close(653)
		
	end subroutine writeObjective
	
	subroutine writeSparseMatrix(fileName,aMat,aCols,aRange,aSize,aDim)
	    implicit none
		
		integer, intent(in) :: aSize, aDim
		character(len=128), intent(in) :: fileName
		real*8, intent(in) :: aMat(aSize)
		integer, intent(in) :: aRange(0:aDim)
		integer, intent(in) :: aCols(aSize)
		
		integer :: i1, i2, i3
		
		open(unit=1400, file=fileName, status='replace', action='write')
		
		write(1400,*) 'matrix:'
		write(1400,*) '    dimension: ', aDim
		write(1400,*) '    nonZeroSize: ', aSize
		write(1400,*) '    elements:'
		write(1400,*) '        # row, column, value'
		do i1 = 1, aDim
		    do i2 = aRange(i1-1)+1, aRange(i1)
			    write(1400,*) '        - [', i1, ', ', aCols(i2), ', ', aMat(i2), ']' 
			enddo
		enddo
		
		close(1400)
	
	end subroutine writeSparseMatrix
	
	subroutine writeLowerTriMatrix(fileName,aMat,aRange,aSize,aDim)
	    implicit none
		
		integer, intent(in) :: aSize, aDim
		character(len=128), intent(in) :: fileName
		real*8, intent(in) :: aMat(aSize)
        integer, intent(in) :: aRange(0:aDim)
		
		integer :: i1, i2, i3
		
		open(unit=1420, file=fileName, status='replace', action='write')
		
		write(1420,*) 'matrix:'
		write(1420,*) '    dimension: ', aDim
		write(1420,*) '    size: ', aSize
		write(1420,*) '    elements:'
		write(1420,*) '        # row, column, value'
		do i1 = 1, aDim
		    i3 = i1 - aRange(i1) + aRange(i1-1) + 1
		    do i2 = aRange(i1-1)+1, aRange(i1)
			    write(1420,*) '        - [', i1, ', ', i3, ', ', aMat(i2), ']'
                i3 = i3 + 1				
			enddo
		enddo
		
		
		close(1420)
		
	end subroutine writeLowerTriMatrix
	
!! Binary
	subroutine writeBinarySolution(stepNum)
	    implicit none
		
		integer, intent(in) :: stepNum
		
		character(len=16) :: stepString
		character(len=32) :: fullFileName
		
		write(stepString,'(I0)') stepNum
        fullFileName = 'SolutionHistory' // stepString // '.out'
	    open(unit=17, file=fullFileName, form='unformatted', status='replace', action='write')
		
		write(17) prevTemp(:)
		write(17) prevTdot(:)
		write(17) prevDisp(:)
		if(intVecSize .ge. 1) then
		    write(17) prevIntDisp(:)
		endif
		write(17) prevVel(:)
		write(17) prevAcc(:)
		
		close(17) 
		
	end subroutine writeBinarySolution
	
!! Debugging
	
	subroutine writeModelData(fileName)
	    implicit none
		
		character(len=128), intent(in) :: fileName
		integer :: i1
		
		open(unit=37, file=fileName, status='replace', action='write')
		
		write(37,*) 'nodeList'
		do i1 = 1, numNodes
		    write(37,*) nodeList(:,i1)
		enddo
		write(37,*) 'elementList'
		do i1 = 1, numEls
		    write(37,*) elementList(:,i1)
		enddo
		write(37,*) 'elementSection'
		write(37,*) elementSection(:)
		write(37,*) 'elementType'
		write(37,*) elementType(:)
		write(37,*) 'elementSurfaces'
		write(37,*) elementSurfaces(:,:)
		write(37,*) 'elementSets'
		write(37,*) elementSets(:)
		write(37,*) 'elSetName'
		write(37,*) elSetName(:)
		write(37,*) 'elSetRange'
		write(37,*) elSetRange(:)
		write(37,*) 'nodeSets'
		write(37,*) nodeSets(:)
		write(37,*) 'ndSetName'
		write(37,*) ndSetName(:)
		write(37,*) 'ndSetRange'
		write(37,*) ndSetRange(:)
		write(37,*) 'sectionType'
		write(37,*) sectionType(:)
		write(37,*) 'sectionMatId'
		write(37,*) sectionMatId(:)
		write(37,*) 'sectionMatName'
		write(37,*) sectionMatName(:)
		write(37,*) 'sectionOrient'
		write(37,*) sectionOrient(:,:)
		write(37,*) 'layupMatName'
		write(37,*) layupMatName(:)
		write(37,*) 'layupMatId'
		write(37,*) layupMatId(:)
		write(37,*) 'layupThickness'
		write(37,*) layupThickness(:)
		write(37,*) 'layupAngle'
		write(37,*) layupAngle(:)
		write(37,*) 'secLayupRange'
		write(37,*) secLayupRange(:)
		write(37,*) 'sectionZOffset'
		write(37,*) sectionZOffset(:)
		write(37,*) 'beamProperties'
		write(37,*) beamProperties(:,:)
		write(37,*) 'beamStiffness'
		write(37,*) beamStiffness(:,:)
		write(37,*) 'beamExpLoadCoef'
		write(37,*) beamExpLoadCoef(:,:)
		write(37,*) 'beamMass'
		write(37,*) beamMass(:,:)
		write(37,*) 'beamThermCond'
		write(37,*) beamThermCond(:,:)
		write(37,*) 'beamSpecHeat'
		write(37,*) beamSpecHeat(:)
		write(37,*) 'materialName'
		write(37,*) materialName(:)
		write(37,*) 'materialDensity'
		write(37,*) materialDensity(:)
		write(37,*) 'materialElastic'
		write(37,*) materialElastic(:,:)
		write(37,*) 'materialStiffMat'
		write(37,*) materialStiffMat(:,:)
		write(37,*) 'materialThermCond'
		write(37,*) materialThermCond(:,:)
		write(37,*) 'materialThermExp'
		write(37,*) materialThermExp(:,:)
		write(37,*) 'materialSpecHeat'
		write(37,*) materialSpecHeat(:)
		write(37,*) 'materialMaxStress'
		write(37,*) materialMaxStress(:,:)
		write(37,*) 'materialMaxStrain'
		write(37,*) materialMaxStrain(:,:)
		write(37,*) 'materialMaxStrnEngy'
		write(37,*) materialMaxStrnEngy(:)
		
		if(numMPC .gt. 0) then
		    write(37,*) 'mpcEqn'
		    write(37,*) mpcEqn(:)
			write(37,*) 'mpcNode'
			write(37,*) mpcNode(:)
			write(37,*) 'mpcDof'
			write(37,*) mpcDof(:)
			write(37,*) 'mpcCoef'
			write(37,*) mpcCoef(:)
			write(37,*) 'mpcRHS'
			write(37,*) mpcRHS(:)
		endif
		
		if(numLds .gt. 0) then
		    write(37,*) 'loadNodes'
		    write(37,*) loadNodes(:)
			write(37,*) 'inputLoads'
	        write(37,*) inputLoads(:,:)
			write(37,*) 'loadsRange'
	        write(37,*) loadsRange(:)
			write(37,*) 'loadsActTime'
	        write(37,*) loadsActTime(:,:)
			write(37,*) 'loadType'
	        write(37,*) loadType(:)
		endif
		
		write(37,*) 'initialDisp'
		write(37,*) initialDisp(:,:)
		write(37,*) 'initialVel'
	    write(37,*) initialVel(:,:)
		write(37,*) 'initalAcc'
	    write(37,*) initialAcc(:,:)
		write(37,*) 'initialTemp'
	    write(37,*) initialTemp(:)
		write(37,*) 'initialTdot'
	    write(37,*) initialTdot(:)
		
		close(37)
	end subroutine writeModelData
	
	subroutine writeDesignVarData(fileName)
	    implicit none
		
		character(len=128), intent(in) :: fileName
		
		open(unit=122, file=fileName, status='replace', action='write')
		
		write(122,*) 'r_dVec'
		write(122,*) r_dVec(:)
		write(122,*) 'dCategory'
		write(122,*) dCategory(:)
		write(122,*) 'dComponent'
		write(122,*) dComponent(:)
		write(122,*) 'dLayer'
		write(122,*) dLayer(:)
		write(122,*) 'dActTime'
		write(122,*) dActTime(:,:)
		write(122,*) 'dElSet'
		write(122,*) dElSet(:)
		write(122,*) 'dNdSet'
		write(122,*) dNdSet(:)
		write(122,*) 'dCoefList'
		write(122,*) dCoefList(:)
		write(122,*) 'dCoefListRange'
		write(122,*) dCoefListRange(:)
		write(122,*) 'elToD'
		write(122,*) elToD(:)
		write(122,*) 'ndToD'
		write(122,*) ndToD(:)
		write(122,*) 'dToNd'
		write(122,*) dToNd(:)
		write(122,*) 'elToDComp'
		write(122,*) elToDComp(:)
		write(122,*) 'dToElComp'
		write(122,*) dToElComp(:)
		write(122,*) 'elToDRange'
		write(122,*) elToDRange(:)
		write(122,*) 'ndToDRange'
		write(122,*) ndToDRange(:)
		write(122,*) 'dToNdRange'
		write(122,*) dToNdRange(:)
		write(122,*) 'elToDCRange'
		write(122,*) elToDCRange(:)
		write(122,*) 'dToElCRange'
		write(122,*) dToElCRange(:)
		write(122,*) 'elToCoef'
		write(122,*) elToCoef(:)
		write(122,*) 'ndToCoef'
		write(122,*) ndToCoef(:)
		write(122,*) 'dToNdCoef'
		write(122,*) dToNdCoef(:)		
		
		close(122)
	end subroutine writeDesignVarData
	
	subroutine writeObjectiveData(fileName)
	    implicit none
		
		character(len=128), intent(in) :: fileName
		
		open(unit=220, file=fileName, status='replace', action='write')
		
		write(220,*) 'objCategory'
		write(220,*) objCategory(:)
		write(220,*) 'objOperator'
	    write(220,*) objOperator(:)
		write(220,*) 'objActTime'
		write(220,*) objActTime(:,:)
		write(220,*) 'objComponent'
		write(220,*) objComponent(:)
		write(220,*) 'objLayer'
		write(220,*) objLayer(:)
		write(220,*) 'objCoef'
		write(220,*) objCoef(:)
		write(220,*) 'objExp'
		write(220,*) objExp(:)
		write(220,*) 'objElSet'
		write(220,*) objElSet(:)
		write(220,*) 'objNdSet'
		write(220,*) objNdSet(:)
		write(220,*) 'objTgtTag'
		write(220,*) objTgtTag(:)
		write(220,*) 'objTgtVal'
		write(220,*) objTgtVal(:)
		write(220,*) 'objTgtRange'
		write(220,*) objTgtRange(:)
		
		close(220)
		
	end subroutine writeObjectiveData
	
	
end module ASenD_output