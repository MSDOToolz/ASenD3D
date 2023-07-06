module ASenD_globalData

!! Model input data

    real*8, allocatable :: nodeList(:,:)
	integer, allocatable :: elementList(:,:)
	integer, allocatable :: elementSection(:)
	integer, allocatable :: elementType(:)
	integer, allocatable :: elementSurfaces(:,:)
	integer, allocatable :: elementSets(:)
	character(len=64), allocatable :: elSetName(:)
	integer, allocatable :: elSetRange(:)
	integer, allocatable :: nodeSets(:)
	character(len=64), allocatable :: ndSetName(:)
	integer, allocatable :: ndSetRange(:)
	character(len=64), allocatable :: sectionType(:)
	integer, allocatable :: sectionMatId(:)
	character(len=64), allocatable :: sectionMatName(:)
	real*8, allocatable :: sectionOrient(:,:)
	character(len=64), allocatable :: layupMatName(:)
	integer, allocatable :: layupMatId(:)
	real*8, allocatable :: layupThickness(:)
	real*8, allocatable :: layupAngle(:)
	integer, allocatable :: secLayupRange(:)
	real*8, allocatable :: sectionZOffset(:)
	real*8, allocatable :: beamProperties(:,:)
	real*8, allocatable :: beamStiffness(:,:)
	real*8, allocatable :: beamExpLoadCoef(:,:)
	real*8, allocatable :: beamMass(:,:)
	real*8, allocatable :: beamThermCond(:,:)
	real*8, allocatable :: beamSpecHeat(:)
	character(len=64), allocatable :: materialName(:)
	real*8, allocatable :: materialDensity(:)
	real*8, allocatable :: materialElastic(:,:)
	real*8, allocatable :: materialStiffMat(:,:)
	real*8, allocatable :: materialThermCond(:,:)
	real*8, allocatable :: materialThermExp(:,:)
	real*8, allocatable :: materialSpecHeat(:)
	real*8, allocatable :: materialMaxStress(:,:)
	real*8, allocatable :: materialMaxStrain(:,:)
	real*8, allocatable :: materialMaxStrnEngy(:)
	integer, allocatable :: mpcEqn(:)
	character(len=64), allocatable :: mpcNode(:)
	integer, allocatable :: mpcDof(:)
	real*8, allocatable :: mpcCoef(:)
	real*8, allocatable :: mpcRHS(:)
	character(len=64), allocatable :: loadNodes(:)
	real*8, allocatable :: inputLoads(:,:)
	integer, allocatable :: loadsRange(:)
	real*8, allocatable :: loadsActTime(:,:)
	character(len=16), allocatable :: loadType(:)
	real*8, allocatable :: initialDisp(:,:)
	real*8, allocatable :: initialVel(:,:)
	real*8, allocatable :: initialAcc(:,:)
	real*8, allocatable :: initialTemp(:)
	real*8, allocatable :: initialTdot(:)
	integer :: numNodes, numEls, numNdSets, ndSetSize, numElSets, elSetSize
	integer :: numSec, layupSize, numMats, numMPC, mpcSize, numLds, sizeLds
	
!! Design varible input data
	
	real*8, allocatable :: r_dVec(:)
	complex*16, allocatable :: c_dVec(:)
	character(len=16), allocatable :: dCategory(:)
	integer, allocatable :: dComponent(:)
	integer, allocatable :: dLayer(:)
	real*8, allocatable :: dActTime(:,:)
	character(len=64), allocatable :: dElSet(:)
	character(len=64), allocatable :: dNdSet(:)
	real*8, allocatable :: dCoefList(:)
	integer, allocatable :: dCoefListRange(:)
	integer :: numDVar, coefListSize
	
	integer, allocatable :: elToD(:)
	integer, allocatable :: dToEl(:)
	integer, allocatable :: ndToD(:)
	integer, allocatable :: dToNd(:)
	integer, allocatable :: elToDComp(:)
	integer, allocatable :: dToElComp(:)
	integer, allocatable :: elToDRange(:)
	integer, allocatable :: dToElRange(:)
	integer, allocatable :: ndToDRange(:)
	integer, allocatable :: dToNdRange(:)
	integer, allocatable :: elToDCRange(:)
	integer, allocatable :: dToElCRange(:)
	real*8, allocatable :: elToCoef(:)
	real*8, allocatable :: dToElCoef(:)
	real*8, allocatable :: ndToCoef(:)
	real*8, allocatable :: dToNdCoef(:)
	integer :: elToDSize, ndToDSize, dToElCSize

!! Objective input data

    character(len=16), allocatable :: objCategory(:)
	character(len=16), allocatable :: objOperator(:)
	real*8, allocatable :: objActTime(:,:)
	integer, allocatable :: objComponent(:)
	integer, allocatable :: objLayer(:)
	real*8, allocatable :: objCoef(:)
	real*8, allocatable :: objExp(:)
	character(len=64), allocatable :: objElSet(:)
	integer, allocatable :: objElSetId(:)
	character(len=64), allocatable :: objNdSet(:)
	integer, allocatable :: objNdSetId(:)
	character(len=32), allocatable :: objTgtTag(:)
	real*8, allocatable :: objTgtVal(:)
	real*8, allocatable :: objTgtExpanded(:)
	integer, allocatable :: objTgtRange(:)
    integer :: numObjTerms, objTgtSize

!! Analysis utilities

    real*8, allocatable :: nodeTemp(:)
	real*8, allocatable :: nodeTdot(:)
    real*8, allocatable :: nodeDisp(:)
	real*8, allocatable :: nodeVel(:)
	real*8, allocatable :: nodeAcc(:)
	real*8, allocatable :: prevTemp(:)
	real*8, allocatable :: prevTdot(:)
    real*8, allocatable :: prevDisp(:)
	real*8, allocatable :: prevVel(:)
	real*8, allocatable :: prevAcc(:)
	real*8, allocatable :: internalDisp(:)
	real*8, allocatable :: prevIntDisp(:)
	real*8, allocatable :: elasticLoad(:)
	real*8, allocatable :: intElasticLoad(:)
	real*8, allocatable :: thermalLoad(:)
	
	real*8, allocatable :: eigenModes(:,:)
	real*8, allocatable :: eigenVals(:)
	real*8, allocatable :: eigenFactors(:)
	real*8, allocatable :: diagMassMat(:)
	
	!! --------------------------
	
	real*8, allocatable :: aFull(:,:)
	
	!! --------------------------
	
	real*8, allocatable :: elasticMat(:)
	integer, allocatable :: elMatCols(:)
	integer, allocatable :: elMatRange(:)
	real*8, allocatable :: elMatLT(:)
	integer, allocatable :: elMatLTRange(:)
	real*8, allocatable :: elMPCMat(:)
	integer, allocatable :: elMPCMatCols(:)
	integer, allocatable :: elMPCMatRange(:)
	real*8, allocatable :: elMPCRHS(:)
	real*8, allocatable :: intElasticMat(:)
	real*8, allocatable :: thermMat(:)
	integer, allocatable :: thermMatCols(:)
	integer, allocatable :: thermMatRange(:)
	real*8, allocatable :: thermMatLT(:)
	integer, allocatable :: thermMatLTRange(:)
	real*8, allocatable :: thermMPCMat(:)
	integer, allocatable :: thermMPCMatCols(:)
	integer, allocatable :: thermMPCMatRange(:)
	real*8, allocatable :: thermMPCRHS(:)
	
	real*8, allocatable :: delDisp(:)
	real*8, allocatable :: swapVec(:)
	real*8, allocatable :: intswapVec(:)
	
    integer, allocatable :: currentRank(:)
	integer, allocatable :: originalRank(:)
	integer, allocatable :: nDofIndex(:,:)
	integer, allocatable :: intMatRange(:)
	integer, allocatable :: intVecRange(:)
	integer, allocatable :: elDepends(:)
	
	real*8, allocatable :: objVal(:)
	real*8, allocatable :: dLdD(:)
	real*8, allocatable :: dLdu(:)
	real*8, allocatable :: dLdv(:)
	real*8, allocatable :: dLda(:)
	real*8, allocatable :: dLdt(:)
	real*8, allocatable :: dLdtdot(:)
	real*8, allocatable :: intdLdu(:)
	real*8, allocatable :: dsumTermdT(:)
	real*8, allocatable :: dsumTermdU(:)
	real*8, allocatable :: intdsumTermdU(:)
	real*8, allocatable :: dsumTermdD(:)
	real*8, allocatable :: dtotVoldD(:)
	real*8, allocatable :: tempAdj(:)
	real*8, allocatable :: tdotAdj(:)
	real*8, allocatable :: dispAdj(:)
	real*8, allocatable :: intDispAdj(:)
	real*8, allocatable :: velAdj(:)
	real*8, allocatable :: accAdj(:)
	real*8, allocatable :: dRudD(:)
	real*8, allocatable :: intdRudD(:)
	real*8, allocatable :: dRtdD(:)
	
	integer :: elMatDim, elMatSize, elMatLTSize, elMPCDim, elMPCSize
    integer :: thermMatDim, thermMatSize, thermMatLTSize, thermMPCDim, thermMPCSize
	integer :: intVecSize, intMatSize

!! Surface utilities

!! Job command input parameters

    character(len=16) :: solverMeth
    integer :: solverBlockDim, solverMaxBW, solveThermal, solveElastic
	integer :: nLGeom, dynamic, writeSolnHist, numTSteps
	real*8 :: nMBeta, nMGamma, delT, simPeriod, rayCoefK, rayCoefM
	real*8 :: loadTime
	integer :: ldRampSteps
	integer :: numEigModes, modalType
	
!! Job utilities

    integer :: lfUnit, abortJob
	

end module ASenD_globalData