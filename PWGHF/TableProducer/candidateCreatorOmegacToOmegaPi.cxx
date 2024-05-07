// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
// 
/// \file 
/// \brief Reconstruction of Omegac0  -> omega pi candidates
/// \author 


#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DCAFitter/DCAFitterN.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/V0.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod;
using namespace o2::aod::cascdata;
using namespace o2::aod::v0data;
using namespace o2::aod::hf_track_index;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;
int toCounter = 0;
// Reconstruction of omegac0 candidates
struct HfCandidateCreatorOmegacToOmegaPi {
  Produces<aod::HfCandOmegaC> rowCandidate;
 
  Configurable<bool> doPvRefit{"doPvRefit", false, "set to true if you do PV refit in trackIndexSkimCreator.cxx"};

  Configurable<bool> propagateToPCA{"propagateToPCA", false, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", true, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> maxDXYIni{"maxDXYIni", 4., "reject (if>0) PCA candidate if tracks DXY exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  Configurable<double> maxChi2{"maxChi2", 100., "discard vertices with chi2/Nprongs > this (or sum{DCAi^2}/Nprongs for abs. distance minimization)"};
  Configurable<bool> refitWithMatCorr{"refitWithMatCorr", true, "when doing propagateTracksToVertex, propagate tracks to vtx with material corrections and rerun minimization"};
  Configurable<bool> rejDiffCollTrack{"rejDiffCollTrack", true, "Reject tracks coming from different collisions"};

  // magnetic field setting from CCDB
  Configurable<bool> isRun2{"isRun2", false, "enable Run 2 or Run 3 GRP objects for magnetic field"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  // cascade cuts
  Configurable<bool> doCascadePreselection{"doCascadePreselection", true, "Use invariant mass and dcaXY cuts to preselect cascade candidates"};
  Configurable<double> massToleranceCascade{"massToleranceCascade", 0.01, "Invariant mass tolerance for cascade"};
  Configurable<float> dcaXYToPVCascadeMax{"dcaXYToPVCascadeMax", 3, "Max cascade DCA to PV in xy plane"};
 //hist
  HistogramRegistry registry{
    "registry",
    {
     
	}
	  };

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  int runNumber;

  using SelectedCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::HfSelCollision>>;
  using MyTracks = soa::Join<aod::TracksWCovDca, aod::HfPvRefitTrack>;
  using FilteredHfTrackAssocSel = soa::Filtered<soa::Join<aod::TrackAssoc, aod::HfSelTrack>>;
  using MyCascTable = soa::Join<aod::CascDatas, aod::CascCovs>; // to use strangeness tracking, use aod::TraCascDatas instead of aod::CascDatas
  using MyV0Table = soa::Join<aod::V0Datas, aod::V0Covs>;
  using MySkimIdx = HfCascLf2Prongs;

  Filter filterSelectCollisions = (aod::hf_sel_collision::whyRejectColl == static_cast<uint16_t>(0)); // filter to use only HF selected collisions
  Filter filterSelectTrackIds = (aod::hf_sel_track::isSelProng > static_cast<uint32_t>(0));

  Preslice<FilteredHfTrackAssocSel> trackIndicesPerCollision = aod::track_association::collisionId; // aod::hf_track_association::collisionId
  Preslice<MyCascTable> cascadesPerCollision = aod::cascdata::collisionId;
  Preslice<MySkimIdx> candidatesPerCollision = hf_track_index::collisionId;
 
  OutputObj<TH1F> hCollision{TH1F("hCollision", "Counter Collison;entries", 500, -2, 2)};
  OutputObj<TH1F> hPionNumber{TH1F("hPionNumber", "Pion Number;entries", 500, 0, 5)};
  OutputObj<TH1F> hPionEta{TH1F("hPionEta", "Pion eta;entries", 500, -2, 2)};
  OutputObj<TH1F> hInvMassCharmBaryon{TH1F("hInvMassOmegac", "Omegac candidate invariant mass;inv mass;entries", 500, 2.2, 3.1)};
  OutputObj<TH1F> hFitterStatus{TH1F("hFitterStatus", "Charm DCAFitter status;status;entries", 3, -0.5, 2.5)};                     // 0 --> vertex(es) found, 1 --> exception found, 2 --> no vertex found (but no exception)
  OutputObj<TH1F> hCandidateCounter{TH1F("hCandidateCounter", "Candidate counter wrt derived data;status;entries", 4, -0.5, 3.5)}; // 0 --> candidates in derived data table, 1 --> candidates passing testbit selection, 2 --> candidates passing fitter step 3 --> candidates filled in new table

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    runNumber = 0;
  };

  void processIdxCombinatorics(SelectedCollisions const& collisions,
                               aod::BCsWithTimestamps const& bcWithTimeStamps,
                               MyTracks const& tracks,
                               FilteredHfTrackAssocSel const& trackIndices,
                               MyCascTable const& cascades,
                               MyV0Table const&,
                               aod::V0sLinked const&){
 

    double massPionFromPDG = MassPiPlus;    // pdg code 211
    double massLambdaFromPDG = MassLambda0; // pdg code 3122
    double massOmegaFromPDG = MassOmegaMinus;
    double massOmegacFromPDG = MassOmegaC0; // pdg code 4332

    // 2-prong vertex fitter to build the omegac vertex
    o2::vertexing::DCAFitterN<2> df;
    df.setPropagateToPCA(propagateToPCA);
    df.setMaxR(maxR);
    df.setMaxDZIni(maxDZIni);
    df.setMaxDXYIni(maxDXYIni);
    df.setMinParamChange(minParamChange);
    df.setMinRelChi2Change(minRelChi2Change);
    df.setMaxChi2(maxChi2);
    df.setUseAbsDCA(useAbsDCA);
    df.setWeightedFinalPCA(useWeightedFinalPCA);

    for (const auto& collision : collisions) {

      // set the magnetic field from CCDB
      auto bc = collision.bc_as<o2::aod::BCsWithTimestamps>();
      initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
      auto magneticField = o2::base::Propagator::Instance()->getNominalBz(); // z component

      df.setBz(magneticField);
      df.setRefitWithMatCorr(refitWithMatCorr);

      // loop over cascades reconstructed by cascadebuilder.cxx
      auto thisCollId = collision.globalIndex();
      auto groupedCascades = cascades.sliceBy(cascadesPerCollision, thisCollId);
      hCollision->Fill(1);
      for (const auto& casc : groupedCascades) {
        
        //----------------accessing particles in the decay chain-------------
        // cascade daughter - charged particle
        auto trackOmegaDauCharged = casc.bachelor_as<MyTracks>(); // kaon <- omega track from MyTracks table //
        // cascade daughter - V0
		int v0index = casc.cascadeId();
        // V0 positive daughter
        auto trackV0Dau0 = casc.posTrack_as<MyTracks>(); // p <- V0 track (positive track) from MyTracks table
        // V0 negative daughter
        auto trackV0Dau1 = casc.negTrack_as<MyTracks>(); // pion <- V0 track (negative track) from MyTracks table
		
        // check that particles come from the same collision
        if (rejDiffCollTrack) {
          if (trackV0Dau0.collisionId() != trackV0Dau1.collisionId()) {
            continue;
          }
          if (trackOmegaDauCharged.collisionId() != trackV0Dau0.collisionId()) {
            continue;
          }
        }

        //-------------------------- V0 info---------------------------
        // pseudorapidity
        double pseudorapV0Dau0 = trackV0Dau0.eta();
        double pseudorapV0Dau1 = trackV0Dau1.eta();

        // pion & p <- V0 tracks
        auto trackParCovV0Dau0 = getTrackParCov(trackV0Dau0);
        auto trackParCovV0Dau1 = getTrackParCov(trackV0Dau1);

        // info from LF table
        std::array<float, 3> pVecV0 = {casc.pxlambda(), casc.pylambda(), casc.pzlambda()}; // pVec stands for vector containing the 3-momentum components
        std::array<float, 3> vertexV0 = {casc.xlambda(), casc.ylambda(), casc.zlambda()};
        std::array<float, 3> pVecV0Dau0 = {casc.pxpos(), casc.pypos(), casc.pzpos()};
        std::array<float, 3> pVecV0Dau1 = {casc.pxneg(), casc.pyneg(), casc.pzneg()};

        //-----------------------------reconstruct cascade track-----------------------------
        // pseudorapidity
        double pseudorapKaFromCas = trackOmegaDauCharged.eta();

        // kaon <- casc track
        auto trackParCovOmegaDauCharged = getTrackParCov(trackOmegaDauCharged);

        // info from LF table
        std::array<float, 3> vertexCasc = {casc.x(), casc.y(), casc.z()};
        std::array<float, 3> pVecCasc = {casc.px(), casc.py(), casc.pz()};
        std::array<float, 21> covCasc = {0.};
        constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
        for (int i = 0; i < 6; i++) {
          covCasc[MomInd[i]] = casc.momentumCovMat()[i];
          covCasc[i] = casc.positionCovMat()[i];
        }
        // create cascade track
        o2::track::TrackParCov trackCasc;
        if (trackOmegaDauCharged.sign() > 0) {
          trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, 1, true);
        } else if (trackOmegaDauCharged.sign() < 0) {
          trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, -1, true);
        } else {
          continue;
        }
        trackCasc.setAbsCharge(1);
        trackCasc.setPID(o2::track::PID::OmegaMinus);

        std::array<float, 3> pVecPionFromCasc = {casc.pxbach(), casc.pybach(), casc.pzbach()};

        //-------------------combining cascade and pion tracks--------------------------
        auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
        for (const auto& trackIndexPion : groupedTrackIndices) {

          // use bachelor selections from HfTrackIndexSkimCreatorTagSelTracks --> bit =2 is CandidateType::CandV0bachelor
          if (!TESTBIT(trackIndexPion.isSelProng(), 2)) {
            continue;
          }

          auto trackPion = trackIndexPion.track_as<MyTracks>();

          if ((rejDiffCollTrack) && (trackOmegaDauCharged.collisionId() != trackPion.collisionId())) {
            continue;
          }

          // ask for opposite sign daughters (charm baryon daughters)
          if (trackPion.sign() * trackOmegaDauCharged.sign() >= 0) {
            continue;
          };
		  
          // pseudorapidity
          double pseudorapPiFromCharmBaryon = trackPion.eta();
	  
          // charm bachelor pion track to be processed with DCAFitter
          auto trackParVarPi = getTrackParCov(trackPion);
	  // statistic pion number 
	  if (toCounter==0){
	     hPionNumber->Fill(4);
	  };
          // reconstruct charm baryon with DCAFitter
	  
          int nVtxFromFitterCharmBaryon = 0;
          try {
            nVtxFromFitterCharmBaryon = df.process(trackCasc, trackParVarPi);
          } catch (...) {
            LOG(error) << "Exception caught in charm DCA fitter process call!";
            hFitterStatus->Fill(1);
	    if (toCounter!=0){
             hPionNumber->Fill(3);
         	 };
            continue;

          }
          if (nVtxFromFitterCharmBaryon == 0) {
            hFitterStatus->Fill(2);
            if (toCounter!=0){
             hPionNumber->Fill(2);
            };
            continue;
          }
          hFitterStatus->Fill(0);
          if (toCounter==0){
             hPionNumber->Fill(1);
          };
          auto vertexCharmBaryonFromFitter = df.getPCACandidate();
          //auto chi2PCACharmBaryon = df.getChi2AtPCACandidate();
          std::array<float, 3> pVecCascAsD;
          std::array<float, 3> pVecPionFromCharmBaryon;
          df.propagateTracksToVertex();
          if (!df.isPropagateTracksToVertexDone()) {
            continue;
          }
          if (toCounter==0){
             hPionNumber->Fill(0);
          };
	      hPionEta->Fill(trackPion.eta());
          df.getTrack(0).getPxPyPzGlo(pVecCascAsD);
          df.getTrack(1).getPxPyPzGlo(pVecPionFromCharmBaryon);
          std::array<float, 3> pVecCharmBaryon = {pVecCascAsD[0] + pVecPionFromCharmBaryon[0], pVecCascAsD[1] + pVecPionFromCharmBaryon[1], pVecCascAsD[2] + pVecPionFromCharmBaryon[2]};

          std::array<float, 3> coordVtxCharmBaryon = df.getPCACandidatePos();
          std::array<float, 6> covVtxCharmBaryon = df.calcPCACovMatrixFlat();

          float dcaxyV0Dau0 = trackV0Dau0.dcaXY();
          float dcaxyV0Dau1 = trackV0Dau1.dcaXY();
          float dcaxyKaFromCasc = trackOmegaDauCharged.dcaXY();

          // DCAz (computed with propagateToDCABxByBz method)
          float dcazV0Dau0 = trackV0Dau0.dcaZ();
          float dcazV0Dau1 = trackV0Dau1.dcaZ();
          float dcazKaFromCasc = trackOmegaDauCharged.dcaZ();

          // primary vertex of the collision
          auto primaryVertex = getPrimaryVertex(collision); // get the associated covariance matrix with auto covMatrixPV = primaryVertex.getCov();
          std::array<float, 3> pvCoord = {collision.posX(), collision.posY(), collision.posZ()};

          if (doPvRefit && ((trackPion.pvRefitSigmaX2() != 1e10f) || (trackPion.pvRefitSigmaY2() != 1e10f) || (trackPion.pvRefitSigmaZ2() != 1e10f))) { // if I asked for PV refit in trackIndexSkimCreator.cxx
            pvCoord[0] = trackPion.pvRefitX();
            pvCoord[1] = trackPion.pvRefitY();
            pvCoord[2] = trackPion.pvRefitZ();

            // o2::dataformats::VertexBase Pvtx;
            primaryVertex.setX(trackPion.pvRefitX());
            primaryVertex.setY(trackPion.pvRefitY());
            primaryVertex.setZ(trackPion.pvRefitZ());
            primaryVertex.setCov(trackPion.pvRefitSigmaX2(), trackPion.pvRefitSigmaXY(), trackPion.pvRefitSigmaY2(), trackPion.pvRefitSigmaXZ(), trackPion.pvRefitSigmaYZ(), trackPion.pvRefitSigmaZ2());

            o2::dataformats::DCA impactParameterV0Dau0;
            o2::dataformats::DCA impactParameterV0Dau1;
            o2::dataformats::DCA impactParameterKaFromCasc;
            o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovV0Dau0, 2.f, matCorr, &impactParameterV0Dau0);
            o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovV0Dau1, 2.f, matCorr, &impactParameterV0Dau1);
            o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovOmegaDauCharged, 2.f, matCorr, &impactParameterKaFromCasc);
            dcaxyV0Dau0 = impactParameterV0Dau0.getY();
            dcaxyV0Dau1 = impactParameterV0Dau1.getY();
            dcaxyKaFromCasc = impactParameterKaFromCasc.getY();
            dcazV0Dau0 = impactParameterV0Dau0.getZ();
            dcazV0Dau1 = impactParameterV0Dau1.getZ();
            dcazKaFromCasc = impactParameterKaFromCasc.getZ();
          }

          // impact parameters
          o2::dataformats::DCA impactParameterCasc;
          o2::dataformats::DCA impactParameterPiFromCharmBaryon;
          o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackCasc, 2.f, matCorr, &impactParameterCasc);
          o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParVarPi, 2.f, matCorr, &impactParameterPiFromCharmBaryon);
          float impactParPiFromCharmBaryonXY = impactParameterPiFromCharmBaryon.getY();
          float impactParPiFromCharmBaryonZ = impactParameterPiFromCharmBaryon.getZ();

          // invariant mass under the hypothesis of particles ID corresponding to the decay chain
          double mLambda = casc.mLambda(); // from LF table, V0 mass under lambda hypothesis
          double mCasc = casc.mOmega();
          const std::array<double, 2> arrMassCharmBaryon = {massOmegaFromPDG, massPionFromPDG};
          double mCharmBaryon = RecoDecay::m(std::array{pVecCascAsD, pVecPionFromCharmBaryon}, arrMassCharmBaryon);

          // computing cosPA
          double cpaV0 = RecoDecay::cpa(vertexCasc, vertexV0, pVecV0);
          double cpaCharmBaryon = RecoDecay::cpa(pvCoord, coordVtxCharmBaryon, pVecCharmBaryon);
          double cpaCasc = RecoDecay::cpa(coordVtxCharmBaryon, vertexCasc, pVecCasc);
          double cpaxyV0 = RecoDecay::cpaXY(vertexCasc, vertexV0, pVecV0);
          double cpaxyCharmBaryon = RecoDecay::cpaXY(pvCoord, coordVtxCharmBaryon, pVecCharmBaryon);
          double cpaxyCasc = RecoDecay::cpaXY(coordVtxCharmBaryon, vertexCasc, pVecCasc);

          // computing decay length and ctau
          double decLenCharmBaryon = RecoDecay::distance(pvCoord, coordVtxCharmBaryon);
          double decLenCascade = RecoDecay::distance(coordVtxCharmBaryon, vertexCasc);
          double decLenV0 = RecoDecay::distance(vertexCasc, vertexV0);

          double phiCharmBaryon, thetaCharmBaryon;
          getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, coordVtxCharmBaryon, phiCharmBaryon, thetaCharmBaryon);
          auto errorDecayLengthCharmBaryon = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), phiCharmBaryon, thetaCharmBaryon) + getRotatedCovMatrixXX(covVtxCharmBaryon, phiCharmBaryon, thetaCharmBaryon));
          auto errorDecayLengthXYCharmBaryon = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), phiCharmBaryon, 0.) + getRotatedCovMatrixXX(covVtxCharmBaryon, phiCharmBaryon, 0.));

          double ctOmegac = RecoDecay::ct(pVecCharmBaryon, decLenCharmBaryon, massOmegacFromPDG);
          double ctCascade = RecoDecay::ct(pVecCasc, decLenCascade, massOmegaFromPDG);
          double ctV0 = RecoDecay::ct(pVecV0, decLenV0, massLambdaFromPDG);

          // computing eta
          double pseudorapCharmBaryon = RecoDecay::eta(pVecCharmBaryon);
          double pseudorapCascade = RecoDecay::eta(pVecCasc);
          double pseudorapV0 = RecoDecay::eta(pVecV0);

          // DCA between daughters
          float dcaCascDau = casc.dcacascdaughters();
          float dcaV0Dau = casc.dcaV0daughters();
          float dcaCharmBaryonDau = std::sqrt(df.getChi2AtPCACandidate());

          // set hfFlag
          int hfFlag = 1 << aod::hf_cand_omegac::DecayType::DecayToOmegaPi;

          // fill test histograms
          hInvMassCharmBaryon->Fill(mCharmBaryon);

          // fill the table
          rowCandidate(collision.globalIndex(),
                       pvCoord[0], pvCoord[1], pvCoord[2],
                       vertexCharmBaryonFromFitter[0], vertexCharmBaryonFromFitter[1], vertexCharmBaryonFromFitter[2],
                       vertexCasc[0], vertexCasc[1], vertexCasc[2],
                       vertexV0[0], vertexV0[1], vertexV0[2],
                       trackOmegaDauCharged.sign(),
                       covVtxCharmBaryon[0], covVtxCharmBaryon[1], covVtxCharmBaryon[2], covVtxCharmBaryon[3], covVtxCharmBaryon[4], covVtxCharmBaryon[5],
					   pVecCharmBaryon[0], pVecCharmBaryon[1], pVecCharmBaryon[2],
                       pVecCasc[0], pVecCasc[1], pVecCasc[2],
                       pVecPionFromCharmBaryon[0], pVecPionFromCharmBaryon[1], pVecPionFromCharmBaryon[2],
                       pVecV0[0], pVecV0[1], pVecV0[2],
                       pVecPionFromCasc[0], pVecPionFromCasc[1], pVecPionFromCasc[2],
                       pVecV0Dau0[0], pVecV0Dau0[1], pVecV0Dau0[2],
                       pVecV0Dau1[0], pVecV0Dau1[1], pVecV0Dau1[2],
                       impactParameterCasc.getY(), impactParPiFromCharmBaryonXY,
                       impactParameterCasc.getZ(), impactParPiFromCharmBaryonZ,
                       std::sqrt(impactParameterCasc.getSigmaY2()), std::sqrt(impactParameterPiFromCharmBaryon.getSigmaY2()),
                       v0index, casc.posTrackId(), casc.negTrackId(),
                       casc.globalIndex(), trackPion.globalIndex(), trackOmegaDauCharged.globalIndex(),
                       mLambda, mCasc, mCharmBaryon,
                       cpaV0, cpaCharmBaryon, cpaCasc, cpaxyV0, cpaxyCharmBaryon, cpaxyCasc,
                       ctOmegac, ctCascade, ctV0, 
                       pseudorapV0Dau0, pseudorapV0Dau1, pseudorapKaFromCas, pseudorapPiFromCharmBaryon,
                       pseudorapCharmBaryon, pseudorapCascade, pseudorapV0,
                       dcaxyV0Dau0, dcaxyV0Dau1, dcaxyKaFromCasc,
                       dcazV0Dau0, dcazV0Dau1, dcazKaFromCasc,
                       dcaCascDau, dcaV0Dau, dcaCharmBaryonDau,
                       decLenCharmBaryon, decLenCascade, decLenV0, errorDecayLengthCharmBaryon, errorDecayLengthXYCharmBaryon,
                       hfFlag);

        } // loop over pions
	toCounter = 1;
      }   // loop over cascades
    }     // close loop collisions
  }       // end of process
  PROCESS_SWITCH(HfCandidateCreatorOmegacToOmegaPi, processIdxCombinatorics, "Do indexes combinatorics", true);

  void processDerivedData(SelectedCollisions const& collisions,
                          aod::BCsWithTimestamps const& bcWithTimeStamps,
                          MyTracks const& tracks,
                          MyCascTable const& cascades,
                          MyV0Table const&,
                          aod::V0sLinked const&,
                          MySkimIdx const& candidates)
  {

    double massPionFromPDG = MassPiPlus;    // pdg code 211
    double massLambdaFromPDG = MassLambda0; // pdg code 3122
	double massOmegaFromPDG = MassOmegaMinus;
    double massOmegacFromPDG = MassOmegaC0; // pdg code 4332

    // 2-prong vertex fitter to build the omegac vertex
    o2::vertexing::DCAFitterN<2> df;
    df.setPropagateToPCA(propagateToPCA);
    df.setMaxR(maxR);
    df.setMaxDZIni(maxDZIni);
    df.setMaxDXYIni(maxDXYIni);
    df.setMinParamChange(minParamChange);
    df.setMinRelChi2Change(minRelChi2Change);
    df.setMaxChi2(maxChi2);
    df.setUseAbsDCA(useAbsDCA);
    df.setWeightedFinalPCA(useWeightedFinalPCA);

    for (const auto& collision : collisions) {

      // set the magnetic field from CCDB
      auto bc = collision.bc_as<o2::aod::BCsWithTimestamps>();
      initCCDB(bc, runNumber, ccdb, isRun2 ? ccdbPathGrp : ccdbPathGrpMag, lut, isRun2);
      auto magneticField = o2::base::Propagator::Instance()->getNominalBz(); // z component

      df.setBz(magneticField);
      df.setRefitWithMatCorr(refitWithMatCorr);

      // loop over cascades reconstructed by cascadebuilder.cxx
      auto thisCollId = collision.globalIndex();
      auto groupedCandidates = candidates.sliceBy(candidatesPerCollision, thisCollId);

      for (const auto& cand : groupedCandidates) {

        hCandidateCounter->Fill(0);

        if (!TESTBIT(cand.hfflag(), aod::hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi)) {
          continue;
        }

        hCandidateCounter->Fill(1);

        auto casc = cand.cascade_as<MyCascTable>();
        auto trackPion = cand.prong0_as<MyTracks>();           // pi <-- charm baryon
        auto trackOmegaDauCharged = casc.bachelor_as<MyTracks>(); // kaon <- Omega track
		int v0index = casc.cascadeId();
        auto trackV0Dau0 = casc.posTrack_as<MyTracks>(); // V0 positive daughter track
        auto trackV0Dau1 = casc.negTrack_as<MyTracks>(); // V0 negative daughter track

        //-------------------------- V0 info---------------------------
        // pseudorapidity
        double pseudorapV0Dau0 = trackV0Dau0.eta();
        double pseudorapV0Dau1 = trackV0Dau1.eta();

        // pion & p <- V0 tracks
        auto trackParCovV0Dau0 = getTrackParCov(trackV0Dau0);
        auto trackParCovV0Dau1 = getTrackParCov(trackV0Dau1);

        // info from LF table
        std::array<float, 3> pVecV0 = {casc.pxlambda(), casc.pylambda(), casc.pzlambda()};
        std::array<float, 3> vertexV0 = {casc.xlambda(), casc.ylambda(), casc.zlambda()};
        std::array<float, 3> pVecV0Dau0 = {casc.pxpos(), casc.pypos(), casc.pzpos()};
        std::array<float, 3> pVecV0Dau1 = {casc.pxneg(), casc.pyneg(), casc.pzneg()};

        //-------------------reconstruct cascade track------------------
        // pseudorapidity
        double pseudorapKaFromCas = trackOmegaDauCharged.eta();

        auto trackParCovOmegaDauCharged = getTrackParCov(trackOmegaDauCharged);

        // info from LF table
        std::array<float, 3> vertexCasc = {casc.x(), casc.y(), casc.z()};
        std::array<float, 3> pVecCasc = {casc.px(), casc.py(), casc.pz()};
        std::array<float, 21> covCasc = {0.};
        constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
        for (int i = 0; i < 6; i++) {
          covCasc[MomInd[i]] = casc.momentumCovMat()[i];
          covCasc[i] = casc.positionCovMat()[i];
        }
        // create cascade track
        o2::track::TrackParCov trackCasc;
        if (trackOmegaDauCharged.sign() > 0) {
          trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, 1, true);
        } else if (trackOmegaDauCharged.sign() < 0) {
          trackCasc = o2::track::TrackParCov(vertexCasc, pVecCasc, covCasc, -1, true);
        } else {
          continue;
        }
        trackCasc.setAbsCharge(1);
        trackCasc.setPID(o2::track::PID::OmegaMinus);

        std::array<float, 3> pVecPionFromCasc = {casc.pxbach(), casc.pybach(), casc.pzbach()};

        //------------reconstruct charm baryon decay vtx---------------
        auto trackParVarPi = getTrackParCov(trackPion); // charm bachelor pion track to be processed with DCAFitter

        // reconstruct charm baryon with DCAFitter
        int nVtxFromFitterCharmBaryon = 0;
        try {
          nVtxFromFitterCharmBaryon = df.process(trackCasc, trackParVarPi);
        } catch (...) {
          LOG(error) << "Exception caught in charm DCA fitter process call!";
          hFitterStatus->Fill(1);
          continue;
        }
        if (nVtxFromFitterCharmBaryon == 0) {
          hFitterStatus->Fill(2);
          continue;
        }
        hFitterStatus->Fill(0);
        hCandidateCounter->Fill(2);
        auto vertexCharmBaryonFromFitter = df.getPCACandidate();
        //auto chi2PCACharmBaryon = df.getChi2AtPCACandidate();
        std::array<float, 3> pVecCascAsD;
        std::array<float, 3> pVecPionFromCharmBaryon;
        df.propagateTracksToVertex();
        if (!df.isPropagateTracksToVertexDone()) {
          continue;
        }
        df.getTrack(0).getPxPyPzGlo(pVecCascAsD);
        df.getTrack(1).getPxPyPzGlo(pVecPionFromCharmBaryon);
        std::array<float, 3> pVecCharmBaryon = {pVecCascAsD[0] + pVecPionFromCharmBaryon[0], pVecCascAsD[1] + pVecPionFromCharmBaryon[1], pVecCascAsD[2] + pVecPionFromCharmBaryon[2]};
        std::array<float, 3> coordVtxCharmBaryon = df.getPCACandidatePos();
        std::array<float, 6> covVtxCharmBaryon = df.calcPCACovMatrixFlat();

        // pseudorapidity
        double pseudorapPiFromCharmBaryon = trackPion.eta();

        // DCAxy (computed with propagateToDCABxByBz method)
        float dcaxyV0Dau0 = trackV0Dau0.dcaXY();
        float dcaxyV0Dau1 = trackV0Dau1.dcaXY();
        float dcaxyKaFromCasc = trackOmegaDauCharged.dcaXY();

        // DCAz (computed with propagateToDCABxByBz method)
        float dcazV0Dau0 = trackV0Dau0.dcaZ();
        float dcazV0Dau1 = trackV0Dau1.dcaZ();
        float dcazKaFromCasc = trackOmegaDauCharged.dcaZ();

        // primary vertex of the collision
        auto primaryVertex = getPrimaryVertex(collision); // get the associated covariance matrix with auto covMatrixPV = primaryVertex.getCov();
        std::array<float, 3> pvCoord = {collision.posX(), collision.posY(), collision.posZ()};

        if (doPvRefit && ((trackPion.pvRefitSigmaX2() != 1e10f) || (trackPion.pvRefitSigmaY2() != 1e10f) || (trackPion.pvRefitSigmaZ2() != 1e10f))) { // if I asked for PV refit in trackIndexSkimCreator.cxx
          pvCoord[0] = trackPion.pvRefitX();
          pvCoord[1] = trackPion.pvRefitY();
          pvCoord[2] = trackPion.pvRefitZ();

          // o2::dataformats::VertexBase Pvtx;
          primaryVertex.setX(trackPion.pvRefitX());
          primaryVertex.setY(trackPion.pvRefitY());
          primaryVertex.setZ(trackPion.pvRefitZ());
          primaryVertex.setCov(trackPion.pvRefitSigmaX2(), trackPion.pvRefitSigmaXY(), trackPion.pvRefitSigmaY2(), trackPion.pvRefitSigmaXZ(), trackPion.pvRefitSigmaYZ(), trackPion.pvRefitSigmaZ2());

          o2::dataformats::DCA impactParameterV0Dau0;
          o2::dataformats::DCA impactParameterV0Dau1;
          o2::dataformats::DCA impactParameterKaFromCasc;
          o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovV0Dau0, 2.f, matCorr, &impactParameterV0Dau0);
          o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovV0Dau1, 2.f, matCorr, &impactParameterV0Dau1);
          o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovOmegaDauCharged, 2.f, matCorr, &impactParameterKaFromCasc);
          dcaxyV0Dau0 = impactParameterV0Dau0.getY();
          dcaxyV0Dau1 = impactParameterV0Dau1.getY();
          dcaxyKaFromCasc = impactParameterKaFromCasc.getY();
          dcazV0Dau0 = impactParameterV0Dau0.getZ();
          dcazV0Dau1 = impactParameterV0Dau1.getZ();
          dcazKaFromCasc = impactParameterKaFromCasc.getZ();
        }

        // impact parameters
        o2::dataformats::DCA impactParameterCasc;
        o2::dataformats::DCA impactParameterPiFromCharmBaryon;
        o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackCasc, 2.f, matCorr, &impactParameterCasc);
        o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParVarPi, 2.f, matCorr, &impactParameterPiFromCharmBaryon);
        float impactParPiFromCharmBaryonXY = impactParameterPiFromCharmBaryon.getY();
        float impactParPiFromCharmBaryonZ = impactParameterPiFromCharmBaryon.getZ();

        // invariant mass under the hypothesis of particles ID corresponding to the decay chain
        double mLambda = casc.mLambda(); // from LF table, V0 mass under lambda hypothesis
        double mCasc = casc.mOmega();
        const std::array<double, 2> arrMassCharmBaryon = {massOmegaFromPDG, massPionFromPDG};
        double mCharmBaryon = RecoDecay::m(std::array{pVecCascAsD, pVecPionFromCharmBaryon}, arrMassCharmBaryon);

        // computing cosPA
        double cpaV0 = RecoDecay::cpa(vertexCasc, vertexV0, pVecV0);
        double cpaCharmBaryon = RecoDecay::cpa(pvCoord, coordVtxCharmBaryon, pVecCharmBaryon);
        double cpaCasc = RecoDecay::cpa(coordVtxCharmBaryon, vertexCasc, pVecCasc);
        double cpaxyV0 = RecoDecay::cpaXY(vertexCasc, vertexV0, pVecV0);
        double cpaxyCharmBaryon = RecoDecay::cpaXY(pvCoord, coordVtxCharmBaryon, pVecCharmBaryon);
        double cpaxyCasc = RecoDecay::cpaXY(coordVtxCharmBaryon, vertexCasc, pVecCasc);

        // computing decay length and ctau
        double decLenCharmBaryon = RecoDecay::distance(pvCoord, coordVtxCharmBaryon);
        double decLenCascade = RecoDecay::distance(coordVtxCharmBaryon, vertexCasc);
        double decLenV0 = RecoDecay::distance(vertexCasc, vertexV0);

        double phiCharmBaryon, thetaCharmBaryon;
        getPointDirection(std::array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, coordVtxCharmBaryon, phiCharmBaryon, thetaCharmBaryon);
        auto errorDecayLengthCharmBaryon = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), phiCharmBaryon, thetaCharmBaryon) + getRotatedCovMatrixXX(covVtxCharmBaryon, phiCharmBaryon, thetaCharmBaryon));
        auto errorDecayLengthXYCharmBaryon = std::sqrt(getRotatedCovMatrixXX(primaryVertex.getCov(), phiCharmBaryon, 0.) + getRotatedCovMatrixXX(covVtxCharmBaryon, phiCharmBaryon, 0.));

        double ctOmegac = RecoDecay::ct(pVecCharmBaryon, decLenCharmBaryon, massOmegacFromPDG);
        double ctCascade = RecoDecay::ct(pVecCasc, decLenCascade, massOmegaFromPDG);
        double ctV0 = RecoDecay::ct(pVecV0, decLenV0, massLambdaFromPDG);

        // computing eta
        double pseudorapCharmBaryon = RecoDecay::eta(pVecCharmBaryon);
        double pseudorapCascade = RecoDecay::eta(pVecCasc);
        double pseudorapV0 = RecoDecay::eta(pVecV0);

        // DCA between daughters
        float dcaCascDau = casc.dcacascdaughters();
        float dcaV0Dau = casc.dcaV0daughters();
        float dcaCharmBaryonDau = std::sqrt(df.getChi2AtPCACandidate());

        // set hfFlag
        int hfFlag = 1 << aod::hf_cand_omegac::DecayType::DecayToOmegaPi;

        // fill test histograms
        hInvMassCharmBaryon->Fill(mCharmBaryon);
        hCandidateCounter->Fill(3);

        // fill the table
        rowCandidate(collision.globalIndex(),
                     pvCoord[0], pvCoord[1], pvCoord[2],
                     vertexCharmBaryonFromFitter[0], vertexCharmBaryonFromFitter[1], vertexCharmBaryonFromFitter[2],
                     vertexCasc[0], vertexCasc[1], vertexCasc[2],
                     vertexV0[0], vertexV0[1], vertexV0[2],
                     trackOmegaDauCharged.sign(),
                     covVtxCharmBaryon[0], covVtxCharmBaryon[1], covVtxCharmBaryon[2], covVtxCharmBaryon[3], covVtxCharmBaryon[4], covVtxCharmBaryon[5],
					 pVecCharmBaryon[0], pVecCharmBaryon[1], pVecCharmBaryon[2],
                     pVecCasc[0], pVecCasc[1], pVecCasc[2],
                     pVecPionFromCharmBaryon[0], pVecPionFromCharmBaryon[1], pVecPionFromCharmBaryon[2],
                     pVecV0[0], pVecV0[1], pVecV0[2],
                     pVecPionFromCasc[0], pVecPionFromCasc[1], pVecPionFromCasc[2],
                     pVecV0Dau0[0], pVecV0Dau0[1], pVecV0Dau0[2],
                     pVecV0Dau1[0], pVecV0Dau1[1], pVecV0Dau1[2],
                     impactParameterCasc.getY(), impactParPiFromCharmBaryonXY,
                     impactParameterCasc.getZ(), impactParPiFromCharmBaryonZ,
                     std::sqrt(impactParameterCasc.getSigmaY2()), std::sqrt(impactParameterPiFromCharmBaryon.getSigmaY2()),
                     v0index, casc.posTrackId(), casc.negTrackId(),
                     casc.globalIndex(), trackPion.globalIndex(), trackOmegaDauCharged.globalIndex(),
                     mLambda, mCasc, mCharmBaryon,
                     cpaV0, cpaCharmBaryon, cpaCasc, cpaxyV0, cpaxyCharmBaryon, cpaxyCasc,
                     ctOmegac, ctCascade, ctV0, 
                     pseudorapV0Dau0, pseudorapV0Dau1, pseudorapKaFromCas, pseudorapPiFromCharmBaryon,
                     pseudorapCharmBaryon, pseudorapCascade, pseudorapV0,
                     dcaxyV0Dau0, dcaxyV0Dau1, dcaxyKaFromCasc,
                     dcazV0Dau0, dcazV0Dau1, dcazKaFromCasc,
                     dcaCascDau, dcaV0Dau, dcaCharmBaryonDau,
                     decLenCharmBaryon, decLenCascade, decLenV0, errorDecayLengthCharmBaryon, errorDecayLengthXYCharmBaryon,
                     hfFlag);

      } // loop over LF Cascade-bachelor candidates
    }   // loop over collisions
  }     // end of process
  PROCESS_SWITCH(HfCandidateCreatorOmegacToOmegaPi, processDerivedData, "Process derived data", false);

}; // end of struct

/// Performs MC matching.
struct HfCandidateCreatorOmegacToOmegaPiMc {
  Produces<aod::HfOmegaCMCRec> rowMCMatchRec;
  Produces<aod::HfOmegaCMCGen> rowMCMatchGen;

  Configurable<bool> matchOmegacMc{"matchOmegacMc", true, "Do MC matching for Omegac0"};

  void init(InitContext const&) {}

  void processDoNoMc(aod::Collisions::iterator const& collision)
  {
    // dummy process function - should not be required in the future
  }
  PROCESS_SWITCH(HfCandidateCreatorOmegacToOmegaPiMc, processDoNoMc, "Do not run MC process function", true);

  void processMc(aod::HfCandOmegaC const& candidates,
                 aod::TracksWMc const& tracks,
                 aod::McParticles const& mcParticles,
				 aod::McCollisionLabels const&)
  {
    float ptCharmBaryonGen = -999.;
	float etaCharmBaryonGen = -999.;
    int indexRec = -1;
    int indexRecCharmBaryon = -1;
    int8_t sign = -9;
    int8_t flag = 0;
    int8_t origin = 0; // to be used for prompt/non prompt
    int8_t debug = 0;
    int8_t debugGenCharmBar = 0;
    int8_t debugGenOmega = 0;
    int8_t debugGenLambda = 0;

    int pdgCodeOmegac0 = Pdg::kOmegaC0;       // 4332            // 3312
	int pdgCodeOmegaMinus = kOmegaMinus;
	int pdgCodeKaMinus = kKMinus;
    int pdgCodeLambda = kLambda0;             // 3122
    int pdgCodePiPlus = kPiPlus;              // 211
    int pdgCodePiMinus = kPiMinus;            // -211
    int pdgCodeProton = kProton;              // 2212
	//bool collisionMatched = false;
    // Match reconstructed candidates.
    for (const auto& candidate : candidates) {
      flag = 0;
      origin = RecoDecay::OriginType::None;
      debug = 0;
	  bool collisionMatched = false;
      auto arrayDaughters = std::array{candidate.piFromCharmBaryon_as<aod::TracksWMc>(), // pi <- charm baryon
                                       candidate.bachelor_as<aod::TracksWMc>(),          // Kaon <- cascade
                                       candidate.posTrack_as<aod::TracksWMc>(),          // p <- lambda
                                       candidate.negTrack_as<aod::TracksWMc>()};         // pi <- lambda
      auto arrayDaughtersCasc = std::array{candidate.bachelor_as<aod::TracksWMc>(),
                                           candidate.posTrack_as<aod::TracksWMc>(),
                                           candidate.negTrack_as<aod::TracksWMc>()};
      auto arrayDaughtersV0 = std::array{candidate.posTrack_as<aod::TracksWMc>(),
                                         candidate.negTrack_as<aod::TracksWMc>()};

      // Omegac matching
      if (matchOmegacMc) {
        // Omegac → pi pi pi p
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, pdgCodeOmegac0, std::array{pdgCodePiPlus, pdgCodeKaMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 3);
        indexRecCharmBaryon = indexRec;
        if (indexRec == -1) {
          debug = 1;
        }
        if (indexRec > -1) {
          // Omega- → kaon pi p
          indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersCasc, pdgCodeOmegaMinus, std::array{pdgCodeKaMinus, pdgCodeProton, pdgCodePiMinus}, true, &sign, 2);
          if (indexRec == -1) {
            debug = 2;
          }
          if (indexRec > -1) {
            // Lambda → p pi
            indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughtersV0, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true, &sign, 1);
            if (indexRec == -1) {
              debug = 3;
            }
            if (indexRec > -1) {
              flag = sign * (1 << aod::hf_cand_omegac::DecayType::OmegaczeroToOmegaPi);
			  collisionMatched = candidate.collision_as<aod::McCollisionLabels>().mcCollisionId() == mcParticles.iteratorAt(indexRecCharmBaryon).mcCollisionId();
            }
          }
        }

        // Check whether the charm baryon is non-prompt (from a b quark).
        if (flag != 0) {
          auto particle = mcParticles.rawIteratorAt(indexRecCharmBaryon);
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, true);
        }
      }
      if (debug == 2 || debug == 3) {
        LOGF(info, "WARNING: Charm baryon decays in the expected final state but the condition on the intermediate states are not fulfilled");
      }
      rowMCMatchRec(flag, debug, origin,collisionMatched);

    } // close loop over candidates

    // Match generated particles.
    for (const auto& particle : mcParticles) {
      ptCharmBaryonGen = -999.;
	  etaCharmBaryonGen = -999.;
      flag = 0;
      sign = -9;
      debugGenCharmBar = 0;
      debugGenOmega = 0;
      debugGenLambda = 0;
      origin = RecoDecay::OriginType::None;
      if (matchOmegacMc) {
        //  Omegac → Omega pi
        if (RecoDecay::isMatchedMCGen(mcParticles, particle, pdgCodeOmegac0, std::array{pdgCodeOmegaMinus, pdgCodePiPlus}, true, &sign)) {
          debugGenCharmBar = 1;
          ptCharmBaryonGen = particle.pt();
		  etaCharmBaryonGen = particle.eta();
          // Omega -> Lambda kaon
          auto cascMC = mcParticles.rawIteratorAt(particle.daughtersIds().front());
          if (RecoDecay::isMatchedMCGen(mcParticles, cascMC, pdgCodeOmegaMinus, std::array{pdgCodeLambda, pdgCodeKaMinus}, true)) {
            debugGenOmega = 1;
            // Lambda -> p pi
            auto v0MC = mcParticles.rawIteratorAt(cascMC.daughtersIds().front());
            if (RecoDecay::isMatchedMCGen(mcParticles, v0MC, pdgCodeLambda, std::array{pdgCodeProton, pdgCodePiMinus}, true)) {
              debugGenLambda = 1;
              flag = sign * (1 << aod::hf_cand_omegac::DecayType::OmegaczeroToOmegaPi);
            }
          }
        }
        // Check whether the charm baryon is non-prompt (from a b quark)
        if (flag != 0) {
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, true);
        }
      }

      rowMCMatchGen(flag, debugGenCharmBar, debugGenOmega, debugGenLambda, ptCharmBaryonGen, etaCharmBaryonGen, origin);
    }
  } // close process
  PROCESS_SWITCH(HfCandidateCreatorOmegacToOmegaPiMc, processMc, "Process MC", true);//false
}; // close struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorOmegacToOmegaPi>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorOmegacToOmegaPiMc>(cfgc)};
}
