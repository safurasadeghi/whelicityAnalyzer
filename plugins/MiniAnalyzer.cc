// system include files
#include <memory>
#include <iostream>
#include <cmath>
#include <vector>
#include <TBranch.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <string>
#include <map>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <TRandom3.h>
#include <TBranch.h>
#include <TClonesArray.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "AnalysisDataFormats/TopObjects/interface/TtEvent.h"
#include "TopQuarkAnalysis/TopKinFitter/interface/TtFullLepKinSolver.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "AnalysisDataFormats/TopObjects/interface/TtDilepEvtSolution.h"
#include "ttamwtsolver.h"




#include "DileptonAnalyticalSolver.h"

using namespace std;
using namespace pat;
using namespace edm;
using namespace reco;
//
// class declaration
//

class MiniAnalyzer : public edm::EDAnalyzer  {
public:
    explicit MiniAnalyzer (const edm::ParameterSet&);
    ~MiniAnalyzer();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;


    // ----------member data ---------------------------
    edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
    edm::EDGetTokenT<pat::MuonCollection> muonToken_;
    edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
    edm::EDGetTokenT<pat::TauCollection> tauToken_;
    edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
    edm::EDGetTokenT<pat::JetCollection>  jetToken_;
    edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
    edm::EDGetTokenT<pat::METCollection> metToken_;
    edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
    edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
    edm::EDGetTokenT<GenEventInfoProduct> genEvtInfo_;
    edm::EDGetTokenT<LHEEventProduct> lheEvtInfo_;
    edm::EDGetTokenT<double> rho_;

    //Boolean variables
    bool isPosMu = false;
    bool isNegMu = false;
    bool isPosEl = false;
    bool isNegEl = false;
    bool isDiMuon = false;
    bool isDiElectron = false;
    bool isElMu = false;
    bool isMuEl = false;
    bool isDiLeptonic = false;


    // function declaration
    bool IsSoftMuon(const pat::Muon&, const reco::Vertex&);

    // histogram declaration
    edm::Service<TFileService> fs;
    const reco::GenParticle* mother;

    TH1F* h_NPV;
    //After Lepton Selection
    TH1F* h_ALS_etaLMuMu;
    TH1F* h_ALS_ptLMuMu;
    TH1F* h_ALS_etaLElEl;
    TH1F* h_ALS_ptLElEl;
    TH1F* h_ALS_etaLElMu;
    TH1F* h_ALS_ptLElMu;
    TH1F* h_ALS_etaLDiLep;
    TH1F* h_ALS_ptLDiLep;
    //After MET Requirements
    TH1F* h_AMS_ptLepMuMu;
    TH1F* h_AMS_ptLepElMu;
    TH1F* h_AMS_ptLepElEl;
    TH1F* h_AMS_ptLepDiLep;
    TH1F* h_AMS_mLepMuMu;
    TH1F* h_AMS_mLepElEl;
    TH1F* h_AMS_mLepElMu;
    TH1F* h_AMS_mLepDiLep;
    TH1F* h_AMS_METMuMu;
    TH1F* h_AMS_METElEl;
    TH1F* h_AMS_METElMu;
    TH1F* h_AMS_METDiLep;
    //After jet Selectin
    TH1F* h_AJS_etaLepMuMu;
    TH1F* h_AJS_ptLepMuMu;
    TH1F* h_AJS_etaLepElEl;
    TH1F* h_AJS_ptLepElEl;
    TH1F* h_AJS_etaLepElMu;
    TH1F* h_AJS_ptLepElMu;
    TH1F* h_AJS_etaLepDiLep;
    TH1F* h_AJS_ptLepDiLep;
    //After Btag requirement
    TH1F* h_ABS_etaLepMuMu;
    TH1F* h_ABS_ptLepMuMu;
    TH1F* h_ABS_etaLepElEl;
    TH1F* h_ABS_ptLepElEl;
    TH1F* h_ABS_etaLepElMu;
    TH1F* h_ABS_ptLepElMu;
    TH1F* h_ABS_etaLepDiLep;
    TH1F* h_ABS_ptLepDiLep;
    TH1F* h_ABS_mLepMuMu;
    TH1F* h_ABS_mLepElEl;
    TH1F* h_ABS_mLepElMu;
    TH1F* h_ABS_mLepDiLep;
    TH1F* h_ABS_METMuMu;
    TH1F* h_ABS_METElEl;
    TH1F* h_ABS_METElMu;
    TH1F* h_ABS_METDiLep;
    TH1F* h_ABS_NJetsMuMu;
    TH1F* h_ABS_NJetsElMu;
    TH1F* h_ABS_NJetsElEl;
    TH1F* h_ABS_NJetsDiLep;
    TH1F* h_ABS_NBJetsMuMu;
    TH1F* h_ABS_NBJetsElEl;
    TH1F* h_ABS_NBJetsElMu;
    TH1F* h_ABS_NBJetsDiLep;
    TH1F* h_ABS_ptLeadingJetMuMu;
    TH1F* h_ABS_ptLeadingJetElEl;
    TH1F* h_ABS_ptLeadingJetElMu;
    TH1F* h_ABS_ptLeadingJetDiLep;
    TH1F* h_ABS_etaLeadingJetMuMu;
    TH1F* h_ABS_etaLeadingJetElEl;
    TH1F* h_ABS_etaLeadingJetElMu;
    TH1F* h_ABS_etaLeadingJetDiLep;

    //After tt reconstruction
    TH1F* h_NBJetsMuMu;
    TH1F* h_NBJetsElMu;
    TH1F* h_NBJetsElEl;
    TH1F* h_NBJetsDiLep;
    TH1F* h_cosMuMu;
    TH1F* h_cosElEl;
    TH1F* h_cosElMu;
    TH1F* h_cosDiLep;
    TH1F* h_cosGenMuMu;
    TH1F* h_cosGenElEl;
    TH1F* h_cosGenElMu;
    TH1F* h_cosGen;
    TH1F* h_mTTbarMuMu;
    TH1F* h_mTTbarElEl;
    TH1F* h_mTTbarElMu;
    TH1F* h_TTbarM;
    TH1F* h_GenTTbarM;
    TH1F* h_ptTMuMu;
    TH1F* h_ptTElEl;
    TH1F* h_ptTElMu;
    TH1F* h_ptTDiLep;
    TH1F* h_yTMuMu;
    TH1F* h_yTElEl;
    TH1F* h_yTElMu;
    TH1F* h_yTDiLep;
    TH1F* h_ptWMuMu;
    TH1F* h_ptWElEl;
    TH1F* h_ptWElMu;
    TH1F* h_ptWDiLep;
    TH1F* h_yWMuMu;
    TH1F* h_yWElEl;
    TH1F* h_yWElMu;
    TH1F* h_yWDiLep;
    TH1F* h_PtMu;
    TH1F* h_etaMu;

    double coriso= 999; // initialise to dummy value
    double coriso2= 999; // initialise to dummy value
    double dileptonMass;
    TLorentzVector Neutrino;
    vector<TLorentzVector> Neutrinos;
    TLorentzVector AntiNeutrino;
    vector<TLorentzVector> AntiNeutrinos;

    TLorentzVector W1;
    TLorentzVector W2;
    TLorentzVector t1;
    TLorentzVector t2;
    TLorentzVector ttbar;
    TLorentzVector lepPos;
    TLorentzVector lepNeg;
    TLorentzVector BJet;
    TLorentzVector BBJet;
    TLorentzVector dilepton;

    //    TtFullLepKinSolver* solver;
    TtAMWTSolver* amwtSolver;
    EffectiveAreas effectiveAreas_;
    edm::EDGetTokenT<TtGenEvent> ttgenEvt_;
    bool isData;

    int NEvent=0;




};


MiniAnalyzer::MiniAnalyzer(const edm::ParameterSet& iConfig):
    vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
    electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
    tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
    photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
    jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
    fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
    metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
    prunedGenToken_(mayConsume<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
    packedGenToken_(mayConsume<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
    genEvtInfo_(mayConsume<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEvtInfo"))),
    lheEvtInfo_(mayConsume<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEvtInfo"))),
    rho_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
    effectiveAreas_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath() ),
    ttgenEvt_( mayConsume<TtGenEvent>(iConfig.getParameter<edm::InputTag>("ttgen"))),
    isData((iConfig.getParameter<bool>("isData")))
{
    TH1::SetDefaultSumw2();

    h_cosMuMu = fs->make<TH1F>("h_cosMuMu",";cos(#theta);",10,-1,1);
    h_cosElEl = fs->make<TH1F>("h_cosElEl",";cos(#theta);",10,-1,1);
    h_cosElMu = fs->make<TH1F>("h_cosElMu",";cos(#theta);",10,-1,1);
    h_cosGenElMu = fs->make<TH1F>("h_cosGenElMu",";cos(#theta);",10,-1,1);
    h_cosGenMuMu = fs->make<TH1F>("h_cosGenMuMu",";cos(#theta);",10,-1,1);
    h_cosGenElEl = fs->make<TH1F>("h_cosGenElEl",";cos(#theta);",10,-1,1);
    h_cosGen = fs->make<TH1F>("h_cosGen",";cos(#theta);",10,-1,1);
    h_GenTTbarM = fs->make<TH1F>("h_GenTTbarM",";M_{TTbar};",100,90,1300);
    h_etaMu = fs->make<TH1F>("h_EtaMu",";#eta_{l};",100,-3,3);
    h_TTbarM = fs->make<TH1F>("h_TTbarM",";M_{TTbar};",100,90,1300);
    h_PtMu = fs->make<TH1F>("h_PtMu",";Mu_{Pt};",100,0.,200.);
    h_NPV = fs->make<TH1F>("h_NPV",";NPV;",100,0,30);
    h_ALS_etaLMuMu = fs->make<TH1F>("h_EtaLepMuMu",";#eta_{l};",100,-3,3);
    h_ALS_etaLElEl = fs->make<TH1F>("h_EtaLepElEl",";#eta_{l};",100,-3,3);
    h_ALS_etaLElMu = fs->make<TH1F>("h_EtaLepElMu",";#eta_{l};",100,-3,3);
    h_ALS_etaLDiLep = fs->make<TH1F>("h_EtaLepDiLep",";#eta_{l};",100,-3,3);
    h_ALS_ptLMuMu = fs->make<TH1F>("h_PtLepMuMu",";Pt_{l};",100,0.,300.);
    h_ALS_ptLElMu = fs->make<TH1F>("h_PtLepElMu",";Pt_{l};",100,0.,300.);
    h_ALS_ptLElEl = fs->make<TH1F>("h_PtLepElEl",";Pt_{l};",100,0.,300.);
    h_ALS_ptLDiLep = fs->make<TH1F>("h_PtLepDiLep",";Pt_{l};",100,0.,300.);

    amwtSolver = new TtAMWTSolver(isData,171.5,173.5,100,80.4,4.8);
}

MiniAnalyzer::~MiniAnalyzer()
{
}

void
MiniAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;


    /////////////FILLINGINPUTS/////////////////////
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);

    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByToken(muonToken_, muons);

    edm::Handle<pat::ElectronCollection> electrons;
    iEvent.getByToken(electronToken_, electrons);

    edm::Handle<pat::JetCollection> jets;
    iEvent.getByToken(jetToken_, jets);

    edm::Handle<pat::METCollection> mets;
    iEvent.getByToken(metToken_, mets);

    Handle<edm::View<reco::GenParticle> > pruned;
    Handle<edm::View<pat::PackedGenParticle> > packed;
    edm::Handle<GenEventInfoProduct> genEvtInfo;
    edm::Handle<LHEEventProduct> lhEvtInfo;
    edm::Handle<TtGenEvent> genEvent;
    if(!isData){
        iEvent.getByToken(prunedGenToken_,pruned);
        iEvent.getByToken(packedGenToken_,packed);
        iEvent.getByToken( genEvtInfo_, genEvtInfo );
        iEvent.getByToken( lheEvtInfo_, lhEvtInfo);
        iEvent.getByToken(ttgenEvt_, genEvent);
    }
    edm::Handle<double> rhoHandle_;
    iEvent.getByToken( rho_,rhoHandle_);
    float rho = *rhoHandle_;

    amwtSolver->SetRho(rho);

    isPosMu = false;
    isNegMu = false;
    isPosEl = false;
    isNegEl = false;
    isDiMuon = false;
    isDiElectron = false;
    isElMu = false;
    isMuEl = false;
    isDiLeptonic = false;

    /////////////////////////////////////////////////
    // make sure we have a good vertex //////////////
    /////////////////////////////////////////////////
    ++NEvent;
    cout << "number of Events " << NEvent << endl;
    if (vertices->empty()) return;
    VertexCollection::const_iterator PV = vertices->end();
    for (VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx, ++PV) {
        if ( !(vtx->isFake())
             && vtx->ndof() >= 4. && vtx->position().Rho() < 2.0
             && fabs(vtx->position().Z()) < 24.0) {
            PV = vtx;
            break;
        }
    }
    if (PV==vertices->end()) return;

    // count how many good vertices we have
    int nGoodVtxs = 0;
    for (VertexCollection::const_iterator vtx = vertices->begin();vtx != vertices->end(); ++vtx) {
        if ( !(vtx->isFake()) && vtx->ndof() >= 4. && vtx->position().Rho() <= 2.0 && fabs(vtx->position().Z()) <= 24.) nGoodVtxs++;
    }
    //////////////EVENT WEIGHT/////////
    /// finding weight so we can
    /// compare different backgrounds
    double theWeight;
    if(isData)
    {
        theWeight = 1;
    }
    else
    {
        theWeight = genEvtInfo->weight();
        int whichWeight = 2;
        theWeight *= lhEvtInfo->weights()[whichWeight].wgt/lhEvtInfo->originalXWGTUP();
    }

    h_NPV->Fill(nGoodVtxs,theWeight);


    //////////////MUON////////////////
    /// finding pos and neg muon with highest pt
    ///
    pat::Muon posMu;
    pat::Muon negMu;



    for (pat::MuonCollection::const_iterator mup = muons->begin();mup != muons->end(); ++mup){
        if( !(mup->charge() > 0)) continue;
        if( !(mup->pt() > 20.0 )) continue;
        if( !(fabs(mup->eta()) < 2.4 )) continue;
        if( !(mup->isPFMuon())) continue;
        if( !(mup->isGlobalMuon() || mup->isTrackerMuon())) continue;
        if( mup->isIsolationValid()){
            reco::MuonPFIsolation pfR04 = mup->pfIsolationR04();
            coriso = pfR04.sumChargedHadronPt + std::max(0., pfR04.sumNeutralHadronEt+pfR04.sumPhotonEt-0.5*pfR04.sumPUPt);
        }
        if (!(coriso/mup->pt() <  0.15)) continue;
        if(mup->pt() > posMu.pt() ) posMu = *mup;
        cout << mup->pt() <<endl;
        cout << posMu.pt() << endl;
        h_PtMu->Fill(mup->pt(),theWeight);
        h_ALS_ptLMuMu->Fill(posMu.pt(),theWeight);
        h_ALS_etaLMuMu->Fill(posMu.eta(),theWeight);
        h_etaMu->Fill(mup->eta(),theWeight);
        isPosMu = true;

    }

    for (pat::MuonCollection::const_iterator mum = muons->begin();mum != muons->end(); ++mum){
        if( !(mum->charge() < 0)) continue;
        if( !(mum->pt() > 20.0 )) continue;
        if( !(fabs(mum->eta()) < 2.4 )) continue;
        if( !(mum->isPFMuon())) continue;
        if( !(mum->isGlobalMuon() || mum->isTrackerMuon())) continue;
        if( mum->isIsolationValid()){
            reco::MuonPFIsolation pfR04 = mum->pfIsolationR04();
            coriso2 = pfR04.sumChargedHadronPt + std::max(0., pfR04.sumNeutralHadronEt+pfR04.sumPhotonEt-0.5*pfR04.sumPUPt);
        }
        if (!(coriso2/mum->pt() <  0.15)) continue;
        if(mum->pt() > negMu.pt() ) negMu = *mum;
        h_PtMu->Fill(mum->pt(),theWeight);
        h_ALS_etaLMuMu->Fill(negMu.eta(),theWeight);
        h_ALS_ptLMuMu->Fill(negMu.pt(),theWeight);
        h_etaMu->Fill(mum->eta(),theWeight);
        isNegMu = true;


    }

    cout << posMu.charge() << " charges " << negMu.charge() << endl;
    ////////////////////ELECTRONS//////////////////////////
    /// electron identification
    ///
    pat::Electron posEl;
    pat::Electron negEl;

    for (pat::ElectronCollection::const_iterator elp = electrons->begin();elp != electrons->end();++elp){
        if( !( elp->charge() > 0 ) ) continue;
        if( !( elp->pt() > 20 ) ) continue;
        if( !( fabs(elp->eta() ) < 2.5)) continue;
        if( !( elp->gsfTrack()->dz() < 0.04 )) continue;
        if( !(elp->passConversionVeto())) continue;
        if( !(elp->gsfTrack()->hitPattern().numberOfLostHits(HitPattern::MISSING_INNER_HITS) <= 0)) continue;
        GsfElectron::PflowIsolationVariables pfIso = elp->pfIsolationVariables();
        static double relCombIsoEA = 999.;
        float abseta = fabs(elp->superCluster()->eta());
        float eA = effectiveAreas_.getEffectiveArea(abseta);
        relCombIsoEA = (( pfIso.sumChargedHadronPt
                          + std::max( 0.0f, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - eA*rho) )
                        / elp->pt() );
        if( !(relCombIsoEA < 0.15)) continue;
        if( !(elp->userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values") > 0.5 )) continue;
        cout << elp->userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values") << "MVA Id "<< endl;
        double dR =0;
        for(pat::MuonCollection::const_iterator dumMu = muons->begin();dumMu != muons->end();++dumMu){
            if(!(dumMu->isGlobalMuon())) continue;
            dR = ROOT::Math::VectorUtil::DeltaR(dumMu->p4(),elp->p4());
            if(!(dR > 0.1)) break;
        }
        if(!(dR > 0.1)) continue;
        if( elp->pt() > posEl.pt() ) posEl = *elp;
        isPosEl = true;

    }

    for (pat::ElectronCollection::const_iterator elm = electrons->begin();elm != electrons->end();++elm){
        if( !( elm->charge() < 0 ) ) continue;
        if( !( elm->pt() > 20 ) ) continue;
        if( !( fabs(elm->eta() ) < 2.5)) continue;
        if( !( elm->gsfTrack()->dz() < 0.04 )) continue;
        if( !(elm->passConversionVeto())) continue;
        if( !(elm->gsfTrack()->hitPattern().numberOfLostHits(HitPattern::MISSING_INNER_HITS) <= 0)) continue;
        GsfElectron::PflowIsolationVariables pfIso = elm->pfIsolationVariables();
        static double relCombIsoEA = 999.;
        float abseta = fabs(elm->superCluster()->eta());
        float eA = effectiveAreas_.getEffectiveArea(abseta);
        relCombIsoEA = (( pfIso.sumChargedHadronPt
                          + std::max( 0.0f, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - eA*rho) )
                        / elm->pt() );
        if( !(relCombIsoEA < 0.15)) continue;
        if( !(elm->userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values") > 0.5 )) continue;
        cout << elm->userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values") << "MVA Id "<< endl;
        double dR =0;
        for(pat::MuonCollection::const_iterator dumMu = muons->begin();dumMu != muons->end();++dumMu){
            if(!(dumMu->isGlobalMuon())) continue;
            dR = ROOT::Math::VectorUtil::DeltaR(dumMu->p4(),elm->p4());
            if(!(dR > 0.1)) break;
        }
        if(!(dR > 0.1)) continue;
        if( elm->pt() > negEl.pt() ) negEl = *elm;
        isNegEl = true;

    }

    //After Lepton Selection
    if(isPosMu && isNegMu) isDiMuon = true;
    if(isPosEl && isNegEl) isDiElectron = true;
    if( isPosEl && isNegMu)  isElMu = true;
    if( isNegEl && isPosMu ) isMuEl = true;
    if(isDiMuon || isDiElectron || isElMu ) isDiLeptonic = true;
    if(isDiMuon)
    {
        h_ALS_etaLMuMu->Fill(posMu.eta(),theWeight);
        h_ALS_etaLMuMu->Fill(negMu.eta(),theWeight);
        h_ALS_ptLMuMu->Fill(posMu.pt(),theWeight);
        h_ALS_ptLMuMu->Fill(negMu.pt(),theWeight);

        h_ALS_etaLDiLep->Fill(posMu.eta(),theWeight);
        h_ALS_etaLDiLep->Fill(negMu.eta(),theWeight);
        h_ALS_ptLDiLep->Fill(posMu.pt(),theWeight);
        h_ALS_ptLDiLep->Fill(negMu.pt(),theWeight);

    }
    if(isDiElectron)
    {
        h_ALS_etaLElEl->Fill(posEl.eta(),theWeight);
        h_ALS_etaLElEl->Fill(negEl.eta(),theWeight);
        h_ALS_ptLElEl->Fill(posEl.pt(),theWeight);
        h_ALS_ptLElEl->Fill(negEl.pt(),theWeight);

        h_ALS_etaLDiLep->Fill(posEl.eta(),theWeight);
        h_ALS_etaLDiLep->Fill(negEl.eta(),theWeight);
        h_ALS_ptLDiLep->Fill(posEl.pt(),theWeight);
        h_ALS_ptLDiLep->Fill(negEl.pt(),theWeight);

    }
    if(isElMu)
    {
        h_ALS_etaLElMu->Fill(posEl.eta(),theWeight);
        h_ALS_etaLElMu->Fill(negMu.eta(),theWeight);
        h_ALS_ptLElMu->Fill(posEl.pt(),theWeight);
        h_ALS_ptLElMu->Fill(negMu.pt(),theWeight);

        h_ALS_etaLDiLep->Fill(posEl.eta(),theWeight);
        h_ALS_etaLDiLep->Fill(negMu.eta(),theWeight);
        h_ALS_ptLDiLep->Fill(posEl.pt(),theWeight);
        h_ALS_ptLDiLep->Fill(negMu.pt(),theWeight);

    }
    if(isMuEl)
    {
        h_ALS_etaLElMu->Fill(posMu.eta(),theWeight);
        h_ALS_etaLElMu->Fill(negEl.eta(),theWeight);
        h_ALS_ptLElMu->Fill(posMu.pt(),theWeight);
        h_ALS_ptLElMu->Fill(negEl.pt(),theWeight);

        h_ALS_etaLDiLep->Fill(posMu.eta(),theWeight);
        h_ALS_etaLDiLep->Fill(negEl.eta(),theWeight);
        h_ALS_ptLDiLep->Fill(posMu.pt(),theWeight);
        h_ALS_ptLDiLep->Fill(negEl.pt(),theWeight);

    }

    ////////////////////METS PF1//////////////////////////
    const pat::MET &met = mets->front();
    if(isDiMuon)
    {
        dileptonMass = posMu.mass() + negMu.mass();
    }
    if(isDiElectron)
    {

    }
    if(isElMu)
    {

    }
    if(isMuEl)
    {

    }

    ///////////////////JETS//////////////////////////////
    vector<pat::Jet> bjets;
    //    pat::Jet j;j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

    for(pat::JetCollection::const_iterator jet_it=jets->begin(); jet_it != jets->end();++jet_it){
        if( !( jet_it->pt() > 30 )) continue;
        if( !( fabs(jet_it->eta()) < 2.4 )) continue;
        if( !(jet_it->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") >  0.460 )) continue;          //loose working point for btaging
        if( !(jet_it->neutralHadronEnergyFraction() < 0.99 && jet_it->neutralEmEnergyFraction() < 0.99 && (jet_it->chargedMultiplicity() + jet_it->neutralMultiplicity())> 1.
              && jet_it->chargedHadronEnergyFraction() > 0. && jet_it->chargedEmEnergyFraction() < 0.99 && jet_it->chargedMultiplicity() > 1.)) continue;
        bjets.push_back(*jet_it);

    }

    cout << "number of bjets"<< bjets.size() << endl;

    cout << posMu.pt() << " Muons " << negMu.pt() << endl;
    if(bjets.size() >= 2) cout << bjets.at(0).p4() << " Jets "<< bjets.at(1).p4() << endl;
    cout << met.pt() << "met"<< endl;

    //////////////////////GEN LEVEL////////////////

    if(!isData){
        TLorentzVector genPosLep;
        TLorentzVector genNegLep;
        //    TLorentzVector genT;
        //    TLorentzVector genTbar;
        TLorentzVector genNu;
        TLorentzVector genNubar;
        TLorentzVector genB;
        TLorentzVector genBbar;
        //    TLorentzVector genWplus;
        //    TLorentzVector genWminus;
        if ( genEvent->isFullLeptonic(true) ){
            genPosLep.SetPtEtaPhiM(genEvent->lepton(true)->pt(),genEvent->lepton(true)->eta(),genEvent->lepton(true)->phi(),genEvent->lepton(true)->mass());
            genNegLep.SetPtEtaPhiM(genEvent->leptonBar(true)->pt(),genEvent->leptonBar(true)->eta(),genEvent->leptonBar(true)->phi(),genEvent->leptonBar(true)->mass());
            //        genT.SetPtEtaPhiM(genEvent->top(true)->pt(),genEvent->top()->eta(),genEvent->top()->phi(),genEvent->top()->mass());
            //        genTbar.SetPtEtaPhiM(genEvent->topBar()->pt(),genEvent->topBar()->eta(),genEvent->topBar()->phi(),genEvent->topBar()->mass());
            genNu.SetPtEtaPhiM(genEvent->neutrino(true)->pt(),genEvent->neutrino(true)->eta(),genEvent->neutrino(true)->phi(),genEvent->neutrino(true)->mass());
            genNubar.SetPtEtaPhiM(genEvent->neutrinoBar(true)->pt(),genEvent->neutrinoBar(true)->eta(),genEvent->neutrinoBar(true)->phi(),genEvent->neutrinoBar(true)->mass());
            genB.SetPtEtaPhiM(genEvent->b()->pt(),genEvent->b()->eta(),genEvent->b()->phi(),genEvent->b()->mass());
            genBbar.SetPtEtaPhiM(genEvent->bBar()->pt(),genEvent->bBar()->eta(),genEvent->bBar()->phi(),genEvent->bBar()->mass());
            //        genWplus.SetPtEtaPhiM(genEvent->wPlus()->pt(),genEvent->wPlus()->eta(),genEvent->wPlus()->phi(),genEvent->wPlus()->mass());
            //        genWminus.SetPtEtaPhiM(genEvent->wMinus()->pt(),genEvent->wMinus()->eta(),genEvent->wMinus()->phi(),genEvent->wMinus()->mass());

            amwtSolver->SetConstraints(genNu.Px()+genNubar.Px(),genNu.Py()+genNubar.Py());
            dilepton = genPosLep + genNegLep;
            if( genPosLep.Pt() > 20 && genNegLep.Pt() > 20 && genB.Pt() > 30 && genBbar.Pt() > 30 && genNu.Pt() > 0.1 && genNubar.Pt() > 0.1 && !(dilepton.M() < 106 || dilepton.M() < 76) ){
                TtAMWTSolver::NeutrinoSolution nuSol =  amwtSolver->NuSolver(genPosLep,genNegLep,genB,genBbar);
                cout << nuSol.neutrino.p4() << " neutrino "<< endl;
                if(nuSol.neutrino.pt()> 0 && nuSol.neutrinoBar.pt() > 0){
                    TLorentzVector nu;
                    TLorentzVector nuBar;
                    nu.SetPtEtaPhiM(nuSol.neutrino.pt(),nuSol.neutrino.eta(),nuSol.neutrino.phi(),nuSol.neutrino.mass());
                    nuBar.SetPtEtaPhiM(nuSol.neutrinoBar.pt(),nuSol.neutrinoBar.eta(),nuSol.neutrinoBar.phi(),nuSol.neutrinoBar.mass());
                    W1 = nu + genPosLep;
                    W2 = nuBar + genNegLep;


                    t1 = W1 + genB;
                    t2 = W2 + genBbar;
                    ttbar = t1 + t2;
                    h_GenTTbarM->Fill(ttbar.M(),theWeight);





                    genPosLep.Boost(-W1.BoostVector());
                    W1.Boost(-t1.BoostVector());
                    float theta1 = (W1.Angle(genPosLep.Vect()));
                    h_cosGen->Fill(TMath::Cos(theta1),theWeight);
                    genNegLep.Boost(-W2.BoostVector());
                    W2.Boost(-t2.BoostVector());
                    float theta2 = (W2.Angle(genNegLep.Vect()));
                    h_cosGen->Fill(TMath::Cos(theta2),theWeight);


                }

            }
        }
        ////gen mu mu
        if ( genEvent->numberOfLeptons(WDecay::kMuon,true) >= 2  ){
            genPosLep.SetPtEtaPhiM(genEvent->muPlus()->pt(),genEvent->muPlus()->eta(),genEvent->muPlus()->phi(),genEvent->muPlus()->mass());
            genNegLep.SetPtEtaPhiM(genEvent->muMinus()->pt(),genEvent->muMinus()->eta(),genEvent->muMinus()->phi(),genEvent->muMinus()->mass());
            //        genT.SetPtEtaPhiM(genEvent->top(true)->pt(),genEvent->top()->eta(),genEvent->top()->phi(),genEvent->top()->mass());
            //        genTbar.SetPtEtaPhiM(genEvent->topBar()->pt(),genEvent->topBar()->eta(),genEvent->topBar()->phi(),genEvent->topBar()->mass());
            genNu.SetPtEtaPhiM(genEvent->neutrino(true)->pt(),genEvent->neutrino(true)->eta(),genEvent->neutrino(true)->phi(),genEvent->neutrino(true)->mass());
            genNubar.SetPtEtaPhiM(genEvent->neutrinoBar(true)->pt(),genEvent->neutrinoBar(true)->eta(),genEvent->neutrinoBar(true)->phi(),genEvent->neutrinoBar(true)->mass());
            genB.SetPtEtaPhiM(genEvent->b()->pt(),genEvent->b()->eta(),genEvent->b()->phi(),genEvent->b()->mass());
            genBbar.SetPtEtaPhiM(genEvent->bBar()->pt(),genEvent->bBar()->eta(),genEvent->bBar()->phi(),genEvent->bBar()->mass());
            //        genWplus.SetPtEtaPhiM(genEvent->wPlus()->pt(),genEvent->wPlus()->eta(),genEvent->wPlus()->phi(),genEvent->wPlus()->mass());
            //        genWminus.SetPtEtaPhiM(genEvent->wMinus()->pt(),genEvent->wMinus()->eta(),genEvent->wMinus()->phi(),genEvent->wMinus()->mass());
            amwtSolver->SetConstraints(genNu.Px()+genNubar.Px(),genNu.Py()+genNubar.Py());
            dilepton = genPosLep + genNegLep;
            if( genPosLep.Pt() > 20 && genNegLep.Pt() > 20 && genB.Pt() > 30 && genBbar.Pt() > 30 && genNu.Pt() > 0.1 && genNubar.Pt() > 0.1 && !(dilepton.M() < 106 || dilepton.M() < 76) ){
                TtAMWTSolver::NeutrinoSolution nuSol =  amwtSolver->NuSolver(genPosLep,genNegLep,genB,genBbar);
                cout << nuSol.neutrino.p4() << " neutrino "<< endl;
                if(nuSol.neutrino.pt()> 0 && nuSol.neutrinoBar.pt() > 0){
                    TLorentzVector nu;
                    TLorentzVector nuBar;
                    nu.SetPtEtaPhiM(nuSol.neutrino.pt(),nuSol.neutrino.eta(),nuSol.neutrino.phi(),nuSol.neutrino.mass());
                    nuBar.SetPtEtaPhiM(nuSol.neutrinoBar.pt(),nuSol.neutrinoBar.eta(),nuSol.neutrinoBar.phi(),nuSol.neutrinoBar.mass());
                    W1 = nu + genPosLep;
                    W2 = nuBar + genNegLep;


                    t1 = W1 + genB;
                    t2 = W2 + genBbar;
                    ttbar = t1 + t2;
                    //                h_GenTTbarM->Fill(ttbar.M(),theWeight);





                    genPosLep.Boost(-W1.BoostVector());
                    W1.Boost(-t1.BoostVector());
                    float theta1 = (W1.Angle(genPosLep.Vect()));
                    h_cosGenMuMu->Fill(TMath::Cos(theta1),theWeight);
                    genNegLep.Boost(-W2.BoostVector());
                    W2.Boost(-t2.BoostVector());
                    float theta2 = (W2.Angle(genNegLep.Vect()));
                    h_cosGenMuMu->Fill(TMath::Cos(theta2),theWeight);


                }

            }
        }
        ////gen electron electron
        if ( genEvent->numberOfLeptons(WDecay::kElec,true) >= 2 ){
            genPosLep.SetPtEtaPhiM(genEvent->ePlus()->pt(),genEvent->ePlus()->eta(),genEvent->ePlus()->phi(),genEvent->ePlus()->mass());
            genNegLep.SetPtEtaPhiM(genEvent->eMinus()->pt(),genEvent->eMinus()->eta(),genEvent->eMinus()->phi(),genEvent->eMinus()->mass());
            //        genT.SetPtEtaPhiM(genEvent->top(true)->pt(),genEvent->top()->eta(),genEvent->top()->phi(),genEvent->top()->mass());
            //        genTbar.SetPtEtaPhiM(genEvent->topBar()->pt(),genEvent->topBar()->eta(),genEvent->topBar()->phi(),genEvent->topBar()->mass());
            genNu.SetPtEtaPhiM(genEvent->neutrino(true)->pt(),genEvent->neutrino(true)->eta(),genEvent->neutrino(true)->phi(),genEvent->neutrino(true)->mass());
            genNubar.SetPtEtaPhiM(genEvent->neutrinoBar(true)->pt(),genEvent->neutrinoBar(true)->eta(),genEvent->neutrinoBar(true)->phi(),genEvent->neutrinoBar(true)->mass());
            genB.SetPtEtaPhiM(genEvent->b()->pt(),genEvent->b()->eta(),genEvent->b()->phi(),genEvent->b()->mass());
            genBbar.SetPtEtaPhiM(genEvent->bBar()->pt(),genEvent->bBar()->eta(),genEvent->bBar()->phi(),genEvent->bBar()->mass());
            //        genWplus.SetPtEtaPhiM(genEvent->wPlus()->pt(),genEvent->wPlus()->eta(),genEvent->wPlus()->phi(),genEvent->wPlus()->mass());
            //        genWminus.SetPtEtaPhiM(genEvent->wMinus()->pt(),genEvent->wMinus()->eta(),genEvent->wMinus()->phi(),genEvent->wMinus()->mass());
            amwtSolver->SetConstraints(genNu.Px()+genNubar.Px(),genNu.Py()+genNubar.Py());
            dilepton = genPosLep + genNegLep;
            if( genPosLep.Pt() > 20 && genNegLep.Pt() > 20 && genB.Pt() > 30 && genBbar.Pt() > 30 && genNu.Pt() > 0.1 && genNubar.Pt() > 0.1 && !(dilepton.M() < 106 || dilepton.M() < 76) ){
                TtAMWTSolver::NeutrinoSolution nuSol =  amwtSolver->NuSolver(genPosLep,genNegLep,genB,genBbar);
                cout << nuSol.neutrino.p4() << " neutrino "<< endl;
                if(nuSol.neutrino.pt()> 0 && nuSol.neutrinoBar.pt() > 0){
                    TLorentzVector nu;
                    TLorentzVector nuBar;
                    nu.SetPtEtaPhiM(nuSol.neutrino.pt(),nuSol.neutrino.eta(),nuSol.neutrino.phi(),nuSol.neutrino.mass());
                    nuBar.SetPtEtaPhiM(nuSol.neutrinoBar.pt(),nuSol.neutrinoBar.eta(),nuSol.neutrinoBar.phi(),nuSol.neutrinoBar.mass());
                    W1 = nu + genPosLep;
                    W2 = nuBar + genNegLep;


                    t1 = W1 + genB;
                    t2 = W2 + genBbar;
                    ttbar = t1 + t2;
                    //                h_GenTTbarM->Fill(ttbar.M(),theWeight);





                    genPosLep.Boost(-W1.BoostVector());
                    W1.Boost(-t1.BoostVector());
                    float theta1 = (W1.Angle(genPosLep.Vect()));
                    h_cosGenElEl->Fill(TMath::Cos(theta1),theWeight);
                    genNegLep.Boost(-W2.BoostVector());
                    W2.Boost(-t2.BoostVector());
                    float theta2 = (W2.Angle(genNegLep.Vect()));
                    h_cosGenElEl->Fill(TMath::Cos(theta2),theWeight);


                }

            }
        }
        ///gen electron muon
        if ( genEvent->numberOfLeptons(WDecay::kElec,true) >= 1 && genEvent->numberOfLeptons(WDecay::kMuon,true) >= 1){
            if( genEvent->ePlus() != nullptr && genEvent->muMinus() != nullptr){
                genPosLep.SetPtEtaPhiM(genEvent->ePlus()->pt(),genEvent->ePlus()->eta(),genEvent->ePlus()->phi(),genEvent->ePlus()->mass());
                genNegLep.SetPtEtaPhiM(genEvent->muMinus()->pt(),genEvent->muMinus()->eta(),genEvent->muMinus()->phi(),genEvent->muMinus()->mass());
                //        genT.SetPtEtaPhiM(genEvent->top(true)->pt(),genEvent->top()->eta(),genEvent->top()->phi(),genEvent->top()->mass());
                //        genTbar.SetPtEtaPhiM(genEvent->topBar()->pt(),genEvent->topBar()->eta(),genEvent->topBar()->phi(),genEvent->topBar()->mass());
                genNu.SetPtEtaPhiM(genEvent->neutrino(true)->pt(),genEvent->neutrino(true)->eta(),genEvent->neutrino(true)->phi(),genEvent->neutrino(true)->mass());
                genNubar.SetPtEtaPhiM(genEvent->neutrinoBar(true)->pt(),genEvent->neutrinoBar(true)->eta(),genEvent->neutrinoBar(true)->phi(),genEvent->neutrinoBar(true)->mass());
                genB.SetPtEtaPhiM(genEvent->b()->pt(),genEvent->b()->eta(),genEvent->b()->phi(),genEvent->b()->mass());
                genBbar.SetPtEtaPhiM(genEvent->bBar()->pt(),genEvent->bBar()->eta(),genEvent->bBar()->phi(),genEvent->bBar()->mass());
                //        genWplus.SetPtEtaPhiM(genEvent->wPlus()->pt(),genEvent->wPlus()->eta(),genEvent->wPlus()->phi(),genEvent->wPlus()->mass());
                //        genWminus.SetPtEtaPhiM(genEvent->wMinus()->pt(),genEvent->wMinus()->eta(),genEvent->wMinus()->phi(),genEvent->wMinus()->mass());
                amwtSolver->SetConstraints(genNu.Px()+genNubar.Px(),genNu.Py()+genNubar.Py());
                dilepton = genPosLep + genNegLep;
                if( genPosLep.Pt() > 20 && genNegLep.Pt() > 20 && genB.Pt() > 30 && genBbar.Pt() > 30 && genNu.Pt() > 0.1 && genNubar.Pt() > 0.1 && !(dilepton.M() < 106 || dilepton.M() < 76) ){
                    TtAMWTSolver::NeutrinoSolution nuSol =  amwtSolver->NuSolver(genPosLep,genNegLep,genB,genBbar);
                    cout << nuSol.neutrino.p4() << " neutrino "<< endl;
                    if(nuSol.neutrino.pt()> 0 && nuSol.neutrinoBar.pt() > 0){
                        TLorentzVector nu;
                        TLorentzVector nuBar;
                        nu.SetPtEtaPhiM(nuSol.neutrino.pt(),nuSol.neutrino.eta(),nuSol.neutrino.phi(),nuSol.neutrino.mass());
                        nuBar.SetPtEtaPhiM(nuSol.neutrinoBar.pt(),nuSol.neutrinoBar.eta(),nuSol.neutrinoBar.phi(),nuSol.neutrinoBar.mass());
                        W1 = nu + genPosLep;
                        W2 = nuBar + genNegLep;


                        t1 = W1 + genB;
                        t2 = W2 + genBbar;
                        ttbar = t1 + t2;
                        //                h_GenTTbarM->Fill(ttbar.M(),theWeight);





                        genPosLep.Boost(-W1.BoostVector());
                        W1.Boost(-t1.BoostVector());
                        float theta1 = (W1.Angle(genPosLep.Vect()));
                        h_cosGenElMu->Fill(TMath::Cos(theta1),theWeight);
                        genNegLep.Boost(-W2.BoostVector());
                        W2.Boost(-t2.BoostVector());
                        float theta2 = (W2.Angle(genNegLep.Vect()));
                        h_cosGenElMu->Fill(TMath::Cos(theta2),theWeight);


                    }

                }
            }
            else if( genEvent->eMinus() != nullptr && genEvent->muPlus() != nullptr)
            {
                genPosLep.SetPtEtaPhiM(genEvent->muPlus()->pt(),genEvent->muPlus()->eta(),genEvent->muPlus()->phi(),genEvent->muPlus()->mass());
                genNegLep.SetPtEtaPhiM(genEvent->eMinus()->pt(),genEvent->eMinus()->eta(),genEvent->eMinus()->phi(),genEvent->eMinus()->mass());
                //        genT.SetPtEtaPhiM(genEvent->top(true)->pt(),genEvent->top()->eta(),genEvent->top()->phi(),genEvent->top()->mass());
                //        genTbar.SetPtEtaPhiM(genEvent->topBar()->pt(),genEvent->topBar()->eta(),genEvent->topBar()->phi(),genEvent->topBar()->mass());
                genNu.SetPtEtaPhiM(genEvent->neutrino(true)->pt(),genEvent->neutrino(true)->eta(),genEvent->neutrino(true)->phi(),genEvent->neutrino(true)->mass());
                genNubar.SetPtEtaPhiM(genEvent->neutrinoBar(true)->pt(),genEvent->neutrinoBar(true)->eta(),genEvent->neutrinoBar(true)->phi(),genEvent->neutrinoBar(true)->mass());
                genB.SetPtEtaPhiM(genEvent->b()->pt(),genEvent->b()->eta(),genEvent->b()->phi(),genEvent->b()->mass());
                genBbar.SetPtEtaPhiM(genEvent->bBar()->pt(),genEvent->bBar()->eta(),genEvent->bBar()->phi(),genEvent->bBar()->mass());
                //        genWplus.SetPtEtaPhiM(genEvent->wPlus()->pt(),genEvent->wPlus()->eta(),genEvent->wPlus()->phi(),genEvent->wPlus()->mass());
                //        genWminus.SetPtEtaPhiM(genEvent->wMinus()->pt(),genEvent->wMinus()->eta(),genEvent->wMinus()->phi(),genEvent->wMinus()->mass());
                amwtSolver->SetConstraints(genNu.Px()+genNubar.Px(),genNu.Py()+genNubar.Py());
                dilepton = genPosLep + genNegLep;
                if( genPosLep.Pt() > 20 && genNegLep.Pt() > 20 && genB.Pt() > 30 && genBbar.Pt() > 30 && genNu.Pt() > 0.1 && genNubar.Pt() > 0.1 && !(dilepton.M() < 106 || dilepton.M() < 76) ){
                    TtAMWTSolver::NeutrinoSolution nuSol =  amwtSolver->NuSolver(genPosLep,genNegLep,genB,genBbar);
                    cout << nuSol.neutrino.p4() << " neutrino "<< endl;
                    if(nuSol.neutrino.pt()> 0 && nuSol.neutrinoBar.pt() > 0){
                        TLorentzVector nu;
                        TLorentzVector nuBar;
                        nu.SetPtEtaPhiM(nuSol.neutrino.pt(),nuSol.neutrino.eta(),nuSol.neutrino.phi(),nuSol.neutrino.mass());
                        nuBar.SetPtEtaPhiM(nuSol.neutrinoBar.pt(),nuSol.neutrinoBar.eta(),nuSol.neutrinoBar.phi(),nuSol.neutrinoBar.mass());
                        W1 = nu + genPosLep;
                        W2 = nuBar + genNegLep;


                        t1 = W1 + genB;
                        t2 = W2 + genBbar;
                        ttbar = t1 + t2;
                        //                h_GenTTbarM->Fill(ttbar.M(),theWeight);





                        genPosLep.Boost(-W1.BoostVector());
                        W1.Boost(-t1.BoostVector());
                        float theta1 = (W1.Angle(genPosLep.Vect()));
                        h_cosGenElMu->Fill(TMath::Cos(theta1),theWeight);
                        genNegLep.Boost(-W2.BoostVector());
                        W2.Boost(-t2.BoostVector());
                        float theta2 = (W2.Angle(genNegLep.Vect()));
                        h_cosGenElMu->Fill(TMath::Cos(theta2),theWeight);


                    }

                }
            }

        }
    }
    //////////////////////TOP RECO/////////////////
    //    TtDilepEvtSolution asol;

    //    asol.setGenEvt(genEvent);
    ////////DiMuon:
    if(isDiMuon && bjets.size() >= 2 && met.pt() > 0
            && !( ROOT::Math::VectorUtil::DeltaR(bjets.at(0).p4(),posMu.p4()) < 0.5 || ROOT::Math::VectorUtil::DeltaR(bjets.at(1).p4(),negMu.p4()) < 0.5 || ROOT::Math::VectorUtil::DeltaR(bjets.at(0).p4(),negMu.p4()) < 0.5 ||  ROOT::Math::VectorUtil::DeltaR(bjets.at(1).p4(),posMu.p4()) < 0.5) ){
        cout << "IM HERE!" << endl;
        lepPos.SetPtEtaPhiE(posMu.pt(),posMu.eta(),posMu.phi(),posMu.energy());
        lepNeg.SetPtEtaPhiE(negMu.pt(),negMu.eta(),negMu.phi(),negMu.energy());
        BJet.SetPtEtaPhiE(bjets.at(0).pt(),bjets.at(0).eta(),bjets.at(0).phi(),bjets.at(0).energy());
        BBJet.SetPtEtaPhiE(bjets.at(1).pt(),bjets.at(1).eta(),bjets.at(1).phi(),bjets.at(1).energy());


        amwtSolver->SetConstraints(met.px(),met.py());
        //        solver->addKinSolInfo(&asol);
        //        solver->useWeightFromMC(true);

        TtAMWTSolver::NeutrinoSolution nuSol =  amwtSolver->NuSolver(lepPos,lepNeg,BJet,BBJet);
        cout << nuSol.neutrino.p4() << " neutrino "<< endl;
        dilepton = lepPos + lepNeg;
        if(nuSol.neutrino.pt()> 0 && nuSol.neutrinoBar.pt() > 0 && !(dilepton.M() < 106 || dilepton.M() < 76)){
            TLorentzVector nu;
            TLorentzVector nuBar;
            nu.SetPtEtaPhiM(nuSol.neutrino.pt(),nuSol.neutrino.eta(),nuSol.neutrino.phi(),nuSol.neutrino.mass());
            nuBar.SetPtEtaPhiM(nuSol.neutrinoBar.pt(),nuSol.neutrinoBar.eta(),nuSol.neutrinoBar.phi(),nuSol.neutrinoBar.mass());
            W1 = nu + lepPos;
            W2 = nuBar + lepNeg;


            t1 = W1 + BJet;
            t2 = W2 + BBJet;
            ttbar = t1 + t2;
            h_TTbarM->Fill(ttbar.M(),theWeight);





            lepPos.Boost(-W1.BoostVector());
            W1.Boost(-t1.BoostVector());
            float theta1 = (W1.Angle(lepPos.Vect()));
            h_cosMuMu->Fill(TMath::Cos(theta1),theWeight);
            lepNeg.Boost(-W2.BoostVector());
            W2.Boost(-t2.BoostVector());
            float theta2 = (W2.Angle(lepNeg.Vect()));
            h_cosMuMu->Fill(TMath::Cos(theta2),theWeight);

        }
    }

    //////DiElectron:
    if(isDiElectron && bjets.size() >= 2 && met.pt() > 0
            && !( ROOT::Math::VectorUtil::DeltaR(bjets.at(0).p4(),posEl.p4()) < 0.5 || ROOT::Math::VectorUtil::DeltaR(bjets.at(1).p4(),negEl.p4()) < 0.5 || ROOT::Math::VectorUtil::DeltaR(bjets.at(0).p4(),negEl.p4()) < 0.5 ||  ROOT::Math::VectorUtil::DeltaR(bjets.at(1).p4(),posEl.p4()) < 0.5) ){
        cout << "IM HERE!" << endl;
        lepPos.SetPtEtaPhiE(posEl.pt(),posEl.eta(),posEl.phi(),posEl.energy());
        lepNeg.SetPtEtaPhiE(negEl.pt(),negEl.eta(),negEl.phi(),negEl.energy());
        BJet.SetPtEtaPhiE(bjets.at(0).pt(),bjets.at(0).eta(),bjets.at(0).phi(),bjets.at(0).energy());
        BBJet.SetPtEtaPhiE(bjets.at(1).pt(),bjets.at(1).eta(),bjets.at(1).phi(),bjets.at(1).energy());



        amwtSolver->SetConstraints(met.px(),met.py());
        TtAMWTSolver::NeutrinoSolution nuSol =  amwtSolver->NuSolver(lepPos,lepNeg,BJet,BBJet);
        cout << nuSol.neutrino.p4() << " neutrino "<< endl;
        dilepton = lepPos + lepNeg;
        if(nuSol.neutrino.pt()> 0 && nuSol.neutrinoBar.pt() > 0 && !(dilepton.M() < 106 || dilepton.M() < 76)){
            TLorentzVector nu;
            TLorentzVector nuBar;
            nu.SetPtEtaPhiM(nuSol.neutrino.pt(),nuSol.neutrino.eta(),nuSol.neutrino.phi(),nuSol.neutrino.mass());
            nuBar.SetPtEtaPhiM(nuSol.neutrinoBar.pt(),nuSol.neutrinoBar.eta(),nuSol.neutrinoBar.phi(),nuSol.neutrinoBar.mass());
            W1 = nu + lepPos;
            W2 = nuBar + lepNeg;


            t1 = W1 + BJet;
            t2 = W2 + BBJet;
            ttbar = t1 + t2;
            h_TTbarM->Fill(ttbar.M(),theWeight);





            lepPos.Boost(-W1.BoostVector());
            W1.Boost(-t1.BoostVector());
            float theta1 = (W1.Angle(lepPos.Vect()));
            h_cosElEl->Fill(TMath::Cos(theta1),theWeight);
            lepNeg.Boost(-W2.BoostVector());
            W2.Boost(-t2.BoostVector());
            float theta2 = (W2.Angle(lepNeg.Vect()));
            h_cosElEl->Fill(TMath::Cos(theta2),theWeight);

        }
    }
    //////ElectronMuon:
    if(isElMu && bjets.size() >= 2 && met.pt() > 0
            && !( ROOT::Math::VectorUtil::DeltaR(bjets.at(0).p4(),posEl.p4()) < 0.5 || ROOT::Math::VectorUtil::DeltaR(bjets.at(1).p4(),negMu.p4()) < 0.5 || ROOT::Math::VectorUtil::DeltaR(bjets.at(0).p4(),negMu.p4()) < 0.5 ||  ROOT::Math::VectorUtil::DeltaR(bjets.at(1).p4(),posEl.p4()) < 0.5) ){
        cout << "IM HERE!" << endl;
        lepPos.SetPtEtaPhiE(posEl.pt(),posEl.eta(),posEl.phi(),posEl.energy());
        lepNeg.SetPtEtaPhiE(negMu.pt(),negMu.eta(),negMu.phi(),negMu.energy());
        BJet.SetPtEtaPhiE(bjets.at(0).pt(),bjets.at(0).eta(),bjets.at(0).phi(),bjets.at(0).energy());
        BBJet.SetPtEtaPhiE(bjets.at(1).pt(),bjets.at(1).eta(),bjets.at(1).phi(),bjets.at(1).energy());



        amwtSolver->SetConstraints(met.px(),met.py());
        TtAMWTSolver::NeutrinoSolution nuSol =  amwtSolver->NuSolver(lepPos,lepNeg,BJet,BBJet);
        cout << nuSol.neutrino.p4() << " neutrino "<< endl;
        dilepton = lepPos + lepNeg;
        if(nuSol.neutrino.pt()> 0 && nuSol.neutrinoBar.pt() > 0 && !(dilepton.M() < 106 || dilepton.M() < 76)){
            TLorentzVector nu;
            TLorentzVector nuBar;
            nu.SetPtEtaPhiM(nuSol.neutrino.pt(),nuSol.neutrino.eta(),nuSol.neutrino.phi(),nuSol.neutrino.mass());
            nuBar.SetPtEtaPhiM(nuSol.neutrinoBar.pt(),nuSol.neutrinoBar.eta(),nuSol.neutrinoBar.phi(),nuSol.neutrinoBar.mass());
            W1 = nu + lepPos;
            W2 = nuBar + lepNeg;


            t1 = W1 + BJet;
            t2 = W2 + BBJet;
            ttbar = t1 + t2;
            h_TTbarM->Fill(ttbar.M(),theWeight);





            lepPos.Boost(-W1.BoostVector());
            W1.Boost(-t1.BoostVector());
            float theta1 = (W1.Angle(lepPos.Vect()));
            h_cosElMu->Fill(TMath::Cos(theta1),theWeight);
            lepNeg.Boost(-W2.BoostVector());
            W2.Boost(-t2.BoostVector());
            float theta2 = (W2.Angle(lepNeg.Vect()));
            h_cosElMu->Fill(TMath::Cos(theta2),theWeight);

        }
    }

    if(isMuEl && bjets.size() >= 2 && met.pt() > 0
            && !( ROOT::Math::VectorUtil::DeltaR(bjets.at(0).p4(),posMu.p4()) < 0.5 || ROOT::Math::VectorUtil::DeltaR(bjets.at(1).p4(),negEl.p4()) < 0.5 || ROOT::Math::VectorUtil::DeltaR(bjets.at(0).p4(),negEl.p4()) < 0.5 ||  ROOT::Math::VectorUtil::DeltaR(bjets.at(1).p4(),posMu.p4()) < 0.5) ){
        cout << "IM HERE!" << endl;
        lepPos.SetPtEtaPhiE(posMu.pt(),posMu.eta(),posMu.phi(),posMu.energy());
        lepNeg.SetPtEtaPhiE(negEl.pt(),negEl.eta(),negEl.phi(),negEl.energy());
        BJet.SetPtEtaPhiE(bjets.at(0).pt(),bjets.at(0).eta(),bjets.at(0).phi(),bjets.at(0).energy());
        BBJet.SetPtEtaPhiE(bjets.at(1).pt(),bjets.at(1).eta(),bjets.at(1).phi(),bjets.at(1).energy());



        amwtSolver->SetConstraints(met.px(),met.py());
        TtAMWTSolver::NeutrinoSolution nuSol =  amwtSolver->NuSolver(lepPos,lepNeg,BJet,BBJet);
        cout << nuSol.neutrino.p4() << " neutrino "<< endl;
        dilepton = lepPos + lepNeg;
        if(nuSol.neutrino.pt()> 0 && nuSol.neutrinoBar.pt() > 0 && !(dilepton.M() < 106 || dilepton.M() < 76)){
            TLorentzVector nu;
            TLorentzVector nuBar;
            nu.SetPtEtaPhiM(nuSol.neutrino.pt(),nuSol.neutrino.eta(),nuSol.neutrino.phi(),nuSol.neutrino.mass());
            nuBar.SetPtEtaPhiM(nuSol.neutrinoBar.pt(),nuSol.neutrinoBar.eta(),nuSol.neutrinoBar.phi(),nuSol.neutrinoBar.mass());
            W1 = nu + lepPos;
            W2 = nuBar + lepNeg;


            t1 = W1 + BJet;
            t2 = W2 + BBJet;
            ttbar = t1 + t2;
            h_TTbarM->Fill(ttbar.M(),theWeight);





            lepPos.Boost(-W1.BoostVector());
            W1.Boost(-t1.BoostVector());
            float theta1 = (W1.Angle(lepPos.Vect()));
            h_cosElMu->Fill(TMath::Cos(theta1),theWeight);
            lepNeg.Boost(-W2.BoostVector());
            W2.Boost(-t2.BoostVector());
            float theta2 = (W2.Angle(lepNeg.Vect()));
            h_cosElMu->Fill(TMath::Cos(theta2),theWeight);

        }
    }





    /////////////////////Print Event Details///////////////////////////////

    //    for (const pat::Muon &mu : *muons) {
    //        if (mu.pt() < 5 || !mu.isLooseMuon()) continue;
    //        printf("muon with pt %4.1f,muon charge %i, dz(PV) %+5.3f, POG loose id %d, tight id %d\n",
    //               mu.pt(),mu.charge(), mu.muonBestTrack()->dz(PV->position()), mu.isLooseMuon(), mu.isTightMuon(*PV));
    //    }


    //    for (const pat::Electron &el : *electrons) {
    //        if (el.pt() < 5) continue;
    //        printf("elec with pt %4.1f, number of lost Hits %i, supercluster eta %+5.3f, sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes),  pass conv veto %d\n",
    //               el.pt(),el.gsfTrack()->numberOfLostHits(), el.superCluster()->eta(), el.sigmaIetaIeta(), el.full5x5_sigmaIetaIeta(),  el.passConversionVeto());
    //    }

    //    edm::Handle<pat::PhotonCollection> photons;
    //    iEvent.getByToken(photonToken_, photons);
    //    for (const pat::Photon &pho : *photons) {
    //        if (pho.pt() < 20 or pho.chargedHadronIso()/pho.pt() > 0.3) continue;
    //        printf("phot with pt %4.1f, supercluster eta %+5.3f, sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes)\n",
    //               pho.pt(), pho.superCluster()->eta(), pho.sigmaIetaIeta(), pho.full5x5_sigmaIetaIeta());
    //    }

    // edm::Handle<pat::TauCollection> taus;
    // iEvent.getByToken(tauToken_, taus);
    // for (const pat::TauCollection &tau : *taus) {
    //     if (tau.pt() < 20) continue;
    //     printf("tau  with pt %4.1f, dxy signif %.1f, ID(byMediumCombinedIsolationDeltaBetaCorr3Hits) %.1f, lead candidate pt %.1f, pdgId %d \n",
    //                 tau.pt(), tau.dxy_Sig(), tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"), tau.leadCand()->pt(), tau.leadCand()->pdgId());
    // }


    //int ijet = 0;
    //    for (const pat::Jet &j : *jets) {
    //        if (j.pt() < 20) continue;
    //        printf("jet  with pt %5.1f (raw pt %5.1f), eta %+4.2f, btag CSV %.3f, CISV %.3f, pileup mva disc %+.2f\n",
    //               j.pt(), j.pt()*j.jecFactor("Uncorrected"), j.eta(), std::max(0.f,j.bDiscriminator("combinedSecondaryVertexBJetTags")), j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"), j.userFloat("pileupJetId:fullDiscriminant"));
    //        // if ((++ijet) == 1) { // for the first jet, let's print the leading constituents
    //        //     std::vector daus(j.daughterPtrVector());
    //        //     std::sort(daus.begin(), daus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); }); // the joys of C++11
    //        //     for (unsigned int i2 = 0, n = daus.size(); i2 < n && i2 <= 3; ++i2) {
    //        //         const pat::PackedCandidate &cand = dynamic_cast<const pat::PackedCandidate &>(*daus[i2]);
    //        //         printf("         constituent %3d: pt %6.2f, dz(pv) %+.3f, pdgId %+3d\n", i2,cand.pt(),cand.dz(PV.position()),cand.pdgId());
    //        //     }
    //        // }
    //    }

    //    edm::Handle<pat::JetCollection> fatjets;
    //    iEvent.getByToken(fatjetToken_, fatjets);
    //    for (const pat::Jet &j : *fatjets) {
    //        printf("AK8j with pt %5.1f (raw pt %5.1f), eta %+4.2f, mass %5.1f ungroomed, %5.1f softdrop, %5.1f pruned, %5.1f trimmed, %5.1f filtered. \n",
    //               j.pt(), j.pt()*j.jecFactor("Uncorrected"), j.eta(), j.mass(), j.userFloat("ak8PFJetsCHSSoftDropMass"), j.userFloat("ak8PFJetsCHSPrunedMass"), j.userFloat("ak8PFJetsCHSTrimmedMass"), j.userFloat("ak8PFJetsCHSFilteredMass"));//, j.userFloat("cmsTopTagPFJetsCHSMassAK8"));

    //        // To get the constituents of the AK8 jets, you have to loop over all of the
    //        // daughters recursively. To save space, the first two constituents are actually
    //        // the Soft Drop SUBJETS, which will then point to their daughters.
    //        // The remaining constituents are those constituents removed by soft drop but
    //        // still in the AK8 jet.
    //        // std::vector<reco::GenParticle> constituents;
    //        //      for ( unsigned ida = 0; ida < j.numberOfDaughters(); ++ida ) {
    //        //   reco::Candidate const * cand = j.daughter(ida);
    //        //   if ( cand->numberOfDaughters() == 0 )
    //        //     constituents.push_back( cand ) ;
    //        //   else {
    //        //     for ( unsigned jda = 0; jda < cand->numberOfDaughters(); ++jda ) {
    //        //       reco::Candidate const * cand2 = cand->daughter(jda);
    //        //       constituents.push_back( cand2 );
    //        //     }
    //        //   }
    //        // }
    //        // std::sort( constituents.begin(), constituents.end(), [] (reco::Candidate const * ida, reco::Candidate const * jda){return ida->pt() > jda->pt();} );

    //        // for ( unsigned int ida = 0; ida < constituents.size(); ++ida ) {
    //        //   const pat::PackedCandidate &cand = dynamic_cast<const pat::PackedCandidate &>(*constituents[ida]);
    //        //   printf("         constituent %3d: pt %6.2f, dz(pv) %+.3f, pdgId %+3d\n", ida,cand.pt(),cand.dz(PV.position()),cand.pdgId());
    //        // }

    //        // auto wSubjets = j.subjets("SoftDrop");
    //        // for ( auto const & iw : wSubjets ) {
    //        //   printf("   w subjet with pt %5.1f (raw pt %5.1f), eta %+4.2f, mass %5.1f ungroomed\n",
    //        //     iw->pt(), iw->pt()*iw->jecFactor("Uncorrected"), iw->eta(), iw->mass() );

    //        // }
    //        // auto tSubjets = j.subjets("CMSTopTag");
    //        // for ( auto const & it : tSubjets ) {
    //        //   printf("   t subjet with pt %5.1f (raw pt %5.1f), eta %+4.2f, mass %5.1f ungroomed\n",
    //        //     it->pt(), it->pt()*it->jecFactor("Uncorrected"), it->eta(), it->mass() );

    //        // }
    //    }



    //    printf("MET: pt %5.1f, phi %+4.2f, sumEt (%.1f). genMET %.1f. MET with JES up/down: %.1f/%.1f\n",
    //           met.pt(), met.phi(), met.sumEt(),
    //           met.genMET()->pt(),
    //           met.shiftedPt(pat::MET::JetEnUp), met.shiftedPt(pat::MET::JetEnDown));

    printf("\n");


}

void MiniAnalyzer::beginJob() {
}
void MiniAnalyzer::endJob() {
}


//____________________________________________________________________________
bool MiniAnalyzer::IsSoftMuon(const pat::Muon& mu , const reco::Vertex& vertex) {
    if (!(muon::isGoodMuon(mu, muon::TMOneStationTight))) return false;
    if (!(mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5)) return false;
    if (!(mu.innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0)) return false;
    if (!(mu.innerTrack()->quality(reco::TrackBase::highPurity))) return false;
    if (!((fabs(mu.innerTrack()->dxy(vertex.position())) < 0.3) && (fabs(mu.innerTrack()->dz(vertex.position())) < 20.))) return false;
    return true;
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
//____________________________________________________________________________
void MiniAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(MiniAnalyzer);
