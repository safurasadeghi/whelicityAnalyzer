#ifndef TTAMWTSOLVER_H
#define TTAMWTSOLVER_H
#include "DataFormats/Candidate/interface/LeafCandidate.h"
//#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include <memory>
#include <vector>
#include <sstream>
#include <iostream>
#include <TRandom3.h>
#include <TH1.h>

#include "TLorentzVector.h"
#include "TMath.h"

#include <LHAPDF/LHAPDF.h>

using namespace std;

class TtAMWTSolver
{

public:

    struct NeutrinoSolution {
        double weight;
        reco::LeafCandidate neutrino;
        reco::LeafCandidate neutrinoBar;
    };

    void SetConstraints(const double xx=0, const double yy=0);
    void SetRho(const double rho_){rho=rho_;}
    TtAMWTSolver(bool isData, const double b, const double e, const double s, const double mW, const double mB);
    ~TtAMWTSolver();
    NeutrinoSolution NuSolver(const TLorentzVector &LV_l, const TLorentzVector &LV_l_, const TLorentzVector &LV_b, const TLorentzVector &LV_b_);
private:

    int nbrJetSmear;
    double get_weight(const TLorentzVector &bquark1, const TLorentzVector &bquark2, const TLorentzVector &lep_p, const TLorentzVector &lep_m, const TLorentzVector &nu1, const TLorentzVector &nu2, const double top_mass) const;
    double get_dalitz_prob(const TLorentzVector &lep, const TLorentzVector &top, double mb, double mw) const;
    void smear_JetMET(const vector<TLorentzVector> &orig_jets, const TVector2 &orig_met, vector<TLorentzVector> &smear_jets, TVector2 &smear_met, TRandom3 *rand3, const TLorentzVector &lep_sum) const;
    void dont_smear_JetMET(const vector<TLorentzVector> &orig_jets, const TVector2 &orig_met, vector<TLorentzVector> &smear_jets, TVector2 &smear_met) const;
    ///
    void FindCoeff(const TLorentzVector& al,
                   const TLorentzVector& l,
                   const TLorentzVector& b_al,
                   const TLorentzVector& b_l,
                   const double mt, const double mat, const double px_miss, const double py_miss,
                   double* koeficienty);
    ///
    void TopRec(const TLorentzVector& al,
                const TLorentzVector& l,
                const TLorentzVector& b_al,
                const TLorentzVector& b_l, const double sol);

    ///
    int quartic(double* koeficienty, double* koreny) const;
    ///
    int cubic(const double* coeffs, double* koreny) const;
    ///
    double sqr(const double x) const {return (x*x);}
    ///
    void SWAP(double& realone, double& realtwo) const;

private:
    ///
    const double topmass_begin;
    ///
    const double topmass_end;
    ///
    const double topmass_step;
    ///
    const double mw;
    ///
    const double mb;
    ///
    double pxmiss_, pymiss_;
    const double e_com;

    double C;
    double D;
    double F;
    double pom;
    double k16;
    double k26;
    double k36;
    double k46;
    double k56;
    double k51;
    double k61;
    double m1;
    double m2;
    double m3;
    double n1;
    double n2;
    double n3;

    ///
    TLorentzVector LV_n, LV_n_, LV_t, LV_t_, LV_tt_t, LV_tt_t_;
    /// provisional
    TLorentzVector genLV_n, genLV_n_;

    TRandom3* rand3;
//    JetResolution* ptResol;
//    JetResolution* etaResol;
//    JetResolution* phiResol;
    double rho;
    JME::JetResolution* ptResol;

    JME::JetResolution* phiResol;
//    JME::JetResolutionScaleFactor* SF;

    double weightmax=0;

    double get_top_pt_prob(const double pt) const;
    double get_2bjet_prob(const TLorentzVector &jet1, const TLorentzVector &jet2, map<double, double> &mapJetPhi2Discr) const;
   string getEnvVar(const string &key) const;
};

#endif // TTAMWTSOLVER_H
