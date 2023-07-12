//
// ********************************************************************
// *                                                                  *
// *                                                                  *
// * This file was obtained from Topas MC Inc under the license       *
// * agreement set forth at http://www.topasmc.org/registration       *
// * Any use of this file constitutes full acceptance of              *
// * this TOPAS MC license agreement.                                 *
// *                                                                  *
// ********************************************************************
//

#ifndef TsActivation_hh
#define TsActivation_hh

#include "TsVNtupleScorer.hh"

class TsActivation : public TsVNtupleScorer
{
public:
    TsActivation(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
    
    virtual ~TsActivation();
    
    
private:
    struct ParticleData {
     ParticleData()
       : fCount(0), fA(0), fZ(0), fEmean(0.), fEmin(0.), fEmax(0.), fTmean(-1.) {}
     ParticleData(G4int count, G4int A, G4int Z, G4double ekin, G4double emin, G4double emax,
                  G4double meanLife)
       : fCount(count), fA(A), fZ(Z), fEmean(ekin), fEmin(emin), fEmax(emax),
         fTmean(meanLife) {}
     G4int     fCount;
     G4int     fA;
     G4int     fZ;
     G4double  fEmean;
     G4double  fEmin;
     G4double  fEmax;
     G4double  fTmean;
    };
    
public:

    G4bool ProcessHits(G4Step*,G4TouchableHistory*);
    void CountParticle(std::map<G4String,ParticleData>& DataMap, G4String name, G4int a, G4int z, G4double Ekin, G4double meanLife);
    void AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer);
    void Merge(std::map<G4String, ParticleData>& destinationMap,
                             const std::map<G4String, ParticleData>& sourceMap);
    void UserHookForEndOfRun();

private:
	// Output variables
    G4String fName;
    G4int fA;
    G4int fZ;
    G4double fEkin;
    G4double fLife;
    G4double fTime;
    G4double fTimeBirth;
    G4bool   fStable;
    G4int fTrackID;
    G4int fTrackID_old;
    const G4ParticleDefinition* fParticle;

    
    
    std::map<G4String,ParticleData> fMetastableIon;
    std::map<G4String,ParticleData> fEmergingParticles;
};
#endif
