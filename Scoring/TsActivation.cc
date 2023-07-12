// Scorer for TsActivation
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

#include "TsActivation.hh"

#include "TsTrackInformation.hh"

#include "G4PSDirectionFlag.hh"
#include "G4VProcess.hh"
#include "G4Track.hh"
#include "G4ParticleTypes.hh"
#include "G4IonTable.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"



TsActivation::TsActivation(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                          G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
                         : TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
    //SetSurfaceScorer();
	pM->SetNeedsSteppingAction();
	
    //Output file with meatsatble ions
    fNtuple->RegisterColumnS(&fName, "Name");
	fNtuple->RegisterColumnI(&fA, "A Mass");
	fNtuple->RegisterColumnI(&fZ, "Z Charge");
    fNtuple->RegisterColumnD(&fEkin, "Ekin", "MeV");
	fNtuple->RegisterColumnD(&fLife, "Mean Life", "s");
    fNtuple->RegisterColumnD(&fTimeBirth, "Time Birth", "s");
    
    fName = "XXX";
    fA = 0;
    fZ = 0;
    fLife = 0.;
    fTime = 0.;
    fEkin = 0.;
    
    fTrackID_old = -1;
}


TsActivation::~TsActivation() {;}


//In LocalRun
G4bool TsActivation::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }

    ResolveSolid(aStep);
    
    G4Track* track = aStep->GetTrack();
    fTrackID = track->GetTrackID();
    const G4ParticleDefinition* particle = track->GetParticleDefinition();
    G4String name     = particle->GetParticleName();
    G4int A        = particle->GetAtomicMass();
    G4int Z        = particle->GetAtomicNumber();
    G4double meanLife = particle->GetPDGLifeTime();
    G4double ekin     = track->GetKineticEnergy();
    G4double time     = track->GetGlobalTime();
    G4bool   stable = particle->GetPDGStable();
    
    //if ((particle->GetPDGStable())&&(ekin == 0.)) fTimeEnd = DBL_MAX;
    
    if(fTrackID != fTrackID_old) //when track is changed
    {
        fName = name;
        fA = A;
        fZ = Z;
        fLife = meanLife;
        fEkin = ekin;
        fTimeBirth = time;
        fStable = stable;
        
        // count population of ions with meanLife != 0.
        if ((G4IonTable::IsIon(particle))&&(meanLife >= 0.))
        {
            CountParticle(fMetastableIon, fName, fA, fZ, fEkin, fLife); //count old track
            fNtuple ->Fill();
        }
        else
             CountParticle(fEmergingParticles, name, A, Z, ekin, fLife); //count old track
        
        
    }
    
    fTrackID_old = fTrackID;
    
    // keep only emerging particles
    //G4StepStatus status = track->GetStep()->GetPostStepPoint()->GetStepStatus();
    //if (status == fWorldBoundary)
    //{
        // histograms: energy flow and activities of emerging particles
     //   cout << "World Boundary\n";
     //   CountParticle(fEmergingParticles, name, A, Z, ekin, time); //count old track
    //}

    return true;
}

void TsActivation::CountParticle(std::map<G4String,ParticleData>& DataMap, G4String name, G4int a, G4int z, G4double Ekin, G4double meanLife)
{
    std::map<G4String, ParticleData>::iterator it = DataMap.find(name);
    if ( it == DataMap.end()) {
      DataMap[name] = ParticleData(1, a, z, Ekin, Ekin, Ekin, meanLife);
    }
    else {
      ParticleData& data = it->second;
      data.fCount++;
      data.fEmean += Ekin;
      //update min max
      G4double emin = data.fEmin;
      if (Ekin < emin) data.fEmin = Ekin;
      G4double emax = data.fEmax;
      if (Ekin > emax) data.fEmax = Ekin;
      data.fTmean = meanLife;
    }

}

//Merge Scoring for MT mode
void TsActivation::AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer) {
    TsVNtupleScorer::AbsorbResultsFromWorkerScorer(workerScorer);

    TsActivation* workerMTScorer = dynamic_cast<TsActivation*>(workerScorer);
    std::map<G4String,ParticleData> fMetastableIon_worker = workerMTScorer->fMetastableIon;
    std::map<G4String,ParticleData> fEmergingParticles_worker = workerMTScorer->fEmergingParticles;
    

    Merge(fMetastableIon, fMetastableIon_worker);
    Merge(fEmergingParticles, fEmergingParticles_worker);

    workerMTScorer->fMetastableIon.clear();
    workerMTScorer->fEmergingParticles.clear();
}

//Merge two datamaps
void TsActivation::Merge(std::map<G4String, ParticleData>& destinationMap,
const std::map<G4String, ParticleData>& sourceMap)
{
    for ( const auto& particleData : sourceMap ) {
      G4String name = particleData.first;
      const ParticleData& localData = particleData.second;
      if ( destinationMap.find(name) == destinationMap.end()) {
        destinationMap[name]
         = ParticleData(localData.fCount,
                        localData.fA,
                        localData.fZ,
                        localData.fEmean,
                        localData.fEmin,
                        localData.fEmax,
                        localData.fTmean);
      }
      else {
        ParticleData& data = destinationMap[name];
        data.fCount += localData.fCount;
        data.fEmean += localData.fEmean;
        G4double emin = localData.fEmin;
        if (emin < data.fEmin) data.fEmin = emin;
        G4double emax = localData.fEmax;
        if (emax > data.fEmax) data.fEmax = emax;
        data.fTmean = localData.fTmean;
      }
    }

}

//Output for EndOfRun. This is the MasterRun, not LocalRun
void TsActivation::UserHookForEndOfRun()
{
    //particles count
     //
     G4cout << "\n List of generated metastable ions:" << G4endl;
        
    for ( const auto& particleData : fMetastableIon ) {
       G4String name = particleData.first;
       ParticleData data = particleData.second;
       G4int count = data.fCount;
       G4double eMean = data.fEmean/count;
       G4double eMin = data.fEmin;
       G4double eMax = data.fEmax;
       G4double meanLife = data.fTmean;
            
       G4cout << "  " << std::setw(13) << name << ": " << std::setw(7) << count
              << "  Emean = " << std::setw(13) << G4BestUnit(eMean, "Energy")
              << "\t( "  << G4BestUnit(eMin, "Energy")
              << " --> " << G4BestUnit(eMax, "Energy") << ")";
       if (meanLife >= 0.)
         G4cout << "\tmean life = " << G4BestUnit(meanLife, "Time")   << G4endl;
       else G4cout << "\tstable" << G4endl;
    }
    
    G4cout << "\n\n List of stable particle produced:" << G4endl;
        
    for ( const auto& particleData : fEmergingParticles ) {
       G4String name = particleData.first;
       ParticleData data = particleData.second;
       G4int count = data.fCount;
       G4double eMean = data.fEmean/count;
       G4double eMin = data.fEmin;
       G4double eMax = data.fEmax;
       G4double meanLife = data.fTmean;
            
       G4cout << "  " << std::setw(13) << name << ": " << std::setw(7) << count
              << "  Emean = " << std::setw(13) << G4BestUnit(eMean, "Energy")
              << "\t( "  << G4BestUnit(eMin, "Energy")
              << " --> " << G4BestUnit(eMax, "Energy") << ")";
       if (meanLife >= 0.)
         G4cout << "\tmean life = " << G4BestUnit(meanLife, "Time")   << G4endl;
       else G4cout << "\tstable" << G4endl;
    }
      

}
