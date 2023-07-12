// Scorer for TsTrackLength
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

#include "TsTrackLength.hh"

#include "TsTrackInformation.hh"

#include "G4PSDirectionFlag.hh"
#include "G4VProcess.hh"

TsTrackLength::TsTrackLength(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                          G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
                         : TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
    //SetSurfaceScorer();
	pM->SetNeedsSteppingAction();
	
	fNtuple->RegisterColumnF(&fTrackLength, "TrackLength", "cm");
	//fNtuple->RegisterColumnF(&fWeight, "Weight", "");
	//fNtuple->RegisterColumnI(&fPType, "PDG");
    //fNtuple->RegisterColumnI(&fTrackID_Current, "TrackID");
    //fNtuple->RegisterColumnI(&fEventID, "EventID");
    
    fTrackID_Current = 1;
    fTrackID_Past = 1;
}


TsTrackLength::~TsTrackLength() {;}


G4bool TsTrackLength::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }

    ResolveSolid(aStep);
    
    fTrackID_Current = aStep->GetTrack()->GetTrackID();
    fPType           = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
    fWeight          = aStep->GetPreStepPoint()->GetWeight();
    fEventID = GetEventID();
    
    G4double StepLength = aStep->GetStepLength();
    
    //if(fTrackID_Current != fTrackID_Past)
    //{
     //   if(fTrackLength != 0)
      //      fNtuple->Fill();
        //else return false;
        
        //fTrackID_Past = fTrackID_Current;
        //fTrackLength = 0.;
    //}
    //else
        fTrackLength += StepLength;

        
    return true;
}

void TsTrackLength::UserHookForEndOfEvent()
{
    fNtuple->Fill();
    fTrackLength = 0.;
    fTrackID_Current = 1;
    fTrackID_Past = 1;
}
