// Scorer for TsEnergyDeposit
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

#include "TsEnergyDeposit.hh"

#include "TsTrackInformation.hh"

#include "G4PSDirectionFlag.hh"
#include "G4VProcess.hh"

TsEnergyDeposit::TsEnergyDeposit(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
		G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
	: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
	//SetSurfaceScorer();
	pM->SetNeedsSteppingAction();

	fNtuple->RegisterColumnF(&fEnergyDeposit, "EnergyDeposit", "MeV");
	//fNtuple->RegisterColumnF(&fWeight, "Weight", "");
	//fNtuple->RegisterColumnI(&fPType, "PDG");
	//fNtuple->RegisterColumnI(&fTrackID_Current, "TrackID");
	//fNtuple->RegisterColumnI(&fEventID, "EventID");

	fTrackID_Current = 1;
	fTrackID_Past = 1;
}


TsEnergyDeposit::~TsEnergyDeposit() {;}


G4bool TsEnergyDeposit::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
	if (!fIsActive) {
		fSkippedWhileInactive++;
		return false;
	}

	ResolveSolid(aStep);

	G4double edep = aStep->GetTotalEnergyDeposit();


	fEnergyDeposit += edep;


	return true;
}

void TsEnergyDeposit::UserHookForEndOfEvent()
{
	fNtuple->Fill();
	fEnergyDeposit = 0.;
}
