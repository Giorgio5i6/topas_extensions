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

#ifndef TsEnergyDeposit_hh
#define TsEnergyDeposit_hh

#include "TsVNtupleScorer.hh"

class TsEnergyDeposit : public TsVNtupleScorer
{
public:
    TsEnergyDeposit(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
    
    virtual ~TsEnergyDeposit();

    G4bool ProcessHits(G4Step*,G4TouchableHistory*);
    void UserHookForEndOfEvent();

private:
	// Output variables
	G4float fEnergyDeposit;
	G4float fWeight;
	G4int fPType;
    G4int fTrackID_Current;
    G4int fTrackID_Past;
    G4int fEventID;
};
#endif
