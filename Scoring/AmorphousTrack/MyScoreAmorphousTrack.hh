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

#ifndef MyScoreAmorphousTrack_hh
#define MyScoreAmorphousTrack_hh

#include "TsTrackerHit_AT.hh"
#include "TsVBinnedScorer.hh"
#include<vector>
#include <unordered_map>

class G4ParticleDefinition;
class MyScoreAmorphousTrack : public TsVBinnedScorer
{
	public:
		MyScoreAmorphousTrack(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
				G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);

		virtual ~MyScoreAmorphousTrack();

		G4bool ProcessHits(G4Step*,G4TouchableHistory*);
		std::vector<G4double> Interpolate(std::vector<std::vector<G4double>>, G4double, G4bool);
		void UserHookForEndOfRun();
		void AccumulateEvent();

		void AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer);
		//void FillHisto(G4double,G4int, G4int, std::vector<double> &);
		void FillHisto(G4double, std::vector<G4double> &);
		void GetLUT(std::unordered_map<std::string, std::vector<std::vector<double>> > & lut, std::string lutfolder, std::string ztype);

	protected:
		void Output();
	private:
		TsTrackerHitsCollection* fHitsCollection;	
		std::unordered_map<int,std::vector<G4double>> fzdDSpectraBin;
		std::unordered_map<int,std::vector<G4double>> fzdDDoseAverzdD;

		std::vector<double> zBinLimit, zBinWidth, zBinCenter;
		G4double zBins, zStart, zEnd, fCountsLowerLimit;

		std::unordered_map<std::string, std::vector<std::vector<G4double>> > flut_zD;
		G4ParticleDefinition* fElectronDefinition;
		G4int fStepCount;
		G4String foutfilename;
		G4String fTypeOfTable;	

};
#endif


